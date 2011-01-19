//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// learning from hypergraphs...
//
// we assume two inputs, one for partition and the other for marginals
//

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <stdexcept>

#include "cicada_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>

#include "lbfgs.h"

path_type forest_path;
path_type intersected_path;

path_type output_path = "-";

int iteration = 100;
bool learn_sgd = false;
bool learn_maxent = false;
bool learn_mira = false;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

int threads = 1;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);


size_t enumerate_forest(const path_type& path);

template <typename Optimizer, typename Generator>
double optimize_online(const size_t& instances,
		       const path_type& forest_path,
		       const path_type& intersected_path,
		       weight_set_type& weights,
		       Generator& generator);
template <typename Optimizer>
double optimize_batch(const size_t& instances,
		      const path_type& forest_path,
		      const path_type& intersected_path,
		      weight_set_type& weights);

struct OptimizeLBFGS;

template <typename Optimizer, typename Generator>
struct OptimizeOnline;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(learn_maxent) + learn_sgd > 1)
      throw std::runtime_error("eitehr learn-{maxent,sgd}");
    if (int(learn_maxent) + learn_sgd == 0)
      learn_maxent = true;
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");

    if (! boost::filesystem::exists(forest_path) || ! boost::filesystem::is_directory(forest_path))
      throw std::runtime_error("no forest?");

    if (! boost::filesystem::exists(intersected_path) || ! boost::filesystem::is_directory(intersected_path))
      throw std::runtime_error("no intersected forest?");
    
    threads = utils::bithack::max(1, threads);

    const size_t instances = enumerate_forest(forest_path);
    
    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    weight_set_type weights;
    double objective = 0.0;
    
    typedef boost::random_number_generator<boost::mt19937> generator_type;
    boost::mt19937 gen;
    gen.seed(time(0) * getpid());
    generator_type generator(gen);

    if (learn_sgd) {
      if (regularize_l1)
	objective = optimize_online<OptimizeOnline<OptimizerSGDL1, generator_type> >(instances, forest_path, intersected_path, weights, generator);
      else
	objective = optimize_online<OptimizeOnline<OptimizerSGDL2, generator_type> >(instances, forest_path, intersected_path, weights, generator);
    } else
      objective = optimize_batch<OptimizeLBFGS>(instances, forest_path, intersected_path, weights);
    
    if (debug)
      std::cerr << "objective: " << objective << std::endl;
    
    utils::compress_ostream os(output_path, 1024 * 1024);
    os.precision(20);
    os << weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Optimizer, typename Generator>
struct OptimizeOnline
{
  typedef Optimizer optimizer_type;
  typedef Generator generator_type;
  typedef std::vector<optimizer_type, std::allocator<optimizer_type> > optimizer_set_type;
  
  OptimizeOnline(const size_t& __instances,
		 const path_type& __forest_path,
		 const path_type& __intersected_path,
		 weight_set_type& __weights,
		 generator_type& __generator)
    : instances(__instances), 
      forest_path(__forest_path),
      intersected_path(__intersected_path),
      weights(__weights),
      generator(__generator) {}
  
  struct Task
  {
    typedef std::pair<path_type, path_type> path_pair_type;

    typedef utils::lockfree_list_queue<path_pair_type, std::allocator<path_pair_type> > queue_type;
    typedef Optimizer optimizer_type;

    typedef typename optimizer_type::weight_type   weight_type;
    typedef typename optimizer_type::gradient_type gradient_type;    
    
    Task(queue_type& __queue,
	 optimizer_type& __optimizer)
      : queue(__queue), optimizer(__optimizer) {}

    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
    
    struct gradients_type
    {
      typedef gradient_type value_type;

      template <typename Index>
      gradient_type& operator[](Index)
      {
	return gradient;
      }

      void clear() { gradient.clear(); }
      
      gradient_type gradient;
    };
    
    struct weight_function
    {
      typedef weight_type value_type;
      
      weight_function(const weight_set_type& __weights, const double& __scale) : weights(__weights), scale(__scale) {}
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e
	return cicada::semiring::traits<value_type>::log(edge.features.dot(weights) * scale);
      }
      
      const weight_set_type& weights;
      const double scale;
    };

    struct feature_function
    {
      typedef gradient_type value_type;

      feature_function(const weight_set_type& __weights, const double& __scale) : weights(__weights), scale(__scale) {}

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e r_e
	gradient_type grad;
	
	const weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights) * scale);
	
	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0)
	    grad[fiter->first] = weight_type(fiter->second) * weight;
	
	return grad;
      }
      
      const weight_set_type& weights;
      const double scale;
    };
    
    void operator()()
    {
      size_t id_forest;
      size_t id_intersected;
      hypergraph_type hypergraph;

      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;

      optimizer.initialize();
      
      path_pair_type paths;
      while (1) {
	queue.pop(paths);
	if (paths.first.empty()) break;
	
	gradients.clear();
	gradients_intersected.clear();
	
	inside.clear();
	inside_intersected.clear();
	
	id_forest = read_forest(paths.first, hypergraph);
	const bool valid_forest = hypergraph.is_valid();
	
	if (valid_forest) {
	  inside.reserve(hypergraph.nodes.size());
	  inside.resize(hypergraph.nodes.size(), weight_type());
	  
	  cicada::inside_outside(hypergraph, inside, gradients,
				 weight_function(optimizer.weights, optimizer.weight_scale),
				 feature_function(optimizer.weights, optimizer.weight_scale));
	}
	
	id_intersected = read_forest(paths.second, hypergraph);
	const bool valid_intersected = hypergraph.is_valid();
	
	if (valid_intersected) {
	  inside_intersected.reserve(hypergraph.nodes.size());
	  inside_intersected.resize(hypergraph.nodes.size(), weight_type());
	  
	  cicada::inside_outside(hypergraph, inside_intersected, gradients_intersected,
				 weight_function(optimizer.weights, optimizer.weight_scale),
				 feature_function(optimizer.weights, optimizer.weight_scale));
	}
	
	if (id_forest != id_intersected)
	  throw std::runtime_error("different segment id?");
	
	if (valid_forest && valid_intersected) {
	  gradient_type& gradient = gradients.gradient;
	  weight_type& Z = inside.back();
	  
	  gradient_type& gradient_intersected = gradients_intersected.gradient;
	  weight_type& Z_intersected = inside_intersected.back();
	  
	  gradient /= Z;
	  gradient_intersected /= Z_intersected;
	  
	  optimizer(gradients_intersected.gradient,
		    gradients.gradient,
		    Z_intersected,
		    Z);
	}
      }
      
      optimizer.finalize();
    }

    size_t read_forest(const path_type& path, hypergraph_type& hypergraph)
    {
      size_t id = 0;
      std::string line;
      
      utils::compress_istream is(path, 1024 * 1024);
      std::getline(is, line);
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! parse_id(id, iter, end))
	throw std::runtime_error("invalid id input");
      
      if (! hypergraph.assign(iter, end))
	throw std::runtime_error("invalid graph format");
      
      if (iter != end)
	throw std::runtime_error("invalid id ||| graph format");
      
      return id;
    }
    
    queue_type&     queue;
    optimizer_type& optimizer;
    
  };
  
  double operator()()
  {
    typedef Task task_type;
    typedef typename task_type::queue_type queue_type;
    typedef typename task_type::path_pair_type path_pair_type;
    
    typedef std::vector<path_pair_type, std::allocator<path_pair_type> > path_pair_set_type;
    
    optimizer_set_type optimizers(threads, optimizer_type(instances, C));
    
    queue_type queue;

    path_pair_set_type paths;
    for (int sample = 0; /**/; ++ sample) {
      const std::string file_name = boost::lexical_cast<std::string>(sample) + ".gz";
      
      const path_type path_forest      = forest_path / file_name;
      const path_type path_intersected = intersected_path / file_name;
      
      if (! boost::filesystem::exists(path_forest)) break;
      if (! boost::filesystem::exists(path_intersected)) continue;
      
      paths.push_back(std::make_pair(path_forest, path_intersected));
    }
    
    weight_set_type weights_mixed;
    double objective = 0.0;
    
    for (int iter = 0; iter < iteration; ++ iter) {
      
      boost::thread_group workers;
      for (int i = 0; i < threads; ++ i)
	workers.add_thread(new boost::thread(task_type(queue, optimizers[i])));
      
      for (size_t pos = 0; pos != paths.size(); ++ pos)
	queue.push(paths[pos]);
      
      for (int i = 0; i < threads; ++ i)
	queue.push(path_pair_type());
      
      std::random_shuffle(paths.begin(), paths.end(), generator);
      
      workers.join_all();
      
      // collect weights from optimizers and perform averaging...
      weights_mixed.clear();
      size_t samples = 0;
      
      typename optimizer_set_type::iterator oiter_end = optimizers.end();
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
	oiter->weights *= oiter->samples;
	
	weights_mixed += oiter->weights;
	samples       += oiter->samples;
      }
      
      weights_mixed *= (1.0 / samples);
      
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter)
	oiter->weights = weights_mixed;
      
      objective = 0.0;
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter)
	objective += oiter->objective;
      if (debug >= 2)
	std::cerr << "objective: " << objective << std::endl;
    }
    
    weights.swap(weights_mixed);
    
    return objective;
  }
  
  size_t instances;
  path_type forest_path;
  path_type intersected_path;
  weight_set_type& weights;
  generator_type& generator;
};

template <typename Optimizer, typename Generator>
double optimize_online(const size_t& instances,
		       const path_type& forest_path,
		       const path_type& intersected_path,
		       weight_set_type& weights,
		       Generator& generator)
{
  return Optimizer(instances, forest_path, intersected_path,  weights, generator)();
}


struct OptimizeLBFGS
{
  
  OptimizeLBFGS(const size_t& instances,
		const path_type& __forest_path,
		const path_type& __intersected_path,
		weight_set_type& __weights)
    : forest_path(__forest_path),
      intersected_path(__intersected_path),
      weights(__weights) {}

  double operator()()
  {
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    
    if (regularize_l1) {
      param.orthantwise_c = C;
      param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    } else
      param.orthantwise_c = 0.0;
    
    param.max_iterations = iteration;
    
    double objective = 0.0;
    weights.clear();
    weights.allocate();
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    return objective;
  }

  
  struct Task
  {
    typedef std::pair<path_type, path_type> path_pair_type;
    typedef utils::lockfree_list_queue<path_pair_type, std::allocator<path_pair_type> > queue_type;

    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_static_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    Task(queue_type&            __queue,
	 const weight_set_type& __weights)
      : queue(__queue), weights(__weights) {}
    
    struct gradients_type
    {
      typedef gradient_type value_type;

      template <typename Index>
      gradient_type& operator[](Index)
      {
	return gradient;
      }

      void clear() { gradient.clear(); }
      
      gradient_type gradient;
    };
    
    struct weight_function
    {
      typedef weight_type value_type;

      weight_function(const weight_set_type& __weights) : weights(__weights) {}
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e
	return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
      }
      
      const weight_set_type& weights;
    };

    struct feature_function
    {
      typedef gradient_type value_type;

      feature_function(const weight_set_type& __weights) : weights(__weights) {}

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e r_e
	gradient_type grad;
	
	const weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights));
	
	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0)
	    grad[fiter->first] = weight_type(fiter->second) * weight;
	
	return grad;
      }
      
      const weight_set_type& weights;
    };


    size_t read_forest(const path_type& path, hypergraph_type& hypergraph)
    {
      size_t id = 0;
      std::string line;
      
      utils::compress_istream is(path, 1024 * 1024);
      std::getline(is, line);
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! parse_id(id, iter, end))
	throw std::runtime_error("invalid id input");
      
      if (! hypergraph.assign(iter, end))
	throw std::runtime_error("invalid graph format");
      
      if (iter != end)
	throw std::runtime_error("invalid id ||| graph format");

      return id;
    }
    
    void operator()()
    {
      path_pair_type paths;

      size_t id_forest;
      size_t id_intersected;
      hypergraph_type hypergraph;
      
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;
      
      gradient_static_type  feature_expectations;

      g.clear();
      objective = 0.0;
      
      while (1) {
	queue.pop(paths);
	if (paths.first.empty()) break;
	
	gradients.clear();
	gradients_intersected.clear();
	
	inside.clear();
	inside_intersected.clear();
	
	id_forest = read_forest(paths.first, hypergraph);
	
	const bool valid_forest = hypergraph.is_valid();
	
	if (valid_forest) {
	  inside.reserve(hypergraph.nodes.size());
	  inside.resize(hypergraph.nodes.size(), weight_type());
	  
	  inside_outside(hypergraph, inside, gradients, weight_function(weights), feature_function(weights));
	}
	
	id_intersected = read_forest(paths.second, hypergraph);
	const bool valid_intersected = hypergraph.is_valid();
	
	if (valid_intersected) {
	  inside_intersected.reserve(hypergraph.nodes.size());
	  inside_intersected.resize(hypergraph.nodes.size(), weight_type());
	  
	  inside_outside(hypergraph, inside_intersected, gradients_intersected, weight_function(weights), feature_function(weights));
	}
	
	if (id_forest != id_intersected)
	  throw std::runtime_error("different segment id?");
	
	if (valid_forest && valid_intersected) {
	  gradient_type& gradient = gradients.gradient;
	  weight_type& Z = inside.back();
	  
	  gradient_type& gradient_intersected = gradients_intersected.gradient;
	  weight_type& Z_intersected = inside_intersected.back();
	  
	  gradient /= Z;
	  gradient_intersected /= Z_intersected;
	  
	  feature_expectations -= gradient_intersected;
	  feature_expectations += gradient;

	  const double margin = log(Z_intersected) - log(Z);
	  
	  objective -= margin;

	  if (debug >= 3)
	    std::cerr << "id: " << id_forest << " margin: " << margin << std::endl;
	}
      }
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(feature_expectations.begin(), feature_expectations.end(), g.begin());
    }

    queue_type&            queue;
    
    const weight_set_type& weights;
    
    double          objective;
    weight_set_type g;
  };

  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    typedef Task                  task_type;
    typedef task_type::queue_type queue_type;
    
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

    OptimizeLBFGS& optimizer = *((OptimizeLBFGS*) instance);
        
    queue_type queue(1024 * threads);

#if 0
    for (feature_type::id_type id = 0; id < optimizer.weights.size(); ++ id)
      if (! feature_type(id).empty() && optimizer.weights[id] != 0.0)
	std::cerr << feature_type(id) << ": " << optimizer.weights[id] << std::endl;
#endif
    
    task_set_type tasks(threads, task_type(queue, optimizer.weights));
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    for (int sample = 0; /**/; ++ sample) {
      const std::string file_name = boost::lexical_cast<std::string>(sample) + ".gz";
      
      const path_type path_forest      = optimizer.forest_path / file_name;
      const path_type path_intersected = optimizer.intersected_path / file_name;
      
      if (! boost::filesystem::exists(path_forest)) break;
      if (! boost::filesystem::exists(path_intersected)) continue;
      
      queue.push(std::make_pair(path_forest, path_intersected));
    }
    
    // collect all the objective and gradients...
    double objective = 0.0;
    std::fill(g, g + n, 0.0);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(std::make_pair(path_type(), path_type()));
    
    workers.join_all();
    
    for (int i = 0; i < threads; ++ i) {
      objective += tasks[i].objective;
      std::transform(tasks[i].g.begin(), tasks[i].g.end(), g, g, std::plus<double>());
    }
        
    // L2...
    if (regularize_l2) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i) {
	g[i] += C * x[i];
	norm += x[i] * x[i];
      }
      objective += 0.5 * C * norm;
    }
    
    if (debug >= 2)
      std::cerr << "objective: " << objective << std::endl;
    
    return objective;
  }
  
  path_type forest_path;
  path_type intersected_path;
  weight_set_type& weights;
};

template <typename Optimizer>
double optimize_batch(const size_t& instances,
		      const path_type& forest_path,
		      const path_type& intersected_path,
		      weight_set_type& weights)
{
  return Optimizer(instances, forest_path, intersected_path, weights)();
}

// I know, it is very stupid...

struct TaskEnumerate
{
  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
  
  TaskEnumerate(queue_type& __queue,
		size_t& __size)
    : queue(__queue),
      size(__size) {}
  
  void operator()()
  {
    path_type       path;
    
    size_t          id;
    hypergraph_type graph;

    std::string line;
    
    while (1) {
      queue.pop_swap(path);
      if (path.empty()) break;
      
      utils::compress_istream is(path);
      std::getline(is, line);

      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();

      if (! parse_id(id, iter, end))
	throw std::runtime_error("invalid id input");
	  
      if (! graph.assign(iter, end))
	throw std::runtime_error("invalid graph format");

      if (iter != end)
	throw std::runtime_error("invalid id ||| graph format");

      ++ size;
    }
  }
  
  queue_type& queue;
  size_t& size;
};


size_t enumerate_forest(const path_type& dir)
{
  typedef TaskEnumerate task_type;
  typedef task_type::queue_type queue_type;
  
  
  queue_type queue;
  std::vector<size_t, std::allocator<size_t> > sizes(threads, 0);
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, sizes[i])));
  
  for (int i = 0; /**/; ++ i) {
    path_type path = dir / (boost::lexical_cast<std::string>(i) + ".gz");
    if (! boost::filesystem::exists(path)) break;
    
    queue.push_swap(path);
  }
  
  for (int i = 0; i < threads; ++ i)
    queue.push(path_type());
  
  workers.join_all();
  
  return std::accumulate(sizes.begin(), sizes.end(), 0);
}



void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("forest",      po::value<path_type>(&forest_path),       "forest path")
    ("intersected", po::value<path_type>(&intersected_path),  "intersected forest path")
    ("output",      po::value<path_type>(&output_path),       "output parameter")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-maxent",  po::bool_switch(&learn_maxent),  "batch LBFGS algorithm")
    ("learn-sgd",     po::bool_switch(&learn_sgd),     "online SGD algorithm")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C"            , po::value<double>(&C),           "regularization constant")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  desc_command.add(opts_command);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options] [operations]\n"
	      << opts_command << std::endl;
    exit(0);
  }
}
