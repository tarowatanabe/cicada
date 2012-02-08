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
#include <deque>

#include "cicada_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>

#include "lbfgs.h"

typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type forest_path;
path_set_type intersected_path;
path_set_type refset_path;

path_type weights_path;
path_set_type weights_history_path;
path_type output_path = "-";
path_type output_objective_path;

path_type bound_lower_file;
path_type bound_upper_file;

int iteration = 100;

bool learn_sgd = false;
bool learn_lbfgs = false;
bool learn_mira = false;
bool learn_xbleu = false;

bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;
double scale = 1.0;
double eta0 = 0.2;
int order = 4;

bool loss_margin = false; // margin by loss, not rank-loss
bool softmax_margin = false;

// scorers
std::string scorer_name = "bleu:order=4,exact=true";
bool scorer_list = false;

bool unite_forest = false;

int threads = 2;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

void read_forest(const path_set_type& forest_path,
		 const path_set_type& intersected_path,
		 hypergraph_set_type& graphs_forest,
		 hypergraph_set_type& graphs_intersected);

template <typename Optimizer>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights);
template <typename Optimizer, typename Generator>
double optimize_online(const hypergraph_set_type& graphs_forest,
		       const hypergraph_set_type& graphs_intersected,
		       weight_set_type& weights,
		       Generator& generator);
struct OptimizeLBFGS;

template <typename Optimizer, typename Generator>
struct OptimizeOnline;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_sgd > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,sgd}");
    if (int(learn_lbfgs) + learn_sgd == 0)
      learn_lbfgs = true;
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));

    if (forest_path.empty())
      throw std::runtime_error("no forest?");
    if (! learn_xbleu && intersected_path.empty())
      throw std::runtime_error("no intersected forest?");
    if (learn_xbleu && refset_path.empty())
      throw std::runtime_error("no reference translations?");

    if (! bound_lower_file.empty())
      if (bound_lower_file != "-" && ! boost::filesystem::exists(bound_lower_file))
	throw std::runtime_error("no lower-bound file? " + bound_lower_file.string());
    
    if (! bound_upper_file.empty())
      if (bound_upper_file != "-" && ! boost::filesystem::exists(bound_upper_file))
	throw std::runtime_error("no upper-bound file? " + bound_upper_file.string());
    
    threads = utils::bithack::max(1, threads);
    
    hypergraph_set_type graphs_forest;
    hypergraph_set_type graphs_intersected;
    
    if (! intersected_path.empty())
      read_forest(forest_path, intersected_path, graphs_forest, graphs_intersected);
    
    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    weight_set_type weights;
    if (! weights_path.empty()) {
      if (! boost::filesystem::exists(weights_path))
	throw std::runtime_error("no path? " + weights_path.string());
      
      utils::compress_istream is(weights_path, 1024 * 1024);
      is >> weights;
    }
    
    weights.allocate();
    
    weight_set_type bounds_lower;
    weight_set_type bounds_upper;
    
    if (! bound_lower_file.empty())
      read_bounds(bound_lower_file, bounds_lower, - std::numeric_limits<double>::infinity());
    
    if (! bound_upper_file.empty())
      read_bounds(bound_upper_file, bounds_upper,   std::numeric_limits<double>::infinity());

    double objective = 0.0;
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());

    if (learn_sgd) {
      if (regularize_l1)
	objective = optimize_online<OptimizeOnline<OptimizerSGDL1, boost::mt19937> >(graphs_forest, graphs_intersected, weights, generator);
      else
	objective = optimize_online<OptimizeOnline<OptimizerSGDL2, boost::mt19937> >(graphs_forest, graphs_intersected, weights, generator);
    } else
      objective = optimize_batch<OptimizeLBFGS>(graphs_forest, graphs_intersected, weights);
    
    if (debug)
      std::cerr << "objective: " << objective << std::endl;

    if (! bounds_lower.empty()) {
      const size_t weights_size = utils::bithack::min(weights.size(), bounds_lower.size());
      
      for (size_t i = 0; i != weights_size; ++ i)
	weights[i] = std::max(weights[i], bounds_lower[i]);
    }
    
    if (! bounds_upper.empty()) {
      const size_t weights_size = utils::bithack::min(weights.size(), bounds_upper.size());
      
      for (size_t i = 0; i != weights_size; ++ i)
	weights[i] = std::min(weights[i], bounds_upper[i]);
    }
    
    utils::compress_ostream os(output_path, 1024 * 1024);
    os.precision(20);
    os << weights;
    
    if (! output_objective_path.empty()) {
      utils::compress_ostream os(output_objective_path, 1024 * 1024);
      os.precision(20);
      os << objective << '\n';
    }
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
  
  OptimizeOnline(const hypergraph_set_type& __graphs_forest,
		 const hypergraph_set_type& __graphs_intersected,
		 weight_set_type& __weights,
		 generator_type& __generator)
    : graphs_forest(__graphs_forest),
      graphs_intersected(__graphs_intersected),
      weights(__weights),
      generator(__generator) {}
  
  struct Task
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    typedef Optimizer optimizer_type;

    typedef typename optimizer_type::weight_type   weight_type;
    typedef typename optimizer_type::gradient_type gradient_type;    
    
    Task(queue_type& __queue,
	 optimizer_type& __optimizer,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected)
      : queue(__queue),
	optimizer(__optimizer),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected)
    {}

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
	return cicada::semiring::traits<value_type>::exp(cicada::dot_product(edge.features, weights) * scale);
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
	
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(edge.features, weights) * scale);
	
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
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;

      optimizer.initialize();
      
      int id;
      while (1) {
	queue.pop(id);
	if (id < 0) break;
	
	gradients.clear();
	gradients_intersected.clear();
	
	inside.clear();
	inside_intersected.clear();

	if (! graphs_forest[id].is_valid() || ! graphs_intersected[id].is_valid()) continue;
	
	inside.reserve(graphs_forest[id].nodes.size());
	inside.resize(graphs_forest[id].nodes.size(), weight_type());
	
	inside_intersected.reserve(graphs_intersected[id].nodes.size());
	inside_intersected.resize(graphs_intersected[id].nodes.size(), weight_type());
	
	cicada::inside_outside(graphs_forest[id], inside, gradients,
			       weight_function(optimizer.weights, optimizer.weight_scale),
			       feature_function(optimizer.weights, optimizer.weight_scale));
	
	cicada::inside_outside(graphs_intersected[id], inside_intersected, gradients_intersected,
			       weight_function(optimizer.weights, optimizer.weight_scale),
			       feature_function(optimizer.weights, optimizer.weight_scale));
	
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
      
      optimizer.finalize();
    }
    
    queue_type&     queue;
    optimizer_type& optimizer;

    const hypergraph_set_type& graphs_forest;
    const hypergraph_set_type& graphs_intersected;
  };
  
  double operator()()
  {
    typedef Task task_type;
    typedef typename task_type::queue_type queue_type;
    typedef std::vector<int, std::allocator<int> > id_set_type;

    id_set_type ids;
    const int id_max = utils::bithack::min(graphs_forest.size(), graphs_intersected.size());
    for (int id = 0; id != id_max; ++ id)
      if (graphs_forest[id].is_valid() && graphs_intersected[id].is_valid())
	ids.push_back(id);
    
    optimizer_set_type optimizers(threads, optimizer_type(ids.size(), C));
    
    queue_type queue;
    
    weight_set_type weights_mixed;
    double objective = 0.0;
    
    for (int iter = 0; iter < iteration; ++ iter) {
      
      boost::thread_group workers;
      for (int i = 0; i < threads; ++ i)
	workers.add_thread(new boost::thread(task_type(queue, optimizers[i], graphs_forest, graphs_intersected)));
      
      for (size_t pos = 0; pos != ids.size(); ++ pos)
	queue.push(ids[pos]);
      
      for (int i = 0; i < threads; ++ i)
	queue.push(-1);
      
      boost::random_number_generator<Generator> gen(generator);
      std::random_shuffle(ids.begin(), ids.end(), gen);
      
      workers.join_all();
      
      // collect weights from optimizers and perform averaging...
      weights_mixed.clear();
      size_t samples = 0;
      
      typename optimizer_set_type::iterator oiter_end = optimizers.end();
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
	oiter->weights *= (oiter->samples + 1);
	
	weights_mixed += oiter->weights;
	samples       += (oiter->samples + 1);
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

  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
  generator_type& generator;
};

template <typename Optimizer, typename Generator>
double optimize_online(const hypergraph_set_type& graphs_forest,
		       const hypergraph_set_type& graphs_intersected,
		       weight_set_type& weights,
		       Generator& generator)
{
  return Optimizer(graphs_forest, graphs_intersected, weights, generator)();
}


struct OptimizeLBFGS
{
  
  OptimizeLBFGS(const hypergraph_set_type& __graphs_forest,
		const hypergraph_set_type& __graphs_intersected,
		weight_set_type& __weights)
    : graphs_forest(__graphs_forest),
      graphs_intersected(__graphs_intersected),
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
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    return objective;
  }

  
  struct Task
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_static_type;
    
    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    Task(queue_type&            __queue,
	 const weight_set_type& __weights,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected,
	 const size_t& __instances)
      : queue(__queue),
	weights(__weights),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected),
	instances(__instances)
    {}
    
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
	return cicada::semiring::traits<value_type>::exp(cicada::dot_product(edge.features, weights));
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
	
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(edge.features, weights));
	
	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0)
	    grad[fiter->first] = weight_type(fiter->second) * weight;
	
	return grad;
      }
      
      const weight_set_type& weights;
    };

    void operator()()
    {
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;
      
      gradient_static_type  feature_expectations;

      g.clear();
      objective = 0.0;
      
      while (1) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
	gradients.clear();
	gradients_intersected.clear();
	
	inside.clear();
	inside_intersected.clear();
	
	inside.reserve(graphs_forest[id].nodes.size());
	inside.resize(graphs_forest[id].nodes.size(), weight_type());
	
	inside_intersected.reserve(graphs_intersected[id].nodes.size());
	inside_intersected.resize(graphs_intersected[id].nodes.size(), weight_type());
	
	cicada::inside_outside(graphs_forest[id], inside, gradients, weight_function(weights), feature_function(weights));
	cicada::inside_outside(graphs_intersected[id], inside_intersected, gradients_intersected, weight_function(weights), feature_function(weights));
	
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
	  std::cerr << "id: " << id << " margin: " << margin << std::endl;
      }
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(feature_expectations.begin(), feature_expectations.end(), g.begin());
      
      // normalize!
      objective /= instances;
      std::transform(g.begin(), g.end(), g.begin(), std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    }

    queue_type&            queue;
    
    const weight_set_type& weights;
    
    const hypergraph_set_type& graphs_forest;
    const hypergraph_set_type& graphs_intersected;
    size_t instances;
    
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
    
    const int id_max = utils::bithack::min(optimizer.graphs_forest.size(), optimizer.graphs_intersected.size());
    size_t instances = 0;
    for (int id = 0; id != id_max; ++ id)
      instances += (optimizer.graphs_forest[id].is_valid() && optimizer.graphs_intersected[id].is_valid());
    
    queue_type queue;
    
    task_set_type tasks(threads, task_type(queue, optimizer.weights, optimizer.graphs_forest, optimizer.graphs_intersected, instances));

    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    for (int id = 0; id != id_max; ++ id)
      if (optimizer.graphs_forest[id].is_valid() && optimizer.graphs_intersected[id].is_valid())
	queue.push(id);
    
    // collect all the objective and gradients...
    double objective = 0.0;
    std::fill(g, g + n, 0.0);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
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
  
  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
};

template <typename Optimizer>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights)
{
  return Optimizer(graphs_forest, graphs_intersected, weights)();
}

void read_forest(const path_set_type& forest_path,
		 const path_set_type& intersected_path,
		 hypergraph_set_type& graphs_forest,
		 hypergraph_set_type& graphs_intersected)
{
  if (unite_forest) {
    size_t id_forest;
    size_t id_intersected;
    
    std::string line;
    
    hypergraph_type graph;
    
    for (path_set_type::const_iterator piter = forest_path.begin(); piter != forest_path.end(); ++ piter) {
      if (debug)
	std::cerr << "reading forest: " << piter->string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_forest = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_forest)) break;
	
	utils::compress_istream is(path_forest);
	std::getline(is, line);
	
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end  = line.end();
      
	if (! parse_id(id_forest, iter, end))
	  throw std::runtime_error("invalid id input: " + path_forest.string());
	if (id_forest != i)
	  throw std::runtime_error("invalid id input: " + path_forest.string());
      
	if (id_forest >= graphs_forest.size())
	  graphs_forest.resize(id_forest + 1);
      
	if (! graph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path_forest.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path_forest.string());
      
	graphs_forest[id_forest].unite(graph);
      }
    }
    
    graphs_intersected.resize(graphs_forest.size());
    
    for (path_set_type::const_iterator piter = intersected_path.begin(); piter != intersected_path.end(); ++ piter) {
    
      if (debug)
	std::cerr << "reading intersected forest: " << piter->string() << std::endl;

      for (size_t i = 0; i < graphs_intersected.size(); ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_intersected = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_intersected)) continue;
	
	utils::compress_istream is(path_intersected);
	std::getline(is, line);
	
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end = line.end();
	
	if (! parse_id(id_intersected, iter, end))
	  throw std::runtime_error("invalid id input" + path_intersected.string());
	if (id_intersected != i)
	  throw std::runtime_error("invalid id input: " + path_intersected.string());
      
	if (! graph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path_intersected.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path_intersected.string());
	
	graphs_intersected[id_intersected].unite(graph);
      }
    }
    
  } else {
    if (forest_path.size() != intersected_path.size())
      throw std::runtime_error("# of forest does not match");
    
    size_t id_forest;
    size_t id_intersected;
    
    std::string line;
    
    for (size_t pos = 0; pos != forest_path.size(); ++ pos) {
      if (debug)
	std::cerr << "reading forest: " << forest_path[pos].string() << " with " << intersected_path[pos].string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_forest      = forest_path[pos] / file_name;
	const path_type path_intersected = intersected_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_forest)) break;
	if (! boost::filesystem::exists(path_intersected)) continue;
	
	{
	  utils::compress_istream is(path_forest);
	  std::getline(is, line);
	  
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator end = line.end();
	  
	  if (! parse_id(id_forest, iter, end))
	    throw std::runtime_error("invalid id input: " + path_forest.string());
	  if (id_forest != i)
	    throw std::runtime_error("invalid id input: " + path_forest.string());
	  
	  graphs_forest.push_back(hypergraph_type());
	  
	  if (! graphs_forest.back().assign(iter, end))
	    throw std::runtime_error("invalid graph format" + path_forest.string());
	  if (iter != end)
	    throw std::runtime_error("invalid id ||| graph format" + path_forest.string());
	}
	
	{
	  utils::compress_istream is(path_intersected);
	  std::getline(is, line);
	  
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator end = line.end();
	  
	  if (! parse_id(id_intersected, iter, end))
	    throw std::runtime_error("invalid id input" + path_intersected.string());
	  if (id_intersected != i)
	    throw std::runtime_error("invalid id input: " + path_intersected.string());
	  
	  graphs_intersected.push_back(hypergraph_type());
	  
	  if (! graphs_intersected.back().assign(iter, end))
	    throw std::runtime_error("invalid graph format" + path_intersected.string());
	  if (iter != end)
	    throw std::runtime_error("invalid id ||| graph format" + path_intersected.string());
	}
      }
    }
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("forest",      po::value<path_set_type>(&forest_path)->multitoken(),      "forest path")
    ("intersected", po::value<path_set_type>(&intersected_path)->multitoken(), "intersected forest path")
    ("oracle",      po::value<path_set_type>(&intersected_path)->multitoken(), "oracle forest path(s) (an alias for --intersected)")
    ("refset",      po::value<path_set_type>(&refset_path)->multitoken(),      "reference translation(s)")
    ("weights",     po::value<path_type>(&weights_path),      "initial parameter")
    ("weights-history", po::value<path_set_type>(&weights_history_path)->multitoken(), "parameter history")
    ("output",      po::value<path_type>(&output_path),       "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("bound-lower", po::value<path_type>(&bound_lower_file),                     "lower bounds definition for feature weights")
    ("bound-upper", po::value<path_type>(&bound_upper_file),                     "upper bounds definition for feature weights")

    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",  po::bool_switch(&learn_lbfgs),  "batch LBFGS algorithm")
    ("learn-sgd",    po::bool_switch(&learn_sgd),    "online SGD algorithm")
    ("learn-xbleu",  po::bool_switch(&learn_xbleu),  "xBLEU algorithm")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    
    ("C",             po::value<double>(&C)->default_value(C),         "regularization constant")
    ("scale",         po::value<double>(&scale)->default_value(scale), "scaling for weight")
    ("eta0",          po::value<double>(&eta0),                        "\\eta_0 for decay")
    ("order",         po::value<int>(&order)->default_value(order),    "ngram order for xBLEU")
    
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")
    
    ("unite",    po::bool_switch(&unite_forest), "unite forest sharing the same id")

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
