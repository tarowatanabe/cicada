// learning from hypergraphs...
//
// we assume two inputs, one for partition and the other for marginals
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"

#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/base64.hpp"
#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

#include "lbfgs.h"

typedef boost::filesystem::path path_type;

typedef cicada::Symbol          symbol_type;
typedef cicada::Vocab           vocab_type;
typedef cicada::Sentence        sentence_type;
typedef cicada::Lattice         lattice_type;
typedef cicada::Rule            rule_type;
typedef cicada::HyperGraph      hypergraph_type;

typedef rule_type::feature_set_type    feature_set_type;
typedef feature_set_type::feature_type feature_type;
typedef cicada::WeightVector<double>   weight_set_type;

typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;


path_type forest_path;
path_type intersected_path;

path_type output_path = "-";

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

bool load = false;

int debug = 0;


void options(int argc, char** argv);
double optimize(const path_type& forest_path,
		const path_type& intersected_path,
		const hypergraph_set_type& graphs_forest,
		const hypergraph_set_type& graphs_intersected,
		weight_set_type& weights);
void enumerate_forest(const path_type& forest_path,
		      const path_type& intersected_path,
		      hypergraph_set_type& graphs_forest,
		      hypergraph_set_type& graphs_intersected);


void bcast_weights(const int rank, weight_set_type& weights);
void send_weights(const weight_set_type& weights);
void reduce_weights(weight_set_type& weights);

enum {
  weights_tag = 1000,
  gradients_tag,
  notify_tag,
  termination_tag,
};

inline
int loop_sleep(bool found, int non_found_iter)
{
  if (! found) {
    boost::thread::yield();
    ++ non_found_iter;
  } else
    non_found_iter = 0;
    
  if (non_found_iter >= 64) {
    struct timespec tm;
    tm.tv_sec = 0;
    tm.tv_nsec = 2000001; // above 2ms
    nanosleep(&tm, NULL);
    
    non_found_iter = 0;
  }
  return non_found_iter;
}

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  try {
    options(argc, argv);
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");

    if (! boost::filesystem::exists(forest_path) || ! boost::filesystem::is_directory(forest_path))
      throw std::runtime_error("no forest?");

    if (! boost::filesystem::exists(intersected_path) || ! boost::filesystem::is_directory(intersected_path))
      throw std::runtime_error("no intersected forest?");

    hypergraph_set_type graphs_forest;
    hypergraph_set_type graphs_intersected;
        
    enumerate_forest(forest_path, intersected_path, graphs_forest, graphs_intersected);

    if (debug && mpi_rank == 0)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    weight_set_type weights;
    
    const double objective = optimize(forest_path, intersected_path, graphs_forest, graphs_intersected, weights);

    if (debug && mpi_rank == 0)
      std::cerr << "objective: " << objective << std::endl;
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_path, 1024 * 1024);
      os.precision(20);
      os << weights;
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    MPI::Abort();
    return 1;
  }
  return 0;
}

struct OptimizeLBFGS
{

  OptimizeLBFGS(const path_type& __forest_path,
		const path_type& __intersected_path,
		const hypergraph_set_type& __graphs_forest,
		const hypergraph_set_type& __graphs_intersected,
		weight_set_type& __weights)
    : forest_path(__forest_path),
      intersected_path(__intersected_path),
      graphs_forest(__graphs_forest),
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
    weights.clear();
    weights.allocate();
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    return objective;
  }

  
  struct Task
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_static_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    Task(const weight_set_type& __weights,
	 const path_type& __dir_forest,
	 const path_type& __dir_intersected,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected)
      : weights(__weights),
	dir_forest(__dir_forest),
	dir_intersected(__dir_intersected),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected) {}
    
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
    
    void operator()()
    {
      const int mpi_rank = MPI::COMM_WORLD.Get_rank();
      const int mpi_size = MPI::COMM_WORLD.Get_size();

      size_t id_forest;
      size_t id_intersected;
      hypergraph_type hypergraph;
      
      std::string line;
      
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;
      
      gradient_static_type  feature_expectations;

      g.clear();
      objective = 0.0;
      
      if (! graphs_forest.empty()) {
	for (int i = 0; i < graphs_forest.size(); ++ i) {
	  const hypergraph_type& hypergraph_forest      = graphs_forest[i];
	  const hypergraph_type& hypergraph_intersected = graphs_intersected[i];

	  if (! hypergraph_forest.is_valid() || ! hypergraph_intersected.is_valid()) continue;
	  
	  gradients.clear();
	  gradients_intersected.clear();
	  
	  inside.clear();
	  inside_intersected.clear();

	  inside.reserve(hypergraph_forest.nodes.size());
	  inside.resize(hypergraph_forest.nodes.size(), weight_type());
	  inside_outside(hypergraph_forest, inside, gradients, weight_function(weights), feature_function(weights));
	  
	  inside_intersected.reserve(hypergraph_intersected.nodes.size());
	  inside_intersected.resize(hypergraph_intersected.nodes.size(), weight_type());
	  inside_outside(hypergraph_intersected, inside_intersected, gradients_intersected, weight_function(weights), feature_function(weights));
	  
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
	    std::cerr << "margin: " << margin << std::endl;
	}
	
      } else {
	for (int i = mpi_rank;/**/; i += mpi_size) {
	  const std::string file_name = boost::lexical_cast<std::string>(i) + ".gz";
	
	  const path_type path_forest      = dir_forest / file_name;
	  const path_type path_intersected = dir_intersected / file_name;
	
	  if (! boost::filesystem::exists(path_forest)) break;
	  if (! boost::filesystem::exists(path_intersected)) continue;
	
	  gradients.clear();
	  gradients_intersected.clear();
	
	  inside.clear();
	  inside_intersected.clear();

	  bool valid_forest = true;
	  bool valid_intersected = true;

	  {
	    utils::compress_istream is(path_forest);
	    std::getline(is, line);

	    std::string::const_iterator iter = line.begin();
	    std::string::const_iterator end = line.end();
	  
	    if (! parse_id(id_forest, iter, end))
	      throw std::runtime_error("invalid id input");
	  
	    if (! hypergraph.assign(iter, end))
	      throw std::runtime_error("invalid graph format");
	  
	    if (iter != end)
	      throw std::runtime_error("invalid id ||| graph format");
	  
	    valid_forest = hypergraph.is_valid();
	  }
	
	  if (valid_forest) {
	    inside.reserve(hypergraph.nodes.size());
	    inside.resize(hypergraph.nodes.size(), weight_type());
	  
	    inside_outside(hypergraph, inside, gradients, weight_function(weights), feature_function(weights));
	  }
	
	  {
	    utils::compress_istream is(path_intersected);
	    std::getline(is, line);

	    std::string::const_iterator iter = line.begin();
	    std::string::const_iterator end = line.end();
	  
	    if (! parse_id(id_intersected, iter, end))
	      throw std::runtime_error("invalid id input");
	  
	    if (! hypergraph.assign(iter, end))
	      throw std::runtime_error("invalid graph format");
	  
	    if (iter != end)
	      throw std::runtime_error("invalid id ||| graph format");
	  
	    valid_intersected = hypergraph.is_valid();
	  }
	
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
      }
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(feature_expectations.begin(), feature_expectations.end(), g.begin());
    }
    
    const weight_set_type& weights;
    const path_type        dir_forest;
    const path_type        dir_intersected;

    const hypergraph_set_type& graphs_forest;
    const hypergraph_set_type& graphs_intersected;
    
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

    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    OptimizeLBFGS& optimizer = *((OptimizeLBFGS*) instance);
    
    // send notification!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
        
    bcast_weights(0, optimizer.weights);
    
    task_type task(optimizer.weights, optimizer.forest_path, optimizer.intersected_path, optimizer.graphs_forest, optimizer.graphs_intersected);
    task();
    
    // collect all the objective and gradients...
    
    reduce_weights(task.g);
    
    std::fill(g, g + n, 0.0);
    std::transform(task.g.begin(), task.g.end(), g, g, std::plus<double>());
    
    double objective = task.objective;
    MPI::COMM_WORLD.Reduce(&task.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);

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
  
  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
};

double optimize(const path_type& forest_path,
		const path_type& intersected_path,
		const hypergraph_set_type& graphs_forest,
		const hypergraph_set_type& graphs_intersected,
		weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    const double objective = OptimizeLBFGS(forest_path, intersected_path, graphs_forest, graphs_intersected, weights)();
    
    // send termination!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, termination_tag);
    
    return objective;
  } else {

    enum {
      NOTIFY = 0,
      TERMINATION,
    };
    
    MPI::Prequest requests[2];

    requests[NOTIFY]      = MPI::COMM_WORLD.Recv_init(0, 0, MPI::INT, 0, notify_tag);
    requests[TERMINATION] = MPI::COMM_WORLD.Recv_init(0, 0, MPI::INT, 0, termination_tag);
    
    for (int i = 0; i < 2; ++ i)
      requests[i].Start();

    while (1) {
      if (MPI::Request::Waitany(2, requests))
	break;
      else {
	typedef OptimizeLBFGS::Task task_type;

	requests[NOTIFY].Start();

	bcast_weights(0, weights);
	
	task_type task(weights, forest_path, intersected_path, graphs_forest, graphs_intersected);
	task();
	
	send_weights(task.g);
	
	double objective = 0.0;
	MPI::COMM_WORLD.Reduce(&task.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    return 0.0;
  }
}

// I know, it is very stupid...

void enumerate_forest(const path_type& forest_path,
		      const path_type& intersected_path,
		      hypergraph_set_type& graphs_forest,
		      hypergraph_set_type& graphs_intersected)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  size_t id_forest;
  size_t id_intersected;
  
  std::string line;
  
  for (int i = mpi_rank; /**/; i += mpi_size) {
    const std::string file_name = boost::lexical_cast<std::string>(i) + ".gz";
    
    const path_type path_forest      = forest_path / file_name;
    const path_type path_intersected = intersected_path / file_name;
    
    if (! boost::filesystem::exists(path_forest)) break;
    if (! boost::filesystem::exists(path_intersected)) continue;
    
    {
      utils::compress_istream is(path_forest);
      std::getline(is, line);
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! parse_id(id_forest, iter, end))
	throw std::runtime_error("invalid id input: " + path_forest.file_string());
      
      graphs_forest.push_back(hypergraph_type());
      
      if (! graphs_forest.back().assign(iter, end))
	throw std::runtime_error("invalid graph format" + path_forest.file_string());
      
      if (! load)
	graphs_forest.clear();
      
      if (iter != end)
	throw std::runtime_error("invalid id ||| graph format" + path_forest.file_string());
    }

    {
      utils::compress_istream is(path_intersected);
      std::getline(is, line);
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! parse_id(id_intersected, iter, end))
	throw std::runtime_error("invalid id input" + path_intersected.file_string());
      
      graphs_intersected.push_back(hypergraph_type());
      
      if (! graphs_intersected.back().assign(iter, end))
	throw std::runtime_error("invalid graph format" + path_intersected.file_string());
      
      if (! load)
	graphs_intersected.clear();

      
      if (iter != end)
	throw std::runtime_error("invalid id ||| graph format" + path_intersected.file_string());
    }
    
    if (id_forest != id_intersected)
      throw std::runtime_error("id do not match");
  }
  
  if (graphs_intersected.size() != graphs_forest.size())
    throw std::runtime_error("# of hypergraphs do not match");

  
  // collect features...
  for (int rank = 0; rank < mpi_size; ++ rank) {
    weight_set_type weights;
    weights.allocate();
    
    for (int id = 0; id < feature_type::allocated(); ++ id)
      if (! feature_type(feature_type::id_type(id)).empty())
	weights[feature_type(id)] = 1.0;
    
    bcast_weights(rank, weights);
  }
}

void reduce_weights(weight_set_type& weights)
{
  typedef utils::mpi_device_source            device_type;
  typedef boost::iostreams::filtering_istream stream_type;

  typedef boost::shared_ptr<device_type> device_ptr_type;
  typedef boost::shared_ptr<stream_type> stream_ptr_type;

  typedef std::vector<device_ptr_type, std::allocator<device_ptr_type> > device_ptr_set_type;
  typedef std::vector<stream_ptr_type, std::allocator<stream_ptr_type> > stream_ptr_set_type;

  typedef boost::tokenizer<utils::space_separator> tokenizer_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  device_ptr_set_type device(mpi_size);
  stream_ptr_set_type stream(mpi_size);

  for (int rank = 1; rank < mpi_size; ++ rank) {
    device[rank].reset(new device_type(rank, weights_tag, 1024 * 1024));
    stream[rank].reset(new stream_type());
    
    stream[rank]->push(boost::iostreams::gzip_decompressor());
    stream[rank]->push(*device[rank]);
  }

  std::string line;
  
  int non_found_iter = 0;
  while (1) {
    bool found = false;
    
    for (int rank = 1; rank < mpi_size; ++ rank)
      while (stream[rank] && device[rank] && device[rank]->test()) {
	if (std::getline(*stream[rank], line)) {
	  tokenizer_type tokenizer(line);
	  
	  tokenizer_type::iterator iter = tokenizer.begin();
	  if (iter == tokenizer.end()) continue;
	  std::string feature = *iter;
	  ++ iter;
	  if (iter == tokenizer.end()) continue;
	  std::string value = *iter;
	  
	  weights[feature] += utils::decode_base64<double>(value);
	} else {
	  stream[rank].reset();
	  device[rank].reset();
	}
	found = true;
      }
    
    if (std::count(device.begin(), device.end(), device_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
}


void send_weights(const weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  boost::iostreams::filtering_ostream os;
  os.push(boost::iostreams::gzip_compressor());
  os.push(utils::mpi_device_sink(0, weights_tag, 1024 * 1024));
  
  for (feature_type::id_type id = 0; id < weights.size(); ++ id)
    if (! feature_type(id).empty() && weights[id] != 0.0)
      os << feature_type(id) << ' ' << utils::encode_base64(weights[id]) << '\n';
}

void bcast_weights(const int rank, weight_set_type& weights)
{
  typedef std::vector<char, std::allocator<char> > buffer_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == rank) {
    boost::iostreams::filtering_ostream os;
    os.push(utils::mpi_device_bcast_sink(rank, 1024));
    
    static const weight_set_type::feature_type __empty;
    
    weight_set_type::const_iterator witer_begin = weights.begin();
    weight_set_type::const_iterator witer_end = weights.end();
    
    for (weight_set_type::const_iterator witer = witer_begin; witer != witer_end; ++ witer)
      if (*witer != 0.0) {
	const weight_set_type::feature_type feature(witer - witer_begin);
	if (feature != __empty)
	  os << feature << ' ' << utils::encode_base64(*witer) << '\n';
      }
  } else {
    weights.clear();
    weights.allocate();
    
    boost::iostreams::filtering_istream is;
    is.push(utils::mpi_device_bcast_source(rank, 1024));
    
    std::string feature;
    std::string value;
    
    while ((is >> feature) && (is >> value))
      weights[feature] = utils::decode_base64<double>(value);
  }
}


void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("forest",      po::value<path_type>(&forest_path),       "forest path")
    ("intersected", po::value<path_type>(&intersected_path),  "intersected forest path")
    ("output",      po::value<path_type>(&output_path),       "output parameter")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C"            , po::value<double>(&C),           "regularization constant")

    ("load", po::bool_switch(&load), "load training data on-memory")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  desc_command.add(opts_command);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {

    if (mpi_rank == 0)
      std::cout << argv[0] << " [options] [operations]\n"
		<< opts_command << std::endl;

    MPI::Finalize();
    exit(0);
  }
}
