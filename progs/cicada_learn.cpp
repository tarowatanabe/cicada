// learning from hypergraphs...
//
// we assume two inputs, one for partition and the other for marginals
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

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

#include "cicada_impl.hpp"

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


path_type forest_path;
path_type intersected_path;

path_type output_path = "-";

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

int threads = 1;

int debug = 0;


void options(int argc, char** argv);
void optimize(const path_type& forest_path,
	      const path_type& intersected_path,
	      weight_set_type& weights);
void enumerate_forest(const path_type& path);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");

    if (! boost::filesystem::exists(forest_path) || ! boost::filesystem::is_directory(forest_path))
      throw std::runtime_error("no forest?");

    if (! boost::filesystem::exists(intersected_path) || ! boost::filesystem::is_directory(intersected_path))
      throw std::runtime_error("no intersected forest?");
    
    threads = utils::bithack::max(1, threads);

    weight_set_type weights;
    
    optimize(forest_path, intersected_path, weights);
    
    utils::compress_ostream os(output_path, 1024 * 1024);
    os << weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct OptimizeLBFGS
{

  OptimizeLBFGS(const path_type& __forest_path,
		const path_type& __intersected_path,
		weight_set_type& __weights)
    : forest_path(__forest_path),
      intersected_path(__intersected_path),
      weights(__weights) {}

  void operator()()
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
  }

  
  struct Task
  {
    
    void operator()()
    {
      
      
      
    }
    
  };

  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    OptimizeLBFGS& optimizer = *((OptimizeLBFGS*) instance);
    
    double objective = 0.0;
    std::fill(g, g + n, 0.0);
    
    static const feature_type featuer_none;
    
    
    
  }
  
  path_type forest_path;
  path_type intersected_path;
  weight_set_type& weights;
};

void optimize(const path_type& forest_path,
	      const path_type& intersected_path,
	      weight_set_type& weights)
{
  OptimizeLBFGS(forest_path, intersected_path, weights)();
}

// I know, it is very stupid...

struct TaskEnumerate
{
  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
  
  TaskEnumerate(queue_type& __queue)
    : queue(__queue) {}
  
  void operator()()
  {
    path_type       path;
    
    int             id_graph;
    std::string     sep;
    hypergraph_type graph;
    
    while (1) {
      queue.pop_swap(path);
      if (path.empty()) break;
      
      utils::compress_istream is(path);
      
      is >> id_graph >> sep >> graph;
      if (sep != "|||")
	throw std::runtime_error("invalid format");
    }
  }

  queue_type& queue;
};


void enumerate_forest(const path_type& dir)
{
  typedef TaskEnumerate task_type;
  typedef task_type::queue_type queue_type;
  
  queue_type queue;
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue)));
  
  for (int i = 0; /**/; ++ i) {
    path_type path = dir / (boost::lexical_cast<std::string>(i) + ".gz");
    if (! boost::filesystem::exists(path)) break;
    
    queue.push_swap(path);
  }
  
  for (int i = 0; i < threads; ++ i)
    queue.push(path_type());
  
  workers.join_all();
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
