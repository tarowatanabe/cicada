//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>

//
// k-best learner
//

#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"
#include "cicada_kbest_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/sgi_hash_set.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>
#include <boost/functional/hash/hash.hpp>

#include "lbfgs.h"

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type kbest_path;
path_set_type oracle_path;
path_type weights_path;
path_type output_path = "-";
path_type output_objective_path;

int iteration = 100;
bool learn_sgd = false;
bool learn_lbfgs = false;
bool learn_mira = false;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

bool unite_kbest = false;

int threads = 1;

int debug = 0;

void options(int argc, char** argv);

void read_kbest(const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles);

template <typename Optimizer>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights);

struct OptimizeLBFGS;

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
    
    if (kbest_path.empty())
      throw std::runtime_error("no kbest?");
    if (oracle_path.empty())
      throw std::runtime_error("no oracke kbest?");
    
    threads = utils::bithack::max(1, threads);
    
    hypothesis_map_type kbests;
    hypothesis_map_type oracles;
    
    read_kbest(kbest_path, oracle_path, kbests, oracles);
    
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
    
    double objective = 0.0;
    
    boost::mt19937 generator;
    generator.seed(time(0) * getpid());
    
    objective = optimize_batch<OptimizeLBFGS>(kbests, oracles, weights);
    
    if (debug)
      std::cerr << "objective: " << objective << std::endl;
    
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

struct OptimizeLBFGS
{
  
  OptimizeLBFGS(const hypothesis_map_type& __kbests,
		const hypothesis_map_type& __oracles,
		weight_set_type& __weights)
    : kbests(__kbests),
      oracles(__oracles),
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
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;

    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    
    Task(queue_type&            __queue,
	 const weight_set_type& __weights,
	 const hypothesis_map_type& __kbests,
	 const hypothesis_map_type& __oracles)
      : queue(__queue),
	weights(__weights),
	kbests(__kbests),
	oracles(__oracles)
    {}
    

    void operator()()
    {
      g.clear();
      objective = 0.0;

      expectation_type  expectations;
      
      expectations.allocate();
      expectations.clear();
      
      while (1) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
	weight_type Z_oracle;
	weight_type Z_kbest;
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter)
	  Z_oracle += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->features.begin(), oiter->features.end(), 0.0));
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter)
	  Z_kbest += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->features.begin(), kiter->features.end(), 0.0));
	
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->features.begin(), oiter->features.end(), 0.0)) / Z_oracle;
	  
	  hypothesis_type::feature_set_type::const_iterator fiter_end = oiter->features.end();
	  for (hypothesis_type::feature_set_type::const_iterator fiter = oiter->features.begin(); fiter != fiter_end; ++ fiter)
	    expectations[fiter->first] -= weight_type(fiter->second) * weight;
	}
	
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->features.begin(), kiter->features.end(), 0.0)) / Z_kbest;
	  
	  hypothesis_type::feature_set_type::const_iterator fiter_end = kiter->features.end();
	  for (hypothesis_type::feature_set_type::const_iterator fiter = kiter->features.begin(); fiter != fiter_end; ++ fiter)
	    expectations[fiter->first] += weight_type(fiter->second) * weight;
	}
	
	const double margin = log(Z_oracle) - log(Z_kbest);
	objective -= margin;
	
	if (debug >= 3)
	  std::cerr << "id: " << id << " margin: " << margin << std::endl;
      }
      
      // transform feature_expectations into g...
      
      g.allocate();
      
      std::copy(expectations.begin(), expectations.end(), g.begin());
    }

    queue_type&            queue;
    
    const weight_set_type& weights;
    
    const hypothesis_map_type& kbests;
    const hypothesis_map_type& oracles;
    
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
        
    queue_type queue;
    
    task_set_type tasks(threads, task_type(queue, optimizer.weights, optimizer.kbests, optimizer.oracles));

    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    const int id_max = utils::bithack::min(optimizer.kbests.size(), optimizer.oracles.size());
    for (int id = 0; id != id_max; ++ id)
      if (! optimizer.kbests[id].empty() && ! optimizer.oracles[id].empty())
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
  
  const hypothesis_map_type& kbests;
  const hypothesis_map_type& oracles;
  
  weight_set_type& weights;
};

template <typename Optimizer>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights)
{
  return Optimizer(kbests, oracles, weights)();
}

struct TaskUnique
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;

  TaskUnique(queue_type& __queue,
	     hypothesis_map_type& __kbests)
    : queue(__queue), kbests(__kbests) {}

#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
				  std::allocator<hypothesis_type> > hypothesis_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
			std::allocator<hypothesis_type> > hypothesis_unique_type;
#endif

  void operator()()
  {
    hypothesis_unique_type hypotheses;

    for (;;) {
      int id = 0;
      queue.pop(id);
      if (id < 0) break;
      
      hypotheses.clear();
      hypotheses.insert(kbests[id].begin(), kbests[id].end());
      
      // clear
      kbests[id].clear();
      hypothesis_set_type(kbests[id]).swap(kbests[id]);
      
      // reserve and assign
      kbests[id].reserve(hypotheses.size());
      kbests[id].insert(kbests[id].end(), hypotheses.begin(), hypotheses.end());
    }
  }
  
  queue_type& queue;
  hypothesis_map_type& kbests;
};

void unique_kbest(hypothesis_map_type& kbests)
{
  typedef TaskUnique task_type;
  typedef task_type::queue_type queue_type;

  queue_type queue;
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, kbests)));
  
  for (size_t id = 0; id != kbests.size(); ++ id)
    if (! kbests[id].empty())
      queue.push(id);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
}

void read_kbest(const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  parser_type parser;
  kbest_feature_type kbest;
  
  if (unite_kbest) {
    size_t id_kbest;
    size_t id_oracle;
    
    for (path_set_type::const_iterator piter = kbest_path.begin(); piter != kbest_path.end(); ++ piter) {
      if (debug)
	std::cerr << "reading kbest: " << piter->string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;

	if (i >= kbests.size())
	  kbests.resize(i + 1);
	
	utils::compress_istream is(path_kbest, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  const size_t& id = boost::fusion::get<0>(kbest);

	  if (id != i)
	    throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(id));
	  	  
	  kbests[i].push_back(hypothesis_type(kbest));
	}
      }
    }
    
    // we will compute unique!
    unique_kbest(kbests);
    
    oracles.resize(kbests.size());
    
    for (path_set_type::const_iterator piter = oracle_path.begin(); piter != oracle_path.end(); ++ piter) {
      if (debug)
	std::cerr << "reading oracles: " << piter->string() << std::endl;
      
      for (size_t i = 0; i < oracles.size(); ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_oracle = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_oracle)) continue;
	
	utils::compress_istream is(path_oracle, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  if (boost::fusion::get<0>(kbest) != i)
	    throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(boost::fusion::get<0>(kbest)));
	  
	  oracles[i].push_back(hypothesis_type(kbest));
	}
      }
    }
    
    // we will compute unique!
    unique_kbest(oracles);
    
  } else {
    // synchronous reading...
    if (kbest_path.size() != oracle_path.size())
      throw std::runtime_error("# of kbests does not match");
    
    for (size_t pos = 0; pos != kbest_path.size(); ++ pos) {
      if (debug)
	std::cerr << "reading kbest: " << kbest_path[pos].string() << " with " << oracle_path[pos].string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest  = kbest_path[pos] / file_name;
	const path_type path_oracle = oracle_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;
	if (! boost::filesystem::exists(path_oracle)) continue;

	kbests.resize(kbests.size() + 1);
	oracles.resize(oracles.size() + 1);

	{
	  utils::compress_istream is(path_kbest, 1024 * 1024);
	  is.unsetf(std::ios::skipws);
	  
	  iter_type iter(is);
	  iter_type iter_end;
	  
	  while (iter != iter_end) {
	    boost::fusion::get<1>(kbest).clear();
	    boost::fusion::get<2>(kbest).clear();
	    
	    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	      if (iter != iter_end)
		throw std::runtime_error("kbest parsing failed");
	    
	    if (boost::fusion::get<0>(kbest) != i)
	      throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(boost::fusion::get<0>(kbest)));
	    
	    kbests.back().push_back(hypothesis_type(kbest));
	  }
	}

	{
	  utils::compress_istream is(path_oracle, 1024 * 1024);
	  is.unsetf(std::ios::skipws);
	  
	  iter_type iter(is);
	  iter_type iter_end;
	  
	  while (iter != iter_end) {
	    boost::fusion::get<1>(kbest).clear();
	    boost::fusion::get<2>(kbest).clear();
	    
	    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	      if (iter != iter_end)
		throw std::runtime_error("kbest parsing failed");

	    if (boost::fusion::get<0>(kbest) != i)
	      throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(boost::fusion::get<0>(kbest)));
	    
	    oracles.back().push_back(hypothesis_type(kbest));
	  }
	}
      }
    }
    
    // uniques...
    unique_kbest(kbests);
    unique_kbest(oracles);
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("kbest",   po::value<path_set_type>(&kbest_path)->multitoken(),  "kbest path")
    ("oracle",  po::value<path_set_type>(&oracle_path)->multitoken(), "oracle kbest path")
    ("weights", po::value<path_type>(&weights_path),                  "initial parameter")
    ("output",  po::value<path_type>(&output_path),                   "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",  po::bool_switch(&learn_lbfgs),  "batch LBFGS algorithm")
    ("learn-sgd",    po::bool_switch(&learn_sgd),    "online SGD algorithm")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C"            , po::value<double>(&C),           "regularization constant")
    
    ("unite",    po::bool_switch(&unite_kbest), "unite kbest sharing the same id")

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
