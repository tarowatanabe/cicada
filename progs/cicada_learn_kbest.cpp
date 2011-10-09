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
#include <boost/numeric/conversion/bounds.hpp>

#include "lbfgs.h"
#include "liblinear/linear.h"

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
bool learn_linear = false;

int linear_solver = L2R_L2LOSS_SVC_DUAL;

bool regularize_l1 = false;
bool regularize_l2 = false;

double C = 1.0;
double eps = std::numeric_limits<double>::infinity();

bool unite_kbest = false;

int threads = 2;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

void read_kbest(const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles);

template <typename Optimizer>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights);
template <typename Optimizer>
double optimize_linear(const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       weight_set_type& weights);

struct OptimizeLinear;
struct OptimizeLBFGS;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_linear > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,linear}");
    if (int(learn_lbfgs) + learn_linear == 0)
      learn_lbfgs = true;
    
    if (learn_lbfgs && regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));
    
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
    
    if (learn_linear)
      objective = optimize_linear<OptimizeLinear>(kbests, oracles, weights);
    else
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

struct OptimizeLinear
{
  typedef struct model        model_type;
  typedef struct parameter    parameter_type;
  typedef struct problem      problem_type;
  typedef struct feature_node feature_node_type;

  typedef std::vector<feature_node_type, std::allocator<feature_node_type> > feature_node_set_type;
  typedef std::vector<feature_node_type*, std::allocator<feature_node_type*> > feature_node_map_type;
  typedef std::vector<int, std::allocator<int> > label_set_type;
  
  typedef size_t size_type;
  typedef size_t offset_type;
  typedef std::vector<offset_type, std::allocator<offset_type> > offset_set_type;

  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
  
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#endif

  static void print_string_stderr(const char *s)
  {
    std::cerr << s << std::flush;
  }

  static void print_string_none(const char *s)
  {
    
  }

  struct Encoder
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;

    Encoder(queue_type& __queue,
	    const hypothesis_map_type& __kbests,
	    const hypothesis_map_type& __oracles)
      : queue(__queue),
	kbests(__kbests),
	oracles(__oracles) {}
    
    queue_type& queue;
    const hypothesis_map_type& kbests;
    const hypothesis_map_type& oracles;
    
    offset_set_type       offsets;
    feature_node_set_type features;
    
    void operator()()
    {
      offsets.clear();
      features.clear();
      
      sentence_unique_type  sentences;
      feature_node_type     feature;
      
      int id = 0;
      
      for (;;) {
	queue.pop(id);
	if (id < 0) break;
	
	sentences.clear();
	for (size_t o = 0; o != oracles[id].size(); ++ o)
	  sentences.insert(oracles[id][o].sentence);
	
	for (size_t o = 0; o != oracles[id].size(); ++ o)
	  for (size_t k = 0; k != kbests[id].size(); ++ k) {
	    const hypothesis_type& oracle = oracles[id][o];
	    const hypothesis_type& kbest  = kbests[id][k];
	    
	    // ignore oracle translations
	    if (sentences.find(kbest.sentence) != sentences.end()) continue;
	    
	    offsets.push_back(features.size());
	    
	    hypothesis_type::feature_set_type::const_iterator oiter = oracle.features.begin();
	    hypothesis_type::feature_set_type::const_iterator oiter_end = oracle.features.end();
	    
	    hypothesis_type::feature_set_type::const_iterator kiter = kbest.features.begin();
	    hypothesis_type::feature_set_type::const_iterator kiter_end = kbest.features.end();
	    
	    while (oiter != oiter_end && kiter != kiter_end) {
	      if (oiter->first < kiter->first) {
		feature.index = oiter->first.id() + 1;
		feature.value = oiter->second;
		features.push_back(feature);
		++ oiter;
	      } else if (kiter->first < oiter->first) {
		feature.index = kiter->first.id() + 1;
		feature.value = - kiter->second;
		features.push_back(feature);
		++ kiter;
	      } else {
		feature.index = oiter->first.id() + 1;
		feature.value = oiter->second - kiter->second;
		if (feature.value != 0.0)
		  features.push_back(feature);
		++ oiter;
		++ kiter;
	      }
	    }
	    
	    for (/**/; oiter != oiter_end; ++ oiter) {
	      feature.index = oiter->first.id() + 1;
	      feature.value = oiter->second;
	      features.push_back(feature);
	    }
	    
	    for (/**/; kiter != kiter_end; ++ kiter) {
	      feature.index = kiter->first.id() + 1;
	      feature.value = - kiter->second;
	      features.push_back(feature);
	    }
	    
	    // termination...
	    // if we have no instances, simply erase the last offsets
	    if (offsets.back() == features.size())
	      offsets.pop_back();
	    else {
	      feature.index = -1;
	      feature.value = 0.0;
	      features.push_back(feature);
	    }
	  }
      }
      
      feature_node_set_type(features).swap(features);
    }
  };
  typedef Encoder encoder_type;
  typedef std::vector<encoder_type, std::allocator<encoder_type> > encoder_set_type;
  
  OptimizeLinear(const hypothesis_map_type& kbests,
		 const hypothesis_map_type& oracles)
    : weights(), objective(0.0)
  {
    // compute unique before processing
    // 

    encoder_type::queue_type queue;
    encoder_set_type encoders(threads, encoder_type(queue, kbests, oracles));
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(encoders[i])));
    
    const size_t id_max = utils::bithack::min(kbests.size(), oracles.size());
    for (size_t id = 0; id != id_max; ++ id)
      if (! kbests[id].empty() && ! oracles[id].empty())
	queue.push(id);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    workers.join_all();
    
    size_t data_size = 0;
    for (int i = 0; i < threads; ++ i)
      data_size += encoders[i].offsets.size();
    
    if (debug)
      std::cerr << "liblinear data size: " << data_size << std::endl;
    
    label_set_type        labels;
    feature_node_map_type features;
    
    labels.reserve(data_size);
    features.reserve(data_size);
    
    labels.resize(data_size, 1);
    for (int i = 0; i < threads; ++ i) {
      for (size_type pos = 0; pos != encoders[i].offsets.size(); ++ pos)
	features.push_back(const_cast<feature_node_type*>(&(*encoders[i].features.begin())) + encoders[i].offsets[pos]);
      
      encoders[i].offsets.clear();
      offset_set_type(encoders[i].offsets).swap(encoders[i].offsets);
    }
    
    problem_type problem;
    
    problem.l = labels.size();
    problem.n = feature_type::allocated();
    problem.y = &(*labels.begin());
    problem.x = &(*features.begin());
    problem.bias = -1;

    parameter_type parameter;
    parameter.solver_type = linear_solver;
    parameter.eps = eps;
    parameter.C = 1.0 / (C * labels.size()); // renormalize!
    parameter.nr_weight    = 0;
    parameter.weight_label = 0;
    parameter.weight       = 0;

    if (parameter.eps == std::numeric_limits<double>::infinity()) {
      if (parameter.solver_type == L2R_LR || parameter.solver_type == L2R_L2LOSS_SVC)
	parameter.eps = 0.01;
      else if (parameter.solver_type == L2R_L2LOSS_SVC_DUAL || parameter.solver_type == L2R_L1LOSS_SVC_DUAL || parameter.solver_type == MCSVM_CS || parameter.solver_type == L2R_LR_DUAL)
	parameter.eps = 0.1;
      else if (parameter.solver_type == L1R_L2LOSS_SVC || parameter.solver_type == L1R_LR)
	parameter.eps = 0.01;
    }

    if (debug >= 2)
      set_print_string_function(print_string_stderr);
    else
      set_print_string_function(print_string_none);
    
    const char* error_message = check_parameter(&problem, &parameter);
    if (error_message)
      throw std::runtime_error(std::string("error: ") + error_message);
    
    static const char* names[] = {"L2R_LR", "L2R_L2LOSS_SVC_DUAL", "L2R_L2LOSS_SVC", "L2R_L1LOSS_SVC_DUAL", "MCSVM_CS",
				  "L1R_L2LOSS_SVC", "L1R_LR", "L2R_LR_DUAL"};
    
    if (debug)
      std::cerr << "solver: " << names[parameter.solver_type] << std::endl;
    
    const model_type* model = train(&problem, &parameter);

    objective = model->objective * C;
    
    // it is an optimization...
    weights.clear();
    for (int j = 0; j != model->nr_feature; ++ j)
      weights[weight_set_type::feature_type(j)] = model->w[j];
    
    free_and_destroy_model(const_cast<model_type**>(&model));
  }
  
public:
  weight_set_type weights;
  double objective;
};

template <typename Optimizer>
double optimize_linear(const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       weight_set_type& weights)
{
  Optimizer optimizer(kbests, oracles);
  
  weights = optimizer.weights;
  
  return optimizer.objective;
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
	 const hypothesis_map_type& __oracles,
	 const size_t& __instances)
      : queue(__queue),
	weights(__weights),
	kbests(__kbests),
	oracles(__oracles),
	instances(__instances)
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
      
      objective /= instances;
      std::transform(g.begin(), g.end(), g.begin(), std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    }

    queue_type&            queue;
    
    const weight_set_type& weights;
    
    const hypothesis_map_type& kbests;
    const hypothesis_map_type& oracles;
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

    const int id_max = utils::bithack::min(optimizer.kbests.size(), optimizer.oracles.size());
    size_t instances = 0;
    for (int id = 0; id != id_max; ++ id)
      instances += (! optimizer.kbests[id].empty() && ! optimizer.oracles[id].empty());
        
    queue_type queue;
    
    task_set_type tasks(threads, task_type(queue, optimizer.weights, optimizer.kbests, optimizer.oracles, instances));

    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    
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

struct TaskReadUnite
{
  typedef std::pair<path_type, path_type> path_pair_type;
  typedef utils::lockfree_list_queue<path_pair_type, std::allocator<path_pair_type> > queue_type;
  
  TaskReadUnite(queue_type& __queue)
    : queue(__queue) {}
  
  void operator()()
  {
    typedef boost::spirit::istream_iterator iter_type;
    typedef kbest_feature_parser<iter_type> parser_type;
    
    parser_type parser;
    kbest_feature_type kbest;

    for (;;) {
      path_pair_type paths;
      queue.pop(paths);
      
      if (paths.first.empty() && paths.second.empty()) break;
      
      const path_type& path           = (! paths.first.empty() ? paths.first : paths.second);
      hypothesis_map_type& hypotheses = (! paths.first.empty() ? kbests      : oracles);
      
      utils::compress_istream is(path, 1024 * 1024);
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
	
	if (id >= hypotheses.size())
	  hypotheses.resize(id + 1);
	
	hypotheses[id].push_back(hypothesis_type(kbest));
      }
    }
  }
  
  queue_type& queue;
  hypothesis_map_type kbests;
  hypothesis_map_type oracles;
};

struct TaskReadSync
{
  typedef std::pair<path_type, path_type> path_pair_type;
  typedef utils::lockfree_list_queue<path_pair_type, std::allocator<path_pair_type> > queue_type;
  
  TaskReadSync(queue_type& __queue)
    : queue(__queue) {}
  
  void operator()()
  {
    typedef boost::spirit::istream_iterator iter_type;
    typedef kbest_feature_parser<iter_type> parser_type;
    
    parser_type parser;
    kbest_feature_type kbest;

    for (;;) {
      path_pair_type paths;
      queue.pop(paths);
      
      if (paths.first.empty()) break;
      
      // we will perform paired reading...
      
      kbests.resize(kbests.size() + 1);
      oracles.resize(oracles.size() + 1);
      
      {
	utils::compress_istream is(paths.first, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	  
	iter_type iter(is);
	iter_type iter_end;
	  
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	    
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  kbests.back().push_back(hypothesis_type(kbest));
	}
      }

      {
	utils::compress_istream is(paths.second, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	  
	iter_type iter(is);
	iter_type iter_end;
	  
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	    
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  oracles.back().push_back(hypothesis_type(kbest));
	}
      }
    }
  }
  
  queue_type& queue;
  hypothesis_map_type kbests;
  hypothesis_map_type oracles;
};

void read_kbest(const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles)
{
  if (unite_kbest) {
    typedef TaskReadUnite task_type;
    typedef task_type::queue_type queue_type;

    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
    
    queue_type queue(threads);
    task_set_type tasks(threads, task_type(queue));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    for (path_set_type::const_iterator piter = kbest_path.begin(); piter != kbest_path.end(); ++ piter) {
      if (debug)
	std::cerr << "reading kbest: " << piter->string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;
	
	queue.push(std::make_pair(path_kbest, path_type()));
      }
    }
    
    for (path_set_type::const_iterator piter = oracle_path.begin(); piter != oracle_path.end(); ++ piter) {
      if (debug)
	std::cerr << "reading oracles: " << piter->string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_oracle = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_oracle)) break;
	
	queue.push(std::make_pair(path_type(), path_oracle));
      }
    }
    
    for (int i = 0; i != threads; ++ i)
      queue.push(std::make_pair(path_type(), path_type()));
    
    workers.join_all();

    size_t kbests_size  = 0;
    size_t oracles_size = 0;
    for (int i = 0; i != threads; ++ i) {
      kbests_size = utils::bithack::max(kbests_size, tasks[i].kbests.size());
      oracles_size = utils::bithack::max(oracles_size, tasks[i].oracles.size());
    }
    
    if (kbests_size != oracles_size)
      throw std::runtime_error("kbest/oracle size do not match");
    
    kbests.reserve(kbests_size);
    kbests.resize(kbests_size);

    oracles.reserve(oracles_size);
    oracles.resize(oracles_size);
    
    for (int i = 0; i != threads; ++ i) {
      for (size_t id = 0; id != tasks[i].kbests.size(); ++ id)
	kbests[id].insert(kbests[id].end(), tasks[i].kbests[id].begin(), tasks[i].kbests[id].end());
      
      tasks[i].kbests.clear();
    }
    
    unique_kbest(kbests);
    
    for (int i = 0; i != threads; ++ i) {
      for (size_t id = 0; id != tasks[i].oracles.size(); ++ id)
	oracles[id].insert(oracles[id].end(), tasks[i].oracles[id].begin(), tasks[i].oracles[id].end());
      
      tasks[i].oracles.clear();
    }
    
    unique_kbest(oracles);
    
  } else {
    typedef TaskReadSync task_type;
    typedef task_type::queue_type queue_type;
    
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
    
    // synchronous reading...
    if (kbest_path.size() != oracle_path.size())
      throw std::runtime_error("# of kbests does not match");
    
    queue_type queue(threads);
    task_set_type tasks(threads, task_type(queue));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    
    for (size_t pos = 0; pos != kbest_path.size(); ++ pos) {
      if (debug)
	std::cerr << "reading kbest: " << kbest_path[pos].string() << " with " << oracle_path[pos].string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest  = kbest_path[pos] / file_name;
	const path_type path_oracle = oracle_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;
	if (! boost::filesystem::exists(path_oracle)) continue;

	queue.push(std::make_pair(path_kbest, path_oracle));
      }
    }
    
    for (int i = 0; i != threads; ++ i)
      queue.push(std::make_pair(path_type(), path_type()));
    
    workers.join_all();

    size_t kbests_size  = 0;
    size_t oracles_size = 0;
    for (int i = 0; i != threads; ++ i) {
      kbests_size  += tasks[i].kbests.size();
      oracles_size += tasks[i].oracles.size();
    }
    
    if (kbests_size != oracles_size)
      throw std::runtime_error("kbest/oracle size do not match");

    kbests.reserve(kbests_size);
    oracles.reserve(oracles_size);
    
    for (int i = 0; i != threads; ++ i) {
      kbests.insert(kbests.end(), tasks[i].kbests.begin(), tasks[i].kbests.end());
      oracles.insert(oracles.end(), tasks[i].oracles.begin(), tasks[i].oracles.end());
      
      tasks[i].kbests.clear();
      tasks[i].oracles.clear();
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
    ("learn-linear", po::bool_switch(&learn_linear), "liblinear algorithm")
    ("solver",       po::value<int>(&linear_solver), "liblinear solver type (default: 1)\n"
     " 0: \tL2-regularized logistic regression (primal)\n"
     " 1: \tL2-regularized L2-loss support vector classification (dual)\n"
     " 2: \tL2-regularized L2-loss support vector classification (primal)\n"
     " 3: \tL2-regularized L1-loss support vector classification (dual)\n"
     " 5: \tL1-regularized L2-loss support vector classification\n"
     " 6: \tL1-regularized logistic regression\n"
     " 7: \tL2-regularized logistic regression (dual)")
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C",             po::value<double>(&C)->default_value(C), "regularization constant")
    ("eps",           po::value<double>(&eps),                 "tolerance for liblinear")
    
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
