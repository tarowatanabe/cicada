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
#include "cicada_text_impl.hpp"

#include "cicada/optimize_qp.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/sgi_hash_set.hpp"
#include "utils/random_seed.hpp"

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

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type kbest_path;
path_set_type oracle_path;
path_type weights_path;
path_type output_path = "-";
path_type output_objective_path;

path_set_type refset_files;

int iteration = 100;
bool learn_sgd = false;
bool learn_lbfgs = false;
bool learn_mira = false;
bool learn_linear = false;
bool learn_svm = false;

int linear_solver = L2R_L2LOSS_SVC_DUAL;

bool regularize_l1 = false;
bool regularize_l2 = false;

double C = 1.0;
double eps = std::numeric_limits<double>::infinity();

bool loss_margin = false; // margin by loss, not rank-loss
bool softmax_margin = false;
bool line_search = false;
bool normalize_vector = false;

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool unite_kbest = false;

int threads = 2;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

void read_kbest(const scorer_document_type& scorers,
		const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles);

void read_refset(const path_set_type& file,
		 scorer_document_type& scorers);

template <typename Optimizer>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights);
template <typename Optimizer>
double optimize_svm(const hypothesis_map_type& kbests,
		    const hypothesis_map_type& oracles,
		    weight_set_type& weights);

struct OptimizeLinear;
struct OptimizeSVM;
struct OptimizeLBFGS;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_linear + learn_svm > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,linear,svm}");
    if (int(learn_lbfgs) + learn_linear + learn_svm == 0)
      learn_lbfgs = true;
    
    if (learn_lbfgs && regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));
    
    if (kbest_path.empty())
      throw std::runtime_error("no kbest?");
    if (oracle_path.empty())
      throw std::runtime_error("no oracke kbest?");
    
    threads = utils::bithack::max(1, threads);

    scorer_document_type scorers(scorer_name);

    if (! refset_files.empty()) {
      read_refset(refset_files, scorers);
      
      if (! unite_kbest && kbest_path.size() > 1) {
	scorer_document_type scorers_iterative(scorer_name);
	scorers_iterative.resize(scorers.size() * kbest_path.size());
	
	for (size_t i = 0; i != kbest_path.size(); ++ i)
	  std::copy(scorers.begin(), scorers.end(), scorers_iterative.begin() + scorers.size() * i);
	
	scorers.swap(scorers_iterative);
      }
    }
    
    hypothesis_map_type kbests;
    hypothesis_map_type oracles;
    
    read_kbest(scorers, kbest_path, oracle_path, kbests, oracles);
    
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
    generator.seed(utils::random_seed());
    
    if (learn_linear)
      objective = optimize_svm<OptimizeLinear>(kbests, oracles, weights);
    else if (learn_svm)
      objective = optimize_svm<OptimizeSVM>(kbests, oracles, weights);
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
    weight_set_type       norms;
    
    void operator()()
    {
      offsets.clear();
      features.clear();
      norms.clear();
      
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

		norms[oiter->first] += feature.value * feature.value;
		
		++ oiter;
	      } else if (kiter->first < oiter->first) {
		feature.index = kiter->first.id() + 1;
		feature.value = - kiter->second;
		features.push_back(feature);
		
		norms[kiter->first] += feature.value * feature.value;

		++ kiter;
	      } else {
		feature.index = oiter->first.id() + 1;
		feature.value = oiter->second - kiter->second;
		if (feature.value != 0.0) {
		  features.push_back(feature);
		  
		  norms[oiter->first] += feature.value * feature.value;
		}
		++ oiter;
		++ kiter;
	      }
	    }
	    
	    for (/**/; oiter != oiter_end; ++ oiter) {
	      feature.index = oiter->first.id() + 1;
	      feature.value = oiter->second;
	      features.push_back(feature);
	      
	      norms[oiter->first] += feature.value * feature.value;
	    }
	    
	    for (/**/; kiter != kiter_end; ++ kiter) {
	      feature.index = kiter->first.id() + 1;
	      feature.value = - kiter->second;
	      features.push_back(feature);
	      
	      norms[kiter->first] += feature.value * feature.value;
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

  struct Normalize
  {
    typedef utils::lockfree_list_queue<size_t, std::allocator<size_t> > queue_type;
    
    Normalize(queue_type& __queue,
	      feature_node_map_type& __features,
	      const weight_set_type& __norms)
      : queue(__queue), features(__features), norms(__norms) {}

    void operator()()
    {
      for (;;) {
	size_t id = 0;
	queue.pop(id);
	if (id == size_t(-1)) break;
	
	for (feature_node_type* feat = features[id]; feat->index != -1; ++ feat)
	  feat->value *= (norms[feat->index - 1] == 0.0 ? 1.0 : norms[feat->index - 1]);
      }
    }
    
    queue_type&    queue;
    
    feature_node_map_type& features;
    const weight_set_type& norms;
  };
  
  struct Gradient
  {
    typedef std::pair<double, double> point_type;
    typedef std::vector<point_type, std::allocator<point_type> > point_set_type;
    
    typedef utils::lockfree_list_queue<size_t, std::allocator<size_t> > queue_type;
    
    Gradient(queue_type& __queue,
	     const feature_node_map_type& __features,
	     const label_set_type&        __labels,
	     const weight_set_type&       __weights,
	     const weight_set_type&       __weights_prev)
      : queue(__queue), features(__features), labels(__labels), weights(__weights), weights_prev(__weights_prev) {}

    void operator()()
    {
      static const double inf = std::numeric_limits<double>::infinity();

      points.clear();
      
      double& grad_pos = grads.first;
      double& grad_neg = grads.second;
      
      grad_pos = 0.0;
      grad_neg = 0.0;
      
      for (;;) {
	size_t id = 0;
	queue.pop(id);
	if (id == size_t(-1)) break;
	
	double margin      = 0.0;
	double margin_prev = 0.0;
	
	for (const feature_node_type* feat = features[id]; feat->index != -1; ++ feat) {
	  margin      += weights[feat->index - 1]      * feat->value;
	  margin_prev += weights_prev[feat->index - 1] * feat->value;
	}
	
	const double bi_pos = margin_prev - margin;
	const double ci_pos = labels[id]  - margin_prev;
	const double ki_pos = (bi_pos != 0.0 ? - ci_pos / bi_pos : - inf);
	
	const double bi_neg = margin_prev + margin;
	const double ci_neg = labels[id]  - margin_prev;
	const double ki_neg = (bi_neg != 0.0 ? - ci_neg / bi_neg : - inf);
	
	if (ki_pos > 0)
	  points.push_back(std::make_pair(ki_pos, bi_pos));
	
	if (ki_neg > 0)
	  points.push_back(std::make_pair(- ki_neg, bi_neg));
	
	grad_pos += bi_pos * ((bi_pos < 0.0 && ki_pos > 0.0) || (bi_pos > 0.0 && ki_pos <= 0.0));
	grad_neg += bi_neg * ((bi_neg < 0.0 && ki_neg > 0.0) || (bi_neg > 0.0 && ki_neg <= 0.0));
      }
      
      std::sort(points.begin(), points.end());
    }

    queue_type&    queue;
    
    const feature_node_map_type& features;
    const label_set_type&        labels;
    const weight_set_type&       weights;
    const weight_set_type&       weights_prev;
    
    point_type     grads;
    point_set_type points;
  };
  
  
  OptimizeLinear(const hypothesis_map_type& kbests,
		 const hypothesis_map_type& oracles,
		 weight_set_type& weights_prev)
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

    weight_set_type norms;
    
    size_t data_size = 0;
    for (int i = 0; i < threads; ++ i) {
      data_size += encoders[i].offsets.size();
      norms += encoders[i].norms;
    }

    
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
    
    if (normalize_vector) {
      typedef Normalize normalizer_type;

      norms.allocate(0.0);
      for (size_t i = 0; i != norms.size(); ++ i) {
	if (norms[i] == 0.0)
	  norms[i] = 1.0;
	else
	  norms[i] = 1.0 / std::sqrt(norms[i] / data_size);
      }
      
      normalizer_type::queue_type queue;
      
      boost::thread_group workers;
      for (int i = 0; i < threads; ++ i)
	workers.add_thread(new boost::thread(normalizer_type(queue, features, norms)));
      
      for (size_t i = 0; i != features.size(); ++ i)
	queue.push(i);
      for (int i = 0; i < threads; ++ i)
	queue.push(size_t(-1));
      
      workers.join_all();
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
    
    // line search...
    
    if (line_search) {
      typedef Gradient gradient_type;
      typedef std::vector<gradient_type, std::allocator<gradient_type> > gradient_set_type;

      typedef gradient_type::point_set_type point_set_type;

      gradient_type::queue_type queue;
      gradient_set_type gradients(threads, gradient_type(queue, features, labels, weights, weights_prev));
      
      boost::thread_group workers;
      for (int i = 0; i < threads; ++ i)
	workers.add_thread(new boost::thread(boost::ref(gradients[i])));
      
      for (size_t i = 0; i != features.size(); ++ i)
	queue.push(i);
      for (int i = 0; i < threads; ++ i)
	queue.push(size_t(-1));
      
      workers.join_all();
      
      // merge points...
      point_set_type points;
      point_set_type points_next;

      double grad_pos = 0.0;
      double grad_neg = 0.0;
      size_t samples = 0;
      for (int i = 0; i < threads; ++ i) {
	grad_pos += gradients[i].grads.first;
	grad_neg += gradients[i].grads.second;
	samples += gradients[i].labels.size();

	if (points.empty())
	  points.swap(gradients[i].points);
	else {
	  points_next.clear();
	  std::merge(points.begin(), points.end(), gradients[i].points.begin(), gradients[i].points.end(), std::back_inserter(points_next));
	  points.swap(points_next);
	}
      }
      
      const double norm_w      = cicada::dot_product(weights, weights);
      const double dot_prod    = cicada::dot_product(weights_prev, weights);
      const double norm_w_prev = cicada::dot_product(weights_prev, weights_prev);
      
      const double a0_pos = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * samples;
      const double b0_pos = (dot_prod - norm_w_prev) * C * samples;
      
      const double a0_neg = (norm_w + 2.0 * dot_prod + norm_w_prev) * C * samples;
      const double b0_neg = (- dot_prod - norm_w_prev) * C * samples;
      
      grad_pos += b0_pos;
      grad_neg += b0_neg;
      
      if (grad_pos < 0.0) {
	double k = 0.0;
	  
	point_set_type::const_iterator piter = std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0));
	point_set_type::const_iterator piter_end = points.end();
	  
	for (/**/; piter != piter_end && grad_pos < 0.0; ++ piter) {
	  const double k_new = piter->first;
	  const double grad_new = grad_pos + std::fabs(piter->second) + a0_pos * (k_new - k);
	    
	  if (grad_new >= 0) {
	    // compute intersection...
	    k = k + grad_pos * (k - k_new) / (grad_new - grad_pos);
	    grad_pos = grad_new;
	    break;
	  } else {
	    k = k_new;
	    grad_pos = grad_new;
	  }
	}
	  
	if (debug >= 3)
	  std::cerr << "grad: " << grad_pos << "  k: " << k << std::endl;

	if (k > 0.0) {
	  if (debug >= 3)
	    std::cerr << "current weights" << std::endl
		      << weights;

	  weights      *= k;
	  weights_prev *= (1.0 - k);
	  weights += weights_prev;

	  if (debug >= 3)
	    std::cerr << "updated weights" << std::endl
		      << weights;
	}
      } else if (grad_neg < 0.0) {
	double k = 0.0;
	  
	point_set_type::const_reverse_iterator piter(std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0)));
	point_set_type::const_reverse_iterator piter_end = points.rend();
	  
	for (/**/; piter != piter_end && grad_neg < 0.0; ++ piter) {
	  const double k_new = - piter->first;
	  const double grad_new = grad_neg + std::fabs(piter->second) + a0_neg * (k_new - k);
	    
	  if (grad_new >= 0) {
	    // compute intersection...
	    k = k + grad_neg * (k - k_new) / (grad_new - grad_neg);
	    grad_neg = grad_new;
	    break;
	  } else {
	    k = k_new;
	    grad_neg = grad_new;
	  }
	}
	  
	if (debug >= 3)
	  std::cerr << "grad: " << grad_neg << "  k: " << - k << std::endl;
	  
	if (k > 0.0) {
	  if (debug >= 3)
	    std::cerr << "current weights" << std::endl
		      << weights;
	  
	  weights      *= - k;
	  weights_prev *= (1.0 + k);
	  weights += weights_prev;

	  if (debug >= 3)
	    std::cerr << "updated weights" << std::endl
		      << weights;
	}
      }
    }
    
    if (normalize_vector)
      for (size_t i = 0; i != weights.size(); ++ i)
	weights[i] *= (norms[i] == 0.0 ? 1.0 : 1.0 / norms[i]);
  }
  
public:
  weight_set_type weights;
  double objective;
};

template <typename Optimizer>
double optimize_svm(const hypothesis_map_type& kbests,
		    const hypothesis_map_type& oracles,
		    weight_set_type& weights)
{
  Optimizer optimizer(kbests, oracles, weights);
  
  weights = optimizer.weights;
  
  return optimizer.objective;
}

struct OptimizeSVM
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef hypothesis_type::feature_value_type feature_value_type;

  struct SampleSet
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    typedef std::vector<size_type, std::allocator<size_type> > offsets_type;

    struct Sample
    {
      typedef features_type::const_iterator const_iterator;

      Sample(const_iterator __first, const_iterator __last) : first(__first), last(__last) {}

      const_iterator begin() const { return first; }
      const_iterator end() const { return last; }
      size_type size() const { return last - first; }
      bool emtpy() const { return first == last; }
      
      const_iterator first;
      const_iterator last;
    };

    typedef Sample sample_type;
    typedef sample_type value_type;
    
    SampleSet() : features(), offsets() { offsets.push_back(0); }
    
    void clear()
    {
      features.clear();
      offsets.clear();
      offsets.push_back(0);
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      features.insert(features.end(), first, last);
      offsets.push_back(features.size());
    }
    
    sample_type operator[](size_type pos) const
    {
      return sample_type(features.begin() + offsets[pos], features.begin() + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }

    void swap(SampleSet& x)
    {
      features.swap(x.features);
      offsets.swap(x.offsets);
    }

    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    features_type features;
    offsets_type  offsets;
  };
  
  typedef SampleSet sample_set_type;

  typedef std::vector<double, std::allocator<double> > loss_set_type;
  typedef std::vector<double, std::allocator<double> > alpha_set_type;
  typedef std::vector<double, std::allocator<double> > f_set_type;
  
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
  
  typedef std::pair<size_type, size_type> pos_pair_type;
  typedef std::vector<pos_pair_type, std::allocator<pos_pair_type> > pos_pair_set_type;

  struct Encoder
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;

    typedef std::pair<double, double> point_type;
    typedef std::vector<point_type, std::allocator<point_type> > point_set_type;
    
    Encoder(queue_type& __queue,
	    const hypothesis_map_type& __kbests,
	    const hypothesis_map_type& __oracles)
      : queue(__queue), kbests(__kbests), oracles(__oracles) {}
    
    void operator()()
    {
      typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
      
      features_type feats;
      sentence_unique_type  sentences;
      
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
	    
	    hypothesis_type::feature_set_type::const_iterator oiter = oracle.features.begin();
	    hypothesis_type::feature_set_type::const_iterator oiter_end = oracle.features.end();
	    
	    hypothesis_type::feature_set_type::const_iterator kiter = kbest.features.begin();
	    hypothesis_type::feature_set_type::const_iterator kiter_end = kbest.features.end();
	    
	    feats.clear();
	
	    while (oiter != oiter_end && kiter != kiter_end) {
	      if (oiter->first < kiter->first) {
		feats.push_back(*oiter);
		++ oiter;
	      } else if (kiter->first < oiter->first) {
		feats.push_back(feature_value_type(kiter->first, - kiter->second));
		++ kiter;
	      } else {
		const double value = oiter->second - kiter->second;
		if (value != 0.0)
		  feats.push_back(feature_value_type(kiter->first, value));
		++ oiter;
		++ kiter;
	      }
	    }
	
	    for (/**/; oiter != oiter_end; ++ oiter)
	      feats.push_back(*oiter);
	
	    for (/**/; kiter != kiter_end; ++ kiter)
	      feats.push_back(feature_value_type(kiter->first, - kiter->second));
	    
	    if (feats.empty()) continue;
	    
	    if (loss_margin) {
	      const double loss = kbest.loss - oracle.loss;
	      
	      // checking...
	      if (loss > 0.0) {
		features.insert(feats.begin(), feats.end());
		losses.push_back(loss);
	      }
	    } else {
	      features.insert(feats.begin(), feats.end());
	      losses.push_back(1.0);
	    }
	  }
      }
    }
    
    queue_type& queue;
    const hypothesis_map_type& kbests;
    const hypothesis_map_type& oracles;
    
    sample_set_type features;
    loss_set_type   losses;

    point_set_type points;
    point_type     grads;
  };

  typedef Encoder encoder_type;
  typedef std::vector<encoder_type, std::allocator<encoder_type> > encoder_set_type;
  
  
  
  struct HMatrix
  {
    HMatrix(const pos_pair_set_type& __positions,
	    const encoder_set_type&   __encoders)
      : positions(__positions), encoders(__encoders) {}

    double operator()(int i, int j) const
    {
      const pos_pair_type& pos_i = positions[i];
      const pos_pair_type& pos_j = positions[j];
      
      return cicada::dot_product(encoders[pos_i.first].features[pos_i.second].begin(), encoders[pos_i.first].features[pos_i.second].end(),
				 encoders[pos_j.first].features[pos_j.second].begin(), encoders[pos_j.first].features[pos_j.second].end(),
				 0.0);
    }
    
    const pos_pair_set_type& positions;
    const encoder_set_type& encoders;
  };
  
  struct MMatrix
  {
    MMatrix(const pos_pair_set_type& __positions,
	    const encoder_set_type&   __encoders)
      : positions(__positions), encoders(__encoders) {}
    
    template <typename W>
    void operator()(W& w, const alpha_set_type& alpha) const
    {
      alpha_set_type::const_iterator aiter = alpha.begin();
      
      const size_type model_size = encoders.size();
      for (size_type i = 0; i != model_size; ++ i) {
	const size_type features_size = encoders[i].features.size();
	
	for (size_type j = 0; j != features_size; ++ j, ++ aiter)
	  if (*aiter > 0.0) {
	    sample_set_type::value_type::const_iterator fiter_end = encoders[i].features[j].end();
	    for (sample_set_type::value_type::const_iterator fiter = encoders[i].features[j].begin(); fiter != fiter_end; ++ fiter) 
	      w[fiter->first] += (*aiter) * fiter->second;
	  }
      }
    }
    
    template <typename W>
    double operator()(const W& w, const size_t& i) const
    {
      const pos_pair_type& pos_i = positions[i];
      
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = encoders[pos_i.first].features[pos_i.second].end();
      for (sample_set_type::value_type::const_iterator fiter = encoders[pos_i.first].features[pos_i.second].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename W>
    void operator()(W& w, const double& update, const size_t& i) const
    {
      const pos_pair_type& pos_i = positions[i];
      
      sample_set_type::value_type::const_iterator fiter_end = encoders[pos_i.first].features[pos_i.second].end();
      for (sample_set_type::value_type::const_iterator fiter = encoders[pos_i.first].features[pos_i.second].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const pos_pair_set_type& positions;
    const encoder_set_type&  encoders;
  };
  
  struct Gradient
  {
    Gradient(encoder_type& __encoder,
	     const weight_set_type& __weights,
	     const weight_set_type& __weights_prev)
      : encoder(__encoder), weights(__weights), weights_prev(__weights_prev) {}

    void operator()()
    {
      static const double inf = std::numeric_limits<double>::infinity();

      encoder.points.clear();
      
      double& grad_pos = encoder.grads.first;
      double& grad_neg = encoder.grads.second;
      
      grad_pos = 0.0;
      grad_neg = 0.0;
      
      for (size_t id = 0; id != encoder.losses.size(); ++ id) {
	const double margin      = cicada::dot_product(weights,      encoder.features[id].begin(), encoder.features[id].end(), 0.0);
	const double margin_prev = cicada::dot_product(weights_prev, encoder.features[id].begin(), encoder.features[id].end(), 0.0);
	
	const double bi_pos = margin_prev - margin;
	const double ci_pos = encoder.losses[id]  - margin_prev;
	const double ki_pos = (bi_pos != 0.0 ? - ci_pos / bi_pos : - inf);
	
	const double bi_neg = margin_prev + margin;
	const double ci_neg = encoder.losses[id]  - margin_prev;
	const double ki_neg = (bi_neg != 0.0 ? - ci_neg / bi_neg : - inf);
	
	if (ki_pos > 0)
	  encoder.points.push_back(std::make_pair(ki_pos, bi_pos));
	
	if (ki_neg > 0)
	  encoder.points.push_back(std::make_pair(- ki_neg, bi_neg));
	
	grad_pos += bi_pos * ((bi_pos < 0.0 && ki_pos > 0.0) || (bi_pos > 0.0 && ki_pos <= 0.0));
	grad_neg += bi_neg * ((bi_neg < 0.0 && ki_neg > 0.0) || (bi_neg > 0.0 && ki_neg <= 0.0));
      }

      std::sort(encoder.points.begin(), encoder.points.end());
    }
    
    encoder_type&          encoder;
    const weight_set_type& weights;
    const weight_set_type& weights_prev;
  };


  OptimizeSVM(const hypothesis_map_type& kbests,
	      const hypothesis_map_type& oracles,
	      weight_set_type& weights_prev)
    : weights(), objective(0.0), tolerance(0.1)
  {
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
    
    // encoding finished!
    
    size_type data_size = 0;
    for (size_type i = 0; i != encoders.size(); ++ i)
      data_size += encoders[i].losses.size();

    if (debug)
      std::cerr << "# of support vectors: " << data_size << std::endl;
    
    pos_pair_set_type positions;
    f_set_type        f;

    positions.reserve(data_size);
    f.reserve(data_size);
    
    for (size_type i = 0; i != encoders.size(); ++ i)
      for (size_type j = 0; j != encoders[i].features.size(); ++ j) {
	positions.push_back(std::make_pair(i, j));
	f.push_back(- encoders[i].losses[j]);
      }
    
    alpha_set_type alpha(data_size, 0.0);
    
    cicada::optimize::QPDCD solver;
    
    HMatrix H(positions, encoders);
    MMatrix M(positions, encoders);
    
    objective = solver(alpha, f, H, M, 1.0 / (C * data_size), tolerance);
    objective *= C;
    
    size_type actives = 0;

    weights.clear();
    alpha_set_type::const_iterator aiter = alpha.begin();
    for (size_type i = 0; i != encoders.size(); ++ i)
      for (size_type j = 0; j != encoders[i].features.size(); ++ j, ++ aiter) 
	if (*aiter > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = encoders[i].features[j].end();
	  for (sample_set_type::value_type::const_iterator fiter = encoders[i].features[j].begin(); fiter != fiter_end; ++ fiter)
	    weights[fiter->first] += (*aiter) * fiter->second; 
	  
	  ++ actives;
	}
    
    if (debug)
      std::cerr << "# of active vectors: " << actives << std::endl;
    
    // line search between the previous solution and the current solution
    if (line_search) {
      typedef encoder_type::point_set_type point_set_type;

      boost::thread_group workers;
      for (int i = 0; i < threads; ++ i)
	workers.add_thread(new boost::thread(Gradient(encoders[i], weights, weights_prev)));
      workers.join_all();
      
      // merge points...
      point_set_type points;
      point_set_type points_next;

      double grad_pos = 0.0;
      double grad_neg = 0.0;
      size_t samples = 0;
      for (int i = 0; i < threads; ++ i) {
	grad_pos += encoders[i].grads.first;
	grad_neg += encoders[i].grads.second;
	samples += encoders[i].losses.size();

	if (points.empty())
	  points.swap(encoders[i].points);
	else {
	  points_next.clear();
	  std::merge(points.begin(), points.end(), encoders[i].points.begin(), encoders[i].points.end(), std::back_inserter(points_next));
	  points.swap(points_next);
	}
      }
      
      const double norm_w      = cicada::dot_product(weights, weights);
      const double dot_prod    = cicada::dot_product(weights_prev, weights);
      const double norm_w_prev = cicada::dot_product(weights_prev, weights_prev);
      
      const double a0_pos = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * samples;
      const double b0_pos = (dot_prod - norm_w_prev) * C * samples;
      
      const double a0_neg = (norm_w + 2.0 * dot_prod + norm_w_prev) * C * samples;
      const double b0_neg = (- dot_prod - norm_w_prev) * C * samples;
      
      grad_pos += b0_pos;
      grad_neg += b0_neg;
      
      if (grad_pos < 0.0) {
	double k = 0.0;
	  
	point_set_type::const_iterator piter = std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0));
	point_set_type::const_iterator piter_end = points.end();
	  
	for (/**/; piter != piter_end && grad_pos < 0.0; ++ piter) {
	  const double k_new = piter->first;
	  const double grad_new = grad_pos + std::fabs(piter->second) + a0_pos * (k_new - k);
	    
	  if (grad_new >= 0) {
	    // compute intersection...
	    k = k + grad_pos * (k - k_new) / (grad_new - grad_pos);
	    grad_pos = grad_new;
	    break;
	  } else {
	    k = k_new;
	    grad_pos = grad_new;
	  }
	}
	  
	if (debug >= 3)
	  std::cerr << "grad: " << grad_pos << "  k: " << k << std::endl;

	if (k > 0.0) {
	  weights      *= k;
	  weights_prev *= (1.0 - k);
	  weights += weights_prev;
	}
      } else if (grad_neg < 0.0) {
	double k = 0.0;
	  
	point_set_type::const_reverse_iterator piter(std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0)));
	point_set_type::const_reverse_iterator piter_end = points.rend();
	  
	for (/**/; piter != piter_end && grad_neg < 0.0; ++ piter) {
	  const double k_new = - piter->first;
	  const double grad_new = grad_neg + std::fabs(piter->second) + a0_neg * (k_new - k);
	    
	  if (grad_new >= 0) {
	    // compute intersection...
	    k = k + grad_neg * (k - k_new) / (grad_new - grad_neg);
	    grad_neg = grad_new;
	    break;
	  } else {
	    k = k_new;
	    grad_neg = grad_new;
	  }
	}
	  
	if (debug >= 3)
	  std::cerr << "grad: " << grad_neg << "  k: " << - k << std::endl;
	  
	if (k > 0.0) {
	  weights      *= - k;
	  weights_prev *= (1.0 + k);
	  weights += weights_prev;
	}
      }
    }
  }
  
public:
  weight_set_type weights;
  double objective;

  const double tolerance;
};


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

      const double cost_factor = (softmax_margin ? 1.0 : 0.0);
      
      while (1) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
	weight_type Z_oracle;
	weight_type Z_kbest;
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter)
	  Z_oracle += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss));
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter)
	  Z_kbest += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss));
	
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss)) / Z_oracle;
	  
	  hypothesis_type::feature_set_type::const_iterator fiter_end = oiter->features.end();
	  for (hypothesis_type::feature_set_type::const_iterator fiter = oiter->features.begin(); fiter != fiter_end; ++ fiter)
	    expectations[fiter->first] -= weight_type(fiter->second) * weight;
	}
	
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss)) / Z_kbest;
	  
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

struct TaskLoss
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  TaskLoss(queue_type& __queue,
	   hypothesis_map_type& __kbests,
	   const scorer_document_type& __scorers)
    : queue(__queue), kbests(__kbests), scorers(__scorers) {}

  void operator()()
  {
    const bool error_metric = scorers.error_metric();
    const double loss_factor = (error_metric ? 1.0 : - 1.0);
    
    for (;;) {
      int id = 0;
      queue.pop(id);
      if (id < 0) break;

      hypothesis_set_type::iterator kiter_end = kbests[id].end();
      for (hypothesis_set_type::iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	kiter->score = scorers[id]->score(sentence_type(sentence_type(kiter->sentence.begin(), kiter->sentence.end())));
	kiter->loss = kiter->score->score() * loss_factor;
      }
    }
  }
  
  queue_type& queue;
  hypothesis_map_type&        kbests;
  const scorer_document_type& scorers;
};


void loss_kbest(hypothesis_map_type& kbests, const scorer_document_type& scorers)
{
  typedef TaskLoss task_type;
  typedef task_type::queue_type queue_type;
  
  queue_type queue;
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, kbests, scorers)));
  
  for (size_t id = 0; id != kbests.size(); ++ id)
    if (! kbests[id].empty())
      queue.push(id);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
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
  typedef boost::fusion::tuple<path_type, path_type, size_t> path_pair_type;
  typedef utils::lockfree_list_queue<path_pair_type, std::allocator<path_pair_type> > queue_type;
  
  TaskReadSync(queue_type& __queue, const scorer_document_type& __scorers)
    : queue(__queue), scorers(__scorers) {}
  
  void operator()()
  {
    typedef boost::spirit::istream_iterator iter_type;
    typedef kbest_feature_parser<iter_type> parser_type;
    
    const bool error_metric = scorers.error_metric();
    const double loss_factor = (error_metric ? 1.0 : - 1.0);
    
    parser_type parser;
    kbest_feature_type kbest;

    for (;;) {
      path_pair_type paths;
      queue.pop(paths);
      
      if (boost::fusion::get<0>(paths).empty()) break;
      
      // we will perform paired reading...
      
      kbests.resize(kbests.size() + 1);
      oracles.resize(oracles.size() + 1);

      const size_t refpos = boost::fusion::get<2>(paths);
      
      if (! scorers.empty())
	if (refpos >= scorers.size())
	  throw std::runtime_error("reference positions outof index");
      
      {
	utils::compress_istream is(boost::fusion::get<0>(paths), 1024 * 1024);
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
	  
	  hypothesis_type& kbest = kbests.back().back();
	  
	  if (! scorers.empty()) {
	    kbest.score = scorers[refpos]->score(sentence_type(kbest.sentence.begin(), kbest.sentence.end()));
	    kbest.loss  = kbest.score->score() * loss_factor;
	  } else
	    kbest.loss = 1;
	}
      }

      {
	utils::compress_istream is(boost::fusion::get<1>(paths), 1024 * 1024);
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

	  hypothesis_type& oracle = oracles.back().back();

	  if (! scorers.empty()) {
	    oracle.score = scorers[refpos]->score(sentence_type(oracle.sentence.begin(), oracle.sentence.end()));
	    oracle.loss  = oracle.score->score() * loss_factor;
	  } else
	    oracle.loss = 0.0;
	}
      }
    }
  }
  
  queue_type& queue;
  const scorer_document_type& scorers;
  
  hypothesis_map_type kbests;
  hypothesis_map_type oracles;
};

void read_kbest(const scorer_document_type& scorers,
		const path_set_type& kbest_path,
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

    if (! scorers.empty()) {
      if (scorers.size() != kbests.size())
	throw std::runtime_error("refset size do not match with kbest size");
      if (scorers.size() != oracles.size())
	throw std::runtime_error("refset size do not match with oracle size");
      
      loss_kbest(kbests, scorers);
      loss_kbest(oracles, scorers);
    } else {
      // fill zero loss to oracles, and one loss to kbests.
      
      for (size_t id = 0; id != kbests.size(); ++ id) {
	hypothesis_set_type::iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter)
	  kiter->loss = 1.0;
      }
      
      for (size_t id = 0; id != oracles.size(); ++ id) {
	hypothesis_set_type::iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter)
	  oiter->loss = 0.0;
      }
    }
    
  } else {
    typedef TaskReadSync task_type;
    
    typedef task_type::queue_type     queue_type;
    typedef task_type::path_pair_type path_pair_type;
    
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
    
    // synchronous reading...
    if (kbest_path.size() != oracle_path.size())
      throw std::runtime_error("# of kbests does not match");
    
    queue_type queue(threads);
    task_set_type tasks(threads, task_type(queue, scorers));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    size_t refpos = 0;
    for (size_t pos = 0; pos != kbest_path.size(); ++ pos) {
      if (debug)
	std::cerr << "reading kbest: " << kbest_path[pos].string() << " with " << oracle_path[pos].string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest  = kbest_path[pos] / file_name;
	const path_type path_oracle = oracle_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;
	if (! boost::filesystem::exists(path_oracle)) continue;
	
	queue.push(path_pair_type(path_kbest, path_oracle, refpos ++));
      }
    }
    
    for (int i = 0; i != threads; ++ i)
      queue.push(path_pair_type(path_type(), path_type(), size_t(-1)));
    
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

void read_refset(const path_set_type& files, scorer_document_type& scorers)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;

  if (files.empty())
    throw std::runtime_error("no reference files?");
    
  scorers.clear();

  parser_type parser;
  id_sentence_type id_sentence;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->string());

    utils::compress_istream is(*fiter, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      id_sentence.second.clear();
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
	if (iter != iter_end)
	  throw std::runtime_error("refset parsing failed");
      
      const int& id = id_sentence.first;
      
      if (id >= static_cast<int>(scorers.size()))
	scorers.resize(id + 1);
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      scorers[id]->insert(id_sentence.second);
    }
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("kbest",   po::value<path_set_type>(&kbest_path)->multitoken(),   "kbest path")
    ("oracle",  po::value<path_set_type>(&oracle_path)->multitoken(),  "oracle kbest path")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    ("weights", po::value<path_type>(&weights_path),                   "initial parameter")
    ("output",  po::value<path_type>(&output_path),                    "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",  po::bool_switch(&learn_lbfgs),  "batch LBFGS algorithm")
    ("learn-linear", po::bool_switch(&learn_linear), "liblinear algorithm")
    ("learn-svm",    po::bool_switch(&learn_svm),    "structural SVM")
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

    ("loss-margin",      po::bool_switch(&loss_margin),      "direct loss margin")
    ("softmax-margin",   po::bool_switch(&softmax_margin),   "softmax margin")
    ("line-search",      po::bool_switch(&line_search),      "perform line search in each iteration")
    ("normalize-vector", po::bool_switch(&normalize_vector), "normalize feature vectors")
    
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")

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
