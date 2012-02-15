//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>

//
// k-best learner
//

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <numeric>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"
#include "cicada_kbest_impl.hpp"
#include "cicada_text_impl.hpp"
#include "cicada_mert_kbest_impl.hpp"

#include "cicada/optimize_qp.hpp"
#include "cicada/optimize.hpp"
#include "cicada/semiring/envelope.hpp"
#include "cicada/feature_vector_compact.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/sgi_hash_set.hpp"
#include "utils/random_seed.hpp"
#include "utils/map_file.hpp"
#include "utils/tempfile.hpp"
#include "utils/mulvector2.hpp"
#include "utils/mathop.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>
#include <boost/functional/hash/hash.hpp>
#include <boost/numeric/conversion/bounds.hpp>

#include "lbfgs.h"
#include "lbfgs_error.hpp"
#include "liblinear/linear.h"

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
typedef std::vector<size_t, std::allocator<size_t> > kbest_map_type;

path_set_type kbest_path;
path_set_type oracle_path;
path_type weights_path;
path_set_type weights_history_path;

path_type output_path = "-";
path_type output_objective_path;

path_type bound_lower_file;
path_type bound_upper_file;

path_set_type refset_files;

int iteration = 100;
bool learn_sgd = false;
bool learn_xbleu = false;
bool learn_lbfgs = false;
bool learn_mira = false;
bool learn_linear = false;
bool learn_svm = false;

int linear_solver = L2R_L2LOSS_SVC_DUAL;

bool regularize_l1 = false;
bool regularize_l2 = false;
bool regularize_entropy = false;

double eps = std::numeric_limits<double>::infinity();
double C = 1.0;
double C2 = 1.0;
double scale = 1.0;
double eta0 = 0.2;
int order = 4;

bool annealing_mode = false;
bool quenching_mode = false;

double temperature = 0.0;
double temperature_start = 1000;
double temperature_end = 0.001;
double temperature_rate = 0.5;

double quench_start = 0.01;
double quench_end = 100;
double quench_rate = 10;

bool loss_margin = false; // margin by loss, not rank-loss
bool softmax_margin = false;
bool line_search = false;
bool mert_search = false;
bool sample_vector = false;
bool direct_loss = false;
bool conservative_loss = false;

bool scale_fixed = false;

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
		hypothesis_map_type& oracles,
		kbest_map_type& kbest_map);

void read_refset(const path_set_type& file,
		 scorer_document_type& scorers);

template <typename Optimizer>
double optimize_xbleu(const hypothesis_map_type& kbests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights);
template <typename Optimizer>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights);
template <typename Optimizer>
double optimize_svm(const hypothesis_map_type& kbests,
		    const hypothesis_map_type& oracles,
		    const weight_set_type& bounds_lower,
		    const weight_set_type& bounds_upper,
		    weight_set_type& weights);

double optimize_mert(const scorer_document_type& scorers,
		     const hypothesis_map_type& kbests,
		     const kbest_map_type& kbest_map,
		     const weight_set_type& weights_prev,
		     weight_set_type& weights);

struct OptimizeLinear;
struct OptimizeSVM;
struct OptimizeLBFGS;
struct OptimizeXBLEU;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_linear + learn_svm + learn_xbleu > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,linear,svm}");
    if (int(learn_lbfgs) + learn_linear + learn_svm + learn_xbleu == 0)
      learn_lbfgs = true;
    
    if (conservative_loss && line_search)
      throw std::runtime_error("we do not allow both conservative-update and line-search");
    
    if (learn_lbfgs && regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));
    
    if (kbest_path.empty())
      throw std::runtime_error("no kbest?");
    if (! learn_xbleu && oracle_path.empty())
      throw std::runtime_error("no oracle kbest?");
    if (learn_xbleu && refset_files.empty())
      throw std::runtime_error("no reference translations?");
    if (learn_xbleu && order <= 0)
      throw std::runtime_error("invalid ngram order");
    
    if (annealing_mode) {
      if (! (temperature_end < temperature_start))
	throw std::runtime_error("temperature should start higher, then decreased");
      if (temperature_rate <= 0.0 || temperature_rate >= 1.0)
	throw std::runtime_error("temperature rate should be 0.0 < rate < 1.0: " + utils::lexical_cast<std::string>(temperature_rate));
    }
    
    if (quenching_mode) {
      if (! (quench_start < quench_end))
	throw std::runtime_error("quenching should start lower, then increased");
      if (quench_rate <= 1.0)
	throw std::runtime_error("quenching rate should be > 1.0: " + utils::lexical_cast<std::string>(quench_rate)); 
    }

    if (! bound_lower_file.empty())
      if (bound_lower_file != "-" && ! boost::filesystem::exists(bound_lower_file))
	throw std::runtime_error("no lower-bound file? " + bound_lower_file.string());
    
    if (! bound_upper_file.empty())
      if (bound_upper_file != "-" && ! boost::filesystem::exists(bound_upper_file))
	throw std::runtime_error("no upper-bound file? " + bound_upper_file.string());
    
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

    if (mert_search && scorers.empty())
      throw std::runtime_error("mert search requires evaluation scores");
    if (sample_vector && scorers.empty())
      throw std::runtime_error("sampling requires evaluation scores");
    
    hypothesis_map_type kbests;
    hypothesis_map_type oracles;
    kbest_map_type      kbest_map;
    
    read_kbest(scorers, kbest_path, oracle_path, kbests, oracles, kbest_map);
    
    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;
    
    weight_set_type weights;
    if (! weights_path.empty()) {
      if (! boost::filesystem::exists(weights_path))
	throw std::runtime_error("no path? " + weights_path.string());
      
      utils::compress_istream is(weights_path, 1024 * 1024);
      is >> weights;
    }
    
    weight_set_type bounds_lower;
    weight_set_type bounds_upper;
    
    if (! bound_lower_file.empty())
      read_bounds(bound_lower_file, bounds_lower, - std::numeric_limits<double>::infinity());
    
    if (! bound_upper_file.empty())
      read_bounds(bound_upper_file, bounds_upper,   std::numeric_limits<double>::infinity());
    
    weights.allocate();
    
    double objective = 0.0;
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());

    const weight_set_type weights_prev = weights;
    
    if (learn_linear)
      objective = optimize_svm<OptimizeLinear>(kbests, oracles, bounds_lower, bounds_upper, weights);
    else if (learn_svm)
      objective = optimize_svm<OptimizeSVM>(kbests, oracles, bounds_lower, bounds_upper, weights);
    else if (learn_xbleu)
      objective = optimize_xbleu<OptimizeXBLEU>(kbests, scorers, weights);
    else
      objective = optimize_batch<OptimizeLBFGS>(kbests, oracles, weights);
    
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

    if (mert_search) {
      const double objective = optimize_mert(scorers, kbests, kbest_map, weights_prev, weights);
      
      if (debug)
	std::cerr << "mert objective: " << objective << std::endl;
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

struct OptimizeLinear
{
  typedef size_t    size_type;
  typedef size_t    offset_type;
  typedef ptrdiff_t difference_type;

  typedef struct model        model_type;
  typedef struct parameter    parameter_type;
  typedef struct problem      problem_type;
  typedef struct feature_node feature_node_type;

  typedef std::vector<feature_node_type, std::allocator<feature_node_type> > feature_node_set_type;
  typedef std::vector<feature_node_type*, std::allocator<feature_node_type*> > feature_node_map_type;
  typedef std::vector<int, std::allocator<int> > label_set_type;
  
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
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef hypothesis_type::feature_value_type feature_value_type;

    struct SampleSet
    {
      typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
      typedef std::vector<size_type, std::allocator<size_type> > offsets_type;

      struct Sample
      {
	typedef const feature_value_type* const_iterator;
      
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
	return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
      }
    
      size_type size() const { return offsets.size() - 1; }
      bool empty() const { return offsets.size() == 1; }
    
      void shrink()
      {
	features_type(features).swap(features);
	offsets_type(offsets).swap(offsets);
      }
    
      void flush()
      {
      
      }
    
      features_type features;
      offsets_type  offsets;
    };
  
    typedef SampleSet sample_set_type;
    typedef std::vector<double, std::allocator<double> > loss_set_type;

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
    
    template <typename Iterator1, typename Iterator2, typename Features>
    void construct_pair(Iterator1 oiter, Iterator1 oiter_end,
			Iterator2 kiter, Iterator2 kiter_end,
			Features& feats)
    {
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
    }
    
    template <typename Iterator>
    void transform_pair(Iterator first, Iterator last)
    {
      feature_node_type feature;
      
      offsets.push_back(features.size());
      
      for (/**/; first != last; ++ first) {
	feature.index = first->first.id() + 1;
	feature.value = first->second;

	features.push_back(feature);
      }
      
      feature.index = -1;
      feature.value = 0.0;
      
      features.push_back(feature);
    }

    struct greater_loss
    {
      greater_loss(const loss_set_type&   __losses) : losses(__losses) {}

      bool operator()(const size_type& x, const size_type& y) const
      {
	return losses[x] > losses[y];
      }
      
      const loss_set_type&   losses;
    };
    
    void operator()()
    {
      typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
      typedef std::vector<size_type, std::allocator<size_type> > pos_set_type;
      
      offsets.clear();
      features.clear();
            
      sentence_unique_type  sentences;
      
      int id = 0;
      
      boost::mt19937 generator;
      generator.seed(utils::random_seed());
      boost::random_number_generator<boost::mt19937> gen(generator);
      
      features_type feats;
      
      pos_set_type    positions;
      sample_set_type features_sample;
      loss_set_type   losses_sample;

      for (;;) {
	queue.pop(id);
	if (id < 0) break;

	if (oracles[id].empty() || kbests[id].empty()) continue;
	
	if (sample_vector) {
	  // first, we collect instances from oracle <-> non-oracle pairs

	  features_sample.clear();
	  losses_sample.clear();
	  
	  sentences.clear();
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    sentences.insert(oracles[id][o].sentence);

	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    for (size_t k = 0; k != kbests[id].size(); ++ k) {
	      const hypothesis_type& oracle = oracles[id][o];
	      const hypothesis_type& kbest  = kbests[id][k];
	    
	      // ignore oracle translations
	      if (sentences.find(kbest.sentence) != sentences.end()) continue;

	      const double loss = kbest.loss - oracle.loss;
	      if (loss <= 0.0) continue;
	      
	      feats.clear();
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	      
	      if (feats.empty()) continue;
	      
	      features_sample.insert(feats.begin(), feats.end());
	      losses_sample.push_back(loss);
	    }
	  
	  // second, collect data from kbests onlly, which is the same as the first examples
	  const size_type sample_size = losses_sample.size();
	  const size_type sample_size_max = sample_size << 2;
	  
	  while (losses_sample.size() < sample_size_max) {
	    const hypothesis_type& hyp1 = kbests[id][gen(kbests[id].size())];
	    const hypothesis_type& hyp2 = kbests[id][gen(kbests[id].size())];
	    
	    const hypothesis_type& kbest  = (hyp1.loss < hyp2.loss ? hyp2 : hyp1);
	    const hypothesis_type& oracle = (hyp1.loss < hyp2.loss ? hyp1 : hyp2);
	    
	    const double loss = kbest.loss - oracle.loss;
	    if (loss <= 1e-4) continue;
	    
	    feats.clear();
	    construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	    
	    if (feats.empty()) continue;
	    
	    features_sample.insert(feats.begin(), feats.end());
	    losses_sample.push_back(loss);
	  }
	  
	  positions.clear();
	  for (size_type i = 0; i != losses_sample.size(); ++ i)
	    positions.push_back(i);
	  
	  std::sort(positions.begin(), positions.end(), greater_loss(losses_sample));
	  
	  for (pos_set_type::const_iterator piter = positions.begin(); piter != positions.begin() + sample_size; ++ piter)
	    transform_pair(features_sample[*piter].begin(), features_sample[*piter].end());
	  
	} else {
	  sentences.clear();
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    sentences.insert(oracles[id][o].sentence);
	  
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    for (size_t k = 0; k != kbests[id].size(); ++ k) {
	      const hypothesis_type& oracle = oracles[id][o];
	      const hypothesis_type& kbest  = kbests[id][k];
	      
	      // ignore oracle translations
	      if (sentences.find(kbest.sentence) != sentences.end()) continue;
	      
	      const double loss = kbest.loss - oracle.loss;
	      if (loss <= 0.0) continue;

	      feats.clear();
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);

	      if (feats.empty()) continue;
	      
	      transform_pair(feats.begin(), feats.end());
	    }
	}
      }
      
      feature_node_set_type(features).swap(features);
    }
  };
  typedef Encoder encoder_type;
  typedef std::vector<encoder_type, std::allocator<encoder_type> > encoder_set_type;

  
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
		 const weight_set_type& bounds_lower,
		 const weight_set_type& bounds_upper,
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
	  const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
	  
	  for (size_t i = 0; i != weights_size; ++ i)
	    weights[i] = k * weights[i] + (1.0 - k) * weights_prev[i];
	  for (size_t i = weights_size; i < weights.size(); ++ i)
	    weights[i] = k * weights[i];
	  for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	    weights[i] = (1.0 - k) * weights_prev[i];
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
	  const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
	  
	  for (size_t i = 0; i != weights_size; ++ i)
	    weights[i] = - k * weights[i] + (1.0 + k) * weights_prev[i];
	  for (size_t i = weights_size; i < weights.size(); ++ i)
	    weights[i] = - k * weights[i];
	  for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	    weights[i] = (1.0 + k) * weights_prev[i];
	}
      }
    }
  }
  
public:
  weight_set_type weights;
  double objective;
};

template <typename Optimizer>
double optimize_svm(const hypothesis_map_type& kbests,
		    const hypothesis_map_type& oracles,
		    const weight_set_type& bounds_lower,
		    const weight_set_type& bounds_upper,
		    weight_set_type& weights)
{
  Optimizer optimizer(kbests, oracles, bounds_lower, bounds_upper, weights);
  
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
      typedef const feature_value_type* const_iterator;
      
      Sample(const_iterator __first, const_iterator __last) : first(__first), last(__last) {}
      
      const_iterator begin() const { return first; }
      const_iterator end() const { return last; }
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    void flush()
    {
      
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
	    const hypothesis_map_type& __oracles,
	    const weight_set_type& __weights)
      : queue(__queue), kbests(__kbests), oracles(__oracles), weights(__weights) {}

    template <typename Iterator1, typename Iterator2, typename Features>
    void construct_pair(Iterator1 oiter, Iterator1 oiter_end,
			Iterator2 kiter, Iterator2 kiter_end,
			Features& feats)
    {
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
    }
    
    struct greater_loss
    {
      greater_loss(const loss_set_type&   __losses) : losses(__losses) {}

      bool operator()(const size_type& x, const size_type& y) const
      {
	return losses[x] > losses[y];
      }
      
      const loss_set_type&   losses;
    };

    void operator()()
    {
      typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
      typedef std::vector<size_type, std::allocator<size_type> > pos_set_type;
      
      features_type feats;
      sentence_unique_type  sentences;
      
      boost::mt19937 generator;
      generator.seed(utils::random_seed());
      boost::random_number_generator<boost::mt19937> gen(generator);
      
      pos_set_type    positions;
      sample_set_type features_sample;
      loss_set_type   losses_sample;
      
      int id = 0;
      
      for (;;) {
	queue.pop(id);
	if (id < 0) break;

	if (oracles[id].empty() || kbests[id].empty()) continue;
	
	if (sample_vector) {
	  features_sample.clear();
	  losses_sample.clear();
	  
	  sentences.clear();
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    sentences.insert(oracles[id][o].sentence);
	  
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    for (size_t k = 0; k != kbests[id].size(); ++ k) {
	      const hypothesis_type& oracle = oracles[id][o];
	      const hypothesis_type& kbest  = kbests[id][k];
	    
	      // ignore oracle translations
	      if (sentences.find(kbest.sentence) != sentences.end()) continue;

	      const double loss = kbest.loss - oracle.loss;
	      if (loss <= 0.0) continue;
	      
	      feats.clear();
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	      
	      if (feats.empty()) continue;
	      
	      features_sample.insert(feats.begin(), feats.end());
	      losses_sample.push_back(loss);
	    }
	  
	  // second, collect data from kbests onlly, which is the same as the first examples
	  const size_type sample_size = losses_sample.size();
	  const size_type sample_size_max = sample_size << 2;
	  
	  while (losses_sample.size() < sample_size_max) {
	    const hypothesis_type& hyp1 = kbests[id][gen(kbests[id].size())];
	    const hypothesis_type& hyp2 = kbests[id][gen(kbests[id].size())];
	    
	    const hypothesis_type& kbest  = (hyp1.loss < hyp2.loss ? hyp2 : hyp1);
	    const hypothesis_type& oracle = (hyp1.loss < hyp2.loss ? hyp1 : hyp2);
	    
	    const double loss = kbest.loss - oracle.loss;
	    if (loss <= 1e-4) continue;
	    
	    feats.clear();
	    construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	    
	    if (feats.empty()) continue;
	    
	    features_sample.insert(feats.begin(), feats.end());
	    losses_sample.push_back(loss);
	  }
	  
	  positions.clear();
	  for (size_type i = 0; i != losses_sample.size(); ++ i)
	    positions.push_back(i);
	  
	  std::sort(positions.begin(), positions.end(), greater_loss(losses_sample));

	  if (conservative_loss) {
	    for (pos_set_type::const_iterator piter = positions.begin(); piter != positions.begin() + sample_size; ++ piter) {
	      const double margin = cicada::dot_product(weights, features_sample[*piter].begin(), features_sample[*piter].end(), 0.0);
	      const double loss = losses_sample[*piter];

	      if (loss_margin) {
		if (loss - margin > 0.0) {
		  features.insert(features_sample[*piter].begin(), features_sample[*piter].end());
		  losses.push_back(loss - margin);
		}
	      } else {
		if (1.0 - margin > 0.0) {
		  features.insert(features_sample[*piter].begin(), features_sample[*piter].end());
		  losses.push_back(1.0 - margin);
		}
	      }
	    }
	  } else {
	    for (pos_set_type::const_iterator piter = positions.begin(); piter != positions.begin() + sample_size; ++ piter) {
	      features.insert(features_sample[*piter].begin(), features_sample[*piter].end());
	      losses.push_back(loss_margin ? losses_sample[*piter] : 1.0);
	    }
	  }
	} else {
	  
	  sentences.clear();
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    sentences.insert(oracles[id][o].sentence);
	
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    for (size_t k = 0; k != kbests[id].size(); ++ k) {
	      const hypothesis_type& oracle = oracles[id][o];
	      const hypothesis_type& kbest  = kbests[id][k];
	    
	      // ignore oracle translations
	      if (sentences.find(kbest.sentence) != sentences.end()) continue;
	      
	      const double loss = kbest.loss - oracle.loss;
	      if (loss <= 0.0) continue;
	      
	      feats.clear();
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	      
	      if (feats.empty()) continue;

	      if (conservative_loss) {
		const double margin = cicada::dot_product(weights, feats.begin(), feats.end(), 0.0);
		
		if (loss_margin) {
		  if (loss - margin > 0.0) {
		    features.insert(feats.begin(), feats.end());
		    losses.push_back(loss - margin);
		  }
		} else {
		  if (1.0 - margin > 0.0) {
		    features.insert(feats.begin(), feats.end());
		    losses.push_back(1.0 - margin);
		  }
		}
	      } else {
		if (loss_margin) {
		  features.insert(feats.begin(), feats.end());
		  losses.push_back(loss);
		} else {
		  features.insert(feats.begin(), feats.end());
		  losses.push_back(1.0);
		}
	      }
	    }
	}
	
	features.flush();
      }
      
      // shrinkt features...
      features.shrink();
    }
    
    queue_type& queue;
    const hypothesis_map_type& kbests;
    const hypothesis_map_type& oracles;
    const weight_set_type&     weights;
    
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
	    const encoder_set_type&  __encoders,
	    const weight_set_type& __bounds_lower,
	    const weight_set_type& __bounds_upper)
      : positions(__positions),
	encoders(__encoders),
	bounds_lower(__bounds_lower),
	bounds_upper(__bounds_upper) {}

    double operator()(int i, int j) const
    {
      const pos_pair_type& pos_i = positions[i];
      const pos_pair_type& pos_j = positions[j];
      
      return cicada::dot_product(encoders[pos_i.first].features[pos_i.second].begin(), encoders[pos_i.first].features[pos_i.second].end(),
				 encoders[pos_j.first].features[pos_j.second].begin(), encoders[pos_j.first].features[pos_j.second].end(),
				 0.0);
    }
    
    const pos_pair_set_type& positions;
    const encoder_set_type&  encoders;
    const weight_set_type&   bounds_lower;
    const weight_set_type&   bounds_upper;
  };
  
  struct MMatrix
  {
    MMatrix(const pos_pair_set_type& __positions,
	    const encoder_set_type&  __encoders,
	    const weight_set_type& __bounds_lower,
	    const weight_set_type& __bounds_upper)
      : positions(__positions),
	encoders(__encoders),
	bounds_lower(__bounds_lower),
	bounds_upper(__bounds_upper){}
    
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
    const weight_set_type&   bounds_lower;
    const weight_set_type&   bounds_upper;
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
	      const weight_set_type& bounds_lower,
	      const weight_set_type& bounds_upper,
	      weight_set_type& weights_prev)
    : weights(), objective(0.0), tolerance(0.1)
  {
    encoder_type::queue_type queue;
    encoder_set_type encoders(threads, encoder_type(queue, kbests, oracles, weights_prev));
    
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
    
    HMatrix H(positions, encoders, bounds_lower, bounds_upper);
    MMatrix M(positions, encoders, bounds_lower, bounds_upper);
    
    objective = solver(alpha, f, H, M, 1.0 / (C * data_size), tolerance, true); // we do not normalize alpha values for compatibility with liblinear
    objective *= C;
    
    size_type actives = 0;
    weights.clear();
    
    if (conservative_loss)
      weights = weights_prev;
    
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
	  const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
	  
	  for (size_t i = 0; i != weights_size; ++ i)
	    weights[i] = k * weights[i] + (1.0 - k) * weights_prev[i];
	  for (size_t i = weights_size; i < weights.size(); ++ i)
	    weights[i] = k * weights[i];
	  for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	    weights[i] = (1.0 - k) * weights_prev[i];
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
	  const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
	  
	  for (size_t i = 0; i != weights_size; ++ i)
	    weights[i] = - k * weights[i] + (1.0 + k) * weights_prev[i];
	  for (size_t i = weights_size; i < weights.size(); ++ i)
	    weights[i] = - k * weights[i];
	  for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	    weights[i] = (1.0 + k) * weights_prev[i];
	}
      }
    }
  }
  
public:
  weight_set_type weights;
  double objective;

  const double tolerance;
};

struct OptimizeXBLEU
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef hypothesis_type::feature_value_type feature_value_type;
  typedef utils::mulvector2<feature_value_type, std::allocator<feature_value_type> > sample_set_type;
  typedef std::vector<sample_set_type, std::allocator<sample_set_type> > sample_map_type;
  
  OptimizeXBLEU(const hypothesis_map_type& __kbests,
		const sample_map_type& __features_kbest,
		const scorer_document_type& __scorers,
		weight_set_type& __weights,
		const double& __lambda,
		const feature_type& __feature_scale)
    : kbests(__kbests),
      features_kbest(__features_kbest),
      scorers(__scorers),
      weights(__weights),
      lambda(__lambda),
      feature_scale(__feature_scale)
  { }
  
  const hypothesis_map_type& kbests;
  const sample_map_type&     features_kbest;
  const scorer_document_type& scorers;
  weight_set_type& weights;
  
  double lambda;
  const feature_type feature_scale;

  double objective_opt;
  weight_set_type weights_opt;
  
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
    
    objective_opt = std::numeric_limits<double>::infinity();
    double objective = 0.0;
    
    const int result = lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeXBLEU::evaluate, 0, this, &param);
    
    if (debug)
      std::cerr << "lbfgs: " << lbfgs_error(result) << std::endl;
    
    // copy from opt weights!
    if (result < 0)
      weights = weights_opt;
    
    return objective;
  }
  
  struct Task
  {
    typedef cicada::semiring::Log<double> weight_type;
    
    static weight_type brevity_penalty(const double x)
    {
      typedef cicada::semiring::traits<weight_type> traits_type;

      // return (std::exp(x) - 1) / (1.0 + std::exp(1000.0 * x)) + 1.0;
      
      return ((traits_type::exp(x) - traits_type::one()) / (traits_type::one() + traits_type::exp(1000.0 * x))) + traits_type::one();
    }
    
    static weight_type derivative_brevity_penalty(const double x)
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
       
      const weight_type expx     = traits_type::exp(x);
      const weight_type expxm1   = expx - traits_type::one();
      const weight_type exp1000x = traits_type::exp(1000.0 * x);
      const weight_type p1exp1000x = traits_type::one() + exp1000x;
      
      return (expx / p1exp1000x) - ((expxm1 * weight_type(1000.0) * exp1000x) / (p1exp1000x * p1exp1000x));
      
      //return expx / (1.0 + exp1000x) - boost::math::expm1(x) * (1000.0 * exp1000x) / ((1.0 + exp1000x) * (1.0 + exp1000x))
    }
    
    static weight_type clip_count(const weight_type& x, const weight_type& clip)
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
      
      //return (x - clip) / (1.0 + std::exp(1000.0 * (x - clip))) + clip;
      return (weight_type(x - clip) / (traits_type::one() + traits_type::exp(1000.0 * (x - clip)))) + weight_type(clip);
    }
    
    static weight_type derivative_clip_count(const weight_type& x, const weight_type& clip)
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
      
      const weight_type exp1000xmc = traits_type::exp(1000.0 * (x - clip));
      const weight_type p1exp1000xmc = exp1000xmc + traits_type::one();
      
      return (traits_type::one() / p1exp1000xmc) - ((weight_type(x - clip) * weight_type(1000.0) * exp1000xmc) / (p1exp1000xmc * p1exp1000xmc));
      
      //return 1.0 / (1.0 + exp1000x) - (x - clip) * (1000.0 * exp1000x) / ((1.0 + exp1000x) * (1.0 + exp1000x));
    }

    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef std::vector<gradient_type, std::allocator<gradient_type> > gradients_type;
    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
    
    typedef std::vector<double, std::allocator<double> > ngram_counts_type;
    typedef std::vector<weight_set_type, std::allocator<weight_set_type> > feature_counts_type;
    
    typedef std::vector<double, std::allocator<double> > margins_type;

    typedef cicada::FeatureVectorUnordered<weight_type, std::allocator<weight_type> > expectation_type;
    
    // queue...
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    
    Task(queue_type& __queue,
	 const hypothesis_map_type& __kbests,
	 const sample_map_type& __features_kbest,
	 const scorer_document_type& __scorers,
	 const weight_set_type& __weights,
	 const feature_type& __feature_scale)
      : queue(__queue),
	kbests(__kbests),
	features_kbest(__features_kbest),
	scorers(__scorers),
	weights(__weights),
	feature_scale(__feature_scale),
	c_matched(order + 1),
	c_hypo(order + 1),
	g_matched(order + 1),
	g_hypo(order + 1),
	g_reference(),
	g_entropy(),
	r(0),
	e(0)
    { }
    
    queue_type& queue;

    const hypothesis_map_type& kbests;
    const sample_map_type& features_kbest;
    const scorer_document_type& scorers;
    const weight_set_type& weights;
    const feature_type feature_scale;
    
    ngram_counts_type   c_matched;
    ngram_counts_type   c_hypo;
    feature_counts_type g_matched;
    feature_counts_type g_hypo;
    weight_set_type     g_reference;
    weight_set_type     g_entropy;
    double r;
    double e;
    
    void operator()()
    {
      weights_type    matched(order + 1);
      weights_type    hypo(order + 1);
      expectation_type expectation;
      
      weights_type    counts_matched(order + 1);
      weights_type    counts_hypo(order + 1);
      
      gradients_type gradients_matched(order + 1);
      gradients_type gradients_hypo(order + 1);
      
      weight_type    reference;
      weight_type    entropy;
      gradient_type  gradient_reference;
      gradient_type  gradient_entropy;
      
      margins_type  margins;

      for (size_t n = 1; n != g_matched.size(); ++ n) {
	gradients_matched[n].allocate();
	gradients_hypo[n].allocate();
	
	g_matched[n].clear();
	g_hypo[n].clear();
      }
      
      gradient_reference.allocate();
      gradient_entropy.allocate();
      g_reference.clear();
      g_entropy.clear();
      
      std::fill(counts_matched.begin(), counts_matched.end(), weight_type());
      std::fill(counts_hypo.begin(), counts_hypo.end(), weight_type());
      std::fill(c_matched.begin(), c_matched.end(), 0.0);
      std::fill(c_hypo.begin(), c_hypo.end(), 0.0);
      r = 0.0;
      e = 0.0;
      
      const double scale = weights[feature_scale];
      
      for (;;) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
	margins.clear();
	  
	std::fill(matched.begin(), matched.end(), weight_type());
	std::fill(hypo.begin(), hypo.end(), weight_type());
	expectation.clear();
	  
	weight_type Z;
	weight_type Z_reference;
	weight_type Z_entropy;
	weight_type dR;
	  
	// first pass... compute margin and Z
	for (size_type k = 0; k != kbests[id].size(); ++ k) {
	  const hypothesis_type& kbest = kbests[id][k];
	  
	  const double margin  = cicada::dot_product(weights, features_kbest[id][k].begin(), features_kbest[id][k].end(), 0.0);
	    
	  margins.push_back(margin);
	  Z += cicada::semiring::traits<weight_type>::exp(margin * scale);
	}

	// second pass... compute sums (counts, etc.)
	for (size_type k = 0; k != kbests[id].size(); ++ k) {
	  const hypothesis_type& kbest = kbests[id][k];
	    
	  const double& margin = margins[k];
	  const weight_type prob = cicada::semiring::traits<weight_type>::exp(margin * scale) / Z;
	    
	  const cicada::eval::Bleu* bleu = dynamic_cast<const cicada::eval::Bleu*>(kbest.score.get());
	  if (! bleu)
	    throw std::runtime_error("no bleu statistics?");
	    
	  // collect scaled bleu stats
	  for (int n = 1; n <= order; ++ n) {
	    if (n - 1 < bleu->ngrams_reference.size())
	      hypo[n] += prob * bleu->ngrams_reference[n - 1];
	    if (n - 1 < bleu->ngrams_hypothesis.size())
	      matched[n] += prob * bleu->ngrams_hypothesis[n - 1];
	  }
	    
	  // collect reference length
	  Z_reference += prob * bleu->length_reference;
	    
	  // collect entropy...
	  Z_entropy -= prob * cicada::semiring::log(prob);
	    
	  // collect expectation
	  sample_set_type::const_reference::const_iterator fiter_end = features_kbest[id][k].end();
	  for (sample_set_type::const_reference::const_iterator fiter = features_kbest[id][k].begin(); fiter != fiter_end; ++ fiter)
	    expectation[fiter->first] += prob * weight_type(fiter->second * scale);
	  
	  expectation[feature_scale] += prob * weight_type(margin);
	    
	  dR += weight_type(1.0 + cicada::semiring::log(prob)) * prob;
	}
	  
	// accumulate
	std::transform(hypo.begin(), hypo.end(), counts_hypo.begin(), counts_hypo.begin(), std::plus<weight_type>());
	std::transform(matched.begin(), matched.end(), counts_matched.begin(), counts_matched.begin(), std::plus<weight_type>());
	  
	reference += Z_reference;
	entropy += Z_entropy;
	  
	expectation_type::const_iterator eiter_end = expectation.end();
	for (expectation_type::const_iterator eiter = expectation.begin(); eiter != eiter_end; ++ eiter) {
	  // collect bleus...
	  for (int n = 1; n <= order; ++ n) {
	    gradients_hypo[n][eiter->first] -= eiter->second * hypo[n];
	    gradients_matched[n][eiter->first] -= eiter->second * matched[n];
	  }
	    
	  // reference lengths
	  gradient_reference[eiter->first] -= eiter->second * Z_reference;
	    
	  // entropy gradient...
	  gradient_entropy[eiter->first] -= - dR * expectation[eiter->first];
	}
	  
	// third pass, collect gradients...
	for (size_type k = 0; k != kbests[id].size(); ++ k) {
	  const hypothesis_type& kbest = kbests[id][k];
	    
	  const double& margin = margins[k];
	  const weight_type prob = cicada::semiring::traits<weight_type>::exp(margin * scale) / Z;
	  const cicada::eval::Bleu* bleu = dynamic_cast<const cicada::eval::Bleu*>(kbest.score.get());
	    
	  // collect feature expectations etc...
	  sample_set_type::const_reference::const_iterator fiter_end = features_kbest[id][k].end();
	  for (sample_set_type::const_reference::const_iterator fiter = features_kbest[id][k].begin(); fiter != fiter_end; ++ fiter) {
	    const weight_type value(fiter->second * scale);
	      
	    // bleu statistics
	    for (int n = 1; n <= order; ++ n) {
	      if (n - 1 < bleu->ngrams_reference.size())
		gradients_hypo[n][fiter->first] += value * prob * bleu->ngrams_reference[n - 1];
		
	      if (n - 1 < bleu->ngrams_hypothesis.size())
		gradients_matched[n][fiter->first] += value * prob * bleu->ngrams_hypothesis[n - 1];
	    }
	      
	    // reference lengths
	    gradient_reference[fiter->first] += value * prob * bleu->length_reference;
	      
	    // entropy: we will collect minus values!
	    gradient_entropy[fiter->first] += - weight_type(1.0 + cicada::semiring::log(prob)) * prob * value;
	  }
	    
	  const weight_type value_scale(margin);
	    
	  for (int n = 1; n <= order; ++ n) {
	    if (n - 1 < bleu->ngrams_reference.size())
	      gradients_hypo[n][feature_scale] += value_scale * prob * bleu->ngrams_reference[n - 1];
	      
	    if (n - 1 < bleu->ngrams_hypothesis.size())
	      gradients_matched[n][feature_scale] += value_scale * prob * bleu->ngrams_hypothesis[n - 1];
	  }
	    
	  gradient_reference[feature_scale] += value_scale * prob * bleu->length_reference;
	  gradient_entropy[feature_scale] += - weight_type(1.0 + cicada::semiring::log(prob)) * prob * value_scale;
	}
      }
      
      // copy from weight space to double space
      std::copy(counts_matched.begin(), counts_matched.end(), c_matched.begin());
      std::copy(counts_hypo.begin(), counts_hypo.end(), c_hypo.begin());
      
      for (size_t n = 1; n != g_matched.size(); ++ n) {
	g_matched[n].allocate();
	g_hypo[n].allocate();
	
	std::copy(gradients_matched[n].begin(), gradients_matched[n].end(), g_matched[n].begin());
	std::copy(gradients_hypo[n].begin(), gradients_hypo[n].end(), g_hypo[n].begin());
      }
      
      g_reference.allocate();
      g_entropy.allocate();
      std::copy(gradient_reference.begin(), gradient_reference.end(), g_reference.begin());
      std::copy(gradient_entropy.begin(), gradient_entropy.end(), g_entropy.begin());

      r = reference;
      e = entropy;
    }
  };
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int size,
				  const lbfgsfloatval_t step)
  {
    typedef Task                  task_type;
    typedef task_type::queue_type queue_type;
    
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

    OptimizeXBLEU& optimizer = *((OptimizeXBLEU*) instance);
    
    queue_type queue;
    task_set_type tasks(threads, task_type(queue,
					   optimizer.kbests,
					   optimizer.features_kbest,
					   optimizer.scorers,
					   optimizer.weights,
					   optimizer.feature_scale));
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    size_type instances = 0;
    for (size_t id = 0; id != optimizer.kbests.size(); ++ id)
      if (! optimizer.kbests[id].empty()) {
	queue.push(id);
	++ instances;
      }
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    // clear g...
    std::fill(g, g + size, 0.0);
    
    workers.join_all();
    
    // collect statistics
    task_type::ngram_counts_type c_matched(order + 1, 0.0);
    task_type::ngram_counts_type c_hypo(order + 1, 0.0);

    task_type::feature_counts_type g_matched(order + 1);
    task_type::feature_counts_type g_hypo(order + 1);
    weight_set_type g_reference;
    weight_set_type g_entropy;
    
    double r(0.0);
    double e(0.0);
      
    for (int i = 0; i < threads; ++ i) {
      std::transform(tasks[i].c_matched.begin(), tasks[i].c_matched.end(), c_matched.begin(), c_matched.begin(), std::plus<double>());
      std::transform(tasks[i].c_hypo.begin(), tasks[i].c_hypo.end(), c_hypo.begin(), c_hypo.begin(), std::plus<double>());
      
      for (int n = 1; n <= order; ++ n) {
	g_matched[n] += tasks[i].g_matched[n];
	g_hypo[n] += tasks[i].g_hypo[n];
      }
      
      g_reference += tasks[i].g_reference;
      g_entropy += tasks[i].g_entropy;
      
      r += tasks[i].r;
      e += tasks[i].e;
    }
    
    // smoothing..
    {
      double smoothing = 1e-40;
      for (int n = 1; n <= order; ++ n) {
	if (c_hypo[n] > 0.0 && c_matched[n] <= 0.0)
	  c_matched[n] = smoothing;
	smoothing *= 0.1;
      }
    }

    // compute P
    double P = 0.0;
    for (int n = 1; n <= order; ++ n)
      if (c_hypo[n] > 0.0)
	P += (1.0 / order) * (utils::mathop::log(c_matched[n]) - utils::mathop::log(c_hypo[n]));

    // compute C and B
    const double C = r / c_hypo[1];
    const double B = task_type::brevity_penalty(1.0 - C);
    
    // for computing g...
    const double exp_P = utils::mathop::exp(P);
    const double C_dC  = C * task_type::derivative_brevity_penalty(1.0 - C);
    
    // xBLEU...
    const double objective_bleu = exp_P * B;
    const double entropy = e / instances;
    
    // entropy
    if (regularize_entropy) {
      // 0.5 * (entropy - C2)^2
      // 
      // thus, derivative is: (entropy - C2) * \nabla entropy
      // we need to consider average of entropy...
      
      std::transform(g_entropy.begin(), g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), (entropy - C2) / instances));
    } else
      std::transform(g_entropy.begin(), g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), - temperature / instances));
    
    for (int n = 1; n <= order; ++ n) 
      if (c_hypo[n] > 0.0) {
	const double factor_matched = - (exp_P * B / order) / c_matched[n];
	const double factor_hypo    = - (exp_P * B / order) / c_hypo[n];
	
	for (size_t i = 0; i != static_cast<size_t>(size); ++ i) {
	  g[i] += factor_matched * g_matched[n][i];
	  g[i] -= factor_hypo * g_hypo[n][i];
	}
      }
    
    if (c_hypo[1] > 0.0) {
      const double factor_ref  = - exp_P * C_dC / r;
      const double factor_hypo = - exp_P * C_dC / c_hypo[1];
      
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i) {
	g[i] -= factor_ref  * g_reference[i];
	g[i] += factor_hypo * g_hypo[1][i];
      }
    }
    
    // we need to minimize negative bleu... + regularized by average entropy...
    double objective = - objective_bleu + (regularize_entropy ? 0.5 * (entropy - C2) * (entropy - C2) : - temperature * entropy);
    
    if (regularize_l2) {
      double norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(size); ++ i) {
	g[i] += optimizer.lambda * x[i] * double(i != optimizer.feature_scale.id());
	norm += x[i] * x[i] * double(i != optimizer.feature_scale.id());
      }
      objective += 0.5 * optimizer.lambda * norm;
    }

    if (scale_fixed)
      g[optimizer.feature_scale.id()] = 0.0;

    if (debug >= 2)
      std::cerr << "objective: " << objective
		<< " xBLEU: " << objective_bleu
		<< " BP: " << B
		<< " entropy: " << entropy
		<< " scale: " << optimizer.weights[optimizer.feature_scale]
		<< std::endl;
    
    
    // keep the best so forth...
    if (objective <= optimizer.objective_opt) {
      optimizer.objective_opt = objective;
      optimizer.weights_opt = optimizer.weights;
    }
    
    return objective;
  }
};


struct OptimizeLBFGS
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
      typedef const feature_value_type* const_iterator;
      
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    features_type features;
    offsets_type  offsets;
  };
  
  typedef SampleSet sample_set_type;
  
  struct sample_pair_type
  {
    typedef std::vector<double, std::allocator<double> > loss_set_type;

    sample_pair_type() : features(), loss(),  offset(0) {}
    sample_pair_type(const hypothesis_set_type& kbests,
		     const hypothesis_set_type& oracles)
      : features(), loss(),  offset(0)
    {
      loss.reserve(kbests.size() + oracles.size());

      hypothesis_set_type::const_iterator oiter_end = oracles.end();
      for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter) {
	features.insert(oiter->features.begin(), oiter->features.end());
	loss.push_back(oiter->loss);
      }
      
      offset = loss.size();
      
      hypothesis_set_type::const_iterator kiter_end = kbests.end();
      for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
	features.insert(kiter->features.begin(), kiter->features.end());
	loss.push_back(kiter->loss);
      }
      
      features.shrink();
    }
    
    size_type oracle_begin() const { return 0; }
    size_type oracle_end() const { return offset; }
    
    size_type kbest_begin() const { return offset; }
    size_type kbest_end() const { return loss.size(); }
    
    size_type size() const { return loss.size(); }
    
    sample_set_type features;
    loss_set_type   loss;
    size_type offset;
  };
  
  typedef std::vector<sample_pair_type, std::allocator<sample_pair_type> > sample_pair_set_type;
  
  OptimizeLBFGS(const hypothesis_map_type& kbests,
		const hypothesis_map_type& oracles,
		weight_set_type& __weights)
    : weights(__weights)
  {
    // transform into sample-pair-set-type
    
    const size_t id_max = utils::bithack::min(kbests.size(), oracles.size());
    
    samples.reserve(id_max);
    for (size_t id = 0; id != id_max; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty())
	samples.push_back(sample_pair_type(kbests[id], oracles[id]));
  }

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
    
    objective_opt = std::numeric_limits<double>::infinity();
    double objective = 0.0;
    
    const int result = lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    if (debug)
      std::cerr << "lbfgs: " << lbfgs_error(result) << std::endl;
    
    // copy from opt weights!
    if (result < 0)
      weights = weights_opt;
    
    return objective;
  }

  
  struct Task
  {
    typedef hypothesis_type::feature_value_type feature_value_type;
    
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;

    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    
    Task(queue_type&            __queue,
	 const weight_set_type& __weights,
	 const sample_pair_set_type& __samples,
	 const size_t& __instances)
      : queue(__queue),
	weights(__weights),
	samples(__samples),
	instances(__instances)
    {}
    

    void operator()()
    {
      typedef std::vector<double, std::allocator<double> > margin_set_type;

      g.clear();
      objective = 0.0;

      expectation_type  expectations;
      
      expectations.allocate();
      expectations.clear();
      
      margin_set_type margins;
      
      const double cost_factor = (softmax_margin ? 1.0 : 0.0);
      
      while (1) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
	weight_type Z_oracle;
	weight_type Z_kbest;
	
	margins.clear();
	margins.resize(samples[id].size());
	
	for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const double loss = samples[id].loss[i];
	  
	  margins[i] = cicada::dot_product(weights, features.begin(), features.end(), cost_factor * loss);
	  
	  Z_oracle += cicada::semiring::traits<weight_type>::exp(margins[i]);
	}
	
	for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const double loss = samples[id].loss[i];
	  
	  margins[i] = cicada::dot_product(weights, features.begin(), features.end(), cost_factor * loss);
	  
	  Z_kbest += cicada::semiring::traits<weight_type>::exp(margins[i]);
	}
	
	for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(margins[i]) / Z_oracle;
	  
	  sample_set_type::value_type::const_iterator fiter_end = features.end();
	  for (sample_set_type::value_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
	    expectations[fiter->first] -= weight_type(fiter->second) * weight;
	}
	
	for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(margins[i]) / Z_kbest;
	  
	  sample_set_type::value_type::const_iterator fiter_end = features.end();
	  for (sample_set_type::value_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
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
    const sample_pair_set_type& samples;
    
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

    const size_t instances = optimizer.samples.size();
        
    queue_type queue;
    
    task_set_type tasks(threads, task_type(queue, optimizer.weights, optimizer.samples, instances));

    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    for (int id = 0; id != instances; ++ id)
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
    
    // keep the best so forth...
    if (objective <= optimizer.objective_opt) {
      optimizer.objective_opt = objective;
      optimizer.weights_opt = optimizer.weights;
    }

    return objective;
  }
  
  sample_pair_set_type samples;
  
  weight_set_type& weights;
  
  double objective_opt;
  weight_set_type weights_opt;
};

template <typename Optimizer>
double optimize_xbleu(const hypothesis_map_type& kbests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights)
{
  const feature_type feature_scale(":feature-scale:");
  
  weights[feature_scale] = scale;
  
  typename Optimizer::sample_map_type features(kbests.size());
  for (size_t id = 0; id != kbests.size(); ++ id)
    if (! kbests[id].empty()) {
      
      for (size_t k = 0; k != kbests[id].size(); ++ k)
	features[id].push_back(kbests[id][k].features.begin(), kbests[id][k].features.end());
    }

  Optimizer optimizer(kbests, features, scorers, weights, C, feature_scale);
  
  double objective = 0.0;
  
  if (annealing_mode) {
    for (temperature = temperature_start; temperature >= temperature_end; temperature *= temperature_rate) {
      if (debug >= 2)
	std::cerr << "temperature: " << temperature << std::endl;
	
      objective = optimizer();
    }
  } else 
    objective = optimizer();
    
  if (quenching_mode) {
    temperature = 0.0;
      
    for (double quench = quench_start; quench <= quench_end; quench *= quench_rate) {
      if (debug >= 2)
	std::cerr << "quench: " << quench << std::endl;
	
      weights[feature_scale] = quench;
	
      objective = optimizer();
    }
  }
    
  return objective;
}


template <typename Optimizer>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights)
{
  return Optimizer(kbests, oracles, weights)();
}

struct EnvelopeTask
{
  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;

  typedef cicada::semiring::Envelope envelope_type;
  typedef std::vector<envelope_type, std::allocator<envelope_type> >  envelope_set_type;

  typedef utils::lockfree_list_queue<int, std::allocator<int> >  queue_type;
  
  EnvelopeTask(queue_type& __queue,
	       segment_document_type&      __segments,
	       const weight_set_type&      __origin,
	       const weight_set_type&      __direction,
	       const hypothesis_map_type&  __kbests,
	       const kbest_map_type&       __kbest_map)
    : queue(__queue),
      segments(__segments),
      origin(__origin),
      direction(__direction),
      kbests(__kbests),
      kbest_map(__kbest_map) {}

  void operator()()
  {
    EnvelopeKBest::line_set_type lines;
    int seg;
    
    EnvelopeKBest envelopes(origin, direction);
    
    while (1) {
      queue.pop(seg);
      if (seg < 0) break;
      
      lines.clear();
      
      for (size_t i = 0; i != kbest_map.size(); ++ i)
	if (kbest_map[i] == seg)
	  envelopes(kbests[i].begin(), kbests[i].end(), std::back_inserter(lines));

      envelopes(lines);
      
      EnvelopeKBest::line_set_type::const_iterator liter_end = lines.end();
      for (EnvelopeKBest::line_set_type::const_iterator liter = lines.begin(); liter != liter_end; ++ liter) {
	const EnvelopeKBest::line_type& line = *liter;
	
	if (debug >= 4)
	  std::cerr << "segment: " << seg << " x: " << line.x << std::endl;
	
	segments[seg].push_back(std::make_pair(line.x, line.hypothesis->score));
      }
    }
  }
  
  queue_type& queue;
  
  segment_document_type& segments;
  
  const weight_set_type& origin;
  const weight_set_type& direction;
  
  const hypothesis_map_type&  kbests;
  const kbest_map_type&       kbest_map;
};

double optimize_mert(const scorer_document_type& scorers,
		     const hypothesis_map_type& kbests,
		     const kbest_map_type& kbest_map,
		     const weight_set_type& weights_prev,
		     weight_set_type& weights)
{
  typedef EnvelopeTask task_type;
  typedef task_type::queue_type queue_type;

  typedef task_type::line_search_type line_search_type;
  
  typedef line_search_type::value_type optimum_type;

  if (kbest_map.empty()) return 0.0;
  
  const weight_set_type& origin = weights_prev;
  weight_set_type direction = weights;
  direction -= weights_prev;
  
  const size_t segment_max = *std::max_element(kbest_map.begin(), kbest_map.end());
  
  task_type::segment_document_type segments(segment_max);
  queue_type queue;
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, segments, origin, direction, kbests, kbest_map)));
  
  for (size_t seg = 0; seg != segment_max; ++ seg)
    queue.push(seg);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
  
  line_search_type line_search;
  
  const optimum_type optimum = line_search(segments, 0.1, 1.1, scorers.error_metric());
  
  const double update = (optimum.lower + optimum.upper) * 0.5;
  
  if (update != 0.0) {
    direction *= update;
    weights = origin;
    weights += direction;
  }

  if (debug >= 2)
    std::cerr << "mert update: " << update << std::endl;
  
  return optimum.objective;
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
  typedef boost::fusion::tuple<path_type, path_type, size_t, size_t> path_pair_type;
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
      const size_t mappos = boost::fusion::get<3>(paths);
      
      if (! scorers.empty())
	if (refpos >= scorers.size())
	  throw std::runtime_error("reference positions outof index");
      
      kbest_map.push_back(mappos);
      
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
  kbest_map_type      kbest_map;
};

void read_kbest(const scorer_document_type& scorers,
		const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles,
		kbest_map_type& kbest_map)
{
  kbest_map.clear();

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
    
    // assign kbest-map
    kbest_map.reserve(kbests_size);
    kbest_map.resize(kbests_size);
    for (size_t seg = 0; seg != kbests_size; ++ seg)
      kbest_map[seg] = seg;

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
	
	queue.push(path_pair_type(path_kbest, path_oracle, refpos ++, i));
      }
    }
    
    for (int i = 0; i != threads; ++ i)
      queue.push(path_pair_type(path_type(), path_type(), size_t(-1), size_t(-1)));
    
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
    kbest_map.reserve(kbests_size);
    
    for (int i = 0; i != threads; ++ i) {
      kbests.insert(kbests.end(), tasks[i].kbests.begin(), tasks[i].kbests.end());
      oracles.insert(oracles.end(), tasks[i].oracles.begin(), tasks[i].oracles.end());
      kbest_map.insert(kbest_map.end(), tasks[i].kbest_map.begin(), tasks[i].kbest_map.end());
      
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
    ("weights-history", po::value<path_set_type>(&weights_history_path)->multitoken(), "parameter history")
    ("output",  po::value<path_type>(&output_path),                    "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")

    ("bound-lower", po::value<path_type>(&bound_lower_file),                     "lower bounds definition for feature weights")
    ("bound-upper", po::value<path_type>(&bound_upper_file),                     "upper bounds definition for feature weights")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",  po::bool_switch(&learn_lbfgs),  "batch LBFGS algorithm")
    ("learn-xbleu",   po::bool_switch(&learn_xbleu),   "xBLEU algorithm")
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
    
    ("regularize-l1",      po::bool_switch(&regularize_l1),      "L1-regularization")
    ("regularize-l2",      po::bool_switch(&regularize_l2),      "L2-regularization")
    ("regularize-entropy", po::bool_switch(&regularize_entropy), " entropy regularization")
    
    ("C",             po::value<double>(&C)->default_value(C), "regularization constant")
    ("C2",            po::value<double>(&C2)->default_value(C2),       "an alternative regularization constant")
    ("scale",         po::value<double>(&scale)->default_value(scale), "scaling for weight")
    ("eta0",          po::value<double>(&eta0),                        "\\eta_0 for decay")
    ("eps",           po::value<double>(&eps),                 "tolerance for liblinear")
    ("order",         po::value<int>(&order)->default_value(order),    "ngram order for xBLEU")

    ("annealing", po::bool_switch(&annealing_mode), "annealing")
    ("quenching", po::bool_switch(&quenching_mode), "quenching")
    
    ("temperature",       po::value<double>(&temperature)->default_value(temperature),             "temperature")
    ("temperature-start", po::value<double>(&temperature_start)->default_value(temperature_start), "start temperature for annealing")
    ("temperature-end",   po::value<double>(&temperature_end)->default_value(temperature_end),     "end temperature for annealing")
    ("temperature-rate",  po::value<double>(&temperature_rate)->default_value(temperature_rate),   "annealing rate")

    ("quench-start", po::value<double>(&quench_start)->default_value(quench_start), "start quench for annealing")
    ("quench-end",   po::value<double>(&quench_end)->default_value(quench_end),     "end quench for annealing")
    ("quench-rate",  po::value<double>(&quench_rate)->default_value(quench_rate),   "quenching rate")

    ("loss-margin",       po::bool_switch(&loss_margin),       "direct loss margin")
    ("softmax-margin",    po::bool_switch(&softmax_margin),    "softmax margin")
    ("line-search",       po::bool_switch(&line_search),       "perform line search in each iteration")
    ("mert-search",       po::bool_switch(&mert_search),       "perform one-dimensional mert")
    ("sample-vector",     po::bool_switch(&sample_vector),     "perform samling")
    ("direct-loss",       po::bool_switch(&direct_loss),       "compute loss by directly treating hypothesis score")
    ("conservative-loss", po::bool_switch(&conservative_loss), "conservative loss")
    
    ("scale-fixed", po::bool_switch(&scale_fixed), "fixed scaling")
    
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
