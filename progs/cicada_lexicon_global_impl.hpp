//
//  Copyright(C) 2009-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <algorithm>
#include <vector>
#include <set>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>

#include <utils/compress_stream.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/repository.hpp>
#include <utils/resource.hpp>
#include <utils/mathop.hpp>
#include <utils/random_seed.hpp>
#include <utils/dense_hash_map.hpp>
#include <utils/dense_hash_set.hpp>

#include "cicada/symbol.hpp"
#include "cicada/vocab.hpp"
#include "cicada/sentence.hpp"

#include "lbfgs.h"

#include "liblinear/linear.h"

typedef boost::filesystem::path path_type;

typedef cicada::Symbol   word_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

struct bitext_type
{
  sentence_type source;
  sentence_type target;

  bitext_type()
    : source(), target() {}
  bitext_type(const sentence_type& __source, const sentence_type& __target)
    : source(__source), target(__target) {}
};

typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
typedef std::vector<word_type, std::allocator<word_type> > word_set_type;

// I know, it is very bad...
extern int    max_iteration;
extern bool   regularize_l1;
extern bool   regularize_l2;
extern double C;

extern int debug;

template <typename Tp>
struct greater_second
{
  bool operator()(const Tp& x, const Tp& y) const
  {
    return x.second > y.second;
  }
};

void read_bitexts(const path_type& path_source,
		  const path_type& path_target,
		  bitext_set_type& bitexts,
		  word_set_type& vocab,
		  const int kbest)
{
  typedef google::dense_hash_map<word_type, size_t, boost::hash<word_type>, std::equal_to<word_type> > word_count_set_type;
  typedef std::vector<std::pair<word_type, size_t>, std::allocator<std::pair<word_type, size_t> > > word_count_sorted_type;
  
  typedef std::set<word_type, std::less<word_type>, std::allocator<word_type> > word_sorted_set_type;

  bitexts.clear();
  
  word_count_set_type  counts;
  word_sorted_set_type uniques_source;

  counts.set_empty_key(word_type());
  
  utils::compress_istream is_src(path_source, 1024 * 1024);
  utils::compress_istream is_trg(path_target, 1024 * 1024);
  
  sentence_type      source;
  sentence_type      target;
  
  while (is_src && is_trg) {
    is_src >> source;
    is_trg >> target;
    
    if (! is_src || ! is_trg) break;
    
    if (source.empty() || target.empty()) continue;
    
    uniques_source.clear();
    uniques_source.insert(source.begin(), source.end());
    
    bitexts.push_back(bitext_type(sentence_type(uniques_source.begin(), uniques_source.end()), target));
    
    sentence_type::const_iterator titer_end = target.end();
    for (sentence_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
      ++ counts[*titer];
  }
  if (is_src || is_trg)
    throw std::runtime_error("# of sentences does not match...");
  
  bitext_set_type(bitexts).swap(bitexts);

  vocab.clear();
  if (kbest <= 0 || static_cast<int>(counts.size()) <= kbest) {
    word_count_set_type::const_iterator citer_end = counts.end();
    for (word_count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      vocab.push_back(citer->first);
  } else {
    word_count_sorted_type sorted(counts.begin(), counts.end());
    
    std::nth_element(sorted.begin(), sorted.begin() + kbest, sorted.end(), greater_second<word_count_sorted_type::value_type>());
    
    const size_t cutoff = sorted[kbest].second;
    
    word_count_sorted_type::const_iterator siter      = sorted.begin();
    word_count_sorted_type::const_iterator siter_end  = sorted.end();
    word_count_sorted_type::const_iterator siter_last = sorted.begin() + kbest;
    
    bool found_equal = false;
    for (/**/; siter != siter_last; ++ siter) {
      vocab.push_back(siter->first);
      found_equal |= (siter->second == cutoff);
    }
    
    if (found_equal)
      for (/**/; siter != siter_end && siter->second == cutoff; ++ siter)
	vocab.push_back(siter->first);
  }
  
  word_set_type(vocab).swap(vocab);
}

struct OptimizerLinearBase
{
  typedef std::vector<double, std::allocator<double> > parameter_set_type;

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
  
  
  OptimizerLinearBase(const bitext_set_type& bitexts,
		      const word_type& word)
  {
    const size_t bias_index = word_type::allocated() + 1;
    offset_set_type offsets;

    bitext_set_type::const_iterator biter_end = bitexts.end();
    for (bitext_set_type::const_iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {

      feature_node_type feat;
      
      sentence_type::const_iterator titer_end = biter->target.end();
      for (sentence_type::const_iterator titer = biter->target.begin(); titer != titer_end; ++ titer) {
	
	labels.push_back(*titer == word ? 1 : -1);
	offsets.push_back(feature_nodes.size());
	
	sentence_type::const_iterator siter_end = biter->source.end();
	for (sentence_type::const_iterator siter = biter->source.begin(); siter != siter_end; ++ siter) {
	  feat.index = siter->id() + 1;
	  feat.value = 1;
	  
	  feature_nodes.push_back(feat);
	}
	
	// bias
	feat.index = bias_index;
	feat.value = 1;
	feature_nodes.push_back(feat);
	
	// termination
	feat.index = -1;
	feat.value = 0;
	feature_nodes.push_back(feat);
      }
    }

    // setup labels and features based on feature_nodes and offsets
    label_set_type(labels).swap(labels);
    feature_node_set_type(feature_nodes).swap(feature_nodes);
    
    features.reserve(offsets.size());
    for (size_type pos = 0; pos != offsets.size(); ++ pos)
      features.push_back(const_cast<feature_node_type*>(&(*feature_nodes.begin())) + offsets[pos]);
  }
  
  static void print_string_stderr(const char *s)
  {
    std::cerr << s << std::flush;
  }
  
  static void print_string_none(const char *s)
  {
    
  }

  label_set_type labels;
  feature_node_map_type features;
  feature_node_set_type feature_nodes;
};

template <int Solver>  
struct OptimizerLinear : public OptimizerLinearBase
{
  typedef OptimizerLinearBase base_type;

  OptimizerLinear(const bitext_set_type& bitexts,
		  const word_type& word)
    : OptimizerLinearBase(bitexts, word) {}
  
  void operator()(parameter_set_type& x)
  {
    problem_type problem;
    
    problem.l = labels.size();
    problem.n = word_type::allocated() + 1; // + 1 for bias
    problem.y = &(*labels.begin());
    problem.x = &(*features.begin());
    problem.bias = 1;
    
    parameter_type parameter;
    parameter.solver_type = Solver;
    parameter.eps = std::numeric_limits<double>::infinity();
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
    
    const model_type* model = train(&problem, &parameter);

    const size_t vocabulary_size = word_type::allocated();
    const word_type::id_type id_bias = vocab_type::EPSILON.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);
    
    // is this correct???
    x[id_bias] = model->w[model->nr_feature];
    
    const double scale = model->label[0];
    for (int j = 0; j != model->nr_feature; ++ j)
      if (model->w[j] != 0.0)
	x[j] = model->w[j] * scale;
    
    free_and_destroy_model(const_cast<model_type**>(&model));
  }
};

struct Optimizer
{
  typedef std::vector<double, std::allocator<double> > parameter_set_type;
  
  const bitext_set_type bitexts;
  const word_type word;
  
  Optimizer(const bitext_set_type& __bitexts,
	    const word_type& __word)
    : bitexts(__bitexts.size()),
      word(__word)
  {
    typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> >  word_set_type;
    
    word_set_type cooc;
    cooc.set_empty_key(word_type());
    
    bitext_set_type::const_iterator biter_end = __bitexts.end();
    for (bitext_set_type::const_iterator biter = __bitexts.begin(); biter != biter_end; ++ biter)
      if (std::find(biter->target.begin(), biter->target.end(), word) != biter->target.end())
	cooc.insert(biter->source.begin(), biter->source.end());
    
    sentence_type source;
    bitext_set_type& bitexts_cooc = const_cast<bitext_set_type&>(bitexts);
    
    bitexts_cooc.clear();
    for (bitext_set_type::const_iterator biter = __bitexts.begin(); biter != biter_end; ++ biter) {
      source.clear();
      sentence_type::const_iterator siter_end = biter->source.end();
      for (sentence_type::const_iterator siter = biter->source.begin(); siter != siter_end; ++ siter)
	if (cooc.find(*siter) != cooc.end())
	  source.push_back(*siter);
      
      bitexts_cooc.push_back(bitext_type(source, biter->target));
    }
  }
};

struct OptimizeAROW : public Optimizer
{
  typedef Optimizer base_type;
  
  OptimizeAROW(const bitext_set_type& __bitexts,
	     const word_type& __word)
    : base_type(__bitexts, __word) {}
  
  void operator()(parameter_set_type& x)
  {
    typedef std::vector<int, std::allocator<int> > position_set_type;

    const size_t vocabulary_size = word_type::allocated();
    
    const word_type::id_type id_bias = vocab_type::EPSILON.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);
    
    parameter_set_type v(vocabulary_size, 1.0);
    
    size_t sample_size = 0;
    position_set_type positions(bitexts.size());
    for (size_t pos = 0; pos < positions.size(); ++ pos) {
      positions[pos] = pos;
      sample_size += bitexts[pos].target.size();
    }

    boost::mt19937 gen;
    gen.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> generator(gen);
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end(), generator);

      double margin_sum = 0.0;
      size_t num_correct = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sentence_type::const_iterator siter_begin = bitext.source.begin();
	sentence_type::const_iterator siter_end   = bitext.source.end();
	
	sentence_type::const_iterator titer_end = bitext.target.end();
	for (sentence_type::const_iterator titer = bitext.target.begin(); titer != titer_end; ++ titer) {
	  const double y = (int(*titer == word) * 2 - 1);
	  
	  double var    = v[id_bias];
	  double margin = x[id_bias];
	  for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	    margin += x[siter->id()];
	    var    += v[siter->id()];
	  }
	  
	  const double score = margin * y;
	  
	  num_correct += bool(score > 1.0);
	  margin_sum += score;
	  
	  if (score < 1.0) {
	    const double beta = 1.0 / (var + C);
	    const double alpha = std::max(0.0, 1.0 - score) * beta;
	    
	    x[id_bias] += alpha * y * v[id_bias];
	    v[id_bias] -= beta * v[id_bias] * v[id_bias];
	    
	    for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	      x[siter->id()] += alpha * y * v[siter->id()];
	      v[siter->id()] -= beta * v[siter->id()] * v[siter->id()];
	    }
	  }
	}
      }
      
      if (debug >= 2)
	std::cerr << "word: " << word
		  << " margin: " << (margin_sum / sample_size)
		  << " ratio: " << (double(num_correct) / sample_size)
		  << std::endl;
      
    }
  }
};

struct OptimizeCW : public Optimizer
{
  typedef Optimizer base_type;
  
  OptimizeCW(const bitext_set_type& __bitexts,
	     const word_type& __word)
    : base_type(__bitexts, __word) {}
  
  void operator()(parameter_set_type& x)
  {
    typedef std::vector<int, std::allocator<int> > position_set_type;

    const size_t vocabulary_size = word_type::allocated();
    
    const word_type::id_type id_bias = vocab_type::EPSILON.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);
    
    parameter_set_type v(vocabulary_size, 1.0);
    
    size_t sample_size = 0;
    position_set_type positions(bitexts.size());
    for (size_t pos = 0; pos < positions.size(); ++ pos) {
      positions[pos] = pos;
      sample_size += bitexts[pos].target.size();
    }
    
    boost::mt19937 gen;
    gen.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> generator(gen);

    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end(), generator);
      
      double margin_sum = 0.0;
      size_t num_correct = 0;

      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sentence_type::const_iterator siter_begin = bitext.source.begin();
	sentence_type::const_iterator siter_end   = bitext.source.end();
	
	sentence_type::const_iterator titer_end = bitext.target.end();
	for (sentence_type::const_iterator titer = bitext.target.begin(); titer != titer_end; ++ titer) {
	  const double y = (int(*titer == word) * 2 - 1);
	  
	  double var    = v[id_bias];
	  double margin = x[id_bias];
	  for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	    margin += x[siter->id()];
	    var    += v[siter->id()];
	  }
	  
	  const double score = margin * y;
	  const double b = 1.0 + 2.0 * C * score;
	  const double alpha = (- b + std::sqrt(b * b - 8.0 * C * (score - C * var))) / (4.0 * C * var);
	  const double beta  = (2.0 * alpha * C) / (1.0 + 2.0 * alpha * C * var);
	  
	  num_correct += bool(score > 1.0);
	  margin_sum += score;
	  
	  if (alpha > 0.0) {
	    x[id_bias] += alpha * y * v[id_bias];
	    v[id_bias] -= beta * v[id_bias] * v[id_bias];
	    
	    for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	      x[siter->id()] += alpha * y * v[siter->id()];
	      v[siter->id()] -= beta * v[siter->id()] * v[siter->id()];
	    }
	  }
	}
      }
     
      if (debug >= 2)
	std::cerr << "word: " << word
		  << " margin: " << (margin_sum / sample_size)
		  << " ratio: " << (double(num_correct) / sample_size)
		  << std::endl;
    }
  }
};


struct OptimizeMIRA : public Optimizer
{
  typedef Optimizer base_type;
  
  OptimizeMIRA(const bitext_set_type& __bitexts,
	       const word_type& __word)
    : base_type(__bitexts, __word) {}
  
  void operator()(parameter_set_type& x)
  {
    typedef std::vector<int, std::allocator<int> > position_set_type;

    const size_t vocabulary_size = word_type::allocated();
    
    const word_type::id_type id_bias = vocab_type::EPSILON.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);

    size_t sample_size = 0;
    position_set_type positions(bitexts.size());
    for (size_t pos = 0; pos < positions.size(); ++ pos) {
      positions[pos] = pos;
      sample_size += bitexts[pos].target.size();
    }
    
    boost::mt19937 gen;
    gen.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> generator(gen);

    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end(), generator);
      
      double margin_sum = 0.0;
      size_t num_correct = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sentence_type::const_iterator siter_begin = bitext.source.begin();
	sentence_type::const_iterator siter_end   = bitext.source.end();
	
	sentence_type::const_iterator titer_end = bitext.target.end();
	for (sentence_type::const_iterator titer = bitext.target.begin(); titer != titer_end; ++ titer) {
	  const double y = (int(*titer == word) * 2 - 1);
	  
	  double margin = x[id_bias];
	  for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	    margin += x[siter->id()];
	  
	  const double score = margin * y;

	  num_correct += bool(score > 1.0);
	  margin_sum += score;
	  
	  if (score <= 1.0) {
	    const double norm = (siter_end - siter_begin) + 1;
	    const double alpha = std::min(1.0 / C, y * (1.0 - score) / norm);
	    
	    x[id_bias] += alpha;
	    for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	      x[siter->id()] += alpha;
	  }
	}
      }
      
      if (debug >= 2)
	std::cerr << "word: " << word
		  << " margin: " << (margin_sum / sample_size)
		  << " ratio: " << (double(num_correct) / sample_size)
		  << std::endl;
    }
  }
};

struct OptimizeSGD : public Optimizer
{
  typedef Optimizer base_type;
  
  OptimizeSGD(const bitext_set_type& __bitexts,
	      const word_type& __word)
    : base_type(__bitexts, __word) {}

  void operator()(parameter_set_type& x)
  {
    if (regularize_l1)
      optimize_l1(x);
    else
      optimize_l2(x);
  }

  
  void apply(double& x, double& penalty, const double& cummulative)
  {
    const double x_half = x;
    if (x > 0.0)
      x = std::max(0.0, x - penalty - cummulative);
    else if (x < 0.0)
      x = std::min(0.0, x - penalty + cummulative);
    penalty += x - x_half;
  }
  
  void optimize_l1(parameter_set_type& x)
  {
    typedef std::vector<int, std::allocator<int> > position_set_type;

    const size_t vocabulary_size = word_type::allocated();
    
    const word_type::id_type id_bias = vocab_type::EPSILON.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);

    size_t sample_size = 0;
    position_set_type positions(bitexts.size());
    for (size_t pos = 0; pos < positions.size(); ++ pos) {
      positions[pos] = pos;
      sample_size += bitexts[pos].target.size();
    }

    size_t epoch = 0;
    const double lambda = C / sample_size;
    
    double penalty = 0.0;
    parameter_set_type penalties(x.size(), 0.0);
    
    boost::mt19937 gen;
    gen.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> generator(gen);

    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end(), generator);
      
      double objective = 0.0;
      double margin_sum = 0.0;
      size_t num_correct_classify = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sentence_type::const_iterator siter_begin = bitext.source.begin();
	sentence_type::const_iterator siter_end   = bitext.source.end();
	
	//const double eta = 1.0 / (1.0 + epoch / sample_size);
	const double eta = 0.2 * std::pow(0.85, double(epoch) / sample_size);
	++ epoch;
	penalty += eta * lambda * bitext.target.size();
	
	double margin = x[id_bias];
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  margin += x[siter->id()];

	const double exp_margin = utils::mathop::exp(margin);
	const double inv_exp_margin = 1.0 / exp_margin;
	
	const double objective_correct   = boost::math::log1p(inv_exp_margin);
	const double objective_incorrect = boost::math::log1p(exp_margin);
	
	const double gradient_correct   =   1.0 / (1.0 + exp_margin);
	const double gradient_incorrect = - 1.0 / (1.0 + inv_exp_margin);
	
	const size_t num_correct   = std::count(bitext.target.begin(), bitext.target.end(), word);
	const size_t num_incorrect = bitext.target.size() - num_correct;
	
	margin_sum           += margin * num_correct - margin * num_incorrect;
	num_correct_classify += (margin > 1.0) * num_correct + (- margin > 1.0) * num_incorrect;
	objective += objective_correct * num_correct + objective_incorrect * num_incorrect;
	
	const double gradient = gradient_correct * num_correct + gradient_incorrect * num_incorrect;
	
	x[id_bias] += eta * gradient;
	apply(x[id_bias], penalties[id_bias], penalty);
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	  x[siter->id()] += eta * gradient;
	  apply(x[siter->id()], penalties[siter->id()], penalty);
	}
      }
      
      if (debug >= 2)
	std::cerr << "word: " << word
		  << " objective: " << objective
		  << " margin: " << (margin_sum / sample_size)
		  << " ratio: " << (double(num_correct_classify) / sample_size)
		  << std::endl;
    }
  }

  void rescale(parameter_set_type& x, double& scale, double& norm, const double& alpha)
  {
    norm *= alpha * alpha;
    if (alpha != 0.0)
      scale *= alpha;
    else {
      scale = 1.0;
      std::fill(x.begin(), x.end(), 0.0);
    }
  }

  void update(double& x, double& scale, double& norm, const double& alpha)
  {
    norm += 2.0 * x * alpha * scale + alpha * alpha;
    x += alpha / scale;
  }
  
  void optimize_l2(parameter_set_type& x)
  {
    // SGD, inspired by Pegasos

    typedef std::vector<int, std::allocator<int> > position_set_type;

    const size_t vocabulary_size = word_type::allocated();
    
    const word_type::id_type id_bias = vocab_type::EPSILON.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);

    double weight_norm = 0.0;
    double weight_scale = 1.0;

    size_t sample_size = 0;
    position_set_type positions(bitexts.size());
    for (size_t pos = 0; pos < positions.size(); ++ pos) {
      positions[pos] = pos;
      sample_size += bitexts[pos].target.size();
    }
    
    size_t epoch = 0;
    const double lambda = C / sample_size;

    boost::mt19937 gen;
    gen.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> generator(gen);

    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end(), generator);
      
      double objective = 0.0;
      double margin_sum = 0.0;
      size_t num_correct_classify = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sentence_type::const_iterator siter_begin = bitext.source.begin();
	sentence_type::const_iterator siter_end   = bitext.source.end();

	//const double eta = 1.0 / (lambda * (epoch + 2));
	const double eta = 0.2 * std::pow(0.85, double(epoch) / sample_size);
	++ epoch;
	
	double margin = x[id_bias];
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  margin += x[siter->id()];
	margin *= weight_scale;

	const double exp_margin = utils::mathop::exp(margin);
	const double inv_exp_margin = 1.0 / exp_margin;

	const double objective_correct   = boost::math::log1p(inv_exp_margin);
	const double objective_incorrect = boost::math::log1p(exp_margin);
	
	const double gradient_correct   =   1.0 / (1.0 + exp_margin);
	const double gradient_incorrect = - 1.0 / (1.0 + inv_exp_margin);
	
	const size_t num_correct   = std::count(bitext.target.begin(), bitext.target.end(), word);
	const size_t num_incorrect = bitext.target.size() - num_correct;
	
	margin_sum           += margin * num_correct - margin * num_incorrect;
	num_correct_classify += (margin > 1.0) * num_correct + (- margin > 1.0) * num_incorrect;
	objective += objective_correct * num_correct + objective_incorrect * num_incorrect;
	
	const double gradient = gradient_correct * num_correct + gradient_incorrect * num_incorrect;
	
	// rescale...
	rescale(x, weight_scale, weight_norm, 1.0 - eta * lambda);

	const double alpha = (eta / bitext.target.size()) * gradient;
	
	update(x[id_bias], weight_scale, weight_norm, alpha);
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  update(x[siter->id()], weight_scale, weight_norm, alpha);
	
	if (weight_norm > 1.0 / lambda)
	  rescale(x, weight_scale, weight_norm, std::sqrt(1.0 / (lambda * weight_norm)));
      }
      
      if (debug >= 2)
	std::cerr << "word: " << word
		  << " objective: " << objective
		  << " margin: " << (margin_sum / sample_size)
		  << " ratio: " << (double(num_correct_classify) / sample_size)
		  << std::endl;
    }
    
    std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::multiplies<double>(), weight_scale));
  }
};




struct OptimizeLBFGS : public Optimizer
{
  typedef Optimizer base_type;
  
  OptimizeLBFGS(const bitext_set_type& __bitexts,
		const word_type& __word)
    : base_type(__bitexts, __word) {}
  
  void operator()(parameter_set_type& x)
  {
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    
    if (regularize_l1) {
      param.orthantwise_c = C;
      param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    } else
      param.orthantwise_c = 0.0;
    
    // maximum iterations...
    param.max_iterations = max_iteration;
    
    double objective = 0.0;
    
    const size_t vocabulary_size = word_type::allocated();
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);
    
    lbfgs(x.size(),
	  &(*x.begin()),
	  &objective,
	  OptimizeLBFGS::evaluate,
	  0,   
	  this,
	  &param);
  }
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    const OptimizeLBFGS& optimizer = *((const OptimizeLBFGS*) instance);
    
    const bitext_set_type& bitexts = optimizer.bitexts;

    double objective = 0.0;
    std::fill(g, g + n, 0.0);
    
    const word_type::id_type id_bias = vocab_type::EPSILON.id();

    double margin_sum = 0.0;
    size_t num_correct_classify = 0;
    size_t sample_size = 0;
    
    bitext_set_type::const_iterator biter_end = bitexts.end();
    for (bitext_set_type::const_iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {

      sentence_type::const_iterator siter_begin = biter->source.begin();
      sentence_type::const_iterator siter_end   = biter->source.end();
      
      double margin = x[id_bias];
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	margin += x[siter->id()];

      const double exp_margin = utils::mathop::exp(margin);
      const double inv_exp_margin = 1.0 / exp_margin;
      
      const double objective_correct   = boost::math::log1p(inv_exp_margin);
      const double objective_incorrect = boost::math::log1p(exp_margin);
      
      const double gradient_correct   = - 1.0 / (1.0 + exp_margin);
      const double gradient_incorrect =   1.0 / (1.0 + inv_exp_margin);
      
      const size_t num_correct   = std::count(biter->target.begin(), biter->target.end(), optimizer.word);
      const size_t num_incorrect = biter->target.size() - num_correct;
      
      margin_sum           += margin * num_correct - margin * num_incorrect;
      num_correct_classify += (margin > 1.0) * num_correct + (- margin > 1.0) * num_incorrect;
      sample_size          += biter->target.size();
      
      objective += objective_correct * num_correct + objective_incorrect * num_incorrect;
      
      const double gradient = gradient_correct * num_correct + gradient_incorrect * num_incorrect;
      g[id_bias] += gradient;
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	g[siter->id()] += gradient;
    }
    
    if (regularize_l2) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i) {
	g[i] += C * x[i];
	norm += x[i] * x[i];
      }
      objective += 0.5 * C * norm;
    }

    if (debug >= 2)
      std::cerr << "word: " << optimizer.word
		<< " objective: " << objective
		<< " margin: " << (margin_sum / sample_size)
		<< " ratio: " << (double(num_correct_classify) / sample_size)
		<< std::endl;
    
    return objective;
  }
};

struct MapReduce
{
  typedef boost::thread thread_type;
  
  typedef utils::lockfree_list_queue<word_type, std::allocator<word_type> > queue_type;
};

template <typename Opt>
struct Mapper
{
  typedef Opt optimizer_type;
  typedef typename optimizer_type::parameter_set_type parameter_set_type;
  
  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::queue_type queue_type;
  
  const bitext_set_type& bitexts;
  queue_type& queue;
  path_type path_lexicon;
  
  Mapper(const bitext_set_type& __bitexts,
	 queue_type&            __queue,
	 const path_type&       __path_lexicon)
    : bitexts(__bitexts),
      queue(__queue),
      path_lexicon(__path_lexicon) {}
  
  void operator()()
  {
    utils::compress_ostream os(path_lexicon, 1024 * 1024);
    os.precision(10);
    
    word_type word;
    parameter_set_type x;
    
    while (1) {
      queue.pop(word);
      if (word == word_type()) break;
      
      optimizer_type optimizer(bitexts, word);
      optimizer(x);
      
      size_t actives = 0;
      for (word_type::id_type id = 0; id < x.size(); ++ id)
	if (x[id] != 0.0) {
	  os << word << ' ' << word_type(id) << ' ' << x[id] << '\n';
	  
	  if (debug >= 3) 
	    std::cerr << "feature: " << word << ' ' << word_type(id) << ' ' << x[id] << std::endl;

	  ++ actives;
	}
      
      if (debug >= 2)
	std::cerr << "# of active features: " << actives << std::endl;
    }
  }
};
