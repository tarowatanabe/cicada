
#include <algorithm>
#include <vector>
#include <set>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>

#include <utils/compress_stream.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/repository.hpp>
#include <utils/resource.hpp>
#include <utils/mathop.hpp>
#include <utils/sgi_hash_set.hpp>

#include "nicttm/Word.hpp"
#include "nicttm/Vocab.hpp"
#include "nicttm/Sentence.hpp"
#include "nicttm/Bitext.hpp"

#include "lbfgs.h"


typedef boost::filesystem::path path_type;

typedef nicttm::Word            word_type;
typedef nicttm::Vocab           vocab_type;
typedef nicttm::Sentence        sent_type;
typedef nicttm::Bitext          bitext_type;

typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
typedef std::vector<word_type, std::allocator<word_type> > word_set_type;

// I know, it is very bad...
extern int    max_iteration;
extern bool   regularize_l1;
extern bool   regularize_l2;
extern double C;

extern int debug;

void read_bitexts(const path_type& path_source,
		  const path_type& path_target,
		  bitext_set_type& bitexts,
		  word_set_type& vocab)
{
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type > > word_unique_set_type;
#else
  typedef sgi::hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type > > word_unique_set_type;
#endif
  
  typedef std::set<word_type, std::less<word_type>, std::allocator<word_type> > word_sorted_set_type;

  bitexts.clear();
  
  word_unique_set_type uniques;
  word_sorted_set_type uniques_source;
  
  utils::compress_istream is_src(path_source, 1024 * 1024);
  utils::compress_istream is_trg(path_target, 1024 * 1024);
  
  sent_type      source;
  sent_type      target;
  
  while (is_src && is_trg) {
    is_src >> source;
    is_trg >> target;
    
    if (! is_src || ! is_trg) break;
    
    if (source.empty() || target.empty()) continue;
    
    uniques_source.clear();
    uniques_source.insert(source.begin(), source.end());
    
    bitexts.push_back(bitext_type(sent_type(uniques_source.begin(), uniques_source.end()), target));
    uniques.insert(target.begin(), target.end());
  }
  if (is_src || is_trg)
    throw std::runtime_error("# of sentences does not match...");
  
  bitext_set_type(bitexts).swap(bitexts);
  
  vocab.clear();
  vocab.reserve(uniques.size());
  vocab.insert(vocab.end(), uniques.begin(), uniques.end());
}

struct Optimizer
{
  typedef std::vector<double, std::allocator<double> > parameter_set_type;
  
  const bitext_set_type& bitexts;
  const word_type word;
  
  Optimizer(const bitext_set_type& __bitexts,
	    const word_type& __word)
    : bitexts(__bitexts),
      word(__word) {}
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
    
    const word_type::id_type id_bias = vocab_type::NONE.id();
    
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
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end());

      double margin_sum = 0.0;
      size_t num_correct = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sent_type::const_iterator siter_begin = bitext.source.begin();
	sent_type::const_iterator siter_end   = bitext.source.end();
	
	sent_type::const_iterator titer_end = bitext.target.end();
	for (sent_type::const_iterator titer = bitext.target.begin(); titer != titer_end; ++ titer) {
	  const double y = (int(*titer == word) * 2 - 1);
	  
	  double var    = v[id_bias];
	  double margin = x[id_bias];
	  for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
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
	    
	    for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
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
    
    const word_type::id_type id_bias = vocab_type::NONE.id();
    
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
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end());
      
      double margin_sum = 0.0;
      size_t num_correct = 0;

      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sent_type::const_iterator siter_begin = bitext.source.begin();
	sent_type::const_iterator siter_end   = bitext.source.end();
	
	sent_type::const_iterator titer_end = bitext.target.end();
	for (sent_type::const_iterator titer = bitext.target.begin(); titer != titer_end; ++ titer) {
	  const double y = (int(*titer == word) * 2 - 1);
	  
	  double var    = v[id_bias];
	  double margin = x[id_bias];
	  for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
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
	    
	    for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
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
    
    const word_type::id_type id_bias = vocab_type::NONE.id();
    
    x.clear();
    x.reserve(vocabulary_size);
    x.resize(vocabulary_size, 0.0);

    size_t sample_size = 0;
    position_set_type positions(bitexts.size());
    for (size_t pos = 0; pos < positions.size(); ++ pos) {
      positions[pos] = pos;
      sample_size += bitexts[pos].target.size();
    }
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end());
      
      double margin_sum = 0.0;
      size_t num_correct = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sent_type::const_iterator siter_begin = bitext.source.begin();
	sent_type::const_iterator siter_end   = bitext.source.end();
	
	sent_type::const_iterator titer_end = bitext.target.end();
	for (sent_type::const_iterator titer = bitext.target.begin(); titer != titer_end; ++ titer) {
	  const double y = (int(*titer == word) * 2 - 1);
	  
	  double margin = x[id_bias];
	  for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	    margin += x[siter->id()];
	  
	  const double score = margin * y;

	  num_correct += bool(score > 1.0);
	  margin_sum += score;
	  
	  if (score <= 1.0) {
	    const double norm = (siter_end - siter_begin) + 1;
	    const double alpha = std::min(1.0 / C, y * (1.0 - score) / norm);
	    
	    x[id_bias] += alpha;
	    for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
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
    
    const word_type::id_type id_bias = vocab_type::NONE.id();
    
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
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end());
      
      double objective = 0.0;
      double margin_sum = 0.0;
      size_t num_correct_classify = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sent_type::const_iterator siter_begin = bitext.source.begin();
	sent_type::const_iterator siter_end   = bitext.source.end();
	
	const double eta = 1.0 / (1.0 + epoch / sample_size);
	++ epoch;
	penalty += eta * lambda * bitext.target.size();
	
	double margin = x[id_bias];
	for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  margin += x[siter->id()];

	const double objective_correct   = boost::math::log1p(utils::mathop::exp(- 1.0 * margin));
	const double objective_incorrect = boost::math::log1p(utils::mathop::exp(  1.0 * margin));
	
	const double gradient_correct   =   1.0 / (1.0 + utils::mathop::exp(  1.0 * margin));
	const double gradient_incorrect = - 1.0 / (1.0 + utils::mathop::exp(- 1.0 * margin));
	
	const size_t num_correct   = std::count(bitext.target.begin(), bitext.target.end(), word);
	const size_t num_incorrect = bitext.target.size() - num_correct;
	
	margin_sum           += margin * num_correct - margin * num_incorrect;
	num_correct_classify += (margin > 1.0) * num_correct + (- margin > 1.0) * num_incorrect;
	objective += objective_correct * num_correct + objective_incorrect * num_incorrect;
	
	const double gradient = gradient_correct * num_correct + gradient_incorrect * num_incorrect;
	
	x[id_bias] += eta * gradient;
	apply(x[id_bias], penalties[id_bias], penalty);
	for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
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
    
    const word_type::id_type id_bias = vocab_type::NONE.id();
    
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

    for (int iter = 0; iter < max_iteration; ++ iter) {
      std::random_shuffle(positions.begin(), positions.end());
      
      double objective = 0.0;
      double margin_sum = 0.0;
      size_t num_correct_classify = 0;
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const bitext_type& bitext = bitexts[*piter];
	
	sent_type::const_iterator siter_begin = bitext.source.begin();
	sent_type::const_iterator siter_end   = bitext.source.end();

	const double eta = 1.0 / (lambda * (epoch + 2));
	++ epoch;
	
	double margin = x[id_bias];
	for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  margin += x[siter->id()];
	margin *= weight_scale;

	const double objective_correct   = boost::math::log1p(utils::mathop::exp(- 1.0 * margin));
	const double objective_incorrect = boost::math::log1p(utils::mathop::exp(  1.0 * margin));
	
	const double gradient_correct   =   1.0 / (1.0 + utils::mathop::exp(  1.0 * margin));
	const double gradient_incorrect = - 1.0 / (1.0 + utils::mathop::exp(- 1.0 * margin));
	
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
	for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
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
    
    const word_type::id_type id_bias = vocab_type::NONE.id();

    double margin_sum = 0.0;
    size_t num_correct_classify = 0;
    size_t sample_size = 0;
    
    bitext_set_type::const_iterator biter_end = bitexts.end();
    for (bitext_set_type::const_iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {

      sent_type::const_iterator siter_begin = biter->source.begin();
      sent_type::const_iterator siter_end   = biter->source.end();
      
      double margin = x[id_bias];
      for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	margin += x[siter->id()];
      
      const double objective_correct   = boost::math::log1p(utils::mathop::exp(- 1.0 * margin));
      const double objective_incorrect = boost::math::log1p(utils::mathop::exp(  1.0 * margin));
      
      const double gradient_correct   = - 1.0 / (1.0 + utils::mathop::exp(  1.0 * margin));
      const double gradient_incorrect =   1.0 / (1.0 + utils::mathop::exp(- 1.0 * margin));
      
      const size_t num_correct   = std::count(biter->target.begin(), biter->target.end(), optimizer.word);
      const size_t num_incorrect = biter->target.size() - num_correct;
      
      margin_sum           += margin * num_correct - margin * num_incorrect;
      num_correct_classify += (margin > 1.0) * num_correct + (- margin > 1.0) * num_incorrect;
      sample_size          += biter->target.size();
      
      objective += objective_correct * num_correct + objective_incorrect * num_incorrect;
      
      const double gradient = gradient_correct * num_correct + gradient_incorrect * num_incorrect;
      g[id_bias] += gradient;
      for (sent_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
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
