//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// pyp-transliteration
//

// this is a simple monotone transliteration model...
// we have only one restaurant + two priors for source/target, resampled from the model.
//
// the lambda sampling is taken from:
// @InProceedings{mochihashi-yamada-ueda:2009:ACLIJCNLP,
//   author    = {Mochihashi, Daichi  and  Yamada, Takeshi  and  Ueda, Naonori},
//   title     = {Bayesian Unsupervised Word Segmentation with Nested Pitman-Yor Language Modeling},
//   booktitle = {Proceedings of the Joint Conference of the 47th Annual Meeting of the ACL and the 4th International Joint Conference on Natural Language Processing of the AFNLP},
//   month     = {August},
//   year      = {2009},
//   address   = {Suntec, Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {100--108},
//   url       = {http://www.aclweb.org/anthology/P/P09/P09-1012}
// }
//


#include <map>
#include <iterator>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>

#include "utils/utf8.hpp"
#include "utils/vector2.hpp"
#include "utils/piece.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/chinese_restaurant_process.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/dense_hash_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

struct PYPTranslit
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef std::string word_type;
  
  typedef utils::piece piece_type;
  typedef utils::piece segment_type;
  
  struct segment_pair_type
  {
    segment_type source;
    segment_type target;
    
    segment_pair_type() : source(), target() {}
    segment_pair_type(const segment_type& __source,
		      const segment_type& __target)
      : source(__source), target(__target) {}
    
    friend
    bool operator==(const segment_pair_type& x, const segment_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    size_t hash_value(segment_pair_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
    }
  };
  
  typedef std::vector<segment_pair_type, std::allocator<segment_pair_type> > derivation_type;
  
  typedef utils::chinese_restaurant_process<segment_pair_type, boost::hash<segment_pair_type>, std::equal_to<segment_pair_type>, std::allocator<segment_pair_type > > table_type;
  
  struct length_base_type
  {
    typedef utils::vector2<double, std::allocator<double> > matrix_type;
    
    length_base_type(const double& __p0_source,
		     const double& __p0_target,
		     const double& __lambda_source,
		     const double& __lambda_target,
		     const double& __strength_shape,
		     const double& __strength_rate)
      : p0_source(__p0_source),
	p0_target(__p0_target),
	lambda_source(__lambda_source),
	lambda_target(__lambda_target),
	strength_shape(__strength_shape),
	strength_rate(__strength_rate)
    {
      initialize(__lambda_source, __lambda_target);
    }
    
    double prob(const size_type source_size, const size_type target_size) const
    {
      if (source_size == 0 || target_size == 0) return 0.0;
      
      if (source_size < 32 && target_size < 32)
	return matrix(source_size, target_size);
      else
	return std::exp(utils::mathop::log_poisson(source_size, lambda_source)
			+ utils::mathop::log_poisson(target_size, lambda_target)
			+ std::log(p0_source) * source_size
			+ std::log(p0_target) * target_size);
    }
    
    void initialize(const double& __lambda_source, const double& __lambda_target)
    {
      lambda_source = __lambda_source;
      lambda_target = __lambda_target;
      
      matrix.clear();
      matrix.resize(32, 32);
      
      for (size_type source_size = 1; source_size != 32; ++ source_size)
	for (size_type target_size = 1; target_size != 32; ++ target_size)
	  matrix(source_size, target_size) = std::exp(utils::mathop::log_poisson(source_size, lambda_source)
						      + utils::mathop::log_poisson(target_size, lambda_target)
						      + std::log(p0_source) * source_size
						      + std::log(p0_target) * target_size);
    }
    
    matrix_type matrix;

    double p0_source;
    double p0_target;
    
    double lambda_source;
    double lambda_target;
    double strength_shape;
    double strength_rate;
  };
  
  PYPTranslit(const double __discount,
	      const double __strength,
	      const double __discount_alpha,
	      const double __discount_beta,
	      const double __strength_shape,
	      const double __strength_rate,
	      const double __p0_source,
	      const double __p0_target,
	      const double __lambda_source,
	      const double __lambda_target,
	      const double __lambda_strength_shape,
	      const double __lambda_strength_rate)
    : table(__discount, __strength, __discount_alpha, __discount_beta, __strength_shape, __strength_rate),
      base(__p0_source, __p0_target, __lambda_source, __lambda_target, __lambda_strength_shape, __lambda_strength_rate)
  { }
  
  
  template <typename Iterator, typename Sampler>
  void increment(Iterator first, Iterator last, Sampler& sampler, const double temperature=1.0)
  {
    for (/**/; first != last; ++ first)
      increment(first->source, first->target, sampler, temperature);
  }

  template <typename Iterator, typename Sampler>
  void decrement(Iterator first, Iterator last, Sampler& sampler)
  {
    for (/**/; first != last; ++ first)
      decrement(first->source, first->target, sampler);
  }
  
  template <typename Sampler>
  bool increment(const piece_type& source, const piece_type& target, Sampler& sampler, const double temperature=1.0)
  {
    return table.increment(segment_pair_type(source, target), base.prob(source.size(), target.size()), sampler, temperature);
  }
  
  template <typename Sampler>
  bool decrement(const piece_type& source, const piece_type& target, Sampler& sampler)
  {
    return table.decrement(segment_pair_type(source, target), sampler);
  }
  
  double prob(const piece_type& source, const piece_type& target) const
  {
    return table.prob(segment_pair_type(source, target), base.prob(source.size(), target.size()));
  }
  
  double log_likelihood() const
  {
    return table.log_likelihood();
  }
  
  double log_likelihood(const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    return table.log_likelihood(discount, strength);
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    double source_size = 0.0;
    double target_size = 0.0;
    typename table_type::const_iterator titer_end = table.end();
    for (typename table_type::const_iterator titer = table.begin(); titer != titer_end; ++ titer) {
      source_size += titer->first.source.size() * titer->second.size_table();
      target_size += titer->first.target.size() * titer->second.size_table();
    }
    
    base.initialize(sampler.gamma(base.strength_shape + source_size, base.strength_rate + table.size_table()),
		    sampler.gamma(base.strength_shape + target_size, base.strength_rate + table.size_table()));
    
    table.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    double source_size = 0.0;
    double target_size = 0.0;
    typename table_type::const_iterator titer_end = table.end();
    for (typename table_type::const_iterator titer = table.begin(); titer != titer_end; ++ titer) {
      source_size += titer->first.source.size() * titer->second.size_table();
      target_size += titer->first.target.size() * titer->second.size_table();
    }
    
    base.initialize(sampler.gamma(base.strength_shape + source_size, base.strength_rate + table.size_table()),
		    sampler.gamma(base.strength_shape + target_size, base.strength_rate + table.size_table()));
    
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
public: 
  table_type        table;
  length_base_type  base;
};


struct PYPGraph
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef PYPTranslit::word_type         word_type;
  typedef PYPTranslit::segment_type      segment_type;
  typedef PYPTranslit::segment_pair_type segment_pair_type;
  typedef PYPTranslit::derivation_type   derivation_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  struct edge_type
  {
    segment_pair_type segment;
    logprob_type      prob;
    size_type         source;
    size_type         target;
    
    edge_type() : segment(), prob(), source(0), target(0) {}
    edge_type(const segment_pair_type& __segment,
	      const logprob_type& __prob,
	      const size_type& __source,
	      const size_type& __target) 
      : segment(__segment), prob(__prob), source(__source), target(__target) {}
  };
  
  typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
  
  typedef utils::vector2<edge_set_type, std::allocator<edge_set_type> > node_set_type;
  typedef utils::vector2<logprob_type, std::allocator<logprob_type> > alpha_type;
  typedef std::vector<size_type, std::allocator<size_type> > position_set_type;

  typedef std::vector<logprob_type, std::allocator<logprob_type> > logprob_set_type;
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  void initialize(const word_type& source, const word_type& target)
  {
    position_source.clear();
    position_target.clear();

    position_source.push_back(0);
    position_target.push_back(0);
    
    word_type::const_iterator siter_begin = source.begin();
    word_type::const_iterator siter_end   = source.end();
    for (word_type::const_iterator siter = source.begin(); siter != siter_end; /**/) {
      siter += utils::utf8_size(*siter);
      
      position_source.push_back(siter - siter_begin);
    }
    
    word_type::const_iterator titer_begin = target.begin();
    word_type::const_iterator titer_end   = target.end();
    for (word_type::const_iterator titer = target.begin(); titer != titer_end; /**/) {
      titer += utils::utf8_size(*titer);
      
      position_target.push_back(titer - titer_begin);
    }

    const size_type src_size = position_source.size();
    const size_type trg_size = position_target.size();
    
    nodes.clear();
    alpha.clear();
    
    nodes.resize(src_size, trg_size);
    alpha.resize(src_size, trg_size);
  }
  
  // forward filter
  logprob_type forward(const segment_type& source, const segment_type& target, const PYPTranslit& translit)
  {
    initialize(source, target);
    
    const size_type src_size = position_source.size();
    const size_type trg_size = position_target.size();
    
    alpha(0,0) = cicada::semiring::traits<logprob_type>::one();
    
    for (size_t src_last = 1; src_last != src_size; ++ src_last)
      for (size_t trg_last = 1; trg_last != trg_size; ++ trg_last) {
	
	edge_set_type& edges = nodes(src_last, trg_last);
	logprob_type& sum = alpha(src_last, trg_last);
	
	for (size_t src_first = 0; src_first != src_last; ++ src_first)
	  for (size_t trg_first = 0; trg_first != trg_last; ++ trg_first) 
	    if (alpha(src_first, trg_first) != cicada::semiring::traits<logprob_type>::zero()) {
	      
	      const segment_pair_type segment(segment_type(source.begin() + position_source[src_first], source.begin() + position_source[src_last]),
					      segment_type(target.begin() + position_target[trg_first], target.begin() + position_target[trg_last]));
	      
	      edges.push_back(edge_type(segment, translit.prob(segment.source, segment.target), src_first, trg_first));
	      
	      sum += edges.back().prob * alpha(src_first, trg_first);
	    }
      }
    
    return alpha(src_size - 1, trg_size - 1);
  }
  
  // backward sampling
  template <typename Sampler>
  logprob_type backward(Sampler& sampler, derivation_type& derivation)
  {
    const size_type src_size = position_source.size();
    const size_type trg_size = position_target.size();
    
    derivation.clear();
    logprob_type prob_derivation = cicada::semiring::traits<logprob_type>::one();

    size_type src_last = src_size - 1;
    size_type trg_last = trg_size - 1;
    
    while (src_last || trg_last) {
      const edge_set_type& edges = nodes(src_last, trg_last);
      
      if (edges.empty())
	throw std::runtime_error("empty edges?");
      
      logprob_type logsum;
      logprobs.clear();
      edge_set_type::const_iterator eiter_end = edges.end();
      for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
	logprobs.push_back(alpha(eiter->source, eiter->target) * eiter->prob);
	logsum += logprobs.back();
      }
      
      probs.clear();
      logprob_set_type::const_iterator liter_end = logprobs.end();
      for (logprob_set_type::const_iterator liter = logprobs.begin(); liter != liter_end; ++ liter)
	probs.push_back(*liter / logsum);
      
      prob_set_type::const_iterator piter = sampler.select(probs.begin(), probs.end());
      
      const size_type pos = piter - probs.begin();
      
      derivation.push_back(edges[pos].segment);
      prob_derivation *= edges[pos].prob;
      
      src_last = edges[pos].source;
      trg_last = edges[pos].target;
    }
   
    return prob_derivation;
  }
  
  node_set_type nodes;
  alpha_type    alpha;
  
  position_set_type position_source;
  position_set_type position_target;
  
  logprob_set_type logprobs;
  prob_set_type    probs;
};

typedef boost::filesystem::path path_type;
typedef utils::sampler<boost::mt19937> sampler_type;

typedef PYPTranslit::size_type       size_type;
typedef PYPTranslit::difference_type difference_type;

typedef PYPTranslit::word_type word_type;
typedef PYPTranslit::segment_type      segment_type;
typedef PYPTranslit::segment_pair_type segment_pair_type;
typedef PYPTranslit::derivation_type   derivation_type;

typedef std::vector<word_type, std::allocator<word_type> >             data_set_type;
typedef std::vector<derivation_type, std::allocator<derivation_type> > derivation_set_type;

path_type train_source_file = "-";
path_type train_target_file = "-";

path_type test_source_file;
path_type test_target_file;

int samples = 30;
int baby_steps = 0;
int anneal_steps = 0;
int resample_rate = 1;
int resample_iterations = 2;
bool slice_sampling = false;

double discount = 0.9;
double strength = 1;

double discount_prior_alpha = 1.0;
double discount_prior_beta  = 1.0;
double strength_prior_shape = 1.0;
double strength_prior_rate  = 1.0;

double lambda_source = 2.0;
double lambda_target = 2.0;
double lambda_shape = 0.2;
double lambda_rate  = 0.1;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);

void read_data(const path_type& path, data_set_type& data);
size_t vocabulary_size(const data_set_type& data);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    if (samples < 0)
      throw std::runtime_error("# of samples must be positive");
            
    if (resample_rate <= 0)
      throw std::runtime_error("resample rate must be >= 1");
    
    if (! slice_sampling && strength < 0.0)
      throw std::runtime_error("negative strength w/o slice sampling is not supported!");

    data_set_type sources;
    data_set_type targets;
    
    read_data(train_source_file, sources);
    read_data(train_target_file, targets);
    
    if (sources.empty())
      throw std::runtime_error("empty source data?");
    
    if (targets.empty())
      throw std::runtime_error("empty target data?");
    
    if (sources.size() != targets.size())
      throw std::runtime_error("the data size do not match");

    const size_t vocab_size_source = vocabulary_size(sources);
    const size_t vocab_size_target = vocabulary_size(targets);
    
    if  (debug >= 2)
      std::cerr << "source charset: " << vocab_size_source << " target charset: " << vocab_size_target << std::endl;
    
    derivation_set_type derivations(sources.size());
    
    sampler_type sampler;
    
    PYPTranslit translit(discount,
			 strength,
			 discount_prior_alpha,
			 discount_prior_beta,
			 strength_prior_shape,
			 strength_prior_rate,
			 1.0 / vocab_size_source,
			 1.0 / vocab_size_target,
			 lambda_source,
			 lambda_target,
			 lambda_shape,
			 lambda_rate);
    

    PYPGraph graph;
    
    // sample parameters, first...
    if (slice_sampling)
      translit.slice_sample_parameters(sampler, resample_iterations);
    else
      translit.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2)
      std::cerr << "discount=" << translit.table.discount()
		<< " strength=" << translit.table.strength()
		<< " lambda-source=" << translit.base.lambda_source
		<< " lambda-target=" << translit.base.lambda_target
		<< std::endl;

    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);

    bool sampling = false;
    int sample_iter = 0;
    
    // then, learn!
    for (size_t iter = 0; sample_iter != samples; ++ iter, sample_iter += sampling) {

      double temperature = 1.0;
      bool anneal_finished = true;
      if (anneal_iter != anneal_last) {
	anneal_finished = false;
	temperature = double(anneal_last - anneal_iter) + 1;
	
	++ anneal_iter;

	if (debug >= 2)
	  std::cerr << "temperature: " << temperature << std::endl;
      }
      
      sampling = anneal_finished;

      if (debug) {
	if (sampling)
	  std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
	else
	  std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
      }
      
      for (size_t i = 0; i != sources.size(); ++ i) 
	if (! sources[i].empty() && ! targets[i].empty()) {

	  translit.decrement(derivations[i].begin(), derivations[i].end(), sampler);

	  // forward filter
	  const PYPGraph::logprob_type logsum = graph.forward(sources[i], targets[i], translit);
	
	  // backward sampling
	  const PYPGraph::logprob_type logderivation = graph.backward(sampler, derivations[i]);
 
	  if (debug >= 3) {
	    std::cerr << "sum=" << logsum << " derivation=" << logderivation << std::endl;
	  
	    derivation_type::const_iterator diter_end = derivations[i].end();
	    for (derivation_type::const_iterator diter = derivations[i].begin(); diter != diter_end; ++ diter)
	      std::cerr << "\tsource=" << diter->source << " target=" << diter->target << std::endl;
	  }
	  
	  translit.increment(derivations[i].begin(), derivations[i].end(), sampler);
	}
      
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  translit.slice_sample_parameters(sampler, resample_iterations);
	else
	  translit.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2)
	  std::cerr << "discount=" << translit.table.discount()
		    << " strength=" << translit.table.strength()
		    << " lambda-source=" << translit.base.lambda_source
		    << " lambda-target=" << translit.base.lambda_target
		    << std::endl;
      }
      
      if (debug)
	std::cerr << "log-likelihood: " << translit.log_likelihood() << std::endl;
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void read_data(const path_type& path, data_set_type& data)
{
  data.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  
  std::string line;
  while (std::getline(is, line))
    data.push_back(line);
}

size_t vocabulary_size(const data_set_type& data)
{
  typedef utils::piece piece_type;
  typedef utils::dense_hash_set<piece_type, boost::hash<piece_type>, std::equal_to<piece_type>, std::allocator<piece_type> >::type vocab_type;
  
  vocab_type vocab;
  vocab.set_empty_key(piece_type());
  
  data_set_type::const_iterator diter_end = data.end();
  for (data_set_type::const_iterator diter = data.begin(); diter != diter_end; ++ diter) {
    const word_type& word = *diter;

    if (word.empty()) continue;
    
    word_type::const_iterator witer_end = word.end();
    for (word_type::const_iterator witer = word.begin(); witer != witer_end; /**/) {
      const size_t char_width = utils::utf8_size(*witer);
      
      vocab.insert(piece_type(witer, witer + char_width));
      
      witer += char_width;
    }
  }
  
  return vocab.size();
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("train-source", po::value<path_type>(&train_source_file), "source train file")
    ("train-target", po::value<path_type>(&train_target_file), "target train file")
    
    ("test-source", po::value<path_type>(&test_source_file), "source test file")
    ("test-target", po::value<path_type>(&test_target_file), "target test file")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    
    ("discount",       po::value<double>(&discount)->default_value(discount),                         "discount ~ Beta(alpha,beta)")
    ("discount-alpha", po::value<double>(&discount_prior_alpha)->default_value(discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("discount-beta",  po::value<double>(&discount_prior_beta)->default_value(discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("strength",       po::value<double>(&strength)->default_value(strength),                         "strength ~ Gamma(shape,rate)")
    ("strength-shape", po::value<double>(&strength_prior_shape)->default_value(strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("strength-rate",  po::value<double>(&strength_prior_rate)->default_value(strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
    ("lambda-source", po::value<double>(&lambda_source)->default_value(lambda_source), "lambda for source")
    ("lambda-target", po::value<double>(&lambda_target)->default_value(lambda_target), "lambda for target")
    ("lambda-shape",  po::value<double>(&lambda_shape)->default_value(lambda_shape),   "lambda ~ Gamma(shape,rate)")
    ("lambda-rate",   po::value<double>(&lambda_rate)->default_value(lambda_rate),     "lambda ~ Gamma(shape,rate)")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}

