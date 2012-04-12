//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-SynAlign based on GHKM work:

// @InProceedings{cohn-blunsom:2009:EMNLP,
//   author    = {Cohn, Trevor  and  Blunsom, Phil},
//   title     = {A {Bayesian} Model of Syntax-Directed Tree to String Grammar Induction},
//   booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//   month     = {August},
//   year      = {2009},
//   address   = {Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {352--361},
//   url       = {http://www.aclweb.org/anthology/D/D09/D09-1037}
// }
//

//
// We do not perform composition, and each node is mapped to corresponding span in the target-side.
//
// How do we handle empty...!

#include <map>
#include <deque>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/tree_rule.hpp>
#include <cicada/tree_rule_compact.hpp>
#include <cicada/rule.hpp>

#include "utils/pyp_parameter.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/chart.hpp"
#include "utils/chunk_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/restaurant.hpp"
#include "utils/restaurant_vector.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/dense_hash_map.hpp"
#include "utils/dense_hash_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/indexed_map.hpp"
#include "utils/indexed_set.hpp"
#include "utils/lockfree_list_queue.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

typedef cicada::Vocab      vocab_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Symbol     symbol_type;
typedef cicada::Symbol     word_type;
typedef cicada::HyperGraph hypergraph_type;

typedef cicada::Rule               rule_type;
typedef rule_type::rule_ptr_type   rule_ptr_type;
typedef rule_type::symbol_set_type symbol_set_type;

//
// syntactic alignment
//
// We assume one side is syntactic, and the other is hiero-like string
//
// Thus, it is possible to induce a synchronous-CFG by mapping non-temrinal symbols to the string-side
// Or, keep them separate like a synchronous-TSG.
//
// our approach: use of hierarchical structure: IBM Model 1 like probability computation when back-off
// 
// TODO:
// How to handle distortion represented by the ordering in [x,3] etc.
//

struct PYP
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
};

struct PYPLexicon
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;

  struct word_pair_type
  {
    word_type source;
    word_type target;
    
    word_pair_type() : source(), target() {}
    word_pair_type(const word_type& __source, const word_type& __target)
      : source(__source), target(__target) {}
    
    friend
    bool operator==(const word_pair_type& x, const word_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    size_t hash_value(word_pair_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.source.id(), x.target.id());
    }
  };
  
  typedef utils::restaurant<word_pair_type, boost::hash<word_pair_type>, std::equal_to<word_pair_type>, std::allocator<word_pair_type > > table_type;

  PYPLexicon(const double& __p0,
	     const parameter_type& parameter)
    : p0(__p0),
      counts0(0),
      table(parameter) {}
  
  template <typename Sampler>
  void increment(const rule_type* source, const symbol_set_type& target, Sampler& sampler, const double temperature=1.0)
  {
    // TODO: we assume penntreebank-style tree...
    
    const word_type src = (source->rhs.size() == 1 && source->rhs.front().is_terminal()
			   ? source->rhs.front()
			   : vocab_type::EPSILON);
    
    size_type num_terminals = 0;
    symbol_set_type::const_iterator titer_end = target.end();
    for (symbol_set_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
      if (titer->is_terminal()) {
	if (table.increment(word_pair_type(src, *titer), p0, sampler, temperature))
	  ++ counts0;
	
	++ num_terminals;
      }
    
    if (! num_terminals)
      if (table.increment(word_pair_type(src, vocab_type::EPSILON), p0, sampler, temperature))
	++ counts0;
  }

  template <typename Sampler>
  void decrement(const rule_type* source, const symbol_set_type& target, Sampler& sampler)
  {
    const word_type src = (source->rhs.size() == 1 && source->rhs.front().is_terminal()
			   ? source->rhs.front()
			   : vocab_type::EPSILON);
    
    size_type num_terminals = 0;
    symbol_set_type::const_iterator titer_end = target.end();
    for (symbol_set_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
      if (titer->is_terminal()) {
	if (table.decrement(word_pair_type(src, *titer), sampler))
	  -- counts0;
	
	++ num_terminals;
      }
    
    if (! num_terminals)
      if (table.decrement(word_pair_type(src, vocab_type::EPSILON), sampler))
	-- counts0;
  }
  
  double prob(const rule_type* source, const symbol_set_type& target) const
  {
    const word_type src = (source->rhs.size() == 1 && source->rhs.front().is_terminal()
			   ? source->rhs.front()
			   : vocab_type::EPSILON);

    double p = 1.0;
    
    size_type num_terminals = 0;
    symbol_set_type::const_iterator titer_end = target.end();
    for (symbol_set_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
      if (titer->is_terminal()) {
	p *= table.prob(word_pair_type(src, *titer), p0);
	++ num_terminals;
      }
    
    if (! num_terminals)
      p *= table.prob(word_pair_type(src, vocab_type::EPSILON), p0);
    
    return p;
  }
  
  double log_likelihood() const
  {
    return table.log_likelihood() + std::log(p0) * counts0;
  }
  
  double log_likelihood(const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    return table.log_likelihood(discount, strength);
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    table.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
  }

  
  double p0;
  size_type counts0;
  table_type table;
};

struct PYPFertility
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  
  typedef utils::restaurant_vector< > table_type;
  typedef utils::unordered_map<symbol_type, table_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			       std::allocator<std::pair<const symbol_type, table_type> > >::type table_set_type;
  
  PYPFertility(const parameter_type& __parameter)
    : p0(1.0 / 10), // heuristic...
      counts0(0),
      fallback(__parameter),
      tables(),
      parameter(__parameter) {}

  template <typename Sampler>
  void increment(const rule_type* source, const symbol_set_type& target, Sampler& sampler, const double temperature=1.0)
  {
    const size_type num_terminals = target.size() - target.arity();
    
    table_set_type::iterator titer = tables.find(source->lhs);
    if (titer == tables.end())
      titer = tables.insert(std::make_pair(source->lhs, table_type(parameter.discount, parameter.strength))).first;
    
    if (titer->second.increment(num_terminals, fallback.prob(num_terminals, p0), sampler, temperature))
      if (fallback.increment(num_terminals, p0, sampler, temperature))
	++ counts0;
  }

  template <typename Sampler>
  void decrement(const rule_type* source, const symbol_set_type& target, Sampler& sampler)
  {
    const size_type num_terminals = target.size() - target.arity();
    
    if (tables[source->lhs].decrement(num_terminals, sampler))
      if (fallback.decrement(num_terminals, sampler))
	-- counts0;
  }
  
  double prob(const rule_type* source, const symbol_set_type& target) const
  {
    const size_type num_terminals = target.size() - target.arity();
    
    const double p = fallback.prob(num_terminals, p0);
    
    table_set_type::const_iterator titer = tables.find(source->lhs);
    
    return (titer != tables.end() ? titer->second.prob(num_terminals, p) : p);
  }
  
  double log_likelihood() const
  {
    double logprob = std::log(p0) * counts0;
    
    logprob += fallback.log_likelihood();
    
    logprob += parameter.log_likelihood();
    table_set_type::const_iterator titer_end = tables.end();
    for (table_set_type::const_iterator titer = tables.begin(); titer != titer_end; ++ titer)
      logprob += titer->second.log_likelihood();
    
    return logprob;
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    fallback.sample_parameters(sampler, num_loop, num_iterations);
    
    for (int iter = 0; iter != num_loop; ++ iter) {
      parameter.discount = sample_discount(tables.begin(), tables.end(), sampler, parameter);
      parameter.strength = sample_strength(tables.begin(), tables.end(), sampler, parameter);
    }
    
    table_set_type::iterator titer_end = tables.end();
    for (table_set_type::iterator titer = tables.begin(); titer != titer_end; ++ titer) {
      titer->second.discount() = parameter.discount;
      titer->second.strength() = parameter.strength;
    }
  }
  
  // we do not support slice-sample-parameters for simplicity...
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    sample_parameters(sampler, num_loop, num_iterations);
  }
  

  template <typename Iterator, typename Sampler>
  double sample_strength(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double x = 0.0;
    double y = 0.0;
    
    for (/**/; first != last; ++ first) {
      x += first->second.sample_log_x(sampler, param.discount, param.strength);
      y += first->second.sample_y(sampler, param.discount, param.strength);
    }
    
    return sampler.gamma(param.strength_shape + y, param.strength_rate - x);
  }
  
  template <typename Iterator, typename Sampler>
  double sample_discount(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double y = 0.0;
    double z = 0.0;
    
    for (/**/; first != last; ++ first) {
      y += first->second.sample_y_inv(sampler, param.discount, param.strength);
      z += first->second.sample_z_inv(sampler, param.discount, param.strength);
    }
    
    return sampler.beta(param.discount_alpha + y, param.discount_beta + z);
  }
  
  double         p0;
  size_type      counts0;
  table_type     fallback;
  table_set_type tables;
  parameter_type parameter;
};

struct PYPDistortion
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;

  struct distortion_type
  {
    symbol_type non_terminal;
    int         distortion;
    
    distortion_type() : non_terminal(), distortion() {}
    distortion_type(const symbol_type& __non_terminal, const int& __distortion)
      : non_terminal(__non_terminal), distortion(__distortion) {}
    
    friend
    bool operator==(const distortion_type& x, const distortion_type& y)
    {
      return x.non_terminal == y.non_terminal && x.distortion == y.distortion;
    }
    
    friend
    size_t hash_value(distortion_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.non_terminal.id(), x.distortion);
    }
  };
  
  typedef utils::restaurant<distortion_type, boost::hash<distortion_type>, std::equal_to<distortion_type>, std::allocator<distortion_type > > table_type;
  typedef utils::unordered_map<symbol_type, table_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			       std::allocator<std::pair<const symbol_type, table_type> > >::type table_set_type;
  
  PYPDistortion(const parameter_type& __parameter)
    : p0(1.0 / 3), // heuristic...
      counts0(0),
      fallback(__parameter),
      tables(),
      parameter(__parameter) {}

  template <typename Sampler>
  void increment(const rule_type* source, const symbol_set_type& target, Sampler& sampler, const double temperature=1.0)
  {
    if (source->rhs.size() == 1 && source->rhs.front().is_terminal()) return;
    
    symbol_set_type::const_iterator titer_begin = target.begin();
    symbol_set_type::const_iterator titer_end   = target.end();
    symbol_set_type::const_iterator titer_pos = titer_begin;
    
    size_type i = 0;
    
    symbol_set_type::const_iterator siter_end = source->rhs.end();
    for (symbol_set_type::const_iterator siter = source->rhs.begin(); siter != siter_end; ++ siter) 
      if (siter->is_non_terminal()) {
	symbol_set_type::const_iterator titer = std::find(titer_begin, titer_end, vocab_type::X.non_terminal(i + 1));
	
	if (titer == titer_end)
	  throw std::runtime_error("invalid distortion!");
	
	const distortion_type distortion(siter->non_terminal(), titer - titer_pos);
	
	table_set_type::iterator iter = tables.find(source->lhs);
	if (iter == tables.end())
	  iter = tables.insert(std::make_pair(source->lhs, table_type(parameter.discount, parameter.strength))).first;
	
	const double p = fallback.prob(distortion, p0);
	
	if (iter->second.increment(distortion, p, sampler, temperature))
	  if (fallback.increment(distortion, p0, sampler, temperature))
	    ++ counts0;
	
	titer_pos = titer + 1;
	++ i;
      }
  }
  
  template <typename Sampler>
  void decrement(const rule_type* source, const symbol_set_type& target, Sampler& sampler)
  {
    if (source->rhs.size() == 1 && source->rhs.front().is_terminal()) return;
    
    symbol_set_type::const_iterator titer_begin = target.begin();
    symbol_set_type::const_iterator titer_end   = target.end();
    symbol_set_type::const_iterator titer_pos = titer_begin;
    
    size_type i = 0;
    
    symbol_set_type::const_iterator siter_end = source->rhs.end();
    for (symbol_set_type::const_iterator siter = source->rhs.begin(); siter != siter_end; ++ siter) 
      if (siter->is_non_terminal()) {
	symbol_set_type::const_iterator titer = std::find(titer_begin, titer_end, vocab_type::X.non_terminal(i + 1));
	
	if (titer == titer_end)
	  throw std::runtime_error("invalid distortion!");
	
	const distortion_type distortion(siter->non_terminal(), titer - titer_pos);
	
	if (tables[source->lhs].decrement(distortion, sampler))
	  if (fallback.decrement(distortion, sampler))
	    -- counts0;
	
	titer_pos = titer + 1;
	++ i;
      }    
  }
  
  double prob(const rule_type* source, const symbol_set_type& target) const
  {
    if (source->rhs.size() == 1 && source->rhs.front().is_terminal()) return 1.0;
    
    double p = 1.0;
    
    symbol_set_type::const_iterator titer_begin = target.begin();
    symbol_set_type::const_iterator titer_end   = target.end();
    symbol_set_type::const_iterator titer_pos = titer_begin;
    
    size_type i = 0;
    symbol_set_type::const_iterator siter_end = source->rhs.end();
    for (symbol_set_type::const_iterator siter = source->rhs.begin(); siter != siter_end; ++ siter) 
      if (siter->is_non_terminal()) {
	symbol_set_type::const_iterator titer = std::find(titer_begin, titer_end, vocab_type::X.non_terminal(i + 1));
	
	if (titer == titer_end)
	  throw std::runtime_error("invalid distortion!");
	
	const distortion_type distortion(siter->non_terminal(), titer - titer_pos);
	
	table_set_type::const_iterator iter = tables.find(source->lhs);
	if (iter == tables.end())
	  p *= fallback.prob(distortion, p0);
	else
	  p *= iter->second.prob(distortion, fallback.prob(distortion, p0));
	
	titer_pos = titer + 1;
	++ i;
      }    
    
    
    return p;
  }
  
  double log_likelihood() const
  {
    double logprob = std::log(p0) * counts0;
    
    logprob += fallback.log_likelihood();
    
    logprob += parameter.log_likelihood();
    table_set_type::const_iterator titer_end = tables.end();
    for (table_set_type::const_iterator titer = tables.begin(); titer != titer_end; ++ titer)
      logprob += titer->second.log_likelihood();
    
    return logprob;
  }

  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    fallback.sample_parameters(sampler, num_loop, num_iterations);
    
    for (int iter = 0; iter != num_loop; ++ iter) {
      parameter.discount = sample_discount(tables.begin(), tables.end(), sampler, parameter);
      parameter.strength = sample_strength(tables.begin(), tables.end(), sampler, parameter);
    }
    
    table_set_type::iterator titer_end = tables.end();
    for (table_set_type::iterator titer = tables.begin(); titer != titer_end; ++ titer) {
      titer->second.discount() = parameter.discount;
      titer->second.strength() = parameter.strength;
    }
  }
  
  // we do not support slice-sample-parameters for simplicity...
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    sample_parameters(sampler, num_loop, num_iterations);
  }
  

  template <typename Iterator, typename Sampler>
  double sample_strength(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double x = 0.0;
    double y = 0.0;
    
    for (/**/; first != last; ++ first) {
      x += first->second.sample_log_x(sampler, param.discount, param.strength);
      y += first->second.sample_y(sampler, param.discount, param.strength);
    }
    
    return sampler.gamma(param.strength_shape + y, param.strength_rate - x);
  }
  
  template <typename Iterator, typename Sampler>
  double sample_discount(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double y = 0.0;
    double z = 0.0;
    
    for (/**/; first != last; ++ first) {
      y += first->second.sample_y_inv(sampler, param.discount, param.strength);
      z += first->second.sample_z_inv(sampler, param.discount, param.strength);
    }
    
    return sampler.beta(param.discount_alpha + y, param.discount_beta + z);
  }

  double         p0;
  size_type      counts0;
  table_type     fallback;
  table_set_type tables;
  parameter_type parameter;
};

struct PYPSynAlign
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  
  struct rule_pair_type
  {
    const rule_type* source;
    symbol_set_type  target;
    
    rule_pair_type() : source(0), target() {}
    rule_pair_type(const rule_type* __source, const symbol_set_type& __target) : source(__source), target(__target) {}
    rule_pair_type(const rule_ptr_type& __source, const symbol_set_type& __target) : source(__source.get()), target(__target) {}
    
    friend
    bool operator==(const rule_pair_type& x, const rule_pair_type& y)
    {
      return ((x.source == y.source) || (x.source && y.source && *(x.source) == *(y.source))) && x.target == y.target;
    }
    
    friend
    size_t hash_value(rule_pair_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.target.begin(), x.target.end(), x.source ? hash_value(*(x.source)) : size_t(0));
    }
  };
  
  typedef utils::restaurant<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>, std::allocator<rule_pair_type > > table_type;
  
  PYPSynAlign(const PYPLexicon& __lexicon,
	      const PYPFertility& __fertility,
	      const PYPDistortion& __distortion,
	      const parameter_type& parameter)
    : lexicon(__lexicon),
      fertility(__fertility),
      distortion(__distortion),
      table(parameter) {}
  
  template <typename Sampler>
  void increment(const rule_type* source, const symbol_set_type& target, Sampler& sampler, const double temperature=1.0)
  {
    increment(rule_pair_type(source, target), sampler, temperature);
  }
  
  template <typename Sampler>
  void increment(const rule_pair_type& rule_pair, Sampler& sampler, const double temperature=1.0)
  {
#if 0
    const double p0 = lexicon.prob(rule_pair.source, rule_pair.target);
    
    
    if (table.increment(rule_pair, p0, sampler, temperature))
      lexicon.increment(rule_pair.source, rule_pair.target, sampler, temperature);
#endif
    
    const double p0 = (lexicon.prob(rule_pair.source, rule_pair.target)
		       * fertility.prob(rule_pair.source, rule_pair.target)
		       * distortion.prob(rule_pair.source, rule_pair.target));
    
    
    if (table.increment(rule_pair, p0, sampler, temperature)) {
      lexicon.increment(rule_pair.source, rule_pair.target, sampler, temperature);
      fertility.increment(rule_pair.source, rule_pair.target, sampler, temperature);
      distortion.increment(rule_pair.source, rule_pair.target, sampler, temperature);
    }
  }
  
  template <typename Sampler>
  void decrement(const rule_type* source, const symbol_set_type& target, Sampler& sampler)
  {
    decrement(rule_pair_type(source, target), sampler);
  }
  
  template <typename Sampler>
  void decrement(const rule_pair_type& rule_pair, Sampler& sampler)
  {
#if 0
    if (table.decrement(rule_pair, sampler))
      lexicon.decrement(rule_pair.source, rule_pair.target, sampler);
#endif
    if (table.decrement(rule_pair, sampler)) {
      lexicon.decrement(rule_pair.source, rule_pair.target, sampler);
      fertility.decrement(rule_pair.source, rule_pair.target, sampler);
      distortion.decrement(rule_pair.source, rule_pair.target, sampler);
    }
  }
  

  double prob(const rule_type* source, const symbol_set_type& target) const
  {
    const double p0 = lexicon.prob(source, target) * fertility.prob(source, target) * distortion.prob(source, target);
    //const double p0 = lexicon.prob(source, target);

    return table.prob(rule_pair_type(source, target), p0);
  }

  double prob(const rule_ptr_type& source, const symbol_set_type& target) const
  {
    return prob(source.get(), target);
  }
  
  double log_likelihood() const
  {
    return lexicon.log_likelihood() + fertility.log_likelihood() + distortion.log_likelihood() + table.log_likelihood();
    //return lexicon.log_likelihood() + table.log_likelihood();
  }
  
  double log_likelihood(const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    return table.log_likelihood(discount, strength);
  }
  
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    lexicon.sample_parameters(sampler, num_loop, num_iterations);

    fertility.sample_parameters(sampler, num_loop, num_iterations);

    distortion.sample_parameters(sampler, num_loop, num_iterations);
    
    table.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    lexicon.slice_sample_parameters(sampler, num_loop, num_iterations);

    fertility.slice_sample_parameters(sampler, num_loop, num_iterations);

    distortion.slice_sample_parameters(sampler, num_loop, num_iterations);
    
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  PYPLexicon    lexicon;
  PYPFertility  fertility;
  PYPDistortion distortion;
  table_type    table;
};

struct PYPGraph
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // for each node in the source-side, we have matching target-span...
  typedef PYPSynAlign::rule_pair_type rule_pair_type;
  typedef std::vector<rule_pair_type, std::allocator<rule_pair_type> > derivation_type;
  
  typedef PYPSynAlign::logprob_type logprob_type;
  typedef PYPSynAlign::prob_type    prob_type;
  
  typedef utils::chart<logprob_type, std::allocator<logprob_type> > beta_type;
  typedef std::vector<beta_type, std::allocator<beta_type> > beta_chart_type;

  typedef utils::chart<double, std::allocator<double> > logprob_target_set_type;
  typedef std::vector<double, std::allocator<double> > logprob_length_set_type;

  struct span_type
  {
    size_type first;
    size_type last;
    
    span_type() : first(0), last(0) {}
    span_type(const size_type& __first, const size_type& __last) : first(__first), last(__last) {}
    
    bool empty() const { return first == last; }
    size_type size() const { return last - first; }

    friend
    bool operator<(const span_type& x, const span_type& y)
    {
      return (x.first < y.first || (!(y.first < x.first) && x.last < y.last));
    }

    friend
    bool operator==(const span_type& x, const span_type& y)
    {
      return x.first == y.first && x.last == y.last;
    }
    
    friend
    size_t hash_value(span_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.last, x.first);
    }
    
    friend
    bool disjoint(const span_type& x, const span_type& y)
    {
      return x.empty() || y.empty() || x.last <= y.first || y.last <= x.first;
    }
  };

  typedef std::vector<span_type, std::allocator<span_type> > span_set_type;
    
  struct span_set_hash : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    size_t operator()(const span_set_type& x) const
    {
      return hasher_type::operator()(x.begin(), x.end(), 0);
    }
  };
  
  typedef utils::unordered_map<span_set_type, symbol_set_type, span_set_hash, std::equal_to<span_set_type>,
			       std::allocator<std::pair<const span_set_type, symbol_set_type> > >::type rule_string_cache_type;
  typedef utils::chart<rule_string_cache_type, std::allocator<rule_string_cache_type> > rule_string_cache_chart_type;
  
  
  struct edge_type
  {
    hypergraph_type::id_type edge;
    symbol_set_type          target;
    span_set_type            antecedent;
    logprob_type             prob;
    
    edge_type() : edge(), target(), antecedent(), prob() {}
    edge_type(const hypergraph_type::id_type& __edge,
	      const symbol_set_type& __target,
	      const span_set_type& __antecedent,
	      const logprob_type& __prob)
      : edge(__edge), target(__target), antecedent(__antecedent), prob(__prob) {}
  };
  
  typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
  
  struct span_edge_map_type
  {
    typedef utils::indexed_set<span_type, boost::hash<span_type>, std::equal_to<span_type>, std::allocator<span_type> > index_type;
    typedef utils::chunk_vector<edge_set_type, 4096/sizeof(edge_set_type), std::allocator<edge_set_type> > mapped_type;
    
    struct value_type
    {
      value_type(const span_type& __first, const edge_set_type& __second) : first(__first), second(__second) {}
      const span_type&     first;
      const edge_set_type& second;
    };
    
    typedef index_type::const_iterator const_iterator;
    
    const_iterator find(const span_type& span) const
    {
      return index.find(span);
    }

    const_iterator begin() const { return index.begin(); }
    const_iterator end() const { return index.end(); }
    
    bool empty() const { return mapped.empty(); }
    size_type size() const { return mapped.size(); }
    
    edge_set_type& operator[](const span_type& span)
    {
      index_type::iterator iter = index.insert(span).first;
      
      if (index.size() > mapped.size())
	mapped.resize(index.size());
      
      return mapped[iter - index.begin()];
    }
    
    value_type operator[](const size_type& pos) const
    {
      return value_type(index[pos], mapped[pos]);
    }

    span_edge_map_type() : index(), mapped() {}
    
    
    index_type  index;
    mapped_type mapped;
  };
  
#if 0
  typedef utils::indexed_map<span_type, edge_set_type, boost::hash<span_type>, std::equal_to<span_type>,
			     std::allocator<std::pair<span_type, edge_set_type> > > span_edge_map_type;
#endif
  typedef std::vector<span_edge_map_type, std::allocator<span_edge_map_type> > span_edge_chart_type;

  typedef std::vector<logprob_type, std::allocator<logprob_type> > logprob_set_type;
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  void initialize(const hypergraph_type& source, const sentence_type& target, const PYPSynAlign& synalign)
  {
    edges.clear();
    edges.reserve(source.nodes.size());
    edges.resize(source.nodes.size(), span_edge_map_type());

    beta.clear();
    beta.reserve(source.nodes.size());
    beta.resize(source.nodes.size(), beta_type(target.size() + 1));

    rule_strings.clear();
    rule_strings.resize(target.size() + 1);
  }
  
  const symbol_set_type& rule_string(const sentence_type& sentence,
				     const span_type& span,
				     const span_set_type& antecedent)
  {
    typedef std::pair<span_type, size_type> span_index_type;
    typedef std::vector<span_index_type, std::allocator<span_index_type> > span_index_set_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
    
    if (span.first >= sentence.size() || span.last > sentence.size() || span.first > span.last)
      throw std::runtime_error("invalid span");
    
    rule_string_cache_type& cache = rule_strings(span.first, span.last);
    
    rule_string_cache_type::iterator iter = cache.find(antecedent);
    if (iter == cache.end()) {
      if (antecedent.empty())
	iter = cache.insert(std::make_pair(antecedent, symbol_set_type(sentence.begin() + span.first, sentence.begin() + span.last))).first;
      else {
	buffer_type buffer;
	span_index_set_type span_index;
	
	span_index.reserve(antecedent.size());
	span_set_type::const_iterator aiter_begin = antecedent.begin();
	span_set_type::const_iterator aiter_end   = antecedent.end();
	for (span_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  span_index.push_back(std::make_pair(*aiter, aiter - aiter_begin));

	std::sort(span_index.begin(), span_index.end());
	
	size_type pos = span.first;
	span_index_set_type::const_iterator siter_end = span_index.end();
	for (span_index_set_type::const_iterator siter = span_index.begin(); siter != siter_end; ++ siter) {
	  if (! siter->first.empty()) {
	    if (siter->first.first < pos)
	      throw std::runtime_error("ivnalid antecedent span");

	    buffer.insert(buffer.end(), sentence.begin() + pos, sentence.begin() + siter->first.first);
	    buffer.push_back(vocab_type::X.non_terminal(siter->second + 1));
	    
	    pos = siter->first.last;
	  } else
	    buffer.push_back(vocab_type::X.non_terminal(siter->second + 1));
	}
	
	if (span.last < pos)
	  throw std::runtime_error("invalid antecedent span");
	
	buffer.insert(buffer.end(), sentence.begin() + pos, sentence.begin() + span.last);
	
	iter = cache.insert(std::make_pair(antecedent, symbol_set_type(buffer.begin(), buffer.end()))).first;
      }
    }
    
    return iter->second;
  }
  
  struct Inside
  {
    Inside(PYPGraph& __graph,
	   const hypergraph_type& __source,
	   const sentence_type&   __target,
	   const PYPSynAlign&     __synalign,
	   const size_type        __terminal_max)
      : graph(__graph), source(__source), target(__target), synalign(__synalign), terminal_max(__terminal_max) {}
    
    template <typename Index>
    void operator()(const hypergraph_type::edge_type& edge, Index& j, const Index& j_ends, span_set_type& antecedent, size_type last)
    {
      if (last == j.size()) {
	const size_type pos = edge.head;
	
	span_type span(target.size(), 0);
	size_type total_hole = 0;
	logprob_type prob_antecedent = cicada::semiring::traits<logprob_type>::one();
	
	for (size_t i = 0; i != edge.tails.size(); ++ i) {
	  antecedent[i] = graph.edges[edge.tails[i]][j[i]].first;
	  
	  const span_type& span_antecedent = antecedent[i];
	  
	  if (! span_antecedent.empty()) {
	    span.first = utils::bithack::min(span.first, span_antecedent.first);
	    span.last  = utils::bithack::max(span.last,  span_antecedent.last);
	    
	    total_hole += span_antecedent.size();
	  }
	  
	  prob_antecedent *= graph.beta[edge.tails[i]](span_antecedent.first, span_antecedent.last);
	}
	
	if (span.first == target.size() && span.last == 0) {
	  // if phrases, we will allow all the lengths..
	  // otherwise, terminal_max!
	      
	  const size_type span_length_min = utils::bithack::branch(pos == source.goal, target.size(), size_type(0));
	  //const size_type span_length_max = utils::bithack::branch(antecedent.empty(), target.size(), utils::bithack::min(target.size(), terminal_max));
	  const size_type span_length_max = utils::bithack::min(target.size(), terminal_max);
	      
	  for (size_type span_length = span_length_min; span_length <= span_length_max; ++ span_length) {
	    const size_type span_first_max = (span_length != 0) * (target.size() - span_length);
	    for (size_type span_first = 0; span_first <= span_first_max; ++ span_first) {
	      const size_type span_last = span_first + span_length;

	      // ITG-like!
	      //if (! antecedent.empty() && span_length) continue;
	      
	      const symbol_set_type& target_string = graph.rule_string(target, span_type(span_first, span_last), antecedent);
		  
	      const logprob_type prob = synalign.prob(edge.rule, target_string);
		  
	      graph.beta[pos](span_first, span_last) = std::max(graph.beta[pos](span_first, span_last), prob * prob_antecedent);
		  
	      graph.edges[pos][span_type(span_first, span_last)].push_back(edge_type(edge.id, target_string, antecedent, prob));
	    }
	  }
	} else {
	  const size_type terminal_hole_size = span.size() - total_hole;
	  const difference_type terminal_span_size = terminal_max - terminal_hole_size;
	  // ITG-like
	  //const difference_type terminal_span_size = 0;
	  
	  const size_type span_length_min = utils::bithack::branch(pos == source.goal, target.size(), span.size());
	  const size_type span_length_max = utils::bithack::min(target.size(), span.size() + terminal_span_size);
	  //
	  // 0 <= span_first
	  // span_last <= target.size()
	  // span_first <= span.first
	  // span.last <= span_last
	  // span_last - span_first < span_length_max
	  //
	      
	  for (size_type span_length = span_length_min; span_length <= span_length_max; ++ span_length) {
	    const size_type span_first_min = utils::bithack::max(difference_type(0), difference_type(span.last) - difference_type(span_length));
	    const size_type span_first_max = utils::bithack::min(span.first, target.size() - span_length);
		
	    for (size_type span_first = span_first_min; span_first <= span_first_max; ++ span_first) {
	      const size_type span_last = span_first + span_length;
		  
	      // ITG-like!
	      //if (span_length - total_hole) continue;
		  
	      const symbol_set_type& target_string = graph.rule_string(target, span_type(span_first, span_last), antecedent);
		  
	      const logprob_type prob = synalign.prob(edge.rule, target_string);
	      
	      graph.beta[pos](span_first, span_last) = std::max(graph.beta[pos](span_first, span_last), prob * prob_antecedent);
	      
	      graph.edges[pos][span_type(span_first, span_last)].push_back(edge_type(edge.id, target_string, antecedent, prob));
	    }
	  }
	}
      } else {
	for (size_type pos = 0; pos != j_ends[last]; ++ pos) {
	  j[last] = pos;
	  
	  const span_type& span = graph.edges[edge.tails[last]][pos].first;
	  
	  // check if valid or not,..
	  bool valid = true;
	  for (size_type prev = 0; prev != last && valid; ++ prev)
	    valid = disjoint(span, graph.edges[edge.tails[prev]][j[prev]].first);
	  
	  if (valid)
	    operator()(edge, j, j_ends, antecedent, last + 1);
	}
      }
    }
    
    PYPGraph& graph;

    const hypergraph_type& source;
    const sentence_type&   target;
    const PYPSynAlign&     synalign;
    const size_type        terminal_max;
  };


  logprob_type inside(const hypergraph_type& source, const sentence_type& target, const PYPSynAlign& synalign, const size_type terminal_max)
  {    
    typedef std::vector<size_type, std::allocator<size_type> > index_set_type;

    initialize(source, target, synalign);

    span_set_type antecedent;
    index_set_type j;
    index_set_type j_ends;

    Inside insider(*this, source, target, synalign, terminal_max);

    hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
      const hypergraph_type::node_type& node = *niter;
      
#if 0
      const size_type pos = node.id;
#endif

      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = source.edges[*eiter];

	antecedent.resize(edge.tails.size());
	j.resize(edge.tails.size());
	j_ends.resize(edge.tails.size());
	
	for (size_type i = 0; i != edge.tails.size(); ++ i) {
	  j_ends[i] = edges[edge.tails[i]].size();
	  j[i]= 0;
	}
	
	insider(edge, j, j_ends, antecedent, 0);
	
#if 0
	for (;;) {
	  span_type span(target.size(), 0);
	  size_type total_hole = 0;
	  logprob_type prob_antecedent = cicada::semiring::traits<logprob_type>::one();
	  
	  bool valid = true;
	  for (size_t i = 0; i != edge.tails.size() && valid; ++ i) {
	    antecedent[i] = edges[edge.tails[i]][j[i]].first;
	    
	    const span_type& span_antecedent = antecedent[i];
	    
	    if (span_antecedent.empty())
	      prob_antecedent *= beta[edge.tails[i]](span_antecedent.first, span_antecedent.last);
	    else {
	      for (size_t prev = 0; prev != i && valid; ++ prev)
		valid = disjoint(span_antecedent, edges[edge.tails[prev]][j[prev]].first);
	      
	      if (! valid) break;
	      
	      span.first = utils::bithack::min(span.first, span_antecedent.first);
	      span.last  = utils::bithack::max(span.last,  span_antecedent.last);
	      
	      total_hole += span_antecedent.size();
	      
	      prob_antecedent *= beta[edge.tails[i]](span_antecedent.first, span_antecedent.last);
	    }
	  }
	  
	  if (valid) {
	    if (span.first == target.size() && span.last == 0) {
	      // if phrases, we will allow all the lengths..
	      // otherwise, terminal_max!
	      
	      const size_type span_length_min = utils::bithack::branch(pos == source.goal, target.size(), size_type(0));
	      //const size_type span_length_max = utils::bithack::branch(antecedent.empty(), target.size(), utils::bithack::min(target.size(), terminal_max));
	      const size_type span_length_max = utils::bithack::min(target.size(), terminal_max);
	      
	      for (size_type span_length = span_length_min; span_length <= span_length_max; ++ span_length) {
		const size_type span_first_max = (span_length != 0) * (target.size() - span_length);
		for (size_type span_first = 0; span_first <= span_first_max; ++ span_first) {
		  const size_type span_last = span_first + span_length;
		  
		  // ITG-like!
		  //if (! antecedent.empty() && span_length) continue;
		  
		  const symbol_set_type& target_string = rule_string(target, span_type(span_first, span_last), antecedent);
		  
		  const logprob_type prob = synalign.prob(edge.rule, target_string);
		  
		  beta[pos](span_first, span_last) = std::max(beta[pos](span_first, span_last), prob * prob_antecedent);
		  
		  edges[pos][span_type(span_first, span_last)].push_back(edge_type(edge.id, target_string, antecedent, prob));
		}
	      }
	    } else {
	      const size_type terminal_hole_size = span.size() - total_hole;
	      const difference_type terminal_span_size = terminal_max - terminal_hole_size;
	      // ITG-like
	      //const difference_type terminal_span_size = 0;
	      
	      const size_type span_length_min = utils::bithack::branch(pos == source.goal, target.size(), span.size());
	      const size_type span_length_max = utils::bithack::min(target.size(), span.size() + terminal_span_size);

	      //
	      // 0 <= span_first
	      // span_last <= target.size()
	      // span_first <= span.first
	      // span.last <= span_last
	      // span_last - span_first < span_length_max
	      //
	      
	      for (size_type span_length = span_length_min; span_length <= span_length_max; ++ span_length) {
		const size_type span_first_min = utils::bithack::max(difference_type(0), difference_type(span.last) - difference_type(span_length));
		const size_type span_first_max = utils::bithack::min(span.first, target.size() - span_length);
		
		for (size_type span_first = span_first_min; span_first <= span_first_max; ++ span_first) {
		  const size_type span_last = span_first + span_length;
		  
		  // ITG-like!
		  //if (span_length - total_hole) continue;
		  
		  const symbol_set_type& target_string = rule_string(target, span_type(span_first, span_last), antecedent);
		  
		  const logprob_type prob = synalign.prob(edge.rule, target_string);
		  
		  beta[pos](span_first, span_last) = std::max(beta[pos](span_first, span_last), prob * prob_antecedent);
		  
		  edges[pos][span_type(span_first, span_last)].push_back(edge_type(edge.id, target_string, antecedent, prob));
		}
	      }
	    }
	  }

	  size_type index = 0;
	  for (/**/; index != edge.tails.size(); ++ index) {
	    ++ j[index];
	    if (j[index] < j_ends[index]) break;
	    j[index] = 0;
	  }
	  
	  // finished!
	  if (index == edge.tails.size()) break;
	}
#endif
      }
    }
    
    return beta[source.goal](0, target.size());
  }
  
  
  // backward sampling
  template <typename Sampler>
  logprob_type outside(const hypergraph_type& source, const sentence_type& target, Sampler& sampler, derivation_type& derivation, const double temperature)
  {
    typedef std::pair<hypergraph_type::id_type, span_type> value_type;
    typedef std::vector<value_type, std::allocator<value_type> > stack_type;

    derivation.clear();
    
    logprob_type prob_derivation = cicada::semiring::traits<logprob_type>::one();

    stack_type stack;
    stack.push_back(std::make_pair(source.goal, span_type(0, target.size())));
    
    while (! stack.empty()) {
      const value_type value = stack.back();
      stack.pop_back();
      
      span_edge_map_type::const_iterator miter = edges[value.first].find(value.second);
      if (miter == edges[value.first].end()) {
	// no derivation!
	derivation.clear();
	return logprob_type();
      }

      //const edge_set_type& edges_mapped = miter->second;
      const edge_set_type& edges_mapped = edges[value.first][miter - edges[value.first].begin()].second;

      logprob_type logsum;
      logprobs.clear();
      edge_set_type::const_iterator eiter_end = edges_mapped.end();
      for (edge_set_type::const_iterator eiter = edges_mapped.begin(); eiter != eiter_end; ++ eiter) {
	logprobs.push_back(eiter->prob);
	logsum += eiter->prob;
      }
      
      probs.clear();
      logprob_set_type::const_iterator liter_end = logprobs.end();
      for (logprob_set_type::const_iterator liter = logprobs.begin(); liter != liter_end; ++ liter)
	probs.push_back(*liter / logsum);
      
      prob_set_type::const_iterator piter = sampler.draw(probs.begin(), probs.end(), temperature);
      const size_type pos = piter - probs.begin();
      
      derivation.push_back(rule_pair_type(source.edges[edges_mapped[pos].edge].rule, edges_mapped[pos].target));
      prob_derivation *= edges_mapped[pos].prob;
      
      if (! source.edges[edges_mapped[pos].edge].tails.empty())
	for (difference_type i = source.edges[edges_mapped[pos].edge].tails.size() - 1; i >= 0; -- i)
	  stack.push_back(std::make_pair(source.edges[edges_mapped[pos].edge].tails[i],
					 edges_mapped[pos].antecedent[i]));
    }
    
    return prob_derivation;
  }


  span_edge_chart_type edges;
  
  beta_chart_type beta;
  rule_string_cache_chart_type rule_strings;

  logprob_set_type logprobs;
  prob_set_type    probs;
};


typedef boost::filesystem::path path_type;
typedef utils::sampler<boost::mt19937> sampler_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
typedef std::vector<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef PYPGraph::size_type size_type;
typedef PYPGraph::derivation_type derivation_type;
typedef std::vector<derivation_type, std::allocator<derivation_type> > derivation_set_type;
typedef std::vector<size_type, std::allocator<size_type> > position_set_type;

struct Task
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;  
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;
  
  Task(queue_type& __mapper,
       queue_type& __reducer,
       const hypergraph_set_type& __sources,
       const sentence_set_type& __targets,
       derivation_set_type& __derivations,
       derivation_set_type& __derivations_prev,
       const PYPSynAlign& __model,
       sampler_type& __sampler,
       const int& __max_terminal)
    : mapper(__mapper),
      reducer(__reducer),
      sources(__sources),
      targets(__targets),
      derivations(__derivations),
      derivations_prev(__derivations_prev),
      model(__model),
      sampler(__sampler),
      max_terminal(__max_terminal) {}

  void operator()()
  {
    size_type pos;
    
    for (;;) {
      mapper.pop(pos);
      
      if (pos == size_type(-1)) break;

      derivations_prev[pos]  = derivations[pos];

      if (! derivations_prev[pos].empty()) {
	derivation_type::const_iterator diter_end = derivations_prev[pos].end();
	for (derivation_type::const_iterator diter = derivations_prev[pos].begin(); diter != diter_end; ++ diter)
	  model.decrement(*diter, sampler);
      }
      
      graph.inside(sources[pos], targets[pos], model, max_terminal);
      
      graph.outside(sources[pos], targets[pos], sampler, derivations[pos], temperature);
      
      derivation_type::const_iterator diter_end = derivations[pos].end();
      for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	model.increment(*diter, sampler, temperature);
      
      reducer.push(pos);
    }
  }
  
  queue_type& mapper;
  queue_type& reducer;
  
  const hypergraph_set_type& sources;
  const sentence_set_type& targets;
  derivation_set_type& derivations;
  derivation_set_type& derivations_prev;
  
  PYPSynAlign  model;
  sampler_type sampler;
  int max_terminal;
  
  double temperature;
  
  PYPGraph graph;
};

struct less_size
{
  less_size(const hypergraph_set_type& __sources,
	    const sentence_set_type& __targets)
    : sources(__sources), targets(__targets) {}

  bool operator()(const size_type& x, const size_type& y) const
  {
    return (sources[x].edges.size() + targets[x].size()) < (sources[y].edges.size() + targets[y].size());
  }
  
  const hypergraph_set_type& sources;
  const sentence_set_type& targets;
};

path_type train_source_file = "-";
path_type train_target_file = "-";

path_type test_source_file;
path_type test_target_file;

int max_terminal = 2;

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

double lexicon_discount = 0.9;
double lexicon_strength = 1;

double lexicon_discount_prior_alpha = 1.0;
double lexicon_discount_prior_beta  = 1.0;
double lexicon_strength_prior_shape = 1.0;
double lexicon_strength_prior_rate  = 1.0;

double fertility_discount = 0.9;
double fertility_strength = 1;

double fertility_discount_prior_alpha = 1.0;
double fertility_discount_prior_beta  = 1.0;
double fertility_strength_prior_shape = 1.0;
double fertility_strength_prior_rate  = 1.0;

double distortion_discount = 0.9;
double distortion_strength = 1;

double distortion_discount_prior_alpha = 1.0;
double distortion_discount_prior_beta  = 1.0;
double distortion_strength_prior_shape = 1.0;
double distortion_strength_prior_rate  = 1.0;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);

size_t read_data(const path_type& path, hypergraph_set_type& graphs);
size_t read_data(const path_type& path, sentence_set_type& sentences);

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

    sentence_set_type   targets;
    const size_type target_vocab_size = read_data(train_target_file, targets);

    if (targets.empty())
      throw std::runtime_error("no target sentence data?");
    
    hypergraph_set_type sources;
    sources.reserve(targets.size());
    
    const size_type source_vocab_size = read_data(train_source_file, sources);
    
    if (sources.size() != targets.size())
      throw std::runtime_error("source/target side do not match!");
    
    PYPLexicon   lexicon((1.0 / source_vocab_size) * (1.0 / target_vocab_size),
			 PYPLexicon::parameter_type(lexicon_discount,
						    lexicon_strength,
						    lexicon_discount_prior_alpha,
						    lexicon_discount_prior_beta,
						    lexicon_strength_prior_shape,
						    lexicon_strength_prior_rate));
    
    PYPFertility   fertility(PYPFertility::parameter_type(fertility_discount,
							  fertility_strength,
							  fertility_discount_prior_alpha,
							  fertility_discount_prior_beta,
							  fertility_strength_prior_shape,
							  fertility_strength_prior_rate));
    
    PYPDistortion   distortion(PYPDistortion::parameter_type(distortion_discount,
							     distortion_strength,
							     distortion_discount_prior_alpha,
							     distortion_discount_prior_beta,
							     distortion_strength_prior_shape,
							     distortion_strength_prior_rate));
    
    PYPSynAlign  synalign(lexicon,
			  fertility,
			  distortion,
			  PYPSynAlign::parameter_type(discount,
						      strength,
						      discount_prior_alpha,
						      discount_prior_beta,
						      strength_prior_shape,
						      strength_prior_rate));

    
    derivation_set_type derivations(sources.size());
    derivation_set_type derivations_prev(sources.size());
    position_set_type positions;
    for (size_t i = 0; i != sources.size(); ++ i)
      if (sources[i].is_valid() && ! targets[i].empty())
	positions.push_back(i);
    position_set_type(positions).swap(positions);
    
    sampler_type sampler;
    
    // sample parameters, first...
    if (slice_sampling)
      synalign.slice_sample_parameters(sampler, resample_iterations);
    else
      synalign.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2)
      std::cerr << "rule: discount=" << synalign.table.discount() << " strength=" << synalign.table.strength() << std::endl
		<< "lexicon: discount=" << synalign.lexicon.table.discount() << " strength=" << synalign.lexicon.table.strength() << std::endl
		<< "fertility: discount=" << synalign.fertility.fallback.discount() << " strength=" << synalign.fertility.fallback.strength() << std::endl
		<< "distortion: discount=" << synalign.distortion.fallback.discount() << " strength=" << synalign.distortion.fallback.strength() << std::endl;
    
    PYPGraph graph;
    
    Task::queue_type queue_mapper;
    Task::queue_type queue_reducer;
    
    std::vector<Task, std::allocator<Task> > tasks(threads, Task(queue_mapper,
								 queue_reducer,
								 sources,
								 targets,
								 derivations,
								 derivations_prev,
								 synalign,
								 sampler,
								 max_terminal));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);

    size_t baby_iter = 0;
    const size_t baby_last = utils::bithack::branch(baby_steps > 0, baby_steps, 0);

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
      
      
      bool baby_finished = true;
      if (baby_iter != baby_last) {
	++ baby_iter;
	baby_finished = false;
      }
      
      sampling = anneal_finished && baby_finished;
      
      if (debug) {
	if (sampling)
	  std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
	else
	  std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
      }

      // assign temperature... and current model
      for (size_type i = 0; i != tasks.size(); ++ i) {
	tasks[i].model = synalign;
	tasks[i].temperature = temperature;
      }
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(positions.begin(), positions.end(), gen);
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_size(sources, targets));
      
      position_set_type::const_iterator piter_end = positions.end();
      position_set_type::const_iterator piter = positions.begin();
      size_type reduced = 0;
      while (piter != piter_end || reduced != positions.size()) {
	
	if (piter != piter_end && queue_mapper.push(*piter, true))
	  ++ piter;
	
	size_type pos = 0;
	if (reduced != positions.size() && queue_reducer.pop(pos, true)) {
	  ++ reduced;
	  
	  if (! derivations_prev[pos].empty()) {
	    derivation_type::const_iterator diter_end = derivations_prev[pos].end();
	    for (derivation_type::const_iterator diter = derivations_prev[pos].begin(); diter != diter_end; ++ diter)
	      synalign.decrement(*diter, sampler);
	  }

	  if (debug >= 3) {
	    std::cerr << "training=" << pos
		      << " nodes=" << sources[pos].nodes.size()
		      << " edges=" << sources[pos].edges.size()
		      << " sentence=" << targets[pos].size() << std::endl;
	    
	    derivation_type::const_iterator diter_end = derivations[pos].end();
	    for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	      std::cerr << "derivation: " << *(diter->source) << " ||| " << diter->target << std::endl;
	  }
	  
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    synalign.increment(*diter, sampler, temperature);
	}
      }

#if 0      
      for (size_t i = 0; i != positions.size(); ++ i) {
	const size_t pos = positions[i];
	
	if (debug >= 3)
	  std::cerr << "training=" << pos << std::endl;
	
	if (! derivations[pos].empty()) {
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    synalign.decrement(*diter, sampler);
	}
	
	const PYPGraph::logprob_type logsum = graph.inside(sources[pos], targets[pos], synalign, max_terminal);
	
	const PYPGraph::logprob_type logderivation = graph.outside(sources[pos], targets[pos], sampler, derivations[pos], temperature);
	
	if (debug >= 3) {
	  std::cerr << "sum=" << logsum << " derivation=" << logderivation << std::endl;
	  
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    std::cerr << "derivation: " << *(diter->source) << " ||| " << diter->target << std::endl;
	}
	
	derivation_type::const_iterator diter_end = derivations[pos].end();
	for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	  synalign.increment(*diter, sampler, temperature);
      }
#endif
      
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  synalign.slice_sample_parameters(sampler, resample_iterations);
	else
	  synalign.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2)
	  std::cerr << "rule: discount=" << synalign.table.discount() << " strength=" << synalign.table.strength() << std::endl
		    << "lexicon: discount=" << synalign.lexicon.table.discount() << " strength=" << synalign.lexicon.table.strength() << std::endl
		    << "fertility: discount=" << synalign.fertility.fallback.discount() << " strength=" << synalign.fertility.fallback.strength() << std::endl
		    << "distortion: discount=" << synalign.distortion.fallback.discount() << " strength=" << synalign.distortion.fallback.strength() << std::endl;
		    
      }
      
      if (debug)
	std::cerr << "log-likelihood: " << synalign.log_likelihood() << std::endl;
	    
    }

    for (int i = 0; i != threads; ++ i)
      queue_mapper.push(size_type(-1));
    
    workers.join_all();

    // dump models...
    // perform testing...?
    
    
    
    
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


size_t read_data(const path_type& path, hypergraph_set_type& graphs)
{
  typedef utils::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_set_type;
  
  word_set_type words;
  words.set_empty_key(word_type());
  
  graphs.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  
  hypergraph_type graph;
  while (is >> graph) {
    graphs.push_back(graph);
    
    hypergraph_type::edge_set_type::const_iterator eiter_end = graphs.back().edges.end();
    for (hypergraph_type::edge_set_type::const_iterator eiter = graphs.back().edges.begin(); eiter != eiter_end; ++ eiter) {
      const hypergraph_type::edge_type& edge = *eiter;

      // assume penn-treeebank style
      if (edge.rule->rhs.size() == 1 && edge.rule->rhs.front().is_terminal())
	words.insert(edge.rule->rhs.front());
    }
  }
  
  return words.size();
}


size_t read_data(const path_type& path, sentence_set_type& sentences)
{
  typedef utils::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_set_type;

  word_set_type words;
  words.set_empty_key(word_type());

  sentences.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  
  sentence_type sentence;
  while (is >> sentence) {
    sentences.push_back(sentence);
    
    words.insert(sentence.begin(), sentence.end());
  }
  
  sentence_set_type(sentences).swap(sentences);

  return words.size();
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

    ("max-terminal", po::value<int>(&max_terminal)->default_value(max_terminal), "max # of terminals in each rule")
    
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
    
    ("lexicon-discount",       po::value<double>(&lexicon_discount)->default_value(lexicon_discount),                         "discount ~ Beta(alpha,beta)")
    ("lexicon-discount-alpha", po::value<double>(&lexicon_discount_prior_alpha)->default_value(lexicon_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("lexicon-discount-beta",  po::value<double>(&lexicon_discount_prior_beta)->default_value(lexicon_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("lexicon-strength",       po::value<double>(&lexicon_strength)->default_value(lexicon_strength),                         "strength ~ Gamma(shape,rate)")
    ("lexicon-strength-shape", po::value<double>(&lexicon_strength_prior_shape)->default_value(lexicon_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("lexicon-strength-rate",  po::value<double>(&lexicon_strength_prior_rate)->default_value(lexicon_strength_prior_rate),   "strength ~ Gamma(shape,rate)")

    ("fertility-discount",       po::value<double>(&fertility_discount)->default_value(fertility_discount),                         "discount ~ Beta(alpha,beta)")
    ("fertility-discount-alpha", po::value<double>(&fertility_discount_prior_alpha)->default_value(fertility_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("fertility-discount-beta",  po::value<double>(&fertility_discount_prior_beta)->default_value(fertility_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("fertility-strength",       po::value<double>(&fertility_strength)->default_value(fertility_strength),                         "strength ~ Gamma(shape,rate)")
    ("fertility-strength-shape", po::value<double>(&fertility_strength_prior_shape)->default_value(fertility_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("fertility-strength-rate",  po::value<double>(&fertility_strength_prior_rate)->default_value(fertility_strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
    ("distortion-discount",       po::value<double>(&distortion_discount)->default_value(distortion_discount),                         "discount ~ Beta(alpha,beta)")
    ("distortion-discount-alpha", po::value<double>(&distortion_discount_prior_alpha)->default_value(distortion_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("distortion-discount-beta",  po::value<double>(&distortion_discount_prior_beta)->default_value(distortion_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("distortion-strength",       po::value<double>(&distortion_strength)->default_value(distortion_strength),                         "strength ~ Gamma(shape,rate)")
    ("distortion-strength-shape", po::value<double>(&distortion_strength_prior_shape)->default_value(distortion_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("distortion-strength-rate",  po::value<double>(&distortion_strength_prior_rate)->default_value(distortion_strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
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


