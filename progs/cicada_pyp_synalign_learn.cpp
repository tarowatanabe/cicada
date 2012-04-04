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
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/dense_hash_map.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/indexed_map.hpp"

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

// unigram model
struct UnigramModel
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef std::vector<double, std::allocator<double> > table_type;
  
  template <typename Iterator>
  UnigramModel(Iterator first, Iterator last) : table(), smooth(-60) {  learn(first, last); }
  UnigramModel() : table(), smooth(-60) { }
  
  template <typename Iterator>
  void learn(Iterator first, Iterator last)
  {
    typedef std::vector<size_type, std::allocator<size_type> > counts_type;

    counts_type counts;
    for (/**/; first != last; ++ first)
      learn(first->begin(), first->end(), counts);
    
    double total = 0.0;
    size_type vocab = 0;
    typename counts_type::const_iterator citer_end = counts.end();
    for (typename counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      total += *citer;
      vocab += (*citer != 0);
    }
    
    smooth = utils::mathop::log(1.0 / vocab);
    
    table.clear();
    table.reserve(counts.size());
    table.resize(counts.size(), smooth);

    // a simple add-one smoothing
    
    const double denom = total + vocab;
    
    table_type::iterator titer = table.end();
    for (typename counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ titer)
      if (*citer)
	*titer = utils::mathop::log((1.0 + *citer) / denom);
  }

  template <typename Iterator, typename Counts>
  void learn(Iterator first, Iterator last, Counts& counts)
  {
    for (/**/; first != last; ++ first) {
      const symbol_type::id_type id = symbol_type(*first).id();
      
      if (id >= counts.size())
	counts.resize(id + 1, 0);
      
      ++ counts[id];
    }
  }
  
  
  double operator()(const symbol_type& x) const
  {
    return (x.id() >= table.size() ? smooth : table[x.id()]);
  }
  
  table_type table;
  double smooth;
};

// PCFG model
struct PCFGModel
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  struct rule_hash
  {
    size_t operator()(const rule_type* x) const
    {
      return (x ? hash_value(*x) : size_t(0));
    }
  };
  
  struct rule_equal
  {
    bool operator()(const rule_type* x, const rule_type* y) const
    {
      return (x == y) || (x && y && *x == *y);
    }
  };
  
  typedef utils::unordered_map<const rule_type*, double, rule_hash, rule_equal,
			       std::allocator<std::pair<const rule_type*, double> > >::type table_type;
  
  template <typename Iterator>
  PCFGModel(Iterator first, Iterator last) : table(), smooth(-60) { learn(first, last); }
  PCFGModel() : table(), smooth(-60) { }
  
  template <typename Iterator>
  void learn(Iterator first, Iterator last)
  {
    typedef utils::unordered_map<const rule_type*, size_type, rule_hash, rule_equal,
				 std::allocator<std::pair<const rule_type*, size_type> > >::type count_type;
    typedef utils::unordered_map<symbol_type, count_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				 std::allocator<std::pair<const symbol_type, count_type> > >::type count_set_type;
    
    count_set_type counts;
    
    for (/**/; first != last; ++ first) {
      const hypergraph_type& graph = *first;
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	const rule_type* rule = &(*edge.rule);
	
	++ counts[rule->lhs][rule];
      }
    }

    table.clear();
    
    // a simple add-one smoothing...
    size_type observed = 0;
    
    typename count_set_type::const_iterator siter_end = counts.end();
    for (typename count_set_type::const_iterator siter = counts.begin(); siter != siter_end; ++ siter) {
      observed += siter->second.size();
      
      size_type total = 0;
      typename count_type::const_iterator citer_end = siter->second.end();
      for (typename count_type::const_iterator citer = siter->second.begin(); citer != citer_end; ++ citer)
	total += citer->second;
      
      const double denom = total + siter->second.size();
      for (typename count_type::const_iterator citer = siter->second.begin(); citer != citer_end; ++ citer)
	table[citer->first] = utils::mathop::log((1.0 + citer->second) / denom);
    }
    
    smooth = utils::mathop::log(1.0 / observed);
  }

  double operator()(const rule_ptr_type& x) const
  {
    return operator()(x.get());
  }
  
  double operator()(const rule_type* x) const
  {
    table_type::const_iterator iter = table.find(x);
    return (iter != table.end() ? iter->second : smooth);
  }
  
  table_type table;
  double smooth;
};

struct LengthModel
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef utils::array_power2<double, 32, std::allocator<double> > cache_type;
  
  LengthModel(const double& __lambda)
    : lambda(__lambda) { initialize(lambda); }

  double operator()(const size_type size) const
  {
    return (size < cache.size() ? cache[size] : utils::mathop::log_geometric(size, lambda));
  }
  
  void initialize(const double& __lambda)
  {
    cache.clear();
    lambda = __lambda;
    
    for (size_type size = 0; size != cache.size(); ++ size)
      cache[size] = utils::mathop::log_geometric(size, lambda);
  }
  
  cache_type cache;
  
  double lambda;
};

struct NonTerminalModel
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef utils::array_power2<double, 32, std::allocator<double> > cache_type;

  NonTerminalModel()
  {
    cache[0] = - std::numeric_limits<double>::infinity();
    
    for (size_type size = 1; size != cache.size(); ++ size)
      cache[size] = utils::mathop::log(1.0 / size);
  }

  double operator()(const size_type size) const
  {
    return (size < cache.size() ? cache[size] : utils::mathop::log(1.0 / size));
  }

  double operator()(const size_type size, const size_type num) const
  {
    double logprob = 0.0;
    for (size_type i = size; i != size + num; ++ i)
      logprob += operator()(i + 1);
    return logprob;
  }
  
  cache_type cache;
};

struct PYPSynAlign
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  struct rule_pair_type
  {
    const rule_type* source;
    symbol_set_type  target;
    
    rule_pair_type() : source(0), target() {}
    rule_pair_type(const rule_type* __source, const symbol_set_type& __target) : source(__source), target(__target) {}
    rule_pair_type(const rule_ptr_type& __source, const symbol_set_type& __target) : source(&(*__source)), target(__target) {}
    
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

  PYPSynAlign(PCFGModel&    __pcfg_source,
	      UnigramModel& __unigram_target,
	      LengthModel&  __length_target,
	      const double __discount,
	      const double __strength,
	      const double __discount_alpha,
	      const double __discount_beta,
	      const double __strength_shape,
	      const double __strength_rate)
    : pcfg_source(__pcfg_source),
      unigram_target(__unigram_target),
      length_target(__length_target),
      non_terminal_target(),
      table(__discount, __strength, __discount_alpha, __discount_beta, __strength_shape, __strength_rate) {}
  
  template <typename Sampler>
  bool increment(const rule_type* source, const symbol_set_type& target, Sampler& sampler, const double temperature=1.0)
  {
    return table.increment(rule_pair_type(source, target), prob_base(source, target), sampler, temperature);
  }

  template <typename Sampler>
  bool increment(const rule_pair_type& rule_pair, Sampler& sampler, const double temperature=1.0)
  {
    return table.increment(rule_pair, prob_base(rule_pair), sampler, temperature);
  }
  
  template <typename Sampler>
  bool decrement(const rule_type* source, const symbol_set_type& target, Sampler& sampler)
  {
    return table.decrement(rule_pair_type(source, target), sampler);
  }
  
  template <typename Sampler>
  bool decrement(const rule_pair_type& rule_pair, Sampler& sampler)
  {
    return table.decrement(rule_pair, sampler);
  }
  
  double prob(const rule_type* source, const symbol_set_type& target, const double base) const
  {
    return table.prob(rule_pair_type(source, target), base);
  }

  double prob(const rule_ptr_type& source, const symbol_set_type& target, const double base) const
  {
    return table.prob(rule_pair_type(&(*source), target), base);
  }

  double logprob(const rule_type* source, const symbol_set_type& target, const double logbase) const
  {
    return cicada::semiring::log(table.prob(rule_pair_type(source, target), cicada::semiring::traits<logprob_type>::exp(logbase)));
  }
  
  double logprob(const rule_ptr_type& source, const symbol_set_type& target, const double logbase) const
  {
    return cicada::semiring::log(table.prob(rule_pair_type(&(*source), target), cicada::semiring::traits<logprob_type>::exp(logbase)));
  }
  
  double prob_base(const rule_pair_type& rule_pair) const
  {
    return prob_base(rule_pair.source, rule_pair.target);
  }
  
  double prob_base(const rule_type* source, const symbol_set_type& target) const
  {
    double logprob = pcfg_source(source);
    
    size_type length = 0;
    symbol_set_type::const_iterator titer_end = target.end();
    for (symbol_set_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) 
      if (titer->is_terminal()) {
	logprob += unigram_target(*titer);
	++ length;
      }
    
    logprob += length_target(length);
    logprob += non_terminal_target(length, target.size() - length);
    
    return utils::mathop::exp(logprob);
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
  void sample_length_parameters(Sampler& sampler)
  {
    // do nothing...!
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    sample_length_parameters(sampler);
    
    table.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    sample_length_parameters(sampler);
    
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  PCFGModel&       pcfg_source;
  UnigramModel&    unigram_target;
  LengthModel&     length_target;
  NonTerminalModel non_terminal_target;
  
  table_type table;
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
    
    edge_type() {}
    edge_type(const hypergraph_type::id_type& __edge,
	      const symbol_set_type& __target,
	      const span_set_type& __antecedent,
	      const logprob_type& __prob)
      : edge(__edge), target(__target), antecedent(__antecedent), prob(__prob) {}
  };
  
  typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
  
  typedef utils::indexed_map<span_type, edge_set_type, boost::hash<span_type>, std::equal_to<span_type>,
			     std::allocator<std::pair<span_type, edge_set_type> > > span_edge_map_type;
  typedef std::vector<span_edge_map_type, std::allocator<span_edge_map_type> > span_edge_chart_type;

  typedef std::vector<logprob_type, std::allocator<logprob_type> > logprob_set_type;
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  void initialize(const hypergraph_type& source, const sentence_type& target, const PYPSynAlign& synalign)
  {
    std::cerr << "source nodes: " << source.nodes.size() << " edges: " << source.edges.size() << std::endl
	      << "targets: " << target.size() << std::endl;
    
    edges.clear();
    edges.resize(source.nodes.size(), span_edge_map_type());

    beta.clear();
    beta.resize(source.nodes.size(), beta_type(target.size() + 1));

    // phrasal log-probabilities
    logprob_targets.clear();
    logprob_targets.resize(target.size() + 1, 0.0);
    
    for (size_type first = 0; first != target.size(); ++ first)
      for (size_type last = first + 1; last <= target.size(); ++ last) {
	std::cerr << "phrase span: " << first << "..." << last << std::endl;
	
	logprob_targets(first, last) = logprob_targets(first, last - 1) + synalign.unigram_target(target[last - 1]);
      }
    
    // length log-probabilities
    logprob_lengths.clear();
    logprob_lengths.resize(target.size() + 1);
    
    for (size_type i = 0; i <= target.size(); ++ i) {
      logprob_lengths[i] = synalign.length_target(i);
      
      std::cerr << "length logprob: " << i << " " << logprob_lengths[i] << std::endl;
    }
    
    // rule-string...
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
    
    rule_string_cache_type& cache = rule_strings(span.first, span.last);
    
    rule_string_cache_type::iterator iter = cache.find(antecedent);
    if (iter == cache.end()) {
      if (antecedent.empty())
	iter = cache.insert(std::make_pair(antecedent, symbol_set_type(sentence.begin() + span.first, sentence.begin() + span.last))).first;
      else {
	buffer_type buffer;
	span_index_set_type span_index;
	
	span_index.reserve(antecedent.size());
	span_set_type::const_iterator aiter_end = antecedent.end();
	for (span_set_type::const_iterator aiter = antecedent.begin(); aiter != aiter_end; ++ aiter)
	  span_index.push_back(std::make_pair(*aiter, aiter - antecedent.begin()));
	
	size_type pos = span.first;
	span_index_set_type::const_iterator siter_end = span_index.end();
	for (span_index_set_type::const_iterator siter = span_index.begin(); siter != siter_end; ++ siter) {
	  if (! siter->first.empty()) {
	    buffer.insert(buffer.end(), sentence.begin() + pos, sentence.begin() + siter->first.first);
	    buffer.push_back(vocab_type::X.non_terminal(siter->second));
	    
	    pos = siter->first.last;
	  } else
	    buffer.push_back(vocab_type::X.non_terminal(siter->second));
	}
	
	buffer.insert(buffer.end(), sentence.begin() + pos, sentence.begin() + span.last);
	
	iter = cache.insert(std::make_pair(antecedent, symbol_set_type(buffer.begin(), buffer.end()))).first;
      }
    }
    
    return iter->second;
  }
  

  logprob_type inside(const hypergraph_type& source, const sentence_type& target, const PYPSynAlign& synalign, const size_type terminal_max)
  {    
    typedef std::vector<int, std::allocator<int> > index_set_type;

    initialize(source, target, synalign);

    span_set_type antecedent;

    hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
      const hypergraph_type::node_type& node = *niter;
      
      const size_type pos = node.id;

      std::cerr << "node pos: " << pos << std::endl;
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = source.edges[*eiter];

	const double logprob_source = synalign.pcfg_source(edge.rule);

	std::cerr << "rule: " << *edge.rule << " logprob: " << logprob_source << std::endl;
	
	index_set_type j_ends(edge.tails.size(), 0);
	index_set_type j(edge.tails.size(), 0);
	
	for (size_type i = 0; i != edge.tails.size(); ++ i)
	  j_ends[i] = edges[edge.tails[i]].size();
	
	for (;;) {
	  span_type span(target.size(), 0);
	  size_type total_hole = 0;
	  logprob_type prob_antecedent = cicada::semiring::traits<logprob_type>::one();
	  double logprob_hole = 0.0;
	  
	  antecedent.clear();
	  bool valid = true;
	  for (size_t i = 0; i != edge.tails.size() && valid; ++ i) {
	    const span_type& span_antecedent = edges[edge.tails[i]][j[i]].first;
	    
	    if (span_antecedent.empty()) continue;
	    
	    for (size_t prev = 0; prev != i && valid; ++ prev)
	      valid = disjoint(span_antecedent, edges[edge.tails[prev]][j[prev]].first);
	    
	    span.first = utils::bithack::min(span.first, span_antecedent.first);
	    span.last  = utils::bithack::max(span.last,  span_antecedent.last);
	    total_hole += span_antecedent.size();
	    logprob_hole += logprob_targets(span_antecedent.first, span_antecedent.last);
	    
	    antecedent.push_back(span_antecedent);
	    
	    prob_antecedent *= beta[edge.tails[i]](span_antecedent.first, span_antecedent.last);
	  }
	  
	  if (valid) {
	    
	    if (span.first == target.size() && span.last == 0) {
	      std::cerr << "phrasal pair or empty antecedents" << std::endl;
	      
	      const size_type span_length_min = utils::bithack::branch(pos == source.goal, target.size(), size_type(0));
	      const size_type span_length_max = utils::bithack::min(target.size(), terminal_max);
	      
	      for (size_type span_length = span_length_min; span_length <= span_length_max; ++ span_length) {
		const size_type span_first_max = (span_length != 0) * (target.size() - span_length);
		for (size_type span_first = 0; span_first <= span_first_max; ++ span_first) {
		  const size_type span_last = span_first + span_length;

		  
		  std::cerr << "phrase span: " << span_first << "..." << span_last << std::endl;

		  const double logprob_terminal = logprob_targets(span_first, span_last);

		  std::cerr << "logprob terminals: " << logprob_terminal << std::endl;
		  
		  const size_type terminal_size = span_length;
		  const double logprob_non_terminal = synalign.non_terminal_target(terminal_size, edge.tails.size());

		  std::cerr << "logprob non-temrinals: " << logprob_non_terminal << std::endl;
		  
		  std::cerr << "logprob length: " << logprob_lengths[terminal_size] << std::endl;
		  
		  const double logprob_prior = logprob_source + logprob_terminal + logprob_non_terminal + logprob_lengths[terminal_size];
		  
		  
		  const symbol_set_type& target_string = rule_string(target, span_type(span_first, span_last), antecedent);
		  
		  const logprob_type prob = (prob_antecedent
					     * cicada::semiring::traits<logprob_type>::exp(synalign.logprob(edge.rule, target_string, logprob_prior)));
		  
		  beta[pos](span_first, span_last) = std::max(beta[pos](span_first, span_last), prob);
		  edges[pos][span_type(span_first, span_last)].push_back(edge_type(edge.id, target_string, antecedent, prob));
		}
	      }
	    } else {
	      std::cerr << "antecedents pair" << std::endl;
	      std::cerr << "span: " << span.first << "..." << span.last << std::endl;

	      const size_type terminal_hole_size = span.size() - total_hole;
	      const difference_type terminal_span_size = terminal_max - terminal_hole_size;

	      const size_type span_length_min = utils::bithack::branch(pos == source.goal, target.size(), span.size());
	      const size_type span_length_max = utils::bithack::min(target.size(), span.size() + terminal_span_size);
	      
	      for (size_type span_length = span_length_min; span_length <= span_length_max; ++ span_length) {
		const size_type span_first_min = utils::bithack::max(difference_type(0), difference_type(span.first) - terminal_span_size);
		const size_type span_first_max = utils::bithack::min(target.size(), span.last + terminal_span_size) - span_length;
		for (size_type span_first = span_first_min; span_first <= span_first_max; ++ span_first) {
		  const size_type span_last = span_first + span_length;
		  const double logprob_terminal = logprob_targets(span_first, span_last) - logprob_hole;
		  
		  const size_type terminal_size = span_length - terminal_hole_size;
		  const double logprob_non_terminal = synalign.non_terminal_target(terminal_size, edge.tails.size());
		  
		  const double logprob_prior = logprob_source + logprob_terminal + logprob_non_terminal + logprob_lengths[terminal_size];
		  
		  const symbol_set_type& target_string = rule_string(target, span_type(span_first, span_last), antecedent);
		  
		  const logprob_type prob = (prob_antecedent
					     * cicada::semiring::traits<logprob_type>::exp(synalign.logprob(edge.rule, target_string, logprob_prior)));
		  
		  beta[pos](span_first, span_last) = std::max(beta[pos](span_first, span_last), prob);
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
      }
    }
    
    return beta[source.goal](0, target.size());
  }

  // backward sampling
  template <typename Sampler>
  logprob_type outside(const hypergraph_type& source, const sentence_type& target, Sampler& sampler, derivation_type& derivation)
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

      const edge_set_type& edges = miter->second;

      logprob_type logsum;
      logprobs.clear();
      edge_set_type::const_iterator eiter_end = edges.end();
      for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
	logprobs.push_back(eiter->prob);
	logsum += eiter->prob;
      }
      
      probs.clear();
      logprob_set_type::const_iterator liter_end = logprobs.end();
      for (logprob_set_type::const_iterator liter = logprobs.begin(); liter != liter_end; ++ liter)
	probs.push_back(*liter / logsum);
      
      prob_set_type::const_iterator piter = sampler.draw(probs.begin(), probs.end());
      const size_type pos = piter - probs.begin();
      
      derivation.push_back(rule_pair_type(source.edges[edges[pos].edge].rule, edges[pos].target));
      prob_derivation *= edges[pos].prob;
      
      if (! source.edges[edges[pos].edge].tails.empty())
	for (difference_type i = source.edges[edges[pos].edge].tails.size() - 1; i >= 0; -- i)
	  stack.push_back(std::make_pair(source.edges[edges[pos].edge].tails[i],
					 edges[pos].antecedent[i]));
    }
    
    return prob_derivation;
  }


  span_edge_chart_type edges;
  
  beta_chart_type beta;
  logprob_target_set_type logprob_targets;
  logprob_length_set_type logprob_lengths;
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

int max_terminal = 4;

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

// target string length...
double lambda = 0.5;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);

void read_data(const path_type& path, hypergraph_set_type& graphs);
void read_data(const path_type& path, sentence_set_type& sentences);

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
    read_data(train_target_file, targets);

    if (targets.empty())
      throw std::runtime_error("no target sentence data?");
    
    hypergraph_set_type sources;
    sources.reserve(targets.size());
    
    read_data(train_source_file, sources);
    
    if (sources.size() != targets.size())
      throw std::runtime_error("source/target side do not match!");
    
    
    PCFGModel    pcfg_source(sources.begin(), sources.end());
    UnigramModel unigram_target(targets.begin(), targets.end());
    LengthModel  length_target(lambda);

    PYPSynAlign  synalign(pcfg_source, unigram_target, length_target,
			  discount,
			  strength,
			  discount_prior_alpha,
			  discount_prior_beta,
			  strength_prior_shape,
			  strength_prior_rate);

    
    derivation_set_type derivations(sources.size());
    position_set_type positions(sources.size());
    for (size_t i = 0; i != sources.size(); ++ i)
      positions[i] = i;
    
    sampler_type sampler;
    
    // sample parameters, first...
    if (slice_sampling)
      synalign.slice_sample_parameters(sampler, resample_iterations);
    else
      synalign.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2)
      std::cerr << "discount=" << synalign.table.discount()
		<< " strength=" << synalign.table.strength()
		<< std::endl;
    
    PYPGraph graph;
    
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
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(positions.begin(), positions.end(), gen);
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_size(sources, targets));
      
      for (size_t i = 0; i != positions.size(); ++ i) {
	const size_t pos = positions[i];
	
	if (! sources[pos].is_valid() || targets[pos].empty()) continue;
	
	if (debug >= 3)
	  std::cerr << "training=" << pos << std::endl;
	
	if (! derivations[pos].empty()) {
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    synalign.decrement(*diter, sampler);
	}
	
	const PYPGraph::logprob_type logsum = graph.inside(sources[pos], targets[pos], synalign, max_terminal);
	
	const PYPGraph::logprob_type logderivation = graph.outside(sources[pos], targets[pos], sampler, derivations[pos]);
	
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
    }

  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct rule_ptr_hash
{
  bool operator()(const rule_ptr_type& x) const
  {
    return hash_value(*x);
  }
};

struct rule_ptr_equal
{
  bool operator()(const rule_ptr_type& x, const rule_ptr_type& y) const
  {
    return *x == *y;
  }
};

void read_data(const path_type& path, hypergraph_set_type& graphs)
{
  typedef utils::unordered_set<rule_ptr_type, rule_ptr_hash, rule_ptr_equal, std::allocator<rule_ptr_type> >::type rule_set_type;
  
  graphs.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  
  rule_set_type rules;
  hypergraph_type graph;
  while (is >> graph) {
#if 0
    hypergraph_type::edge_set_type::iterator eiter_end = graph.edges.end();
    for (hypergraph_type::edge_set_type::iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter)
      eiter->rule = *(rules.insert(eiter->rule).first);
#endif
    
    graphs.push_back(graph);
  }
  
  // uniquify rules...
}


void read_data(const path_type& path, sentence_set_type& sentences)
{
  sentences.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  
  sentence_type sentence;
  while (is >> sentence)
    sentences.push_back(sentence);
  
  sentence_set_type(sentences).swap(sentences);
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
    
    ("lambda",        po::value<double>(&lambda)->default_value(lambda),               "lambda for target terminals (Geometric distribution)")
    
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


