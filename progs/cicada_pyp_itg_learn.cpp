//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-PIALIGN!
// based on

// @InProceedings{neubig-EtAl:2011:ACL-HLT2011,
//   author    = {Neubig, Graham  and  Watanabe, Taro  and  Sumita, Eiichiro  and  Mori, Shinsuke  and  Kawahara, Tatsuya},
//   title     = {An Unsupervised Model for Joint Phrase Alignment and Extraction},
//   booktitle = {Proceedings of the 49th Annual Meeting of the Association for Computational Linguistics: Human Language Technologies},
//   month     = {June},
//   year      = {2011},
//   address   = {Portland, Oregon, USA},
//   publisher = {Association for Computational Linguistics},
//   pages     = {632--641},
//   url       = {http://www.aclweb.org/anthology/P11-1064}
// }


// parsing by
//
// @InProceedings{saers-nivre-wu:2009:IWPT09,
//   author    = {Saers, Markus  and  Nivre, Joakim  and  Wu, Dekai},
//   title     = {Learning Stochastic Bracketing Inversion Transduction Grammars with a Cubic Time Biparsing Algorithm},
//   booktitle = {Proceedings of the 11th International Conference on Parsing Technologies (IWPT'09)},
//   month     = {October},
//   year      = {2009},
//   address   = {Paris, France},
//   publisher = {Association for Computational Linguistics},
//   pages     = {29--32},
//   url       = {http://www.aclweb.org/anthology/W09-3804}
// }

// an extension: we employ infinite-ITG by memorizing the permutation
//
// we use std::string (or, symbol?) as an indiction of "permutation" flag
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/progress.hpp>

#include <map>
#include <deque>
#include <iterator>
#include <numeric>
#include <queue>
#include <set>

#include <cicada/sentence.hpp>
#include <cicada/alignment.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring/logprob.hpp>

#include "utils/vector2.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/pyp_parameter.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/chart.hpp"
#include "utils/bichart.hpp"
#include "utils/chunk_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/restaurant_sync.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/compact_map.hpp"
#include "utils/compact_set.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/symbol_set.hpp"
#include "utils/indexed_set.hpp"
#include "utils/rwticket.hpp"
#include "utils/atomicop.hpp"
#include "utils/hashmurmur3.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/array.hpp>

typedef cicada::Vocab      vocab_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Alignment  alignment_type;
typedef cicada::Symbol     symbol_type;
typedef cicada::Symbol     word_type;
typedef cicada::HyperGraph hypergraph_type;

struct PYP
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef uint32_t  id_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

  typedef enum {
    TERMINAL = 0,
    STRAIGHT,
    INVERTED,
  } itg_type;

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
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.source.id(), x.target.id());
    }
  };

  struct word_pair_unassigned : public utils::unassigned<word_type>
  {
    typedef utils::unassigned<word_type> unassigned_type;

    word_pair_type operator()() const
    {
      return word_pair_type(unassigned_type::operator()(),
			    unassigned_type::operator()());
    }
  };

  struct span_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    size_type first;
    size_type last;
    
    span_type() : first(0), last(0) {}
    span_type(const size_type& __first, const size_type& __last)
      : first(__first), last(__last) { }

    bool empty() const { return first == last; }
    size_type size() const { return last - first; }
    
    friend
    bool operator==(const span_type& x, const span_type& y)
    {
      return x.first == y.first && x.last == y.last;
    }
    
    friend
    size_t hash_value(span_type const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.first, x.last);
    }

    friend
    bool disjoint(const span_type& x, const span_type& y)
    {
      return x.last <= y.first || y.last <= x.first;
    }

    friend
    bool adjacent(const span_type& x, const span_type& y)
    {
      return x.last == y.first || y.last == x.first;
    }
    
  };

  struct span_pair_type
  {
    typedef span_type::size_type      size_type;
    typedef span_type::difference_type difference_type;

    span_type source;
    span_type target;
    
    span_pair_type() : source(), target() {}
    span_pair_type(const span_type& __source, const span_type& __target)
      : source(__source), target(__target) {}
    span_pair_type(size_type s1, size_type s2, size_type t1, size_type t2)
      : source(s1, s2), target(t1, t2) {}
    
    bool empty() const { return source.empty() && target.empty(); }
    size_type size() const { return source.size() + target.size(); }
    
    friend
    bool operator==(const span_pair_type& x, const span_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    size_t hash_value(span_pair_type const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x);
    }
    
    friend
    bool disjoint(const span_pair_type& x, const span_pair_type& y)
    {
      return disjoint(x.source, y.source) && disjoint(x.target, y.target);
    }

    friend
    bool adjacent(const span_pair_type& x, const span_pair_type& y)
    {
      return adjacent(x.source, y.source) && adjacent(x.target, y.target);
    }
    
  };

  // How to differentiate base or generative?
  struct rule_type
  {
    span_pair_type span;
    span_pair_type left;
    span_pair_type right;
    itg_type       itg;

    id_type        word_pair;

    rule_type()
      : span(), left(), right(), itg(), word_pair(id_type(-1)) {}
    rule_type(const span_pair_type& __span, const itg_type& __itg)
      : span(__span), left(), right(), itg(__itg), word_pair(id_type(-1)) {}
    rule_type(const span_pair_type& __span, const span_pair_type& __left, const span_pair_type& __right, const itg_type& __itg)
      : span(__span), left(__left), right(__right), itg(__itg), word_pair(id_type(-1)) {}
    
    bool is_terminal() const { return left.empty() && right.empty(); }
    bool is_straight() const { return ! is_terminal() && left.target.last  == right.target.first; }
    bool is_inverted() const { return ! is_terminal() && left.target.first == right.target.last; }
  };
};

struct PYPTerminal
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  typedef PYP::id_type         id_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef PYP::word_pair_type   word_pair_type;

  typedef utils::pyp_parameter parameter_type;

  // do we use pialign style indexing...?
  typedef utils::indexed_set<word_pair_type, boost::hash<word_pair_type>, std::equal_to<word_pair_type>,
			     std::allocator<word_pair_type> > word_pair_set_type;

  typedef utils::restaurant_sync<> table_type;
  
  PYPTerminal(const parameter_type& parameter,
	      const double __p0_source,
	      const double __p0_target,
	      const double __p0_epsilon)
    : table(parameter),
      p0_source(__p0_source),
      p0_target(__p0_target),
      p0_epsilon(__p0_epsilon)
  { }

  id_type word_pair_id(const word_type& source, const word_type& target)
  {
    word_pair_set_type::iterator iter = word_pairs.insert(word_pair_type(source, target)).first;
    return iter - word_pairs.begin();
  }
  
  std::pair<id_type, bool> find_word_pair(const word_type& source, const word_type& target) const
  {
    word_pair_set_type::const_iterator iter = word_pairs.find(word_pair_type(source, target));
    return std::make_pair(iter - word_pairs.begin(), iter != word_pairs.end());
  }

  template <typename Iterator, typename Sampler>
  void increment(Iterator first, Iterator last, Sampler& sampler, const double temperature)
  {
    for (/**/; first != last; ++ first)
      table.increment(*first, prior(*first), sampler, temperature);
  }

  template <typename Iterator, typename Sampler>
  void decrement(Iterator first, Iterator last, Sampler& sampler)
  {
    table.decrement(first, last, sampler);
  }
  
  double prob(const id_type id) const
  {
    return table.prob(id, prior(id));
  }
  
  double prob(const word_type& source, const word_type& target) const
  {
    const std::pair<id_type, bool> result = find_word_pair(source, target);
    
    if (result.second)
      return table.prob(result.first, prior(source, target));
    else
      return table.prob(prior(source, target));
  }
  
  double prior(const id_type id) const
  {
    return prior(word_pairs[id].source, word_pairs[id].target);
  }
  
  double prior(const word_type& source, const word_type& target) const
  {
    return ((source != vocab_type::EPSILON ? (1.0 - p0_epsilon) * p0_source : p0_epsilon)
	    * (target != vocab_type::EPSILON ? (1.0 - p0_epsilon) * p0_target : p0_epsilon));
  }
  
  double log_likelihood() const
  {
    double ll = table.log_likelihood();
    
    size_type counts = 0;
    size_type counts_epsilon_target = 0;
    size_type counts_source_epsilon = 0;
    
    for (PYP::id_type id = 0; id != table.size(); ++ id)
      if (! table[id].empty()) {
	const word_type& source = word_pairs[id].source;
	const word_type& target = word_pairs[id].target;
	
	counts                += (source != vocab_type::EPSILON && target != vocab_type::EPSILON) * table[id].size_table();
	counts_epsilon_target += (source == vocab_type::EPSILON) * table[id].size_table();
	counts_source_epsilon += (target == vocab_type::EPSILON) * table[id].size_table();
      }
    
    ll += std::log((1.0 - p0_epsilon) * (1.0 - p0_epsilon) * p0_source * p0_target) * counts;
    ll += std::log((1.0 - p0_epsilon) * p0_source * p0_epsilon) * counts_source_epsilon;
    ll += std::log((1.0 - p0_epsilon) * p0_target * p0_epsilon) * counts_epsilon_target;
    
    return ll;
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
  
  word_pair_set_type word_pairs;
  
  table_type table;
  
  double     p0_source;
  double     p0_target;
  double     p0_epsilon;
};


struct PYPRule
{
  // rule selection probabilities
  // terminal, stright, inversion
  
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef PYP::rule_type rule_type;
  
  typedef PYP::itg_type itg_type;

  typedef uint32_t id_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  typedef utils::restaurant_sync<> table_type;

  PYPRule(const double& __p0_terminal, const parameter_type& parameter)
    : p0_terminal(__p0_terminal), p0((1.0 - __p0_terminal) * 0.5), counts0_terminal(0), counts0(0), table(parameter)
  {
    table.reserve(3);
    table.resize(3);
  }
  
  template <typename Sampler>
  void increment(const rule_type& rule, Sampler& sampler, const double temperature=1.0)
  {
    const itg_type itg = (rule.is_terminal() ? PYP::TERMINAL : (rule.is_straight() ? PYP::STRAIGHT : PYP::INVERTED));
    
    table.increment(itg, itg == PYP::TERMINAL ? p0_terminal : p0, sampler, temperature);
  }

  
  template <typename Sampler>
  void decrement(const rule_type& rule, Sampler& sampler)
  {
    const itg_type itg = (rule.is_terminal() ? PYP::TERMINAL : (rule.is_straight() ? PYP::STRAIGHT : PYP::INVERTED));
    
    table.decrement(itg, sampler);
  }
  
  double prob(const rule_type& rule) const
  {
    const itg_type itg = (rule.is_terminal() ? PYP::TERMINAL : (rule.is_straight() ? PYP::STRAIGHT : PYP::INVERTED));
    
    return table.prob(itg, itg == PYP::TERMINAL ? p0_terminal : p0);
  }
  
  double prob_terminal() const
  {
    return table.prob(PYP::TERMINAL, p0_terminal);
  }

  double prob_straight() const
  {
    return table.prob(PYP::STRAIGHT, p0);
  }
  
  double prob_inverted() const
  {
    return table.prob(PYP::INVERTED, p0);
  }
  
  double log_likelihood() const
  {
    return table.log_likelihood() + std::log(p0_terminal) * counts0_terminal + std::log(p0) * counts0;
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    table.sample_parameters(sampler, num_loop, num_iterations);
    
    counts0_terminal = table[PYP::TERMINAL].size_table();
    counts0          = table.size_table() - table[PYP::TERMINAL].size_table();
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    table.slice_sample_parameters(sampler, num_loop, num_iterations);

    counts0_terminal = table[PYP::TERMINAL].size_table();
    counts0          = table.size_table() - table[PYP::TERMINAL].size_table();
  }
  
  double     p0_terminal;
  double     p0;
  size_type  counts0_terminal;
  size_type  counts0;
  table_type table;
};

struct PYPTree
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::rule_type rule_type;
  typedef PYP::itg_type  itg_type;

  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  typedef utils::restaurant_sync<> table_type;
  
  
};

struct PYPITG
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::rule_type rule_type;
  typedef PYP::itg_type  itg_type;

  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

  PYPITG(const PYPRule&   __rule,
	 const PYPTerminal& __terminal)
    : rule(__rule), terminal(__terminal) {}


  double log_likelihood() const
  {
    return rule.log_likelihood() + terminal.log_likelihood();
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    rule.sample_parameters(sampler, num_loop, num_iterations);

    terminal.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    rule.slice_sample_parameters(sampler, num_loop, num_iterations);
    
    terminal.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  PYPRule    rule;
  PYPTerminal terminal;
};

struct PYPGraph
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::id_type id_type;

  typedef PYP::span_type      span_type;
  typedef PYP::span_pair_type span_pair_type;
  
  typedef PYP::rule_type rule_type;
  
  typedef std::vector<rule_type, std::allocator<rule_type> > derivation_type;

  typedef PYPTerminal::logprob_type logprob_type;
  typedef PYPTerminal::prob_type    prob_type;
  
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  typedef utils::bichart<logprob_type, std::allocator<logprob_type> > chart_type;
  
  struct edge_type
  {
    rule_type    rule;
    logprob_type prob;
    
    edge_type() : rule(), prob() {}
    edge_type(const rule_type& __rule, const logprob_type& __prob)
      : rule(__rule), prob(__prob) {}
  };
  
  typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
  typedef utils::bichart<edge_set_type, std::allocator<edge_set_type> > edge_chart_type;
  
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  typedef std::vector<span_pair_set_type, std::allocator<span_pair_set_type> > agenda_type;
  
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > stack_type;

  typedef std::pair<span_pair_type, span_pair_type> span_pairs_type;
  
  struct span_pairs_unassigned
  {
    span_pairs_type operator()() const
    {
      return span_pairs_type(span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)),
			     span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)));
    }
  };
  
  typedef utils::compact_set<span_pairs_type,
			     span_pairs_unassigned, span_pairs_unassigned,
			     utils::hashmurmur3<size_t>, std::equal_to<span_pairs_type>,
			     std::allocator<span_pairs_type> > span_pairs_unique_type;

  typedef utils::chart<logprob_type, std::allocator<logprob_type> > chart_mono_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > alpha_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > beta_type;

  typedef utils::vector2<logprob_type, std::allocator<logprob_type> > matrix_type;

  typedef std::pair<logprob_type, span_pair_type> score_span_pair_type;
  typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;

  
  void initialize(const sentence_type& source,
		  const sentence_type& target,
		  const PYPITG& model)
  {
    //std::cerr << "initialize" << std::endl;
    matrix.clear();
    matrix.reserve(source.size() + 1, target.size() + 1);
    matrix.resize(source.size() + 1, target.size() + 1);
    
    logprob_term = model.rule.prob_terminal();
    logprob_str  = model.rule.prob_straight();
    logprob_inv  = model.rule.prob_inverted();

    for (size_type src = 0; src <= source.size(); ++ src)
      for (size_type trg = (src == 0); trg <= target.size(); ++ trg)
	matrix(src, trg) = model.terminal.prob(src == 0 ? vocab_type::EPSILON : source[src - 1],
					       trg == 0 ? vocab_type::EPSILON : target[trg - 1]);
  }
  
  void forward_backward(const sentence_type& sentence, const chart_mono_type& chart, alpha_type& alpha, beta_type& beta)
  {
    const logprob_type logprob_rule = std::max(logprob_str, logprob_inv);

    // forward...
    alpha[0] = 1.0;
    for (size_type last = 1; last <= sentence.size(); ++ last)
      for (size_type first = 0; first != last; ++ first)
	alpha[last] = std::max(alpha[last], alpha[first] * chart(first, last) * (first != 0 ? logprob_rule : logprob_type(1.0)));
    
    // backward...
    beta[sentence.size()] = 1.0;
    for (difference_type first = sentence.size() - 1; first >= 0; -- first)
      for (size_type last = first + 1; last <= sentence.size(); ++ last)
	beta[first] = std::max(beta[first], chart(first, last) * beta[last] * (last != sentence.size() ? logprob_rule : logprob_type(1.0)));
  }
  
  // sort by less so that we can pop from a greater item
  struct heap_compare
  {
    bool operator()(const score_span_pair_type& x, const score_span_pair_type& y) const
    {
      return x.first < y.first;
    }
  };
  
  // forward filtering
  std::pair<logprob_type, bool> forward(const sentence_type& source,
					const sentence_type& target,
					const logprob_type beam)
  {
    //std::cerr << "forward" << std::endl;

    
    // initialize...
    chart.clear();
    chart.reserve(source.size() + 1, target.size() + 1);
    chart.resize(source.size() + 1, target.size() + 1);
    
    edges.clear();
    edges.reserve(source.size() + 1, target.size() + 1);
    edges.resize(source.size() + 1, target.size() + 1);
    
    agenda.clear();
    agenda.reserve(source.size() + target.size() + 1);
    agenda.resize(source.size() + target.size() + 1);

    chart_source.clear();
    chart_source.reserve(source.size() + 1);
    chart_source.resize(source.size() + 1);

    chart_target.clear();
    chart_target.reserve(target.size() + 1);
    chart_target.resize(target.size() + 1);

    alpha_source.clear();
    alpha_target.clear();
    beta_source.clear();
    beta_target.clear();
    
    alpha_source.resize(source.size() + 1);
    alpha_target.resize(target.size() + 1);
    beta_source.resize(source.size() + 1);
    beta_target.resize(target.size() + 1);
    
    for (size_type src = 0; src <= source.size(); ++ src)
      for (size_type trg = 0; trg <= target.size(); ++ trg) {
	const size_type source_first = src;
	const size_type target_first = trg;

	if (src < source.size() && trg < target.size()) {
	  const logprob_type prob_term = matrix(source_first + 1, target_first + 1) * logprob_term;
	  
	  const size_type source_last = source_first + 1;
	  const size_type target_last = target_first + 1;
	    
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	    
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::TERMINAL), prob_term));
	    
	  chart(source_first, source_last, target_first, target_last) = prob_term;
	    
	  chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), prob_term);
	  chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), prob_term);
	    
	  agenda[span_pair.size()].push_back(span_pair);
	}
	
	// epsilon-to-many
	if (src < source.size()) {
	  const size_type source_last = source_first + 1;
	  const size_type target_last = target_first;
	  
	  const logprob_type prob_term = logprob_term * matrix(source_last, 0);
	  
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::TERMINAL), prob_term));
	    
	  chart(source_first, source_last, target_first, target_last) = prob_term;
	    
	  chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), prob_term);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
	
	// epsilon-to-many
	if (trg < target.size()) {
	  const size_type target_last = target_first + 1;
	  const size_type source_last = source_first;
	    
	  const logprob_type prob_term = logprob_term * matrix(0, target_last);
	  
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::TERMINAL), prob_term));
	  
	  chart(source_first, source_last, target_first, target_last) = prob_term;
	    
	  chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), prob_term);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
      }

    // forward-backward to compute estiamtes...
    forward_backward(source, chart_source, alpha_source, beta_source);
    forward_backward(target, chart_target, alpha_target, beta_target);
      
    
    // start parsing...
    
    span_pairs_unique_type spans_unique;
     
    // traverse agenda, smallest first...
    const size_type length_max = source.size() + target.size();

    for (size_type length = 1; length != length_max; ++ length) 
      if (! agenda[length].empty()) {
	span_pair_set_type& spans = agenda[length];

	// construct heap...
	heap.clear();
	heap.reserve(spans.size());
	span_pair_set_type::const_iterator siter_end = spans.end();
	for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)  {
	  const logprob_type score = (chart(siter->source.first, siter->source.last, siter->target.first, siter->target.last)
				      * std::min(alpha_source[siter->source.first] * beta_source[siter->source.last],
						 alpha_target[siter->target.first] * beta_target[siter->target.last]));
	  
	  heap.push_back(score_span_pair_type(score, *siter));
	  std::push_heap(heap.begin(), heap.end(), heap_compare());
	}

	heap_type::iterator hiter_begin = heap.begin();
	heap_type::iterator hiter       = heap.end();
	heap_type::iterator hiter_end   = heap.end();
	  
	const logprob_type logprob_threshold = hiter_begin->first * beam;
	
	for (/**/; hiter_begin != hiter && hiter_begin->first > logprob_threshold; -- hiter)
	  std::pop_heap(hiter_begin, hiter, heap_compare());
	
	// we will process from hiter to hiter_end...
	spans_unique.clear();
	
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const span_pair_type& span_pair = iter->second;
	  
	  // we borrow the notation...
	  
	  const difference_type l = length;
	  const difference_type s = span_pair.source.first;
	  const difference_type t = span_pair.source.last;
	  const difference_type u = span_pair.target.first;
	  const difference_type v = span_pair.target.last;
	  
	  const difference_type T = source.size();
	  const difference_type V = target.size();

	  // remember, we have processed only upto length l. thus, do not try to combine with spans greather than l!
	  // also, keep used pair of spans in order to remove duplicates.
	  
	  for (difference_type S = utils::bithack::max(s - l, difference_type(0)); S <= s; ++ S) {
	    const difference_type L = l - (s - S);
	    
	    // straight
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == s); ++ U) {
	      // parent span: StUv
	      // span1: SsUu
	      // span2: stuv

	      if (edges(S, s, U, u).empty()) continue;
	      
	      const span_pair_type  span1(S, s, U, u);
	      const span_pair_type& span2(span_pair);

	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_str
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));

	      const span_pair_type span_head(S, t, U, v);
	      
	      chart(S, t, U, v) += logprob;
	      edges(S, t, U, v).push_back(edge_type(rule_type(span_head, span1, span2, PYP::STRAIGHT), logprob));
	      
	      if (edges(S, t, U, v).size() == 1)
		agenda[span_head.size()].push_back(span_head);
	    }
	    
	    // inversion
	    for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: StuU
	      // span1: SsvU
	      // span2: stuv
	      
	      if (edges(S, s, v, U).empty()) continue;
	      
	      const span_pair_type  span1(S, s, v, U);
	      const span_pair_type& span2(span_pair);

	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_inv
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
	      
	      const span_pair_type span_head(S, t, u, U);

	      chart(S, t, u, U) += logprob;
	      edges(S, t, u, U).push_back(edge_type(rule_type(span_head, span1, span2, PYP::INVERTED), logprob));

	      if (edges(S, t, u, U).size() == 1)
		agenda[span_head.size()].push_back(span_head);
	    }
	  }
	  
	  for (difference_type S = t; S <= utils::bithack::min(t + l, T); ++ S) {
	    const difference_type L = l - (S - t);
	    
	    // inversion
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == t); ++ U) {
	      // parent span: sSUv
	      // span1: stuv
	      // span2: tSUu
	      
	      if (edges(t, S, U, u).empty()) continue;
	      
	      const span_pair_type& span1(span_pair);
	      const span_pair_type  span2(t, S, U, u);

	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_inv
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));

	      const span_pair_type span_head(s, S, U, v);
	      
	      chart(s, S, U, v) += logprob;
	      edges(s, S, U, v).push_back(edge_type(rule_type(span_head, span1, span2, PYP::INVERTED), logprob));

	      if (edges(s, S, U, v).size() == 1)
		agenda[span_head.size()].push_back(span_head);
	    }
	    
	    // straight
	    for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: sSuU
	      // span1: stuv
	      // span2: tSvU
	      
	      if (edges(t, S, v, U).empty()) continue;
	      
	      const span_pair_type& span1(span_pair);
	      const span_pair_type  span2(t, S, v, U);
	      
	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_str
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
	      
	      const span_pair_type span_head(s, S, u, U);
	      
	      chart(s, S, u, U) += logprob;
	      edges(s, S, u, U).push_back(edge_type(rule_type(span_head, span1, span2, PYP::STRAIGHT), logprob));
	      
	      if (edges(s, S, u, U).size() == 1)
		agenda[span_head.size()].push_back(span_head);
	    }
	  }
	}
      }
    
    return std::make_pair(chart(0, source.size(), 0, target.size()), ! edges(0, source.size(), 0, target.size()).empty());
  }
  
  // backward sampling
  template <typename Sampler>
  logprob_type backward(const sentence_type& source,
			const sentence_type& target,
			derivation_type& derivation,
			Sampler& sampler,
			const double temperature)
  {
    //std::cerr << "backward" << std::endl;

    derivation.clear();
    logprob_type prob_derivation = cicada::semiring::traits<logprob_type>::one(); 
    
    // top-down sampling of the hypergraph...
    
    stack.clear();
    stack.push_back(span_pair_type(0, source.size(), 0, target.size()));
    
    while (! stack.empty()) {
      const span_pair_type span_pair = stack.back();
      stack.pop_back();
      
      const edge_set_type& edges_span = edges(span_pair.source.first, span_pair.source.last, span_pair.target.first, span_pair.target.last);
      
      if (edges_span.empty()) {
	// no derivation???
	derivation.clear();
	return logprob_type();
      }
      
      logprob_type logsum;
      edge_set_type::const_iterator eiter_end = edges_span.end();
      for (edge_set_type::const_iterator eiter = edges_span.begin(); eiter != eiter_end; ++ eiter)
	logsum += eiter->prob;
      
      probs.clear();
      for (edge_set_type::const_iterator eiter = edges_span.begin(); eiter != eiter_end; ++ eiter)
	probs.push_back(eiter->prob / logsum);
      
      const size_type pos_sampled = sampler.draw(probs.begin(), probs.end(), temperature) - probs.begin();
      
      // we will push in a right first manner, so that when popped, we will derive left-to-right traversal.
      
      const rule_type& rule = edges_span[pos_sampled].rule;
      logprob_type prob = edges_span[pos_sampled].prob;
      
      if (! rule.right.empty()) {
	prob /= chart(rule.right.source.first, rule.right.source.last, rule.right.target.first, rule.right.target.last);
	
	stack.push_back(rule.right);
      }
      
      if (! rule.left.empty()) {
	prob /= chart(rule.left.source.first, rule.left.source.last, rule.left.target.first, rule.left.target.last);
	
	stack.push_back(rule.left);
      }
      
      prob_derivation *= prob;
      
      derivation.push_back(rule);
    }

    return prob_derivation;
  }

  logprob_type logprob_term;
  logprob_type logprob_str;
  logprob_type logprob_inv;
  matrix_type  matrix;

  chart_type      chart;
  edge_chart_type edges;
  
  agenda_type     agenda;
  stack_type      stack;
  heap_type       heap;
  
  chart_mono_type chart_source;
  chart_mono_type chart_target;
  
  alpha_type alpha_source;
  alpha_type alpha_target;
  beta_type  beta_source;
  beta_type  beta_target;
  
  prob_set_type probs;
};

struct PYPViterbi
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::span_type      span_type;
  typedef PYP::span_pair_type span_pair_type;
  
  typedef PYP::rule_type rule_type;
  
  typedef std::vector<rule_type, std::allocator<rule_type> > derivation_type;

  typedef PYPTerminal::logprob_type logprob_type;
  typedef PYPTerminal::prob_type    prob_type;
  
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  typedef utils::bichart<logprob_type, std::allocator<logprob_type> > chart_type;
  
  struct edge_type
  {
    rule_type    rule;
    logprob_type prob;
    
    edge_type() : rule(), prob() {}
    edge_type(const rule_type& __rule, const logprob_type& __prob)
      : rule(__rule), prob(__prob) {}
  };
  
  typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
  typedef utils::bichart<edge_set_type, std::allocator<edge_set_type> > edge_chart_type;
  
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  typedef std::vector<span_pair_set_type, std::allocator<span_pair_set_type> > agenda_type;
  
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > stack_type;

  typedef std::pair<span_pair_type, span_pair_type> span_pairs_type;
  
  struct span_pairs_unassigned
  {
    span_pairs_type operator()() const
    {
      return span_pairs_type(span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)),
			     span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)));
    }
  };
  
  typedef utils::compact_set<span_pairs_type,
			     span_pairs_unassigned, span_pairs_unassigned,
			     utils::hashmurmur3<size_t>, std::equal_to<span_pairs_type>,
			     std::allocator<span_pairs_type> > span_pairs_unique_type;
  
  typedef utils::chart<logprob_type, std::allocator<logprob_type> > chart_mono_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > alpha_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > beta_type;

  typedef utils::vector2<logprob_type, std::allocator<logprob_type> > matrix_type;

  typedef std::pair<logprob_type, span_pair_type> score_span_pair_type;
  typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;

  
  void initialize(const sentence_type& source,
		  const sentence_type& target,
		  const PYPITG& model)
  {
    //std::cerr << "initialize" << std::endl;
    matrix.clear();
    matrix.reserve(source.size() + 1, target.size() + 1);
    matrix.resize(source.size() + 1, target.size() + 1);
        
    logprob_term = model.rule.prob_terminal();
    logprob_str  = model.rule.prob_straight();
    logprob_inv  = model.rule.prob_inverted();

    for (size_type src = 0; src <= source.size(); ++ src)
      for (size_type trg = (src == 0); trg <= target.size(); ++ trg)
	matrix(src, trg) = model.terminal.prob(src == 0 ? vocab_type::EPSILON : source[src - 1],
					       trg == 0 ? vocab_type::EPSILON : target[trg - 1]);
    
  }
  
  void forward_backward(const sentence_type& sentence, const chart_mono_type& chart, alpha_type& alpha, beta_type& beta)
  {
    const logprob_type logprob_rule = std::max(logprob_str, logprob_inv);

    // forward...
    alpha[0] = 1.0;
    for (size_type last = 1; last <= sentence.size(); ++ last)
      for (size_type first = 0; first != last; ++ first)
	alpha[last] = std::max(alpha[last], alpha[first] * chart(first, last) * (first != 0 ? logprob_rule : logprob_type(1.0)));
    
    // backward...
    beta[sentence.size()] = 1.0;
    for (difference_type first = sentence.size() - 1; first >= 0; -- first)
      for (size_type last = first + 1; last <= sentence.size(); ++ last)
	beta[first] = std::max(beta[first], chart(first, last) * beta[last] * (last != sentence.size() ? logprob_rule : logprob_type(1.0)));
  }
  
  // sort by less so that we can pop from a greater item
  struct heap_compare
  {
    bool operator()(const score_span_pair_type& x, const score_span_pair_type& y) const
    {
      return x.first < y.first;
    }
  };
  
  // forward filtering
  std::pair<logprob_type, bool> forward(const sentence_type& source,
					const sentence_type& target,
					const logprob_type beam)
  {
    //std::cerr << "forward" << std::endl;

    
    // initialize...
    chart.clear();
    chart.reserve(source.size() + 1, target.size() + 1);
    chart.resize(source.size() + 1, target.size() + 1);
    
    edges.clear();
    edges.reserve(source.size() + 1, target.size() + 1);
    edges.resize(source.size() + 1, target.size() + 1);
    
    agenda.clear();
    agenda.reserve(source.size() + target.size() + 1);
    agenda.resize(source.size() + target.size() + 1);

    chart_source.clear();
    chart_source.reserve(source.size() + 1);
    chart_source.resize(source.size() + 1);

    chart_target.clear();
    chart_target.reserve(target.size() + 1);
    chart_target.resize(target.size() + 1);

    alpha_source.clear();
    alpha_target.clear();
    beta_source.clear();
    beta_target.clear();
    
    alpha_source.resize(source.size() + 1);
    alpha_target.resize(target.size() + 1);
    beta_source.resize(source.size() + 1);
    beta_target.resize(target.size() + 1);
    
    for (size_type src = 0; src <= source.size(); ++ src)
      for (size_type trg = 0; trg <= target.size(); ++ trg) {
	const size_type source_first = src;
	const size_type target_first = trg;

	if (src < source.size() && trg < target.size()) {
	  const logprob_type prob_term = matrix(source_first + 1, target_first + 1) * logprob_term;
	  
	  const size_type source_last = source_first + 1;
	  const size_type target_last = target_first + 1;
	  
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::TERMINAL), prob_term));
	  
	  chart(source_first, source_last, target_first, target_last) = prob_term;
	  
	  chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), prob_term);
	  chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), prob_term);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
	
	// epsilon-to-many
	if (src < source.size()) {
	  const size_type source_last = source_first + 1;
	  const size_type target_last = target_first;
	  
	  const logprob_type prob_term = logprob_term * matrix(source_last, 0);
	  
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	    
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::TERMINAL), prob_term));
	    
	  chart(source_first, source_last, target_first, target_last) = prob_term;
	    
	  chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), prob_term);
	    
	  agenda[span_pair.size()].push_back(span_pair);
	}
	
	// epsilon-to-many
	if (trg < target.size()) {
	  const size_type target_last = target_first + 1;
	  const size_type source_last = source_first;
	    
	  const logprob_type prob_term = logprob_term * matrix(0, target_last);
	  
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::TERMINAL), prob_term));
	  
	  chart(source_first, source_last, target_first, target_last) = prob_term;
	  
	  chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), prob_term);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
      }

    // forward-backward to compute estiamtes...
    forward_backward(source, chart_source, alpha_source, beta_source);
    forward_backward(target, chart_target, alpha_target, beta_target);
      
    
    // start parsing...
    
    span_pairs_unique_type spans_unique;
    
    // traverse agenda, smallest first...
    const size_type length_max = source.size() + target.size();

    for (size_type length = 1; length != length_max; ++ length) 
      if (! agenda[length].empty()) {
	span_pair_set_type& spans = agenda[length];

	// construct heap...
	heap.clear();
	heap.reserve(spans.size());
	span_pair_set_type::const_iterator siter_end = spans.end();
	for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)  {
	  const logprob_type score = (chart(siter->source.first, siter->source.last, siter->target.first, siter->target.last)
				      * std::min(alpha_source[siter->source.first] * beta_source[siter->source.last],
						 alpha_target[siter->target.first] * beta_target[siter->target.last]));
	  
	  heap.push_back(score_span_pair_type(score, *siter));
	  std::push_heap(heap.begin(), heap.end(), heap_compare());
	}

	heap_type::iterator hiter_begin = heap.begin();
	heap_type::iterator hiter       = heap.end();
	heap_type::iterator hiter_end   = heap.end();
	
	const logprob_type logprob_threshold = hiter_begin->first * beam;
	for (/**/; hiter_begin != hiter && hiter_begin->first > logprob_threshold; -- hiter)
	  std::pop_heap(hiter_begin, hiter, heap_compare());
	
	// we will process from hiter to hiter_end...
	spans_unique.clear();
	
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const span_pair_type& span_pair = iter->second;
	  
	  // we borrow the notation...
	  
	  const difference_type l = length;
	  const difference_type s = span_pair.source.first;
	  const difference_type t = span_pair.source.last;
	  const difference_type u = span_pair.target.first;
	  const difference_type v = span_pair.target.last;
	  
	  const difference_type T = source.size();
	  const difference_type V = target.size();
	  
	  // remember, we have processed only upto length l. thus, do not try to combine with spans greather than l!
	  // also, keep used pair of spans in order to remove duplicates.
	  
	  for (difference_type S = utils::bithack::max(s - l, difference_type(0)); S <= s; ++ S) {
	    const difference_type L = l - (s - S);
	    
	    // straight
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == s); ++ U) {
	      // parent span: StUv
	      // span1: SsUu
	      // span2: stuv

	      if (edges(S, s, U, u).empty()) continue;
	      
	      const span_pair_type  span1(S, s, U, u);
	      const span_pair_type& span2(span_pair);

	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_str
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
	      
	      if (logprob > chart(S, t, U, v)) {
		const span_pair_type span_head(S, t, U, v);
		
		chart(S, t, U, v) = logprob;
		
		if (edges(S, t, U, v).empty()) {
		  edges(S, t, U, v).push_back(edge_type(rule_type(span_head, span1, span2, PYP::STRAIGHT), logprob));
		  
		  agenda[span_head.size()].push_back(span_head);
		} else
		  edges(S, t, U, v).back() = edge_type(rule_type(span_head, span1, span2, PYP::STRAIGHT), logprob);
	      }
	    }
	    
	    // inversion
	    for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: StuU
	      // span1: SsvU
	      // span2: stuv
	      
	      if (edges(S, s, v, U).empty()) continue;
	      
	      const span_pair_type  span1(S, s, v, U);
	      const span_pair_type& span2(span_pair);

	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_inv
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));

	      if (logprob > chart(S, t, u, U)) {
		const span_pair_type span_head(S, t, u, U);
		
		chart(S, t, u, U) = logprob;
		
		if (edges(S, t, u, U).empty()) {
		  edges(S, t, u, U).push_back(edge_type(rule_type(span_head, span1, span2, PYP::INVERTED), logprob));
		  
		  agenda[span_head.size()].push_back(span_head);
		} else
		  edges(S, t, u, U).back() = edge_type(rule_type(span_head, span1, span2, PYP::INVERTED), logprob);
	      }
	    }
	  }
	  
	  for (difference_type S = t; S <= utils::bithack::min(t + l, T); ++ S) {
	    const difference_type L = l - (S - t);
	    
	    // inversion
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == t); ++ U) {
	      // parent span: sSUv
	      // span1: stuv
	      // span2: tSUu
	      
	      if (edges(t, S, U, u).empty()) continue;
	      
	      const span_pair_type& span1(span_pair);
	      const span_pair_type  span2(t, S, U, u);

	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_inv
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));

	      if (logprob > chart(s, S, U, v)) {
		const span_pair_type span_head(s, S, U, v);
		
		chart(s, S, U, v) = logprob;

		if (edges(s, S, U, v).empty()) {
		  edges(s, S, U, v).push_back(edge_type(rule_type(span_head, span1, span2, PYP::INVERTED), logprob));
		  
		  agenda[span_head.size()].push_back(span_head);
		} else
		  edges(s, S, U, v).back() = edge_type(rule_type(span_head, span1, span2, PYP::INVERTED), logprob);
	      }
	    }
	    
	    // straight
	    for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: sSuU
	      // span1: stuv
	      // span2: tSvU
	      
	      if (edges(t, S, v, U).empty()) continue;
	      
	      const span_pair_type& span1(span_pair);
	      const span_pair_type  span2(t, S, v, U);
	      
	      if (! spans_unique.insert(std::make_pair(span1, span2)).second) continue;
	      
	      const logprob_type logprob = (logprob_str
					    * chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    * chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
	      
	      if (logprob > chart(s, S, u, U)) {
		const span_pair_type span_head(s, S, u, U);
		
		chart(s, S, u, U) = logprob;
		
		if (edges(s, S, u, U).empty()) {
		  edges(s, S, u, U).push_back(edge_type(rule_type(span_head, span1, span2, PYP::STRAIGHT), logprob));
		  
		  agenda[span_head.size()].push_back(span_head);
		} else
		  edges(s, S, u, U).back() = edge_type(rule_type(span_head, span1, span2, PYP::STRAIGHT), logprob);
	      }
	    }
	  }
	}
      }
    
    return std::make_pair(chart(0, source.size(), 0, target.size()), ! edges(0, source.size(), 0, target.size()).empty());
  }
  
  // backward Viterbi
  logprob_type backward(const sentence_type& source,
			const sentence_type& target,
			derivation_type& derivation)
  {
    derivation.clear();
    logprob_type prob_derivation = cicada::semiring::traits<logprob_type>::one(); 
    
    stack.clear();
    stack.push_back(span_pair_type(0, source.size(), 0, target.size()));
    
    while (! stack.empty()) {
      const span_pair_type span_pair = stack.back();
      stack.pop_back();
      
      const edge_set_type& edges_span = edges(span_pair.source.first, span_pair.source.last, span_pair.target.first, span_pair.target.last);
      
      if (edges_span.empty()) {
	// no derivation???
	derivation.clear();
	return logprob_type();
      }
      
      // we will push in a right first manner, so that when popped, we will derive left-to-right traversal.
      
      const rule_type& rule = edges_span.front().rule;
      logprob_type prob = edges_span.front().prob;
      
      if (! rule.right.empty()) {
	prob /= chart(rule.right.source.first, rule.right.source.last, rule.right.target.first, rule.right.target.last);
	
	stack.push_back(rule.right);
      }
      
      if (! rule.left.empty()) {
	prob /= chart(rule.left.source.first, rule.left.source.last, rule.left.target.first, rule.left.target.last);
	
	stack.push_back(rule.left);
      }
      
      prob_derivation *= prob;
      
      derivation.push_back(rule);
    }
    
    return prob_derivation;
  }

  logprob_type logprob_term;
  logprob_type logprob_str;
  logprob_type logprob_inv;
  matrix_type  matrix;

  chart_type      chart;
  edge_chart_type edges;
  
  agenda_type     agenda;
  stack_type      stack;
  heap_type       heap;
  
  chart_mono_type chart_source;
  chart_mono_type chart_target;
  
  alpha_type alpha_source;
  alpha_type alpha_target;
  beta_type  beta_source;
  beta_type  beta_target;
  
  prob_set_type probs;
};


typedef boost::filesystem::path path_type;
typedef utils::sampler<boost::mt19937> sampler_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;

typedef PYP::size_type size_type;
typedef PYPGraph::derivation_type derivation_type;
typedef std::vector<derivation_type, std::allocator<derivation_type> > derivation_set_type;
typedef std::vector<size_type, std::allocator<size_type> > position_set_type;

struct Time
{
  double initialize;
  double increment;
  double decrement;
  double forward;
  double backward;
  
  Time() : initialize(0), increment(0), decrement(0), forward(0), backward(0) {}

  Time& operator+=(const Time& x)
  {
    initialize += x.initialize;
    increment  += x.increment;
    decrement  += x.decrement;
    forward    += x.forward;
    backward   += x.backward;

    return *this;
  }
};

struct Counter
{
  Counter() : counter(0) {}
  
  void increment()
  {
    utils::atomicop::fetch_and_add(counter, size_type(1));
  }
  
  void wait(size_type target)
  {
    for (;;) {
      for (int i = 0; i < 64; ++ i) {
	if (counter == target)
	  return;
	else
	  boost::thread::yield();
      }
      
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    }
  }

  void clear() { counter = 0; }
  
  volatile size_type counter;
};
typedef Counter counter_type;

struct Task
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;  
  typedef PYP::id_type         id_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;
  
  typedef PYPITG::logprob_type logprob_type;
  typedef PYPITG::prob_type    prob_type;

  typedef std::vector<size_type, std::allocator<size_type> > count_set_type;
  typedef  std::vector<id_type, std::allocator<id_type> > id_set_type;

  count_set_type counts;
  id_set_type    ids;

  Task(queue_type& __mapper,
       queue_type& __reducer_derivation,
       counter_type& __reducer,
       const sentence_set_type& __sources,
       const sentence_set_type& __targets,
       derivation_set_type& __derivations,
       PYPITG& __model,
       const sampler_type& __sampler,
       const logprob_type& __beam)
    : mapper(__mapper),
      reducer_derivation(__reducer_derivation),
      reducer(__reducer),
      sources(__sources),
      targets(__targets),
      derivations(__derivations),
      model(__model),
      sampler(__sampler),
      beam(__beam) {}

  void operator()()
  {
    PYPGraph graph;
    
    size_type pos;
    
    for (;;) {
      mapper.pop(pos);
      
      if (pos == size_type(-1)) break;

      utils::resource res1;
      
      if (! derivations[pos].empty()) {
	counts.clear();
	counts.resize(3, size_type(0));
	ids.clear();
	
	derivation_type::const_iterator diter_end = derivations[pos].end();
	for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter) {
	  ++ counts[diter->is_terminal() ? PYP::TERMINAL : (diter->is_straight() ? PYP::STRAIGHT : PYP::INVERTED)];
	  
	  if (diter->is_terminal())
	    ids.push_back(diter->word_pair);
	}
	
	for (size_type i = 0; i != counts.size(); ++ i)
	  if (counts[i])
	    model.rule.table.decrement(i, counts[i], sampler);
	
	std::sort(ids.begin(), ids.end(), std::greater<id_type>());

	model.terminal.decrement(ids.begin(), ids.end(), sampler);
      }

      utils::resource res2;
      
      graph.initialize(sources[pos], targets[pos], model);

      utils::resource res3;
      
      graph.forward(sources[pos], targets[pos], beam);
      
      utils::resource res4;
      
      graph.backward(sources[pos], targets[pos], derivations[pos], sampler, temperature);
      
      utils::resource res5;
      
      {
	counts.clear();
	counts.resize(3, size_type(0));
	ids.clear();
	
	derivation_type::iterator diter_end = derivations[pos].end();
	for (derivation_type::iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter) {
	  
	  ++ counts[diter->is_terminal() ? PYP::TERMINAL : (diter->is_straight() ? PYP::STRAIGHT : PYP::INVERTED)];
	  
	  if (diter->is_terminal()) {
	    diter->word_pair = model.terminal.word_pair_id(diter->span.source.empty() ? vocab_type::EPSILON : sources[pos][diter->span.source.first],
							   diter->span.target.empty() ? vocab_type::EPSILON : targets[pos][diter->span.target.first]);
	    
	    ids.push_back(diter->word_pair);
	  }
	}
	
	for (size_type i = 0; i != counts.size(); ++ i)
	  if (counts[i])
	    model.rule.table.increment(i, counts[i], i == PYP::TERMINAL ? model.rule.p0_terminal : model.rule.p0, sampler, temperature);
	
	std::sort(ids.begin(), ids.end(), std::greater<id_type>());
	
	model.terminal.increment(ids.begin(), ids.end(), sampler, temperature);
      }
      
      utils::resource res6;
      
      time.decrement  += res2.thread_time() - res1.thread_time();
      time.initialize += res3.thread_time() - res2.thread_time();
      time.forward    += res4.thread_time() - res3.thread_time();
      time.backward   += res5.thread_time() - res4.thread_time();
      time.increment  += res6.thread_time() - res5.thread_time();

      reducer_derivation.push(pos);
      
      reducer.increment();
    }
  }
  
  queue_type&   mapper;
  queue_type&   reducer_derivation;
  counter_type& reducer;
  
  const sentence_set_type& sources;
  const sentence_set_type& targets;
  derivation_set_type& derivations;
  
  PYPITG&  model;
  sampler_type sampler;
  
  logprob_type beam;
  
  double temperature;
  Time   time;
};

struct less_size
{
  less_size(const sentence_set_type& __sources,
	    const sentence_set_type& __targets)
    : sources(__sources), targets(__targets) {}

  bool operator()(const size_type& x, const size_type& y) const
  {
    return (sources[x].size() + targets[x].size()) < (sources[y].size() + targets[y].size());
  }
  
  const sentence_set_type& sources;
  const sentence_set_type& targets;
};

inline
path_type add_suffix(const path_type& path, const std::string& suffix)
{
  bool has_suffix_gz  = false;
  bool has_suffix_bz2 = false;
  
  path_type path_added = path;
  
  if (path.extension() == ".gz") {
    path_added = path.parent_path() / path.stem();
    has_suffix_gz = true;
  } else if (path.extension() == ".bz2") {
    path_added = path.parent_path() / path.stem();
    has_suffix_bz2 = true;
  }
  
  path_added = path_added.string() + suffix;
  
  if (has_suffix_gz)
    path_added = path_added.string() + ".gz";
  else if (has_suffix_bz2)
    path_added = path_added.string() + ".bz2";
  
  return path_added;
}

path_type train_source_file = "-";
path_type train_target_file = "-";

path_type test_source_file;
path_type test_target_file;

path_type output_file;      // dump model file
path_type output_test_file; // dump test output file

path_type output_sample_file;
path_type output_model_file;

double beam = 1e-4;

int samples = 1;
int burns = 10;
int baby_steps = 1;
int anneal_steps = 1;
int resample_rate = 1;
int resample_iterations = 2;
bool slice_sampling = false;

bool sample_hypergraph = false;
bool sample_alignment = false;

double epsilon_prior = 1e-4;

double rule_prior_terminal = 1.0 / 3;

double rule_discount_alpha = 1.0;
double rule_discount_beta  = 1.0;
double rule_strength_shape = 1.0;
double rule_strength_rate  = 1.0;

double terminal_discount_alpha = 1.0;
double terminal_discount_beta  = 1.0;
double terminal_strength_shape = 10.0;
double terminal_strength_rate  = 0.1;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);

void prepare(const sentence_set_type& sources,
	     const sentence_set_type& targets,
	     PYPITG& model);
size_t read_data(const path_type& path,
		 sentence_set_type& sentences);

void viterbi(const path_type& output_file,
	     const path_type& source_file,
	     const path_type& target_file,
	     const PYPITG& model);

template <typename Dumper>
struct Reducer
{
  typedef typename Task::size_type  size_type;
  typedef typename Task::queue_type queue_type;
  
  Reducer(const path_type& path,
	  queue_type& queue,
	  const sentence_set_type& sources,
	  const sentence_set_type& targets,
	  const derivation_set_type& derivations)
    : path_(path),
      queue_(queue),
      sources_(sources),
      targets_(targets),
      derivations_(derivations) {}
  

  void operator()()
  {
    if (path_.empty()) {
      size_type pos;
      
      for (;;) {
	queue_.pop(pos);
	if (pos == size_type(-1)) break;
      }
    } else {
      typedef std::vector<size_type, std::allocator<size_type> > heap_type;
      
      utils::compress_ostream os(path_, 1024 * 1024);
      
      size_type pos;
      
      size_type id = 0;
      heap_type heap;
      
      for (;;) {
	queue_.pop(pos);
	if (pos == size_type(-1)) break;
	
	if (id == pos) {
	  if (! derivations_[id].empty())
	    dumper_(os, sources_[id], targets_[id], derivations_[id]);
	  os << '\n';
	  ++ id;
	} else {
	  heap.push_back(pos);
	  std::push_heap(heap.begin(), heap.end(), std::greater<size_type>());
	}

	while (! heap.empty() && id == heap.front()) {
	  if (! derivations_[id].empty())
	    dumper_(os, sources_[id], targets_[id], derivations_[id]);
	  os << '\n';
	  ++ id;
	  
	  std::pop_heap(heap.begin(), heap.end(), std::greater<size_type>());
	  heap.pop_back();
	}
      }
      
      while (! heap.empty() && id == heap.front()) {
	if (! derivations_[id].empty())
	  dumper_(os, sources_[id], targets_[id], derivations_[id]);
	os << '\n';
	++ id;
	
	std::pop_heap(heap.begin(), heap.end(), std::greater<size_type>());
	heap.pop_back();
      }
      
      if (! heap.empty())
	throw std::runtime_error("error while dumping...?");
    }
  }
  
  const path_type path_;
  
  queue_type& queue_;
  
  const sentence_set_type&   sources_;
  const sentence_set_type&   targets_;
  const derivation_set_type& derivations_;
  
  Dumper dumper_;
};

struct DumpNone
{
  std::ostream& operator()(std::ostream& os, const sentence_type& source, const sentence_type& target, const derivation_type& derivation)
  {
    return os;
  }
};

struct DumpDerivation
{
  typedef std::vector<std::string, std::allocator<std::string> > stack_type;
  
  std::ostream& operator()(std::ostream& os, const sentence_type& source, const sentence_type& target, const derivation_type& derivation)
  {
    stack.clear();
    derivation_type::const_iterator diter_end = derivation.end();
    for (derivation_type::const_iterator diter = derivation.begin(); diter != diter_end; ++ diter) {
      if (diter->is_terminal()) {
	
	const word_type& word_source = (! diter->span.source.empty() ? source[diter->span.source.first] : vocab_type::EPSILON);
	const word_type& word_target = (! diter->span.target.empty() ? target[diter->span.target.first] : vocab_type::EPSILON);
	
	os << "((( " << word_source << " ||| " << word_target << " )))";
	
	while (! stack.empty() && stack.back() != " ") {
	  os << stack.back();
	  stack.pop_back();
	}
	
	if (! stack.empty() && stack.back() == " ") {
	  os << stack.back();
	  stack.pop_back();
	}
	
      } else if (diter->is_straight()) {
	os << "[ ";
	stack.push_back(" ]");
	stack.push_back(" ");
      } else {
	os << "< ";
	stack.push_back(" >");
	stack.push_back(" ");
      }
    }
    
    return os;
  }
  
  stack_type stack;
};

struct DumpAlignment
{
  std::ostream& operator()(std::ostream& os, const sentence_type& source, const sentence_type& target, const derivation_type& derivation)
  {
    alignment.clear();
    
    derivation_type::const_iterator diter_end = derivation.end();
    for (derivation_type::const_iterator diter = derivation.begin(); diter != diter_end; ++ diter) {
      if (diter->is_terminal() && ! diter->span.source.empty() && ! diter->span.target.empty()) {
	
	for (size_type src = diter->span.source.first; src != diter->span.source.last; ++ src)
	  for (size_type trg = diter->span.target.first; trg != diter->span.target.last; ++ trg)
	    alignment.push_back(std::make_pair(src, trg));
      }
      
    }
    
    std::sort(alignment.begin(), alignment.end());
    
    os << alignment;
    
    return os;
  }
  
  alignment_type alignment;
};

struct DumpHypergraph
{
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > stack_type;
  typedef hypergraph_type::rule_type     rule_type;
  typedef hypergraph_type::rule_ptr_type rule_ptr_type;

  DumpHypergraph()
  {
    rule_type::symbol_set_type straight(2);
    rule_type::symbol_set_type inverted(2);
    
    straight[0] = vocab_type::X1;
    straight[1] = vocab_type::X2;
    
    inverted[0] = vocab_type::X2;
    inverted[1] = vocab_type::X1;
	  
    rule_straight = rule_type::create(rule_type(vocab_type::X, straight));
    rule_inverted = rule_type::create(rule_type(vocab_type::X, inverted));
  }

  std::ostream& operator()(std::ostream& os, const sentence_type& source, const sentence_type& target, const derivation_type& derivation)
  {
    hypergraph_type::edge_type::node_set_type tails(2);
    
    graph_source.clear();
    graph_target.clear();
    
    graph_source.goal = graph_source.add_node().id;
    graph_target.goal = graph_target.add_node().id;
    
    stack.clear();
    stack.push_back(0);
	      
    derivation_type::const_iterator diter_end = derivation.end();
    for (derivation_type::const_iterator diter = derivation.begin(); diter != diter_end; ++ diter) {
      const hypergraph_type::id_type head = stack.back();
      stack.pop_back();
		
      if (diter->is_terminal()) {
	hypergraph_type::edge_type& edge_source = graph_source.add_edge();
	hypergraph_type::edge_type& edge_target = graph_target.add_edge();
		  

	if (diter->span.source.empty())
	  edge_source.rule = rule_type::create(rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1));
	else
	  edge_source.rule = rule_type::create(rule_type(vocab_type::X,
							 source.begin() + diter->span.source.first,
							 source.begin() + diter->span.source.last));

	if (diter->span.target.empty())
	  edge_target.rule = rule_type::create(rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1));
	else
	  edge_target.rule = rule_type::create(rule_type(vocab_type::X,
							 target.begin() + diter->span.target.first,
							 target.begin() + diter->span.target.last));
		  
	graph_source.connect_edge(edge_source.id, head);
	graph_target.connect_edge(edge_target.id, head);
		  
      } else {
	hypergraph_type::node_type& node_source1 = graph_source.add_node();
	hypergraph_type::node_type& node_source2 = graph_source.add_node();

	hypergraph_type::node_type& node_target1 = graph_target.add_node();
	hypergraph_type::node_type& node_target2 = graph_target.add_node();
		  
	// push in reverse order!
	stack.push_back(node_source2.id);
	stack.push_back(node_source1.id);
		  
	tails[0] = node_source1.id;
	tails[1] = node_source2.id;
		  
	hypergraph_type::edge_type& edge_source = graph_source.add_edge(tails.begin(), tails.end());
	hypergraph_type::edge_type& edge_target = graph_target.add_edge(tails.begin(), tails.end());

	if (diter->is_straight()) {
	  edge_source.rule = rule_straight;
	  edge_target.rule = rule_straight;
	} else {
	  edge_source.rule = rule_inverted;
	  edge_target.rule = rule_inverted;
	}
	
	graph_source.connect_edge(edge_source.id, head);
	graph_target.connect_edge(edge_target.id, head);
      }
    }
    
    graph_source.topologically_sort();
    graph_target.topologically_sort();
    
    os << graph_source << " ||| " << graph_target;
    
    return os;
  }
  
  stack_type stack;
  
  hypergraph_type graph_source;
  hypergraph_type graph_target;

  rule_ptr_type rule_straight;
  rule_ptr_type rule_inverted;
};

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    if (samples < 0)
      throw std::runtime_error("# of samples must be positive");
    
    if (resample_rate <= 0)
      throw std::runtime_error("resample rate must be >= 1");
    
    if (beam < 0.0 || beam > 1.0)
      throw std::runtime_error("invalid beam width");
    
    sentence_set_type   sources;
    sentence_set_type   targets;
    
    const size_type source_vocab_size = read_data(train_source_file, sources);
    const size_type target_vocab_size = read_data(train_target_file, targets);

    if (sources.empty())
      throw std::runtime_error("no source sentence data?");

    if (targets.empty())
      throw std::runtime_error("no target sentence data?");
    
    if (sources.size() != targets.size())
      throw std::runtime_error("source/target side do not match!");

    PYPRule model_rule(rule_prior_terminal,
		       PYPRule::parameter_type(rule_discount_alpha,
					       rule_discount_beta,
					       rule_strength_shape,
					       rule_strength_rate));

    PYPTerminal model_terminal(PYPTerminal::parameter_type(terminal_discount_alpha,
							terminal_discount_beta,
							terminal_strength_shape,
							terminal_strength_rate),
			     1.0 / double(source_vocab_size),
			     1.0 / double(target_vocab_size),
			     epsilon_prior);
    
    PYPITG model(model_rule, model_terminal);

    // prepare model!
    prepare(sources, targets, model);
    
    derivation_set_type derivations(sources.size());
    position_set_type positions;
    for (size_t i = 0; i != sources.size(); ++ i)
      if (! sources[i].empty() && ! targets[i].empty())
	positions.push_back(i);
    position_set_type(positions).swap(positions);
    
    sampler_type sampler;
    
    // sample parameters, first...
    if (slice_sampling)
      model.slice_sample_parameters(sampler, resample_iterations);
    else
      model.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2)
      std::cerr << "rule: discount=" << model.rule.table.discount() << " strength=" << model.rule.table.strength() << std::endl
		<< "terminal=" << model.rule.prob_terminal() << " straight=" << model.rule.prob_straight() << " inverted=" << model.rule.prob_inverted() << std::endl
		<< "terminal: discount=" << model.terminal.table.discount() << " strength=" << model.terminal.table.strength() << std::endl;
    
    Task::queue_type mapper(threads);
    Task::queue_type reducer_derivation;
    Counter          reducer;
    
    std::vector<Task, std::allocator<Task> > tasks(threads, Task(mapper,
								 reducer_derivation,
								 reducer,
								 sources,
								 targets,
								 derivations,
								 model,
								 sampler,
								 beam));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);
    
    size_t baby_iter = 0;
    const size_t baby_last = utils::bithack::branch(baby_steps > 0, baby_steps, 0);
    
    size_t burn_iter = 0;
    const size_t burn_last = utils::bithack::branch(burns > 0, burns, 0);
    
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
      
      bool burn_finished = true;
      if (burn_iter != burn_last) {
	++ burn_iter;
	burn_finished = false;
      }
      
      sampling = anneal_finished && baby_finished && burn_finished;
      
      if (debug) {
	if (sampling)
	  std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
	else
	  std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
      }
      
      std::auto_ptr<boost::thread> dumper;
      
      if (sampling && ! output_sample_file.empty()) {
	// assign path...
	const path_type path = add_suffix(output_sample_file, "." + utils::lexical_cast<std::string>(sample_iter + 1));

	if (sample_hypergraph)
	  dumper.reset(new boost::thread(Reducer<DumpHypergraph>(path,
								 reducer_derivation,
								 sources,
								 targets,
								 derivations)));
	else if (sample_alignment)
	  dumper.reset(new boost::thread(Reducer<DumpAlignment>(path,
								reducer_derivation,
								sources,
								targets,
								derivations)));
	else
	  dumper.reset(new boost::thread(Reducer<DumpDerivation>(path,
								 reducer_derivation,
								 sources,
								 targets,
								 derivations)));
      } else
	dumper.reset(new boost::thread(Reducer<DumpNone>(path_type(),
							 reducer_derivation,
							 sources,
							 targets,
							 derivations)));
      
      // assign temperature...
      Time time;
      
      for (size_type i = 0; i != tasks.size(); ++ i) {
	tasks[i].temperature = temperature;
	
	time += tasks[i].time;
      }
      
      // shuffle
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(positions.begin(), positions.end(), gen);
      
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_size(sources, targets));
      
      std::auto_ptr<boost::progress_display> progress(debug
						      ? new boost::progress_display(positions.size(), std::cerr, "", "", "")
						      : 0);
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	mapper.push(*piter);
	
	if (debug)
	  ++ (*progress);
      }
      
      reducer.wait(positions.size());
      reducer.clear();
            
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  model.slice_sample_parameters(sampler, resample_iterations);
	else
	  model.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2)
	  std::cerr << "rule: discount=" << model.rule.table.discount() << " strength=" << model.rule.table.strength() << std::endl
		    << "terminal=" << model.rule.prob_terminal() << " straight=" << model.rule.prob_straight() << " inverted=" << model.rule.prob_inverted() << std::endl
		    << "terminal: discount=" << model.terminal.table.discount() << " strength=" << model.terminal.table.strength() << std::endl;
      }
      
      if (debug >= 2) {
	Time time_end;
	for (size_type i = 0; i != tasks.size(); ++ i)
	  time_end += tasks[i].time;
	
	std::cerr << "initialize: " << (time_end.initialize - time.initialize) / tasks.size() << " seconds" << std::endl
		  << "forward: " << (time_end.forward - time.forward) / tasks.size() << " seconds" << std::endl
		  << "backward: " << (time_end.backward - time.backward) / tasks.size() << " seconds" << std::endl
		  << "increment: " << (time_end.increment - time.increment) / tasks.size() << " seconds" << std::endl
		  << "decrement: " << (time_end.decrement - time.decrement) / tasks.size() << " seconds" << std::endl;
      }
      
      if (debug)
	std::cerr << "log-likelihood: " << model.log_likelihood() << std::endl;

      if (sampling && ! output_model_file.empty()) {
	const path_type path = add_suffix(output_model_file, "." + utils::lexical_cast<std::string>(sample_iter + 1));
	
	utils::compress_ostream os(path, 1024 * 1024);
	os.precision(20);
	
	for (PYP::id_type id = 0; id != model.terminal.table.size(); ++ id)
	  if (! model.terminal.table[id].empty())
	    os << model.terminal.word_pairs[id].source
	       << ' ' << model.terminal.word_pairs[id].target
	       << ' ' << model.terminal.prob(id) << '\n';
      }
      
      reducer_derivation.push(size_type(-1));
      dumper->join();
    }
    
    for (int i = 0; i != threads; ++ i)
      mapper.push(size_type(-1));
    
    workers.join_all();
    
    // compute viterbi with test...
    if (! output_test_file.empty()) {
      if (test_source_file.empty() || ! boost::filesystem::exists(test_source_file))
	throw std::runtime_error("test file output specified, but no source data?");
      if (test_target_file.empty() || ! boost::filesystem::exists(test_target_file))
	throw std::runtime_error("test file output specified, but no target data?");
      
      viterbi(output_test_file, test_source_file, test_target_file, model);
    }
    
    // TODO: output final model file...
    if (! output_file.empty()) {
      
      
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Tp>
struct greater_psecond
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

struct PrepareMapper
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef PYP::word_pair_type word_pair_type;

  typedef utils::compact_map<word_pair_type, size_t,
			     PYP::word_pair_unassigned, PYP::word_pair_unassigned,
			     boost::hash<word_pair_type>, std::equal_to<word_pair_type>,
			     std::allocator<std::pair<const word_pair_type, size_t> >  > count_set_type;
  
  PrepareMapper(const sentence_set_type& __sources,
		const sentence_set_type& __targets)
    : sources(__sources), targets(__targets) {  }
  
  void operator()()
  {
    for (size_type i = first; i != last; ++ i)
      if (! sources[i].empty() && ! targets[i].empty())
	for (size_t src = 0; src <= sources[i].size(); ++ src)
	  for (size_t trg = (src == 0); trg <= targets[i].size(); ++ trg)
	    ++ counts[word_pair_type(src == 0 ? vocab_type::EPSILON : sources[i][src - 1],
				     trg == 0 ? vocab_type::EPSILON : targets[i][trg - 1])];
  }
  
  count_set_type counts;
  
  const sentence_set_type& sources;
  const sentence_set_type& targets;
  
  size_type first;
  size_type last;
};

void prepare(const sentence_set_type& sources,
	     const sentence_set_type& targets,
	     PYPITG& model)
{
  typedef PYP::word_pair_type word_pair_type;
  typedef PrepareMapper::count_set_type count_set_type;
  typedef std::vector<const count_set_type::value_type*, std::allocator<const count_set_type::value_type*> > sorted_type;

  std::vector<PrepareMapper, std::allocator<PrepareMapper> > mapper(threads, PrepareMapper(sources, targets));
  
  const size_type interval = (sources.size() + threads - 1) / threads;
  size_type first = 0;
  
  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i) {
    const size_type last = utils::bithack::min(first + interval, sources.size());
    
    mapper[i].first = first;
    mapper[i].last  = last;
    
    workers.add_thread(new boost::thread(boost::ref(mapper[i])));
    
    first = last;
  }
  workers.join_all();

  count_set_type counts;
  
  for (int i = 0; i != threads; ++ i) {
    if (counts.empty())
      counts.swap(mapper[i].counts);
    else {
      count_set_type::const_iterator citer_end = mapper[i].counts.end();
      for (count_set_type::const_iterator citer = mapper[i].counts.begin(); citer != citer_end; ++ citer)
	counts[citer->first] += citer->second;
    }
    
    mapper[i].counts.clear();
  }
  
  sorted_type sorted;
  sorted.reserve(counts.size());
  
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
    sorted.push_back(&(*citer));
  
  std::sort(sorted.begin(), sorted.end(), greater_psecond<const count_set_type::value_type>());

  sorted_type::const_iterator siter_end = sorted.end();
  for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
    model.terminal.word_pair_id((*siter)->first.source, (*siter)->first.target);
  
  model.terminal.table.reserve(model.terminal.word_pairs.size());
  model.terminal.table.resize(model.terminal.word_pairs.size());
  model.rule.table.reserve(3);
  model.rule.table.resize(3);
}

struct ViterbiMapReduce
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  typedef PYP::id_type         id_type;

  struct bitext_type
  {
    size_type       id;
    sentence_type   source;
    sentence_type   target;
    derivation_type derivation;
    
    bitext_type() : id(size_type(-1)), source(), target(), derivation() {}

    friend
    bool operator<(const bitext_type& x, const bitext_type& y)
    {
      return x.id < y.id;
    }
    
    void swap(bitext_type& x)
    {
      std::swap(id, x.id);
      source.swap(x.source);
      target.swap(x.target);
      derivation.swap(x.derivation);
    }
  };
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type > > queue_type;
};

namespace std
{
  inline
  void swap(ViterbiMapReduce::bitext_type& x, ViterbiMapReduce::bitext_type& y)
  {
    x.swap(y);
  }
};

struct ViterbiMapper
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  typedef PYP::id_type         id_type;

  typedef PYP::logprob_type logprob_type;
  typedef PYP::prob_type    prob_type;

  typedef ViterbiMapReduce::bitext_type bitext_type;
  typedef ViterbiMapReduce::queue_type  queue_type;
  
  ViterbiMapper(queue_type& __mapper,
		queue_type& __reducer,
		const PYPITG& __model,
		const logprob_type& __beam)
    : mapper(__mapper), reducer(__reducer), model(__model), beam(__beam) {}
  
  void operator()()
  {
    PYPViterbi viterbi;
    
    bitext_type bitext;
    
    for (;;) {
      mapper.pop_swap(bitext);
      
      if (bitext.id == size_type(-1)) break;
      
      bitext.derivation.clear();
      
      if (! bitext.source.empty() && ! bitext.target.empty()) {
	viterbi.initialize(bitext.source, bitext.target, model);
	
	viterbi.forward(bitext.source, bitext.target, beam);
	
	viterbi.backward(bitext.source, bitext.target, bitext.derivation);
      }
      
      reducer.push_swap(bitext);
    }
  }
  
  queue_type&   mapper;
  queue_type&   reducer;
  const PYPITG& model;
  const logprob_type beam;
};

template <typename Dumper>
struct ViterbiReducer
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  typedef PYP::id_type         id_type;

  typedef PYP::logprob_type logprob_type;
  typedef PYP::prob_type    prob_type;

  typedef ViterbiMapReduce::bitext_type bitext_type;
  typedef ViterbiMapReduce::queue_type  queue_type;  
  
  ViterbiReducer(queue_type& __reducer,
		 std::ostream& __os,
		 const int __debug)
    : reducer(__reducer),
      os(__os),
      debug(__debug) {}

  void operator()()
  {
    typedef std::set<bitext_type, std::less<bitext_type>, std::allocator<bitext_type> > bitext_set_type;

    bitext_type     bitext;
    bitext_set_type bitexts;

    Dumper dumper;
    
    size_type id = 0;
    for (;;) {
      reducer.pop_swap(bitext);
      
      if (bitext.id == size_type(-1)) break;
      
      if (bitext.id != id)
	bitexts.insert(bitext);
      else {
	if (! bitext.derivation.empty())
	  dumper(os, bitext.source, bitext.target, bitext.derivation);
	os << '\n';
	++ id;
	
	if (debug) {
	  if (id % 10000 == 0)
	    std::cerr << '.';
	  if (id % 1000000 == 0)
	    std::cerr << '\n';
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id == id) {
	const bitext_type& bitext = *bitexts.begin();
	
	if (! bitext.derivation.empty())
	  dumper(os, bitext.source, bitext.target, bitext.derivation);
	os << '\n';
	++ id;
	
	if (debug) {
	  if (id % 10000 == 0)
	    std::cerr << '.';
	  if (id % 1000000 == 0)
	    std::cerr << '\n';
	}

	bitexts.erase(bitexts.begin());
      }
    }

    while (! bitexts.empty() && bitexts.begin()->id == id) {
      const bitext_type& bitext = *bitexts.begin();
      
      if (! bitext.derivation.empty())
	dumper(os, bitext.source, bitext.target, bitext.derivation);
      os << '\n';
      ++ id;

      if (debug) {
	if (id % 10000 == 0)
	  std::cerr << '.';
	if (id % 1000000 == 0)
	  std::cerr << '\n';
      }
      
      bitexts.erase(bitexts.begin());
    }
    
    if (! bitexts.empty())
      throw std::runtime_error("bitexts stil renam???");
    
    if (debug && id >= 10000 && id % 1000000 != 0)
      std::cerr << std::endl;
  }
  
  queue_type& reducer;
  std::ostream& os;
  int debug;
};

void viterbi(const path_type& output_file,
	     const path_type& source_file,
	     const path_type& target_file,
	     const PYPITG& model)
{
  typedef ViterbiMapReduce::bitext_type bitext_type;
  typedef ViterbiMapReduce::queue_type  queue_type;  
  
  const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
  
  utils::compress_ostream os(output_file, 1024 * 1024* (! flush_output));
  
  queue_type mapper(1024 * threads);
  queue_type reducer;
  
  boost::thread_group workers;
  boost::thread_group dumper;
  
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(ViterbiMapper(mapper, reducer, model, beam)));
  
  if (sample_hypergraph)
    dumper.add_thread(new boost::thread(ViterbiReducer<DumpHypergraph>(reducer, os, debug)));
  else if (sample_alignment)
    dumper.add_thread(new boost::thread(ViterbiReducer<DumpAlignment>(reducer, os, debug)));
  else
    dumper.add_thread(new boost::thread(ViterbiReducer<DumpDerivation>(reducer, os, debug)));
  
  utils::compress_istream is_src(source_file, 1024 * 1024);
  utils::compress_istream is_trg(target_file, 1024 * 1024);
  
  bitext_type bitext;
  bitext.id = 0;
  
  for (;;) {
    is_src >> bitext.source;
    is_trg >> bitext.target;
    
    if (! is_src || ! is_trg) break;
    
    mapper.push(bitext);
    
    ++ bitext.id;
  }
  
  if (is_src || is_trg)
    throw std::runtime_error("# of lines do not match...");
  
  // join all the workers...
  for (int i = 0; i != threads; ++ i)
    mapper.push(bitext_type());
  workers.join_all();
  
  // terminate dumper...
  reducer.push(bitext_type());
  dumper.join_all();
}


size_t read_data(const path_type& path, sentence_set_type& sentences)
{
  typedef utils::compact_set<word_type,
			     utils::unassigned<word_type>, utils::unassigned<word_type>,
			     boost::hash<word_type>, std::equal_to<word_type>,
			     std::allocator<word_type> > word_set_type;

  word_set_type words;

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
    
    ("output",      po::value<path_type>(&output_file),      "output file for model")
    ("output-test", po::value<path_type>(&output_test_file), "output file for the inferred test")
    
    ("output-sample", po::value<path_type>(&output_sample_file), "output derivation for each sample file")
    ("output-model",  po::value<path_type>(&output_model_file),  "output model for each sample file (or phrase table)")
    
    ("beam",                po::value<double>(&beam)->default_value(beam),                            "beam threshold")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("burns",               po::value<int>(&burns)->default_value(burns),                             "# of burn-ins")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    ("hypergraph",          po::bool_switch(&sample_hypergraph),                                      "dump sampled derivation in hypergraph")
    ("alignment",           po::bool_switch(&sample_alignment),                                       "dump sampled derivation in alignment")

    ("epsilon-prior",       po::value<double>(&epsilon_prior)->default_value(epsilon_prior),             "prior for epsilon")
    
    ("rule-prior-terminal", po::value<double>(&rule_prior_terminal)->default_value(rule_prior_terminal), "prior for terminal")
    
    ("rule-discount-alpha", po::value<double>(&rule_discount_alpha)->default_value(rule_discount_alpha), "discount ~ Beta(alpha,beta)")
    ("rule-discount-beta",  po::value<double>(&rule_discount_beta)->default_value(rule_discount_beta),   "discount ~ Beta(alpha,beta)")

    ("rule-strength-shape", po::value<double>(&rule_strength_shape)->default_value(rule_strength_shape), "strength ~ Gamma(shape,rate)")
    ("rule-strength-rate",  po::value<double>(&rule_strength_rate)->default_value(rule_strength_rate),   "strength ~ Gamma(shape,rate)")

    
    ("terminal-discount-alpha", po::value<double>(&terminal_discount_alpha)->default_value(terminal_discount_alpha), "discount ~ Beta(alpha,beta)")
    ("terminal-discount-beta",  po::value<double>(&terminal_discount_beta)->default_value(terminal_discount_beta),   "discount ~ Beta(alpha,beta)")

    ("terminal-strength-shape", po::value<double>(&terminal_strength_shape)->default_value(terminal_strength_shape), "strength ~ Gamma(shape,rate)")
    ("terminal-strength-rate",  po::value<double>(&terminal_strength_rate)->default_value(terminal_strength_rate),   "strength ~ Gamma(shape,rate)")

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


