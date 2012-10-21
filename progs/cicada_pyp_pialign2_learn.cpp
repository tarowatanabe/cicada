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

// a simplified pialign:
//
// we do not use any base-measure, but use the 1.0 / (vocabulary for source * vocab size for target)
// and use epsilon-prior, as in cicada_pyp_itg_learn
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <map>
#include <deque>
#include <iterator>
#include <numeric>
#include <queue>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring/logprob.hpp>
#include <cicada/discounter.hpp>

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
#include "utils/restaurant.hpp"
#include "utils/restaurant_vector.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/dense_hash_map.hpp"
#include "utils/dense_hash_set.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/symbol_set.hpp"
#include "utils/unique_set.hpp"
#include "utils/rwticket.hpp"
#include "utils/spinlock.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/array.hpp>

typedef cicada::Vocab      vocab_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Symbol     symbol_type;
typedef cicada::Symbol     word_type;
typedef cicada::HyperGraph hypergraph_type;

struct PYP
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef enum {
    TERMINAL = 0,
    STRAIGHT,
    INVERTED,
    GENERATIVE,
    BASE,
  } itg_type;

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
      typedef utils::hashmurmur<size_t> hasher_type;
      
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
      typedef utils::hashmurmur<size_t> hasher_type;
      
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

    rule_type()
      : span(), left(), right(), itg() {}
    rule_type(const span_pair_type& __span, const itg_type& __itg)
      : span(__span), left(), right(), itg(__itg) {}
    rule_type(const span_pair_type& __span, const span_pair_type& __left, const span_pair_type& __right, const itg_type& __itg)
      : span(__span), left(__left), right(__right), itg(__itg) {}
    
    bool is_terminal() const { return left.empty() && right.empty(); }
    bool is_straight() const { return ! is_terminal() && left.target.last  == right.target.first; }
    bool is_inverted() const { return ! is_terminal() && left.target.first == right.target.last; }
  };

  struct phrase_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef const word_type* iterator;
    typedef const word_type* const_iterator;
    
    phrase_type() : first(0), last(0) {}
    phrase_type(const_iterator __first, const_iterator __last)
      : first(__first), last(__last) { }
    template <typename Iterator>
    phrase_type(Iterator __first, Iterator __last)
      : first(&(*__first)), last(&(*__last)) { }
    
    bool empty() const { return first == last; }
    size_type size() const { return last - first; }
    
    const_iterator begin() const { return first; }
    const_iterator end() const { return last; }
    
    friend
    bool operator==(const phrase_type& x, const phrase_type& y)
    {
      return (x.size() == y.size() && std::equal(x.begin(), x.end(), y.begin()));
    }
    
    friend
    size_t hash_value(phrase_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.begin(), x.end(), 0);
    }
    
    friend
    std::ostream& operator<<(std::ostream& os, const phrase_type& x)
    {
      if (! x.empty()) {
	std::copy(x.begin(), x.end() - 1, std::ostream_iterator<word_type>(os, " "));
	os << *(x.end() - 1);
      }
      return os;
    }
    
    const_iterator first;
    const_iterator last;
  };

  struct phrase_pair_type
  {
    phrase_type source;
    phrase_type target;
    
    phrase_pair_type() : source(), target() {}
    phrase_pair_type(const phrase_type& __source, const phrase_type& __target)
      : source(__source), target(__target) {}
    template <typename IteratorSource, typename IteratorTarget>
    phrase_pair_type(IteratorSource source_first, IteratorSource source_last,
		     IteratorTarget target_first, IteratorTarget target_last)
      : source(source_first, source_last), target(target_first, target_last) {}
    
    friend
    bool operator==(const phrase_pair_type& x, const phrase_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    size_t hash_value(phrase_pair_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
    }
  };
  
  struct non_terminal_trie_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint32_t  id_type;
    
    typedef std::vector<id_type, std::allocator<id_type> > id_set_type;
    
    //
    // multiple of 4: L, R, l, r
    //
    
    void initialize(const size_type depth)
    {
      parent_.clear();
      child_.clear();
      last_.clear();
      
      // 4, 16, 64, 256, 1024 ...
      
      // compute parent_ and last_
      parent_.resize(4, id_type(-1));
      last_.push_back(parent_.size());
      id_type prev = 0;
      for (size_type order = 0; order != depth; ++ order) {
	id_type last = parent_.size();
	for (id_type id = prev; id != last; ++ id) {
	  child_.push_back(parent_.size());
	  parent_.resize(parent_.size() + 4, id);
	}
	
	prev = last;
	last_.push_back(parent_.size());
      }
    }
    
    id_set_type parent_;
    id_set_type child_;
    id_set_type last_;
  };
};


struct PYPRule
{
  // rule selection probabilities
  // terminal, stright, inversion
  
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef PYP::phrase_type      phrase_type;
  typedef PYP::phrase_pair_type phrase_pair_type;

  typedef PYP::rule_type rule_type;
  
  typedef PYP::itg_type itg_type;

  typedef uint32_t id_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  typedef utils::restaurant_vector<> table_type;

  typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;

  PYPRule(const double& __p0_terminal, const parameter_type& parameter)
    : p0_terminal(__p0_terminal), p0((1.0 - __p0_terminal) * 0.5), counts0_terminal(0), counts0(0), table(parameter) {}
  
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
    
    if (PYP::TERMINAL < table.size()) {
      counts0_terminal = table[PYP::TERMINAL].size_table();
      counts0          = table.size_table() - table[PYP::TERMINAL].size_table();
    }
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
    
    if (PYP::TERMINAL < table.size()) {
      counts0_terminal = table[PYP::TERMINAL].size_table();
      counts0          = table.size_table() - table[PYP::TERMINAL].size_table();
    }
  }
  
  double     p0_terminal;
  double     p0;
  size_type  counts0_terminal;
  size_type  counts0;
  table_type table;
};

struct PYPPhrase
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::phrase_type      phrase_type;
  typedef PYP::phrase_pair_type phrase_pair_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  
  typedef utils::symbol_set<phrase_type, boost::hash<phrase_type>, std::equal_to<phrase_type>, std::allocator<phrase_type > > phrase_set_type;
  
  typedef phrase_set_type::index_type id_type;
  typedef std::pair<id_type, id_type> id_pair_type;
  
  typedef utils::symbol_set<id_pair_type, utils::hashmurmur<size_t>, std::equal_to<id_pair_type>, std::allocator<id_pair_type> > phrase_pair_set_type;
  
  
  //typedef utils::restaurant<id_pair_type, utils::hashmurmur<size_t>, std::equal_to<id_pair_type>, std::allocator<id_pair_type > > table_type;
  typedef utils::restaurant_vector<> table_type;
  
  PYPPhrase(const double& __p0,
	    const parameter_type& parameter)
    : p0(__p0),
      counts0(0),
      table(parameter),
      phrases(),
      phrase_pairs() {}
  
  std::pair<id_type, bool> phrase_find(const phrase_type& phrase) const
  {
    phrase_set_type::const_iterator iter = phrases.find(phrase);
    return std::make_pair(iter - phrases.begin(), iter != phrases.end());
  }
  
  
  id_type phrase_id(const phrase_type& phrase)
  {
    phrase_set_type::iterator iter = phrases.insert(phrase).first;
    return iter - phrases.begin();
  }
  
  std::pair<id_type, bool> phrase_pair_find(const id_type& source, const id_type& target) const
  {
    phrase_pair_set_type::const_iterator iter = phrase_pairs.find(std::make_pair(source, target));
    return std::make_pair(iter - phrase_pairs.begin(), iter != phrase_pairs.end());
  }
  
  
  id_type phrase_pair_id(const id_type& source, const id_type& target)
  {
    phrase_pair_set_type::iterator iter = phrase_pairs.insert(std::make_pair(source, target)).first;
    return iter - phrase_pairs.begin();
  }
  
  template <typename Sampler>
  void increment_existing(const phrase_pair_type& phrase_pair, const bool base, Sampler& sampler, const double temperature=1.0)
  {
    const id_type id_source = phrase_id(phrase_pair.source);
    const id_type id_target = phrase_id(phrase_pair.target);
    
    const id_type id_pair = phrase_pair_id(id_source, id_target);
    
    table.increment_existing(id_pair, sampler);

    counts0 += (base && phrase_pair.source.size() <= 1 && phrase_pair.target.size() <= 1);
  }

  template <typename Sampler>
  void increment_new(const phrase_pair_type& phrase_pair, const bool base, Sampler& sampler, const double temperature=1.0)
  {
    const id_type id_source = phrase_id(phrase_pair.source);
    const id_type id_target = phrase_id(phrase_pair.target);
    
    const id_type id_pair = phrase_pair_id(id_source, id_target);
    
    table.increment_new(id_pair, sampler);

    counts0 += (base && phrase_pair.source.size() <= 1 && phrase_pair.target.size() <= 1);
  }
  
  template <typename Sampler>
  void decrement(const phrase_pair_type& phrase_pair, const bool base, Sampler& sampler)
  {
    const id_type id_source = phrase_id(phrase_pair.source);
    const id_type id_target = phrase_id(phrase_pair.target);
    
    const id_type id_pair = phrase_pair_id(id_source, id_target);
    
    table.decrement(id_pair, sampler);

    counts0 -= (base && phrase_pair.source.size() <= 1 && phrase_pair.target.size() <= 1);
  }
  
  logprob_type logprob_fallback() const
  {
    return table.prob(logprob_type(1.0));
  }

  logprob_type logprob_base() const
  {
    return table.prob(p0);
  }
  
  std::pair<logprob_type, bool> prob(const phrase_type& source, const phrase_type& target) const
  {
    phrase_set_type::const_iterator siter = phrases.find(source);
    phrase_set_type::const_iterator titer = phrases.find(target);
    
    if (siter == phrases.end() || titer == phrases.end())
      return std::make_pair(table.prob(p0), false);
    else {
      phrase_pair_set_type::const_iterator piter = phrase_pairs.find(std::make_pair(siter - phrases.begin(),
										    titer - phrases.begin()));
      
      if (piter == phrase_pairs.end())
	return std::make_pair(table.prob(p0), false);
      else
	return table.prob_model(piter - phrase_pairs.begin(), p0);
    }
  }
  
  double log_likelihood() const
  {
    return table.log_likelihood() + cicada::semiring::log(p0) * counts0;
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

  void prune(bool debug)
  {
    table.prune();
    
    // erase unused phrase entry in phrases...
    std::vector<bool, std::allocator<bool> > inserted(phrases.size());
    
    size_type phrase_pair = 0;
    size_type phrase_pair_erased = 0;
    
    for (id_type id = 0; id != table.size(); ++ id) {
      if (! table[id].empty()) {
	const id_pair_type& pair = phrase_pairs[id];
	
	inserted[pair.first] = true;
	inserted[pair.second] = true;
	
	++ phrase_pair;
      } else {
	phrase_pairs.erase(id);
	
	++ phrase_pair_erased;
      }
    }
    
    size_type phrase = 0;
    size_type erased = 0;
    for (id_type id = 0; id != inserted.size(); ++ id) {
      if (! inserted[id]) {
	phrases.erase(id);
	++ erased;
      } else
	++ phrase;
    }
    
    if (debug)
      std::cerr << "# of phrases: " << phrase << " erased: " << erased << std::endl
		<< "# of phrase pairs: " << phrase_pair << " erased: " << phrase_pair_erased << std::endl;
  }
  
  logprob_type p0;
  size_type    counts0;
  
  table_type           table;
  phrase_set_type      phrases;
  phrase_pair_set_type phrase_pairs;
};

struct PYPEpsilon
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::phrase_type      phrase_type;
  typedef PYP::phrase_pair_type phrase_pair_type;
  
  typedef PYP::rule_type rule_type;
  typedef PYP::itg_type  itg_type;

  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

  PYPEpsilon(const double& __p0_epsilon)
    : p0(1.0 - __p0_epsilon),
      p0_epsilon(__p0_epsilon),
      counts0(0),
      counts0_epsilon(0) {}
  
  template <typename Sampler>
  void increment(const phrase_pair_type& phrase_pair, Sampler& sampler, const double temperature=1.0)
  {
    if (phrase_pair.source.size() <= 1 && phrase_pair.target.size() <= 1) {
      counts0 += ! phrase_pair.source.empty();
      counts0 += ! phrase_pair.target.empty();
      counts0_epsilon += phrase_pair.source.empty();
      counts0_epsilon += phrase_pair.target.empty();
    }
  }
  
  template <typename Sampler>
  void decrement(const phrase_pair_type& phrase_pair, Sampler& sampler)
  {
    if (phrase_pair.source.size() <= 1 && phrase_pair.target.size() <= 1) {
      counts0 -= ! phrase_pair.source.empty();
      counts0 -= ! phrase_pair.target.empty();      
      counts0_epsilon -= phrase_pair.source.empty();
      counts0_epsilon -= phrase_pair.target.empty();
    }
  }
  
  double log_likelihood() const
  {
    return cicada::semiring::log(p0) * counts0 + cicada::semiring::log(p0_epsilon) * counts0_epsilon;
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    
  }

  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    
  }
  
  logprob_type p0;
  logprob_type p0_epsilon;
  size_type    counts0;
  size_type    counts0_epsilon;
};

struct PYPPiAlign
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::phrase_type      phrase_type;
  typedef PYP::phrase_pair_type phrase_pair_type;
  
  typedef PYP::rule_type rule_type;
  typedef PYP::itg_type  itg_type;

  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

  typedef utils::rwticket mutex_type;
  
  PYPPiAlign(const PYPRule&   __rule,
	     const PYPPhrase& __phrase,
	     const PYPEpsilon& __epsilon)
    : rule(__rule), phrase(__phrase), epsilon(__epsilon) {}
  
  template <typename Sampler>
  void increment(const sentence_type& source, const sentence_type& target, const rule_type& r, Sampler& sampler, const double temperature)
  {
    const phrase_pair_type phrase_pair(source.begin() + r.span.source.first, source.begin() + r.span.source.last,
				       target.begin() + r.span.target.first, target.begin() + r.span.target.last);
    
    rule.increment(r, sampler, temperature);
    
    if (r.itg == PYP::GENERATIVE)
      phrase.increment_existing(phrase_pair, r.itg != PYP::GENERATIVE, sampler, temperature);
    else
      phrase.increment_new(phrase_pair, r.itg != PYP::GENERATIVE, sampler, temperature);
    
    if (r.itg != PYP::GENERATIVE)
      epsilon.increment(phrase_pair, sampler, temperature);
  }
  
  template <typename Sampler>
  void decrement(const sentence_type& source, const sentence_type& target, const rule_type& r, Sampler& sampler)
  {
    const phrase_pair_type phrase_pair(source.begin() + r.span.source.first, source.begin() + r.span.source.last,
				       target.begin() + r.span.target.first, target.begin() + r.span.target.last);
    
    rule.decrement(r, sampler);
    
    phrase.decrement(phrase_pair, r.itg != PYP::GENERATIVE, sampler);
    
    if (r.itg != PYP::GENERATIVE)
      epsilon.decrement(phrase_pair, sampler);
  }

  double log_likelihood() const
  {
    return rule.log_likelihood() + phrase.log_likelihood() + epsilon.log_likelihood();
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    rule.sample_parameters(sampler, num_loop, num_iterations);

    phrase.sample_parameters(sampler, num_loop, num_iterations);
    
    epsilon.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    rule.slice_sample_parameters(sampler, num_loop, num_iterations);
    
    phrase.slice_sample_parameters(sampler, num_loop, num_iterations);

    epsilon.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  PYPRule    rule;
  PYPPhrase  phrase;
  PYPEpsilon epsilon;
  
  mutex_type mutex;
};

struct PYPGraph
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::phrase_type      phrase_type;
  typedef PYP::phrase_pair_type phrase_pair_type;
  
  typedef PYP::span_type      span_type;
  typedef PYP::span_pair_type span_pair_type;
  
  typedef PYP::rule_type rule_type;
  
  typedef std::vector<rule_type, std::allocator<rule_type> > derivation_type;

  typedef PYPPhrase::logprob_type logprob_type;
  typedef PYPPhrase::prob_type    prob_type;
  
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  typedef utils::bichart<logprob_type, std::allocator<logprob_type> > chart_type;
  typedef utils::bichart<logprob_type, std::allocator<logprob_type> > base_type;
  
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
  typedef utils::dense_hash_set<span_pairs_type, utils::hashmurmur<size_t>, std::equal_to<span_pairs_type>,
				std::allocator<span_pairs_type> >::type span_pairs_unique_type;

  typedef utils::chart<logprob_type, std::allocator<logprob_type> > unigram_type;
  
  typedef utils::chart<logprob_type, std::allocator<logprob_type> > model1_type;
  typedef std::vector<model1_type, std::allocator<model1_type> > model1_chart_type;
  
  typedef utils::chart<logprob_type, std::allocator<logprob_type> > chart_mono_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > alpha_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > beta_type;

  typedef std::pair<logprob_type, span_pair_type> score_span_pair_type;
  typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;

  static
  logprob_type geometric_mean(const logprob_type& x)
  {
    return cicada::semiring::traits<logprob_type>::exp(cicada::semiring::log(x) * 0.5);
  }
  
  void initialize(const sentence_type& source,
		  const sentence_type& target,
		  const PYPPiAlign& model)
  {
    //std::cerr << "initialize" << std::endl;

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
    
    
    //std::cerr << "initialize chart" << std::endl;
    
    
    logprob_base    = model.phrase.logprob_base() * model.epsilon.p0 * model.epsilon.p0;
    logprob_epsilon = model.phrase.logprob_base() * model.epsilon.p0 * model.epsilon.p0_epsilon;
    
    logprob_fallback = model.phrase.logprob_fallback();
    logprob_term = model.rule.prob_terminal() * logprob_fallback;
    logprob_str  = model.rule.prob_straight() * logprob_fallback;
    logprob_inv  = model.rule.prob_inverted() * logprob_fallback;
    
    // initialization for bases
    for (size_type source_first = 0; source_first <= source.size(); ++ source_first)
      for (size_type target_first = 0; target_first <= target.size(); ++ target_first) {

	// epsilons...
	if (source_first < source.size()) {
	  const size_type source_last = source_first + 1;
	  const size_type target_last = target_first;
	  
	  const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	  const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);
	  
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);

	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::BASE), logprob_epsilon));
	  
	  chart(source_first, source_last, target_first, target_last) += logprob_epsilon;
	  
	  chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_epsilon);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
	
	// epsilons...
	if (target_first < target.size()) {
	  const size_type target_last = target_first + 1;
	  const size_type source_last = source_first;
	    
	  const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	  const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);

	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::BASE), logprob_epsilon));
	  
	  chart(source_first, source_last, target_first, target_last) += logprob_epsilon;
	  
	  chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_epsilon);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
	
	// word-pair
	
	if (source_first < source.size() && target_first < target.size()) {
	  const size_type source_last = source_first + 1;
	  const size_type target_last = target_first + 1;
	  
	  const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	  const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);

	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  
	  edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::BASE), logprob_base));
	  
	  chart(source_first, source_last, target_first, target_last) += logprob_base;
	  
	  chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_base);
	  chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_base);
	  
	  agenda[span_pair.size()].push_back(span_pair);
	}
      }
    
    
    // actually, we should be carefull with the chart assignment, since we need to consult the base of phrase-model + rule-model

    // epsilons.. 
    const std::pair<PYPPhrase::id_type, bool> empty_id = model.phrase.phrase_find(phrase_type(source.begin(), source.begin()));

    for (size_type source_first = 0; source_first <= source.size(); ++ source_first)
      for (size_type target_first = 0; target_first <= target.size(); ++ target_first) {
		
	if (empty_id.second && source_first < source.size())
	  for (size_type source_last = source_first + 1; source_last <= source.size(); ++ source_last) {
	    const size_type target_last = target_first;
	    
	    const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	    const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);

	    const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	    
	    const std::pair<PYPPhrase::id_type, bool> source_id = model.phrase.phrase_find(phrase_source);
	    
	    if (! source_id.second) continue;
	    
	    const std::pair<PYPPhrase::id_type, bool> pair_id = model.phrase.phrase_pair_find(source_id.first, empty_id.first);

	    if (! pair_id.second) continue;
	    
	    const std::pair<logprob_type, bool> logprob_gen = model.phrase.table.prob_model(pair_id.first, logprob_type(0.0));
	    
	    if (! logprob_gen.second) continue;
	    
	    if (edges(source_first, source_last, target_first, target_last).empty())
	      agenda[span_pair.size()].push_back(span_pair);
	    
	    edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::GENERATIVE), logprob_gen.first));
	    
	    chart(source_first, source_last, target_first, target_last) += logprob_gen.first;
	    
	    chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_gen.first);
	  }
	
	if (empty_id.second && target_first < target.size())
	  for (size_type target_last = target_first + 1; target_last <= target.size(); ++ target_last) {
	    const size_type source_last = source_first;
	    
	    const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	    const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);

	    const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	    	    
	    const std::pair<PYPPhrase::id_type, bool> target_id = model.phrase.phrase_find(phrase_target);
	    
	    if (! target_id.second) continue;

	    const std::pair<PYPPhrase::id_type, bool> pair_id = model.phrase.phrase_pair_find(empty_id.first, target_id.first);
	    
	    if (! pair_id.second) continue;
	    
	    const std::pair<logprob_type, bool> logprob_gen = model.phrase.table.prob_model(pair_id.first, logprob_type(0.0));
	    
	    if (! logprob_gen.second) continue;
	    
	    if (edges(source_first, source_last, target_first, target_last).empty())
	      agenda[span_pair.size()].push_back(span_pair);
	    
	    edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::GENERATIVE), logprob_gen.first));
	    
	    chart(source_first, source_last, target_first, target_last) += logprob_gen.first;
	    
	    chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_gen.first);
	  }
	
	// phrases... is it correct?
	
	if (source_first < source.size() && target_first < target.size())
	  for (size_type source_last = source_first + 1; source_last <= source.size(); ++ source_last) {
	    const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	    
	    const std::pair<PYPPhrase::id_type, bool> source_id = model.phrase.phrase_find(phrase_source);
	  
	    if (source_id.second)
	      for (size_type target_last = target_first + 1; target_last <= target.size(); ++ target_last) {
		const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);
		
		const span_pair_type span_pair(source_first, source_last, target_first, target_last);
		
		// generative model
		const std::pair<PYPPhrase::id_type, bool> target_id = model.phrase.phrase_find(phrase_target);
		
		if (! target_id.second) continue;
		
		const std::pair<PYPPhrase::id_type, bool> pair_id = model.phrase.phrase_pair_find(source_id.first, target_id.first);
		
		if (! pair_id.second) continue;
		
		const std::pair<logprob_type, bool> logprob_gen = model.phrase.table.prob_model(pair_id.first, logprob_type(0.0));
		
		if (! logprob_gen.second) continue;
		
		if (edges(source_first, source_last, target_first, target_last).empty())
		  agenda[span_pair.size()].push_back(span_pair);
		
		edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::GENERATIVE), logprob_gen.first));
		
		chart(source_first, source_last, target_first, target_last) += logprob_gen.first;
		
		chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_gen.first);
		chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_gen.first);
	      }
	  }
      }
    
    // forward/backward to compute alpha_{source,target} and beta_{source,target}
    forward_backward(source, chart_source, alpha_source, beta_source);
    forward_backward(target, chart_target, alpha_target, beta_target);
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

    span_pairs_unique_type spans_unique;
    spans_unique.set_empty_key(span_pairs_type(span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)),
					       span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1))));
    
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

    // we will transform the pre-order stack operation into post-order operation via stack_derivation...
    // HOW?
    stack_derivation.clear();
    
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

  logprob_type logprob_fallback;
  logprob_type logprob_base;
  logprob_type logprob_epsilon;
  logprob_type logprob_term;
  logprob_type logprob_str;
  logprob_type logprob_inv;

  chart_type      chart;
  edge_chart_type edges;
  
  agenda_type     agenda;
  stack_type      stack;
  derivation_type stack_derivation;
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
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;
  
  typedef PYPPiAlign::logprob_type logprob_type;
  typedef PYPPiAlign::prob_type    prob_type;

  Task(queue_type& __mapper,
       counter_type& __reducer,
       const sentence_set_type& __sources,
       const sentence_set_type& __targets,
       derivation_set_type& __derivations,
       PYPPiAlign& __model,
       const sampler_type& __sampler,
       const logprob_type& __beam)
    : mapper(__mapper),
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
      
      {
	utils::resource start;
	
	graph.initialize(sources[pos], targets[pos], model);
	
	utils::resource end;
	
	time.initialize += end.thread_time() - start.thread_time();
      }
      
      utils::resource start;
      
      graph.forward(sources[pos], targets[pos], beam);
      
      utils::resource middle;
      
      graph.backward(sources[pos], targets[pos], derivations[pos], sampler, temperature);
      
      utils::resource end;
      
      time.forward  += middle.thread_time() - start.thread_time();
      time.backward += end.thread_time() - middle.thread_time();      
      
      reducer.increment();
    }
  }
  
  queue_type&   mapper;
  counter_type& reducer;
  
  const sentence_set_type& sources;
  const sentence_set_type& targets;
  derivation_set_type& derivations;
  
  PYPPiAlign&  model;
  sampler_type sampler;
  
  logprob_type beam;
  
  double temperature;
  
  Time time;
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

path_type output_sample_file;
path_type output_model_file;

int max_sentence_length = 40;
double beam = 1e-4;

int samples = 1;
int burns = 10;
int baby_steps = 1;
int anneal_steps = 1;
int resample_rate = 1;
int resample_iterations = 2;
bool slice_sampling = false;
bool sample_hypergraph = false;

double epsilon_prior = 1e-2;

double rule_prior_terminal = 0.1;

double rule_discount_alpha = 1.0;
double rule_discount_beta  = 1.0;
double rule_strength_shape = 1.0;
double rule_strength_rate  = 1.0;

double phrase_discount_alpha = 1.0;
double phrase_discount_beta  = 1.0;
double phrase_strength_shape = 1.0;
double phrase_strength_rate  = 1.0;

int blocks  = 64;
int threads = 1;
int debug = 0;

void options(int argc, char** argv);

size_t read_data(const path_type& path, sentence_set_type& sentences);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);

    blocks  = utils::bithack::max(blocks, 1);
    threads = utils::bithack::max(threads, 1);
    
    if (samples < 0)
      throw std::runtime_error("# of samples must be positive");
    
    if (resample_rate <= 0)
      throw std::runtime_error("resample rate must be >= 1");
    
    if (beam < 0.0 || beam > 1.0)
      throw std::runtime_error("invalid beam width");
    
    if (rule_prior_terminal <= 0.0 || rule_prior_terminal >= 1.0)
      throw std::runtime_error("invalid rule-prior for terminal");
    
    
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

    
    PYPPhrase model_phrase(1.0 / (double(source_vocab_size) * double(target_vocab_size)),
			   PYPPhrase::parameter_type(phrase_discount_alpha,
						     phrase_discount_beta,
						     phrase_strength_shape,
						     phrase_strength_rate));
    
    PYPPiAlign model(model_rule, model_phrase, PYPEpsilon(epsilon_prior));
    
    derivation_set_type derivations(sources.size());
    position_set_type positions;
    for (size_t i = 0; i != sources.size(); ++ i)
      if (! sources[i].empty() && ! targets[i].empty()
	  && (max_sentence_length <= 0
	      || (static_cast<int>(sources[i].size()) <= max_sentence_length
		  && static_cast<int>(targets[i].size()) <= max_sentence_length)))
	positions.push_back(i);
    position_set_type(positions).swap(positions);
    
    sampler_type sampler;

    // sample parameters, first...
    if (slice_sampling)
      model.slice_sample_parameters(sampler, resample_iterations);
    else
      model.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2)
      std::cerr << "rule: discount=" << model.rule.table.discount()
		<< " strength=" << model.rule.table.strength() << std::endl
		<< "fallback=" << model.phrase.table.prob(1.0) << " terminal=" << model.rule.prob_terminal() << " straight=" << model.rule.prob_straight() << " inverted=" << model.rule.prob_inverted() << std::endl
		<< "phrase: discount=" << model.phrase.table.discount()
		<< " strength=" << model.phrase.table.strength() << std::endl;
    
    Task::queue_type queue_mapper;
    Counter reducer;
    
    std::vector<Task, std::allocator<Task> > tasks(threads, Task(queue_mapper,
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
      
      position_set_type::const_iterator piter_end = positions.end();
      position_set_type::const_iterator piter = positions.begin();
      
      Time time_model;
      position_set_type mapped;
      
      size_type reduced_total = 0;
      while (piter != piter_end) {
	mapped.clear();
	
	position_set_type::const_iterator piter_last = std::min(piter + blocks, piter_end);
	for (/**/; piter != piter_last; ++ piter) {
	  const size_type pos = *piter;
	  
	  if (! derivations[pos].empty()) {
	    utils::resource start;
	    
	    derivation_type::const_reverse_iterator diter_end = derivations[pos].rend();
	    for (derivation_type::const_reverse_iterator diter = derivations[pos].rbegin(); diter != diter_end; ++ diter)
	      model.decrement(sources[pos], targets[pos], *diter, sampler);
	    
	    utils::resource end;
	    
	    time_model.decrement += end.thread_time() - start.thread_time();
	  }
	  
	  mapped.push_back(pos);
	}

	reducer.clear();
	
	position_set_type::const_iterator miter_end = mapped.end();
	for (position_set_type::const_iterator miter = mapped.begin(); miter != miter_end; ++ miter)
	  queue_mapper.push(*miter);

	reducer.wait(mapped.size());
	
	for (size_type reduced = 0; reduced != mapped.size(); ++ reduced, ++ reduced_total) {
	  if (debug) {
	    if ((reduced_total + 1) % 10000 == 0)
	      std::cerr << '.';
	    if ((reduced_total + 1) % 1000000 == 0)
	      std::cerr << '\n';
	  }
	}
	
	for (position_set_type::const_iterator miter = mapped.begin(); miter != miter_end; ++ miter) {
	  const size_type pos = *miter;

	  utils::resource start;
	  
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    model.increment(sources[pos], targets[pos], *diter, sampler, temperature);

	  utils::resource end;
	  
	  time_model.increment += end.thread_time() - start.thread_time();
	}
      }
      
      if (debug && positions.size() >= 10000 && positions.size() % 1000000 != 0)
	std::cerr << std::endl;
      
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  model.slice_sample_parameters(sampler, resample_iterations);
	else
	  model.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2)
	  std::cerr << "rule: discount=" << model.rule.table.discount()
		    << " strength=" << model.rule.table.strength() << std::endl
		    << "fallback=" << model.phrase.table.prob(1.0) << " terminal=" << model.rule.prob_terminal() << " straight=" << model.rule.prob_straight() << " inverted=" << model.rule.prob_inverted() << std::endl
		    << "phrase: discount=" << model.phrase.table.discount()
		    << " strength=" << model.phrase.table.strength() << std::endl;
      }

      if (debug >= 2) {
	Time time_end;
	for (size_type i = 0; i != tasks.size(); ++ i)
	  time_end += tasks[i].time;
	
	std::cerr << "initialize: " << (time_end.initialize - time.initialize) / tasks.size() << " seconds" << std::endl
		  << "forward: " << (time_end.forward - time.forward) / tasks.size() << " seconds" << std::endl
		  << "backward: " << (time_end.backward - time.backward) / tasks.size() << " seconds" << std::endl
		  << "increment: " << time_model.increment << " seconds" << std::endl
		  << "decrement: " << time_model.decrement << " seconds" << std::endl;
      }
      
      if (debug)
	std::cerr << "log-likelihood: " << model.log_likelihood() << std::endl;
      
      // perform model pruning..
      model.phrase.prune(debug >= 2);
      
      if (sampling && ! output_sample_file.empty()) {
	// dump derivations..!
	
	const path_type path = add_suffix(output_sample_file, "." + utils::lexical_cast<std::string>(sample_iter + 1));
	
	utils::compress_ostream os(path, 1024 * 1024);

	if (sample_hypergraph) {
	  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > stack_type;
	  typedef hypergraph_type::rule_type     rule_type;
	  typedef hypergraph_type::rule_ptr_type rule_ptr_type;

	  stack_type stack;
	  
	  hypergraph_type graph_source;
	  hypergraph_type graph_target;
	  
	  rule_type::symbol_set_type straight(2);
	  rule_type::symbol_set_type inverted(2);
	  
	  straight[0] = vocab_type::X1;
	  straight[1] = vocab_type::X2;
	  
	  inverted[0] = vocab_type::X2;
	  inverted[1] = vocab_type::X1;

	  hypergraph_type::edge_type::node_set_type tails(2);
	  
	  const rule_ptr_type rule_straight(rule_type::create(rule_type(vocab_type::X, straight)));
	  const rule_ptr_type rule_inverted(rule_type::create(rule_type(vocab_type::X, inverted)));
	  
	  for (size_type pos = 0; pos != derivations.size(); ++ pos) {
	    graph_source.clear();
	    graph_target.clear();
	    
	    if (! derivations[pos].empty()) {
	      // construct pair of hypergrpahs
	      
	      graph_source.goal = graph_source.add_node().id;
	      graph_target.goal = graph_target.add_node().id;
	      
	      stack.clear();
	      stack.push_back(0);
	      
	      derivation_type::const_iterator diter_end = derivations[pos].end();
	      for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter) {
		const hypergraph_type::id_type head = stack.back();
		stack.pop_back();
		
		if (diter->is_terminal()) {
		  hypergraph_type::edge_type& edge_source = graph_source.add_edge();
		  hypergraph_type::edge_type& edge_target = graph_target.add_edge();
		  

		  if (diter->span.source.empty())
		    edge_source.rule = rule_type::create(rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1));
		  else
		    edge_source.rule = rule_type::create(rule_type(vocab_type::X,
								   sources[pos].begin() + diter->span.source.first,
								   sources[pos].begin() + diter->span.source.last));

		  if (diter->span.target.empty())
		    edge_target.rule = rule_type::create(rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1));
		  else
		    edge_target.rule = rule_type::create(rule_type(vocab_type::X,
								   targets[pos].begin() + diter->span.target.first,
								   targets[pos].begin() + diter->span.target.last));
		  
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
	    }
	    
	    os << graph_source << " ||| " << graph_target << '\n';
	  }
	} else {
	  typedef std::vector<std::string, std::allocator<std::string> > stack_type;
	  
	  stack_type stack;
	  
	  for (size_type pos = 0; pos != derivations.size(); ++ pos) {
	    if (! derivations[pos].empty()) {
	      // we need to transform the stack-structure into tree-struct... HOW?
	    
	      stack.clear();
	      derivation_type::const_iterator diter_end = derivations[pos].end();
	      for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter) {
		if (diter->is_terminal()) {

		  if (diter->itg == PYP::GENERATIVE)
		    os << "((( "
		       << PYP::phrase_type(sources[pos].begin() + diter->span.source.first, sources[pos].begin() + diter->span.source.last)
		       << " ||| "
		       << PYP::phrase_type(targets[pos].begin() + diter->span.target.first, targets[pos].begin() + diter->span.target.last)
		       << " )))";
		  else
		    os << "{{{ "
		       << PYP::phrase_type(sources[pos].begin() + diter->span.source.first, sources[pos].begin() + diter->span.source.last)
		       << " ||| "
		       << PYP::phrase_type(targets[pos].begin() + diter->span.target.first, targets[pos].begin() + diter->span.target.last)
		       << " }}}";
		
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
	    }
	    os << '\n';
	  }
	}
      }
      
      if (sampling && ! output_model_file.empty()) {
	// output phrase-table...
	// How to generate lexicalized-reordering table...?
	//
	
	//
	// we will use the derivations to compute "alignment point", then extract lexicalized reordering score...
	//
	
	typedef PYP::phrase_type      phrase_type;
	typedef PYP::phrase_pair_type phrase_pair_type;
	typedef utils::dense_hash_map<phrase_pair_type, double, boost::hash<phrase_pair_type>, std::equal_to<phrase_pair_type>,
				      std::allocator<std::pair<const phrase_pair_type, double> > >::type phrase_pair_set_type;
	
	typedef utils::dense_hash_map<phrase_type, double, boost::hash<phrase_type>, std::equal_to<phrase_type>,
				      std::allocator<std::pair<const phrase_type, double> > >::type phrase_set_type;
	
	typedef boost::array<size_type, 5> reordering_type;
	typedef utils::dense_hash_map<phrase_pair_type, reordering_type, boost::hash<phrase_pair_type>, std::equal_to<phrase_pair_type>,
				      std::allocator<std::pair<const phrase_pair_type, reordering_type> > >::type reordering_pair_set_type;
	
	typedef utils::vector2<bool, std::allocator<bool> > matrix_type;
	
	matrix_type matrix;
	reordering_pair_set_type reorderings;
	reorderings.set_empty_key(phrase_pair_type());
	
	for (size_type pos = 0; pos != derivations.size(); ++ pos)
	  if (! derivations[pos].empty()) {
	    const size_type source_size = sources[pos].size();
	    const size_type target_size = targets[pos].size();
	    
	    matrix.clear();
	    matrix.resize(source_size + 2, target_size + 2);
	    
	    matrix(0, 0) = true;  // BOS
	    matrix(source_size + 1, target_size + 1) = true; // EOS
	    
	    derivation_type::const_iterator diter_end = derivations[pos].end();
	    for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	      if (! diter->span.source.empty() && ! diter->span.target.empty()) {
		matrix(diter->span.source.first + 1, diter->span.target.first + 1) = true;
		matrix(diter->span.source.first + 1, diter->span.target.last) = true;
		matrix(diter->span.source.last, diter->span.target.first + 1) = true;
		matrix(diter->span.source.last, diter->span.target.last) = true;
	      }
	    
	    for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	      if (! diter->span.source.empty() && ! diter->span.target.empty()) {
		const PYP::span_pair_type& span = diter->span;

		const bool connected_left_top     = matrix(span.source.first,    span.target.first);
		const bool connected_right_top    = matrix(span.source.last + 1, span.target.first);
		const bool connected_left_bottom  = matrix(span.source.first,    span.target.last + 1);
		const bool connected_right_bottom = matrix(span.source.last + 1, span.target.last + 1);
		
		const phrase_pair_type phrase_pair(sources[pos].begin() + span.source.first, sources[pos].begin() + span.source.last,
						   targets[pos].begin() + span.target.first, targets[pos].begin() + span.target.last);
		
		reordering_type& reorder = reorderings[phrase_pair];
		
		reorder[0] += 1;
		reorder[1] += (  connected_left_top && ! connected_right_top);
		reorder[2] += (! connected_left_top &&   connected_right_top);
		reorder[3] += (  connected_left_bottom && ! connected_right_bottom);
		reorder[4] += (! connected_left_bottom &&   connected_right_bottom);
	      }
	  }
	
	phrase_pair_set_type phrases;
	phrase_set_type      phrases_source;
	phrase_set_type      phrases_target;

	phrases.set_empty_key(phrase_pair_type());
	phrases_source.set_empty_key(phrase_type());
	phrases_target.set_empty_key(phrase_type());

	for (PYPPhrase::id_type id = 0; id != model.phrase.table.size(); ++ id)
	  if (! model.phrase.table[id].empty()) {
	    const PYPPhrase::id_pair_type& pair = model.phrase.phrase_pairs[id];
	    
	    const phrase_type& phrase_source = model.phrase.phrases[pair.first];
	    const phrase_type& phrase_target = model.phrase.phrases[pair.second];
	    
	    if (! phrase_source.empty() && ! phrase_target.empty()) {
	      const double prob = model.phrase.table.prob(id, 0.0);
	      
	      phrases[phrase_pair_type(phrase_source, phrase_target)] = prob;
	      phrases_source[phrase_source] += prob;
	      phrases_target[phrase_target] += prob;
	    }
	  }
	
	const path_type path = add_suffix(output_model_file, "." + utils::lexical_cast<std::string>(sample_iter + 1));
	
	utils::compress_ostream os(path, 1024 * 1024);
	os.precision(10);

	const double penalty = std::exp(1);
	
	phrase_pair_set_type::const_iterator piter_end = phrases.end();
	for (phrase_pair_set_type::const_iterator piter = phrases.begin(); piter != piter_end; ++ piter) {
	  const reordering_type& reorder = reorderings[piter->first];
	  
	  const double count = reorder[0];
	  const double count_prev_mono   = reorder[1];
	  const double count_prev_swap   = reorder[2];
	  const double count_prev_others = count - (count_prev_mono + count_prev_swap);
	  const double count_next_mono   = reorder[3];
	  const double count_next_swap   = reorder[4];
	  const double count_next_others = count - (count_next_mono + count_next_swap);
	  
	  const double prob_prev_mono   = (0.5 + count_prev_mono)   / (0.5 * 3 + count);
	  const double prob_prev_swap   = (0.5 + count_prev_swap)   / (0.5 * 3 + count);
	  const double prob_prev_others = (0.5 + count_prev_others) / (0.5 * 3 + count);
	  
	  const double prob_next_mono   = (0.5 + count_next_mono)   / (0.5 * 3 + count);
	  const double prob_next_swap   = (0.5 + count_next_swap)   / (0.5 * 3 + count);
	  const double prob_next_others = (0.5 + count_next_others) / (0.5 * 3 + count);
	  
	  os << piter->first.source << " ||| " << piter->first.target
	     << " |||"
	     << " " << (piter->second / phrases_source[piter->first.source])
	     << " " << (piter->second / phrases_target[piter->first.target])
	     << " " << piter->second
	     << " " << penalty
	     << " |||"
	     << " " << prob_prev_mono
	     << " " << prob_prev_swap
	     << " " << prob_prev_others
	     << " " << prob_next_mono
	     << " " << prob_next_swap
	     << " " << prob_next_others
	     << '\n';
	}
      }
    }
    
    for (int i = 0; i != threads; ++ i)
      queue_mapper.push(size_type(-1));
    
    workers.join_all();
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
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
    
    ("output-sample", po::value<path_type>(&output_sample_file), "output sample file")
    ("output-model",  po::value<path_type>(&output_model_file),  "output model file (or phrase table)")
    
    
    ("max-sentence-length", po::value<int>(&max_sentence_length)->default_value(max_sentence_length), "max sentence length")
    ("beam",                po::value<double>(&beam)->default_value(beam),                            "beam threshold")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("burns",               po::value<int>(&burns)->default_value(burns),                             "# of burn-ins")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    ("hypergraph",          po::bool_switch(&sample_hypergraph),                                      "dump sampled derivation in hypergraph")

    ("epsilon-prior",       po::value<double>(&epsilon_prior)->default_value(epsilon_prior),             "prior for epsilon")
    ("rule-prior-terminal", po::value<double>(&rule_prior_terminal)->default_value(rule_prior_terminal), "prior for terminal")
    
    ("rule-discount-alpha", po::value<double>(&rule_discount_alpha)->default_value(rule_discount_alpha), "discount ~ Beta(alpha,beta)")
    ("rule-discount-beta",  po::value<double>(&rule_discount_beta)->default_value(rule_discount_beta),   "discount ~ Beta(alpha,beta)")

    ("rule-strength-shape", po::value<double>(&rule_strength_shape)->default_value(rule_strength_shape), "strength ~ Gamma(shape,rate)")
    ("rule-strength-rate",  po::value<double>(&rule_strength_rate)->default_value(rule_strength_rate),   "strength ~ Gamma(shape,rate)")

    ("phrase-discount-alpha", po::value<double>(&phrase_discount_alpha)->default_value(phrase_discount_alpha), "discount ~ Beta(alpha,beta)")
    ("phrase-discount-beta",  po::value<double>(&phrase_discount_beta)->default_value(phrase_discount_beta),   "discount ~ Beta(alpha,beta)")

    ("phrase-strength-shape", po::value<double>(&phrase_strength_shape)->default_value(phrase_strength_shape), "strength ~ Gamma(shape,rate)")
    ("phrase-strength-rate",  po::value<double>(&phrase_strength_rate)->default_value(phrase_strength_rate),   "strength ~ Gamma(shape,rate)")
    
    ("blocks",  po::value<int>(&blocks),  "# of blocks")
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


