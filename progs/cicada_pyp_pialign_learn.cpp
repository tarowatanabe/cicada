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
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/indexed_map.hpp"
#include "utils/indexed_set.hpp"
#include "utils/unique_set.hpp"

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

};

// a base mearure for PYPLexicon
// we will revise this and split into two: p(target|source) and p(source|target)
//

struct LexiconModel
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

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
  
  typedef utils::dense_hash_map<word_pair_type, double, boost::hash<word_pair_type>, std::equal_to<word_pair_type>,
				std::allocator<std::pair<const word_pair_type, double> > >::type table_type;
  
  typedef boost::filesystem::path path_type;
  
  LexiconModel(const double __smooth=1e-7)
    : table(), smooth(__smooth)
  {
    table.set_empty_key(word_pair_type());
  }
  
  LexiconModel(const path_type& path)
    : table(), smooth()
  {
    table.set_empty_key(word_pair_type());
    
    open(path);
  }
  
  void open(const path_type& path)
  {
    typedef utils::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_set_type;
        
    typedef boost::fusion::tuple<std::string, std::string, double > lexicon_parsed_type;
    typedef boost::spirit::istream_iterator iterator_type;
    
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    qi::rule<iterator_type, std::string(), standard::blank_type>         word;
    qi::rule<iterator_type, lexicon_parsed_type(), standard::blank_type> parser; 
    
    word   %= qi::lexeme[+(standard::char_ - standard::space)];
    parser %= word >> word >> qi::double_ >> (qi::eol | qi::eoi);
    
    word_set_type words;
    words.set_empty_key(word_type());
    table.clear();
    
    utils::compress_istream is(path, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    lexicon_parsed_type lexicon_parsed;

    iterator_type iter(is);
    iterator_type iter_end;

    while (iter != iter_end) {
      boost::fusion::get<0>(lexicon_parsed).clear();
      boost::fusion::get<1>(lexicon_parsed).clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, lexicon_parsed))
	if (iter != iter_end)
	  throw std::runtime_error("global lexicon parsing failed");
      
      const word_type target(boost::fusion::get<0>(lexicon_parsed));
      const word_type source(boost::fusion::get<1>(lexicon_parsed));
      const double&   prob(boost::fusion::get<2>(lexicon_parsed));
      
      table[word_pair_type(source, target)] = prob;
      
      words.insert(target);
    }
    
    smooth = 1.0 / words.size();
  }
  
  double operator()(const word_type& source, const word_type& target) const
  {
    table_type::const_iterator iter = table.find(word_pair_type(source, target));
    
    return (iter != table.end() ? iter->second : smooth);
  }
  
  table_type table;
  double smooth;
};

struct PYPLexicon
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::phrase_type phrase_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  
  typedef utils::restaurant<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type > > table_type;
  typedef utils::alloc_vector<table_type, std::allocator<table_type> > table_set_type;

  typedef LexiconModel lexicon_type;
  
  PYPLexicon(const lexicon_type& __lexicon_source_target,
	     const lexicon_type& __lexicon_target_source,
	     const parameter_type& __parameter)
    : lexicon_source_target(&__lexicon_source_target),
      lexicon_target_source(&__lexicon_target_source),
      tables_source_target(),
      tables_target_source(),
      parameter_source_target(__parameter),
      parameter_target_source(__parameter) {}
  
  template <typename Sampler>
  void increment(const phrase_type& source, const phrase_type& target, Sampler& sampler, const double temperature=1.0)
  {
    if (source.empty() && target.empty())
      throw std::runtime_error("invalid phrase pair");
    
    increment(source, target, tables_source_target, *lexicon_source_target, parameter_source_target, sampler, temperature);
    increment(target, source, tables_target_source, *lexicon_target_source, parameter_target_source, sampler, temperature);
  }
  
  template <typename Sampler>
  void increment(const phrase_type& source,
		 const phrase_type& target,
		 table_set_type& tables,
		 const lexicon_type& lexicon,
		 const parameter_type& parameter,
		 Sampler& sampler,
		 const double temperature=1.0)
  {
    if (source.empty()) {
      if (! tables.exists(vocab_type::EPSILON.id()))
	tables[vocab_type::EPSILON.id()] = table_type(parameter.discount, parameter.strength);
      
      table_type& table = tables[vocab_type::EPSILON.id()];
      
      phrase_type::const_iterator titer_end = target.end();
      for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	table.increment(*titer, lexicon(vocab_type::EPSILON, *titer), sampler, temperature);
    } else if (! target.empty()) {
      phrase_type::const_iterator siter_end = source.end();
      for (phrase_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter) {
	if (! tables.exists(siter->id()))
	  tables[siter->id()] = table_type(parameter.discount, parameter.strength);
	
	table_type& table = tables[siter->id()];
	
	phrase_type::const_iterator titer_end = target.end();
	for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	  table.increment(*titer, lexicon(*siter, *titer), sampler, temperature);
      }
    }
  }

  template <typename Sampler>
  void decrement(const phrase_type& source, const phrase_type& target, Sampler& sampler)
  {
    if (source.empty() && target.empty())
      throw std::runtime_error("invalid phrase pair");
    
    decrement(source, target, tables_source_target, sampler);
    decrement(target, source, tables_target_source, sampler);
  }
  
  template <typename Sampler>
  void decrement(const phrase_type& source, const phrase_type& target, table_set_type& tables, Sampler& sampler)
  {
    if (source.empty()) {
      table_type& table = tables[vocab_type::EPSILON.id()];
      
      phrase_type::const_iterator titer_end = target.end();
      for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	table.decrement(*titer, sampler);
    } else if (! target.empty()) {
      phrase_type::const_iterator siter_end = source.end();
      for (phrase_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter) {
	table_type& table = tables[siter->id()];
	
	phrase_type::const_iterator titer_end = target.end();
	for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	  table.decrement(*titer, sampler);
      }
    }
  }
  
  double prob_source_target(const word_type& source, const word_type& target) const
  {
    return prob(source, target, tables_source_target, *lexicon_source_target);
  }
  
  double prob_target_source(const word_type& target, const word_type& source) const
  {
    return prob(target, source, tables_target_source, *lexicon_target_source);
  }
  
  double prob(const word_type& source, const word_type& target, const table_set_type& tables, const lexicon_type& lexicon) const
  {
    if (! tables.exists(source.id()))
      return lexicon(source, target);
    else
      return tables[source.id()].prob(target, lexicon(source, target));
  }
  
  double prob(const phrase_type& source, const phrase_type& target) const
  {
    if (source.empty() && target.empty())
      throw std::runtime_error("invalid phrase pair");
    
    const double logprob_source_target = logprob(source, target, tables_source_target, *lexicon_source_target);
    const double logprob_target_source = logprob(target, source, tables_target_source, *lexicon_target_source);
    
    if (source.empty() || target.empty())
      return std::exp(logprob_source_target + logprob_target_source);
    else
      return std::exp((logprob_source_target + logprob_target_source) * 0.5);
  }

  double prob_source_target(const phrase_type& source, const phrase_type& target) const
  {
    return std::exp(logprob(source, target, tables_source_target, *lexicon_source_target));
  }
  
  double prob_target_source(const phrase_type& target, const phrase_type& source) const
  {
    return std::exp(logprob(target, source, tables_target_source, *lexicon_target_source));
  }

  double logprob(const phrase_type& source, const phrase_type& target, const table_set_type& tables, const lexicon_type& lexicon) const
  {
    double lp = 0.0;

    if (source.empty()) {
      if (! tables.exists(vocab_type::EPSILON.id())) {
	phrase_type::const_iterator titer_end = target.end();
	for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	  lp += std::log(lexicon(vocab_type::EPSILON, *titer));
      } else {
	const table_type& table = tables[vocab_type::EPSILON.id()];
	
	phrase_type::const_iterator titer_end = target.end();
	for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	  lp += std::log(table.prob(*titer, lexicon(vocab_type::EPSILON, *titer)));
      }
    } else if (! target.empty()) {
      phrase_type::const_iterator titer_end = target.end();
      for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) {
	double sum = 0.0;
	phrase_type::const_iterator siter_end = source.end();
	for (phrase_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter) {
	  if (! tables.exists(siter->id()))
	    sum += lexicon(*siter, *titer);
	  else
	    sum += tables[siter->id()].prob(*titer, lexicon(*siter, *titer));
	}
	
	lp += std::log(sum);
      }
    }

    return lp;
  }
  
  double log_likelihood() const
  {
    double logprob = parameter_source_target.log_likelihood() + parameter_target_source.log_likelihood();
    
    table_set_type::const_iterator siter_end = tables_source_target.end();
    for (table_set_type::const_iterator siter = tables_source_target.begin(); siter != siter_end; ++ siter)
      if (*siter)
	logprob += (*siter)->log_likelihood();
    
    table_set_type::const_iterator titer_end = tables_target_source.end();
    for (table_set_type::const_iterator titer = tables_target_source.begin(); titer != titer_end; ++ titer)
      if (*titer)
	logprob += (*titer)->log_likelihood();
    
    return logprob;
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    for (int iter = 0; iter != num_loop; ++ iter) {
      parameter_source_target.discount = sample_discount(tables_source_target.begin(), tables_source_target.end(), sampler, parameter_source_target);
      parameter_source_target.strength = sample_strength(tables_source_target.begin(), tables_source_target.end(), sampler, parameter_source_target);

      parameter_target_source.discount = sample_discount(tables_target_source.begin(), tables_target_source.end(), sampler, parameter_target_source);
      parameter_target_source.strength = sample_strength(tables_target_source.begin(), tables_target_source.end(), sampler, parameter_target_source);
    }
    
    assign_parameters(tables_source_target.begin(), tables_source_target.end(), parameter_source_target);
    assign_parameters(tables_target_source.begin(), tables_target_source.end(), parameter_target_source);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    // currently, we do not support slice sampling...
    sample_parameters(sampler, num_loop, num_iterations);
  }

  template <typename Iterator, typename Sampler>
  double sample_strength(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double x = 0.0;
    double y = 0.0;
    
    for (/**/; first != last; ++ first) 
      if (*first) {
	x += (*first)->sample_log_x(sampler, param.discount, param.strength);
	y += (*first)->sample_y(sampler, param.discount, param.strength);
      }
    
    return sampler.gamma(param.strength_shape + y, param.strength_rate - x);
  }
  
  template <typename Iterator, typename Sampler>
  double sample_discount(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double y = 0.0;
    double z = 0.0;
    
    for (/**/; first != last; ++ first) 
      if (*first) {
	y += (*first)->sample_y_inv(sampler, param.discount, param.strength);
	z += (*first)->sample_z_inv(sampler, param.discount, param.strength);
      }
    
    return sampler.beta(param.discount_alpha + y, param.discount_beta + z);
  }

  template <typename Iterator>
  void assign_parameters(Iterator first, Iterator last, const parameter_type& param)
  {
    for (/**/; first != last; ++ first)
      if (*first) {
	(*first)->discount() = param.discount;
	(*first)->strength() = param.strength;
      }
  }
  
  const lexicon_type* lexicon_source_target;
  const lexicon_type* lexicon_target_source;
  
  table_set_type tables_source_target;
  table_set_type tables_target_source;
  parameter_type parameter_source_target;
  parameter_type parameter_target_source;
};

struct LengthModel
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::array_power2<double, 32, std::allocator<double> > cache_type;
  
  LengthModel(const double& __lambda,
	      const double& __strength_shape,
	      const double& __strength_rate)
    : lambda(__lambda),
      strength_shape(__strength_shape),
      strength_rate(__strength_rate),
      cache()
  {
    assign(__lambda);
  }

  double prob(const size_type size) const
  {
    return std::exp(logprob(size));
  }
  
  double logprob(const size_type size) const
  {
    if (size < cache.size())
      return cache[size];
    else
      return utils::mathop::log_poisson(size, lambda);
  }
  
  void assign(const double __lambda)
  {
    lambda = __lambda;
    
    cache[0] = 0.0; // this is wrong... but counter the epsilon problem.
    for (size_type size = 1; size != cache.size(); ++ size)
      cache[size] = utils::mathop::log_poisson(size, lambda);
  }
  
  template <typename Sampler>
  void sample_parameters(const double& alpha, const double& beta, Sampler& sampler)
  {
    assign(sampler.gamma(alpha + strength_shape, beta + strength_rate));
  }
  
  double lambda;
  double strength_shape;
  double strength_rate;

  cache_type cache;
};

// length prior...
struct PYPLength
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::phrase_type phrase_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef std::vector<size_type, std::allocator<size_type> > count_set_type;

  PYPLength(const LengthModel& __length_source,
	    const LengthModel& __length_target)
    : length_source(__length_source),
      length_target(__length_target),
      counts_source(),
      counts_target()
  {}

  template <typename Sampler>
  void increment(const phrase_type& source, const phrase_type& target, Sampler& sampler, const double temperature=1.0)
  {
    if (! source.empty()) {
      if (source.size() >= counts_source.size())
	counts_source.resize(source.size() + 1, 0);
      
      ++ counts_source[source.size()];
    }

    if (! target.empty()) {
      if (target.size() >= counts_target.size())
	counts_target.resize(target.size() + 1, 0);
      
      ++ counts_target[target.size()];
    }
  }

  template <typename Sampler>
  void decrement(const phrase_type& source, const phrase_type& target, Sampler& sampler)
  {
    if (! source.empty()) {
      if (source.size() >= counts_source.size() || ! counts_source[source.size()])
	throw std::runtime_error("invalid decrment");
      
      -- counts_source[source.size()];
    }

    if (! target.empty()) {
      if (target.size() >= counts_target.size() || ! counts_target[target.size()])
	throw std::runtime_error("invalid decrment");
      
      -- counts_target[target.size()];
    }
  }
  
  double prob(const phrase_type& source, const phrase_type& target) const
  {
    if (source.empty() && target.empty())
      throw std::runtime_error("invalid phrase pair");
    
    return std::exp(length_source.logprob(source.size()) + length_target.logprob(target.size()));
  }
  
  logprob_type logprob(const phrase_type& source, const phrase_type& target) const
  {
    if (source.empty() && target.empty())
      throw std::runtime_error("invalid phrase pair");
    
    return cicada::semiring::traits<logprob_type>::exp(length_source.logprob(source.size()) + length_target.logprob(target.size()));
  }

  double log_likelihood() const
  {
    double logprob = 0.0;
    
    for (size_type source = 1; source < counts_source.size(); ++ source)
      logprob += length_source.logprob(source) * counts_source[source];
    
    for (size_type target = 1; target < counts_target.size(); ++ target)
      logprob += length_target.logprob(target) * counts_target[target];
    
    return logprob;
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler)
  {
    double source_alpha = 0.0;
    double source_beta = 0.0;
    double target_alpha = 0.0;
    double target_beta = 0.0;
    
    for (size_type source = 1; source < counts_source.size(); ++ source) {
      source_alpha += source * counts_source[source];
      source_beta += counts_source[source];
    }

    for (size_type target = 1; target < counts_target.size(); ++ target) {
      target_alpha += target * counts_target[target];
      target_beta += counts_target[target];
    }

    length_source.sample_parameters(source_alpha, source_beta, sampler);
    length_target.sample_parameters(target_alpha, target_beta, sampler);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler)
  {
    sample_parameters(sampler);
  }
  
  LengthModel length_source;
  LengthModel length_target;

  count_set_type counts_source;
  count_set_type counts_target;
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
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::pyp_parameter parameter_type;
  typedef utils::restaurant_vector<> table_type;
    
  PYPRule(const parameter_type& parameter)
    : p0(1.0 / 3), counts0(0), table(parameter) {}

  template <typename Sampler>
  void increment(const rule_type& rule, Sampler& sampler, const double temperature=1.0)
  {
    const itg_type itg = (rule.is_terminal() ? PYP::TERMINAL : (rule.is_straight() ? PYP::STRAIGHT : PYP::INVERTED));
    
    if (table.increment(itg, p0, sampler, temperature))
      ++ counts0;
  }
  
  template <typename Sampler>
  void decrement(const rule_type& rule, Sampler& sampler)
  {
    const itg_type itg = (rule.is_terminal() ? PYP::TERMINAL : (rule.is_straight() ? PYP::STRAIGHT : PYP::INVERTED));
    
    if (table.decrement(itg, sampler))
      -- counts0;
  }
  
  double prob(const rule_type& rule) const
  {
    const itg_type itg = (rule.is_terminal() ? PYP::TERMINAL : (rule.is_straight() ? PYP::STRAIGHT : PYP::INVERTED));
    
    return table.prob(itg, p0);
  }
  
  double prob_terminal() const
  {
    return table.prob(PYP::TERMINAL, p0);
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
    return table.log_likelihood() + std::log(p0) * counts0;
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
  
  double     p0;
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
  
  typedef utils::restaurant<phrase_pair_type, boost::hash<phrase_pair_type>, std::equal_to<phrase_pair_type>, std::allocator<phrase_pair_type > > table_type;
  
  PYPPhrase(const PYPLexicon& __lexicon,
	    const PYPLength&  __length,
	    const parameter_type& parameter)
    : lexicon(__lexicon),
      length(__length),
      table(parameter) {}

  template <typename Sampler>
  void increment_existing(const phrase_pair_type& phrase_pair, const bool leaf, Sampler& sampler, const double temperature=1.0)
  {
    if (table.increment_existing(phrase_pair, sampler)) {
      length.increment(phrase_pair.source, phrase_pair.target, sampler, temperature);

      lexicon.increment(phrase_pair.source, phrase_pair.target, sampler, temperature);
    }
  }

  template <typename Sampler>
  void increment_new(const phrase_pair_type& phrase_pair, const bool leaf, Sampler& sampler, const double temperature=1.0)
  {
    if (table.increment_new(phrase_pair, sampler)) {
      length.increment(phrase_pair.source, phrase_pair.target, sampler, temperature);
      
      lexicon.increment(phrase_pair.source, phrase_pair.target, sampler, temperature);
    }
  }
  
  
  template <typename Sampler>
  void decrement(const phrase_pair_type& phrase_pair, const bool leaf, Sampler& sampler)
  {
    if (table.decrement(phrase_pair, sampler)) {
      length.decrement(phrase_pair.source, phrase_pair.target, sampler);
      
      lexicon.decrement(phrase_pair.source, phrase_pair.target, sampler);
    }
  }
  
  template <typename LogProb>
  LogProb prob(const LogProb& base) const
  {
    return table.prob(base);
  }

  template <typename LogProb>
  std::pair<LogProb, bool> prob(const phrase_type& source, const phrase_type& target, const LogProb& base) const
  {
    return table.prob_model(phrase_pair_type(source, target), base);
  }
  
  double log_likelihood() const
  {
    //return lexicon.log_likelihood() + length.log_likelihood() + table.log_likelihood();
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
    length.sample_parameters(sampler);
    
    lexicon.sample_parameters(sampler, num_loop, num_iterations);
    
    table.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    length.slice_sample_parameters(sampler);
    
    lexicon.slice_sample_parameters(sampler, num_loop, num_iterations);
    
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  PYPLexicon lexicon;
  PYPLength  length;
  
  table_type table;
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
  
  PYPPiAlign(const PYPRule&   __rule,
	     const PYPPhrase& __phrase)
    : rule(__rule), phrase(__phrase) {}

  template <typename Sampler>
  void increment(const sentence_type& source, const sentence_type& target, const rule_type& r, Sampler& sampler, const double temperature)
  {
    rule.increment(r, sampler, temperature);

    const phrase_pair_type phrase_pair(source.begin() + r.span.source.first, source.begin() + r.span.source.last,
				       target.begin() + r.span.target.first, target.begin() + r.span.target.last);
    
    if (r.itg == PYP::GENERATIVE)
      phrase.increment_existing(phrase_pair, r.is_terminal(), sampler, temperature);
    else
      phrase.increment_new(phrase_pair, r.is_terminal(), sampler, temperature);
  }

  template <typename Sampler>
  void decrement(const sentence_type& source, const sentence_type& target, const rule_type& r, Sampler& sampler)
  {
    rule.decrement(r, sampler);
    
    phrase.decrement(phrase_pair_type(source.begin() + r.span.source.first, source.begin() + r.span.source.last,
				      target.begin() + r.span.target.first, target.begin() + r.span.target.last),
		     r.is_terminal(), 
		     sampler);
  }

  double log_likelihood() const
  {
    return rule.log_likelihood() + phrase.log_likelihood();
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    rule.sample_parameters(sampler, num_loop, num_iterations);

    phrase.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    rule.slice_sample_parameters(sampler, num_loop, num_iterations);
    
    phrase.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  PYPRule   rule;
  PYPPhrase phrase;
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
  
  typedef utils::chart<logprob_type, std::allocator<logprob_type> > model1_type;
  typedef std::vector<model1_type, std::allocator<model1_type> > model1_chart_type;
  
  typedef utils::chart<logprob_type, std::allocator<logprob_type> > chart_mono_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > alpha_type;
  typedef std::vector<logprob_type, std::allocator<logprob_type> > beta_type;

  typedef std::pair<logprob_type, span_pair_type> score_span_pair_type;
  typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;
  
  void initialize(const sentence_type& source,
		  const sentence_type& target,
		  const PYPPiAlign& model,
		  const size_type max_length)
  {
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
    
    base.clear();
    base.reserve(source.size() + 1, target.size() + 1);
    base.resize(source.size() + 1, target.size() + 1, logprob_type(1));

    base_source.clear();
    base_source.reserve(source.size() + 1, target.size() + 1);
    base_source.resize(source.size() + 1, target.size() + 1, logprob_type(1));

    base_target.clear();
    base_target.reserve(source.size() + 1, target.size() + 1);
    base_target.resize(source.size() + 1, target.size() + 1, logprob_type(1));
    
    model1_source.clear();
    model1_target.clear();
    model1_source.reserve(target.size());
    model1_target.reserve(source.size());
    model1_source.resize(target.size(), model1_type(source.size() + 1));
    model1_target.resize(source.size(), model1_type(target.size() + 1));
    
    epsilon_source.clear();
    epsilon_target.clear();
    epsilon_source.reserve(target.size());
    epsilon_target.reserve(source.size());
    epsilon_source.resize(target.size());
    epsilon_target.resize(source.size());

    // initialize model1 probabilities...
    //std::cerr << "model1 probabilities" << std::endl;
    
    for (size_type target_pos = 0; target_pos != target.size(); ++ target_pos) {
      epsilon_source[target_pos] = model.phrase.lexicon.prob_source_target(vocab_type::EPSILON, target[target_pos]);
      
      for (size_type source_first = 0; source_first != source.size(); ++ source_first) {
	double sum = 0.0;
	for (size_type source_last = source_first + 1; source_last <= source.size(); ++ source_last) {
	  sum += model.phrase.lexicon.prob_source_target(source[source_last - 1], target[target_pos]);
	  model1_source[target_pos](source_first, source_last) = sum;
	}
      }
    }
    
    for (size_type source_pos = 0; source_pos != source.size(); ++ source_pos) {
      epsilon_target[source_pos] = model.phrase.lexicon.prob_target_source(vocab_type::EPSILON, source[source_pos]);
      
      for (size_type target_first = 0; target_first != target.size(); ++ target_first) {
	double sum = 0.0;
	for (size_type target_last = target_first + 1; target_last <= target.size(); ++ target_last) {
	  sum += model.phrase.lexicon.prob_target_source(target[target_last - 1], source[source_pos]);
	  model1_target[source_pos](target_first, target_last) = sum;
	}
      }
    }

    //std::cerr << "initialize chart" << std::endl;
    
    const logprob_type logprob_fallback = model.phrase.prob(logprob_type(1.0));
    const logprob_type logprob_term = model.rule.prob_terminal() * logprob_fallback;
    
    // initialize base...
    // actually, we should be carefull with the chart assignment, since we need to consult the base of phrase-model + rule-model
    for (size_type source_first = 0; source_first != source.size(); ++ source_first)
      for (size_type target_first = 0; target_first != target.size(); ++ target_first) {
		
	// epsilons.. 
	for (size_type source_last = source_first + 1; source_last <= source.size(); ++ source_last) {
	  const size_type target_last = target_first;
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	  const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);

	  const logprob_type loglength = model.phrase.length.logprob(phrase_source, phrase_target);
	  
	  logprob_type& logbase = base(source_first, source_last, target_first, target_last);
	  logbase = base(source_first, source_last - 1, target_first, target_last) * epsilon_target[source_last - 1];
	  
	  const std::pair<logprob_type, bool> logprob_gen = model.phrase.prob(phrase_source, phrase_target, logprob_type(0.0));
	  const logprob_type logprob_base = logprob_term * logbase * loglength;

	  logprob_type& logprob_chart = chart(source_first, source_last, target_first, target_last);
	  
	  if (logprob_gen.second) {
	    edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::GENERATIVE), logprob_gen.first));
	    logprob_chart += logprob_gen.first;

	    chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_gen.first);
	  }
	  
	  if (source_last - source_first <= max_length) {
	    edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::BASE), logprob_base));
	    logprob_chart += logprob_base;
	    
	    chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_base);
	  }
	  
	  if (! edges(source_first, source_last, target_first, target_last).empty())
	    agenda[span_pair.size()].push_back(span_pair);
	}
	
	for (size_type target_last = target_first + 1; target_last <= target.size(); ++ target_last) {
	  const size_type source_last = source_first;
	  const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	  const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	  const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);
	  
	  const logprob_type loglength = model.phrase.length.logprob(phrase_source, phrase_target);
	  
	  logprob_type& logbase = base(source_first, source_last, target_first, target_last);
	  logbase = base(source_first, source_last, target_first, target_last - 1) * epsilon_source[target_last - 1];
	  
	  const std::pair<logprob_type, bool> logprob_gen = model.phrase.prob(phrase_source, phrase_target, logprob_type(0.0));
	  const logprob_type logprob_base = logprob_term * logbase * loglength;

	  logprob_type& logprob_chart = chart(source_first, source_last, target_first, target_last);
	  
	  if (logprob_gen.second) {
	    edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::GENERATIVE), logprob_gen.first));
	    logprob_chart += logprob_gen.first;
	    
	    chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_gen.first);
	  }
	  
	  if (target_last - target_first <= max_length) {
	    edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::BASE), logprob_base));
	    logprob_chart += logprob_base;
	    
	    chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_base);
	  }
	  
	  if (! edges(source_first, source_last, target_first, target_last).empty())
	    agenda[span_pair.size()].push_back(span_pair);
	}
	
	// phrases... is it correct?
	for (size_type source_last = source_first + 1; source_last <= source.size(); ++ source_last)
	  for (size_type target_last = target_first + 1; target_last <= target.size(); ++ target_last) {
	    const span_pair_type span_pair(source_first, source_last, target_first, target_last);
	    const phrase_type phrase_source(source.begin() + source_first, source.begin() + source_last);
	    const phrase_type phrase_target(target.begin() + target_first, target.begin() + target_last);
	    
	    logprob_type& logbase_source = base_source(source_first, source_last, target_first, target_last);
	    logprob_type& logbase_target = base_target(source_first, source_last, target_first, target_last);
	    
	    logbase_source = (base_source(source_first, source_last, target_first, target_last - 1)
			      * model1_source[target_last - 1](source_first, source_last));
	    
	    logbase_target = (base_target(source_first, source_last - 1, target_first, target_last)
			      * model1_target[source_last - 1](target_first, target_last));
	    
	    logprob_type& logbase = base(source_first, source_last, target_first, target_last);
	    logbase = cicada::semiring::traits<logprob_type>::exp(cicada::semiring::log(logbase_source) * 0.5
								  + cicada::semiring::log(logbase_target) * 0.5);
	    
	    const logprob_type loglength = model.phrase.length.logprob(phrase_source, phrase_target);
	    
	    const std::pair<logprob_type, bool> logprob_gen = model.phrase.prob(phrase_source, phrase_target, logprob_type(0.0));
	    const logprob_type logprob_base = logprob_term * logbase * loglength;

	    logprob_type& logprob_chart = chart(source_first, source_last, target_first, target_last);
	    
	    if (logprob_gen.second) {
	      edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::GENERATIVE), logprob_gen.first));
	      logprob_chart += logprob_gen.first;
	      
	      chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_gen.first);
	      chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_gen.first);
	    }
	    
	    if (source_last - source_first <= max_length && target_last - target_first <= max_length) {
	      edges(source_first, source_last, target_first, target_last).push_back(edge_type(rule_type(span_pair, PYP::BASE), logprob_base));
	      logprob_chart += logprob_base;
	      
	      chart_source(source_first, source_last) = std::max(chart_source(source_first, source_last), logprob_base);
	      chart_target(target_first, target_last) = std::max(chart_target(target_first, target_last), logprob_base);
	    }
	    
	    if (! edges(source_first, source_last, target_first, target_last).empty())
	      agenda[span_pair.size()].push_back(span_pair);
	  }
      }
    
    // forward/backward to compute alpha_{source,target} and beta_{source,target}
    forward_backward(source, chart_source, alpha_source, beta_source);
    forward_backward(target, chart_target, alpha_target, beta_target);
  }
  
  void forward_backward(const sentence_type& sentence, const chart_mono_type& chart, alpha_type& alpha, beta_type& beta)
  {
    // forward...
    alpha[0] = 1.0;
    for (size_type last = 1; last <= sentence.size(); ++ last)
      for (size_type first = 0; first != last; ++ first)
	alpha[last] = std::max(alpha[last], alpha[first] * chart(first, last));
    
    // backward...
    beta[sentence.size()] = 1.0;
    for (difference_type first = sentence.size() - 1; first >= 0; -- first)
      for (size_type last = first + 1; last <= sentence.size(); ++ last)
	beta[first] = std::max(beta[first], chart(first, last) * beta[last]);
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
  logprob_type forward(const sentence_type& source,
		       const sentence_type& target,
		       const PYPPiAlign& model,
		       const logprob_type beam,
		       const size_type max_length)
  {
    //std::cerr << "initialize" << std::endl;
    
    initialize(source, target, model, max_length);
    
    //std::cerr << "forward" << std::endl;

    span_pairs_unique_type spans_unique;
    spans_unique.set_empty_key(span_pairs_type(span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)),
					       span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1))));

    const logprob_type logprob_fallback = model.phrase.prob(logprob_type(1.0));
    const logprob_type logprob_str = model.rule.prob_straight() * logprob_fallback;
    const logprob_type logprob_inv = model.rule.prob_inverted() * logprob_fallback;
    
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
	
	// erase spans in hiter_begin to hiter
	for (heap_type::iterator iter = hiter_begin; iter != hiter; ++ iter) {
	  const span_pair_type& span_pair = iter->second;
	  
	  edges(span_pair.source.first, span_pair.source.last, span_pair.target.first, span_pair.target.last).clear();
	}
	
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

    return chart(0, source.size(), 0, target.size());
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

  chart_type chart;
  edge_chart_type edges;
  
  agenda_type agenda;
  stack_type stack;
  heap_type  heap;

  chart_mono_type chart_source;
  chart_mono_type chart_target;
  
  alpha_type alpha_source;
  alpha_type alpha_target;
  beta_type  beta_source;
  beta_type  beta_target;
  
  base_type base;
  base_type base_source;
  base_type base_target;
  
  model1_chart_type model1_source;
  model1_chart_type model1_target;

  prob_set_type epsilon_source;
  prob_set_type epsilon_target;
  
  prob_set_type probs;
};


typedef boost::filesystem::path path_type;
typedef utils::sampler<boost::mt19937> sampler_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;

typedef PYP::size_type size_type;
typedef PYPGraph::derivation_type derivation_type;
typedef std::vector<derivation_type, std::allocator<derivation_type> > derivation_set_type;
typedef std::vector<size_type, std::allocator<size_type> > position_set_type;

struct Task
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;  
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;
  
  typedef PYPPiAlign::logprob_type logprob_type;
  typedef PYPPiAlign::prob_type    prob_type;

  Task(queue_type& __mapper,
       queue_type& __reducer,
       const sentence_set_type& __sources,
       const sentence_set_type& __targets,
       derivation_set_type& __derivations,
       derivation_set_type& __derivations_prev,
       const PYPPiAlign& __model,
       const sampler_type& __sampler,
       const logprob_type& __beam,
       const int& __max_length)
    : mapper(__mapper),
      reducer(__reducer),
      sources(__sources),
      targets(__targets),
      derivations(__derivations),
      derivations_prev(__derivations_prev),
      model(__model),
      sampler(__sampler),
      beam(__beam),
      max_length(__max_length) {}

  void operator()()
  {
    PYPGraph graph;
    
    size_type pos;
    
    for (;;) {
      mapper.pop(pos);
      
      if (pos == size_type(-1)) break;
	
      derivations_prev[pos] = derivations[pos];

      // decrement model...
      if (! derivations_prev[pos].empty()) {
	derivation_type::const_iterator diter_end = derivations_prev[pos].end();
	for (derivation_type::const_iterator diter = derivations_prev[pos].begin(); diter != diter_end; ++ diter)
	  model.decrement(sources[pos], targets[pos], *diter, sampler);
      }
      
      graph.forward(sources[pos], targets[pos], model, beam, max_length);
      
      graph.backward(sources[pos], targets[pos], derivations[pos], sampler, temperature);
      
      // increment model...
      derivation_type::const_iterator diter_end = derivations[pos].end();
      for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	model.increment(sources[pos], targets[pos], *diter, sampler, temperature);
      
      reducer.push(pos);
    }
  }
  
  queue_type& mapper;
  queue_type& reducer;
  
  const sentence_set_type& sources;
  const sentence_set_type& targets;
  derivation_set_type& derivations;
  derivation_set_type& derivations_prev;
  
  PYPPiAlign  model;
  sampler_type sampler;
  
  logprob_type beam;
  int max_length;
  
  double temperature;
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

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

int max_phrase_length = 7;
int max_sentence_length = 40;
double beam = 1e-7;

int samples = 1;
int burns = 10;
int baby_steps = 1;
int anneal_steps = 1;
int resample_rate = 1;
int resample_iterations = 2;
bool slice_sampling = false;
bool sample_hypergraph = false;

double rule_discount = 0.9;
double rule_strength = 1;

double rule_discount_prior_alpha = 1.0;
double rule_discount_prior_beta  = 1.0;
double rule_strength_prior_shape = 1.0;
double rule_strength_prior_rate  = 1.0;

double phrase_discount = 0.9;
double phrase_strength = 1;

double phrase_discount_prior_alpha = 1.0;
double phrase_discount_prior_beta  = 1.0;
double phrase_strength_prior_shape = 1.0;
double phrase_strength_prior_rate  = 1.0;

double lexicon_discount = 0.9;
double lexicon_strength = 1;

double lexicon_discount_prior_alpha = 1.0;
double lexicon_discount_prior_beta  = 1.0;
double lexicon_strength_prior_shape = 1.0;
double lexicon_strength_prior_rate  = 1.0;

double lambda_source = 2.0;
double lambda_target = 2.0;
double lambda_shape = 1e-2;
double lambda_rate  = 1e+7;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);

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

    LexiconModel lexicon_source_target(1.0 / target_vocab_size);
    LexiconModel lexicon_target_source(1.0 / source_vocab_size);

    {
      boost::thread_group workers_open;
      
      if (! lexicon_source_target_file.empty())
	workers_open.add_thread(new boost::thread(boost::bind(&LexiconModel::open,
							      boost::ref(lexicon_source_target),
							      boost::cref(lexicon_source_target_file))));
      if (! lexicon_target_source_file.empty())
	workers_open.add_thread(new boost::thread(boost::bind(&LexiconModel::open,
							      boost::ref(lexicon_target_source),
							      boost::cref(lexicon_target_source_file))));    
      
      workers_open.join_all();
    }

    PYPRule model_rule(PYPRule::parameter_type(rule_discount,
					       rule_strength,
					       rule_discount_prior_alpha,
					       rule_discount_prior_beta,
					       rule_strength_prior_shape,
					       rule_strength_prior_rate));

    PYPLexicon model_lexicon(lexicon_source_target,
			     lexicon_target_source,
			     PYPLexicon::parameter_type(lexicon_discount,
							lexicon_strength,
							lexicon_discount_prior_alpha,
							lexicon_discount_prior_beta,
							lexicon_strength_prior_shape,
							lexicon_strength_prior_rate));

    PYPLength model_length(LengthModel(lambda_source, lambda_shape, lambda_rate),
			   LengthModel(lambda_target, lambda_shape, lambda_rate));
    
    PYPPhrase model_phrase(model_lexicon,
			   model_length,
			   PYPPhrase::parameter_type(phrase_discount,
						     phrase_strength,
						     phrase_discount_prior_alpha,
						     phrase_discount_prior_beta,
						     phrase_strength_prior_shape,
						     phrase_strength_prior_rate));
    
    PYPPiAlign model(model_rule, model_phrase);
    
    derivation_set_type derivations(sources.size());
    derivation_set_type derivations_prev(sources.size());
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
		<< "terminal=" << model.rule.prob_terminal() << " straight=" << model.rule.prob_straight() << " inverted=" << model.rule.prob_inverted() << std::endl
		<< "phrase: discount=" << model.phrase.table.discount()
		<< " strength=" << model.phrase.table.strength() << std::endl
		<< "lexicon p(target|source): discount=" << model.phrase.lexicon.parameter_source_target.discount
		<< " strength=" << model.phrase.lexicon.parameter_source_target.strength << std::endl
		<< "lexicon p(source|target): discount=" << model.phrase.lexicon.parameter_target_source.discount
		<< " strength=" << model.phrase.lexicon.parameter_target_source.strength << std::endl
		<< "length source: lambda=" << model.phrase.length.length_source.lambda << std::endl
		<< "length target: lambda=" << model.phrase.length.length_target.lambda << std::endl;
    
    Task::queue_type queue_mapper;
    Task::queue_type queue_reducer;
    
    std::vector<Task, std::allocator<Task> > tasks(threads, Task(queue_mapper,
								 queue_reducer,
								 sources,
								 targets,
								 derivations,
								 derivations_prev,
								 model,
								 sampler,
								 beam, 
								 max_phrase_length));
    
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

      // assign temperature and model...
      for (size_type i = 0; i != tasks.size(); ++ i) {
	tasks[i].model = model;
	tasks[i].temperature = temperature;
      }
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(positions.begin(), positions.end(), gen);
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_size(sources, targets));
      
      position_set_type::const_iterator piter_end = positions.end();
      position_set_type::const_iterator piter = positions.begin();
      
      size_type invalid = 0;
      size_type reduced = 0;
      while (piter != piter_end || reduced != positions.size()) {
	
	for (int i = 0; i != threads && piter != piter_end && queue_mapper.push(*piter, true); ++ i)
	  ++ piter;
	
	size_type pos = 0;
	if (reduced != positions.size() && queue_reducer.pop(pos, piter != piter_end)) {
	  ++ reduced;
	  
	  if (! derivations_prev[pos].empty()) {
	    derivation_type::const_iterator diter_end = derivations_prev[pos].end();
	    for (derivation_type::const_iterator diter = derivations_prev[pos].begin(); diter != diter_end; ++ diter)
	      model.decrement(sources[pos], targets[pos], *diter, sampler);
	  }

	  if (debug >= 3) {
	    std::cerr << "training=" << pos << std::endl
		      << "source=" << sources[pos] << std::endl
		      << "target=" << targets[pos] << std::endl;
	    
	    
	    derivation_type::const_iterator diter_end = derivations[pos].end();
	    for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter) {
	      
	      std::cerr << "derivation: ";
	      switch (diter->itg) {
	      case PYP::TERMINAL:   std::cerr << "ter"; break;
	      case PYP::STRAIGHT:   std::cerr << "str"; break;
	      case PYP::INVERTED:   std::cerr << "inv"; break;
	      case PYP::GENERATIVE: std::cerr << "gen"; break;
	      case PYP::BASE:       std::cerr << "bas"; break;
	      default: std::cerr << "UNK";
	      }
	      
	      std::cerr << " source: " << diter->span.source.first << "..." << diter->span.source.last
			<< " target: " << diter->span.target.first << "..." << diter->span.target.last;
	      
	      if (diter->itg == PYP::GENERATIVE || diter->itg == PYP::TERMINAL || diter->itg == PYP::BASE)
		std::cerr << " pair: "
			  << PYP::phrase_type(sources[pos].begin() + diter->span.source.first, sources[pos].begin() + diter->span.source.last)
			  << " ||| "
			  << PYP::phrase_type(targets[pos].begin() + diter->span.target.first, targets[pos].begin() + diter->span.target.last);
	      
	      std::cerr << std::endl;
	    }
	  }
	  
	  invalid += derivations[pos].empty();

	  // increment model...
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    model.increment(sources[pos], targets[pos], *diter, sampler, temperature);
	  
	  if (debug) {
	    if ((reduced + 1) % 10000 == 0)
	      std::cerr << '.';
	    if ((reduced + 1) % 1000000 == 0)
	      std::cerr << '\n';
	  }
	}
      }
      
      if (debug && (reduced + 1) >= 10000 && (reduced + 1) % 1000000 != 0)
	std::cerr << std::endl;

      if (debug)
	std::cerr << "training: " << positions.size() << " empty derivations: " << invalid << std::endl;
      
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  model.slice_sample_parameters(sampler, resample_iterations);
	else
	  model.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2)
	  std::cerr << "rule: discount=" << model.rule.table.discount()
		    << " strength=" << model.rule.table.strength() << std::endl
		    << "terminal=" << model.rule.prob_terminal() << " straight=" << model.rule.prob_straight() << " inverted=" << model.rule.prob_inverted() << std::endl
		    << "phrase: discount=" << model.phrase.table.discount()
		    << " strength=" << model.phrase.table.strength() << std::endl
		    << "lexicon p(target|source): discount=" << model.phrase.lexicon.parameter_source_target.discount
		    << " strength=" << model.phrase.lexicon.parameter_source_target.strength << std::endl
		    << "lexicon p(source|target): discount=" << model.phrase.lexicon.parameter_target_source.discount
		    << " strength=" << model.phrase.lexicon.parameter_target_source.strength << std::endl
		    << "length source: lambda=" << model.phrase.length.length_source.lambda << std::endl
		    << "length target: lambda=" << model.phrase.length.length_target.lambda << std::endl;
      }
      
      if (debug)
	std::cerr << "log-likelihood: " << model.log_likelihood() << std::endl;
      
      if (sampling && ! output_sample_file.empty()) {
	typedef std::vector<std::string, std::allocator<std::string> > stack_type;
	
	// dump derivations..!
	
	const path_type path = add_suffix(output_sample_file, "." + utils::lexical_cast<std::string>(sample_iter + 1));
	
	utils::compress_ostream os(path, 1024 * 1024);

	if (sample_hypergraph) {
	  typedef hypergraph_type::rule_type     rule_type;
	  typedef hypergraph_type::rule_ptr_type rule_ptr_type;
	  
	  hypergraph_type graph_source;
	  hypergraph_type graph_target;
	  
	  for (size_type pos = 0; pos != derivations.size(); ++ pos) {
	    graph_source.clear();
	    graph_target.clear();
	    
	    if (! derivations[pos].empty()) {
	      // construct pair of hypergrpahs
	      
	      
	      
	      
	      
	    }
	    
	    os << graph_source << " ||| " << graph_target << '\n';
	  }
	} else {
	  stack_type stack;
	  
	  for (size_type pos = 0; pos != derivations.size(); ++ pos) {
	    if (! derivations[pos].empty()) {
	      // we need to transform the stack-structure into tree-struct... HOW?
	    
	      stack.clear();
	      derivation_type::const_iterator diter_end = derivations[pos].end();
	      for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter) {
		if (diter->is_terminal()) {
		  os << "((( "
		     << PYP::phrase_type(sources[pos].begin() + diter->span.source.first, sources[pos].begin() + diter->span.source.last)
		     << " ||| "
		     << PYP::phrase_type(targets[pos].begin() + diter->span.target.first, targets[pos].begin() + diter->span.target.last)
		     << " )))";
		
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
				      std::allocator<std::pair<const phrase_pair_type, reordering_type> > >::type reordering_set_type;
	
	typedef utils::vector2<bool, std::allocator<bool> > matrix_type;
	
	matrix_type matrix;
	reordering_set_type reorderings;
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
	
	PYPPhrase::table_type::const_iterator titer_end = model.phrase.table.end();
	for (PYPPhrase::table_type::const_iterator titer = model.phrase.table.begin(); titer != titer_end; ++ titer) 
	  if (! titer->first.source.empty() && ! titer->first.target.empty()) {
	    const double prob = model.phrase.table.prob(titer->first, 0.0);
	    
	    phrases[titer->first] = prob;
	    phrases_source[titer->first.source] += prob;
	    phrases_target[titer->first.target] += prob;
	  }
	
	const path_type path = add_suffix(output_model_file, "." + utils::lexical_cast<std::string>(sample_iter + 1));
	
	utils::compress_ostream os(path, 1024 * 1024);

	const double penalty = std::exp(1);
	
	phrase_pair_set_type::const_iterator piter_end = phrases.end();
	for (phrase_pair_set_type::const_iterator piter = phrases.begin(); piter != piter_end; ++ piter) {
	  const reordering_type& reorder = reorderings[piter->first];
	  
	  os << piter->first.source << " ||| " << piter->first.target
	     << " |||"
	     << " " << (piter->second / phrases_source[piter->first.source])
	     << " " << (piter->second / phrases_target[piter->first.target])
	     << " " << piter->second
	     << " " << model.phrase.lexicon.prob_source_target(piter->first.source, piter->first.target)
	     << " " << model.phrase.lexicon.prob_target_source(piter->first.target, piter->first.source)
	     << " " << penalty
	     << " |||"
	     << " " << reorder[0]
	     << " " << reorder[1]
	     << " " << reorder[2]
	     << " " << reorder[3]
	     << " " << reorder[4]
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

    ("lexicon-source-target", po::value<path_type>(&lexicon_source_target_file), "lexicon file for p(target | source)")
    ("lexicon-target-source", po::value<path_type>(&lexicon_target_source_file), "lexicon file for p(source | target)")

    ("max-phrase-length",   po::value<int>(&max_phrase_length)->default_value(max_phrase_length),     "max phrase length")
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
    
    ("rule-discount",       po::value<double>(&rule_discount)->default_value(rule_discount),                         "discount ~ Beta(alpha,beta)")
    ("rule-discount-alpha", po::value<double>(&rule_discount_prior_alpha)->default_value(rule_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("rule-discount-beta",  po::value<double>(&rule_discount_prior_beta)->default_value(rule_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("rule-strength",       po::value<double>(&rule_strength)->default_value(rule_strength),                         "strength ~ Gamma(shape,rate)")
    ("rule-strength-shape", po::value<double>(&rule_strength_prior_shape)->default_value(rule_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("rule-strength-rate",  po::value<double>(&rule_strength_prior_rate)->default_value(rule_strength_prior_rate),   "strength ~ Gamma(shape,rate)")

    ("phrase-discount",       po::value<double>(&phrase_discount)->default_value(phrase_discount),                         "discount ~ Beta(alpha,beta)")
    ("phrase-discount-alpha", po::value<double>(&phrase_discount_prior_alpha)->default_value(phrase_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("phrase-discount-beta",  po::value<double>(&phrase_discount_prior_beta)->default_value(phrase_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("phrase-strength",       po::value<double>(&phrase_strength)->default_value(phrase_strength),                         "strength ~ Gamma(shape,rate)")
    ("phrase-strength-shape", po::value<double>(&phrase_strength_prior_shape)->default_value(phrase_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("phrase-strength-rate",  po::value<double>(&phrase_strength_prior_rate)->default_value(phrase_strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
    ("lexicon-discount",       po::value<double>(&lexicon_discount)->default_value(lexicon_discount),                         "discount ~ Beta(alpha,beta)")
    ("lexicon-discount-alpha", po::value<double>(&lexicon_discount_prior_alpha)->default_value(lexicon_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("lexicon-discount-beta",  po::value<double>(&lexicon_discount_prior_beta)->default_value(lexicon_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("lexicon-strength",       po::value<double>(&lexicon_strength)->default_value(lexicon_strength),                         "strength ~ Gamma(shape,rate)")
    ("lexicon-strength-shape", po::value<double>(&lexicon_strength_prior_shape)->default_value(lexicon_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("lexicon-strength-rate",  po::value<double>(&lexicon_strength_prior_rate)->default_value(lexicon_strength_prior_rate),   "strength ~ Gamma(shape,rate)")

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


