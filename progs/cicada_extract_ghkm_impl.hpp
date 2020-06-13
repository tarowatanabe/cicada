//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXTRACT_GHKM_IMPL__HPP__
#define __CICADA__EXTRACT_GHKM_IMPL__HPP__ 1

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>


#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>
#include <boost/array.hpp>

#include <string>
#include <vector>
#include <set>
#include <deque>

#include "cicada/sentence.hpp"
#include "cicada/alignment.hpp"
#include "cicada/vocab.hpp"
#include "cicada/tree_rule.hpp"
#include "cicada/hypergraph.hpp"
#include "cicada/operation/functional.hpp"
#include "cicada/semiring.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/span_edge.hpp"

#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/chart.hpp"
#include "utils/simple_vector.hpp"
#include "utils/small_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/indexed_set.hpp"
#include "utils/map_file_allocator.hpp"

#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/vector_set.hpp>
#include <utils/compact_map.hpp>
#include <utils/compact_set.hpp>
#include <utils/hashmurmur3.hpp>
#include <utils/getline.hpp>

struct Bitext
{
  typedef cicada::HyperGraph hypergraph_type;
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Alignment  alignment_type;
  
  hypergraph_type source;
  sentence_type   target;
  alignment_type  alignment;
  
  Bitext() : source(), target(), alignment() {}
  Bitext(const hypergraph_type& __source,
	 const sentence_type&   __target,
	 const alignment_type&  __alignment)
    : source(__source), target(__target), alignment(__alignment) {}
  
  void swap(Bitext& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    alignment.swap(x.alignment);
  }

  void clear()
  {
    source.clear();
    target.clear();
    alignment.clear();
  }

  void assign(const std::string& line)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end  = line.end();
    
    if ((! source.assign(iter, end))
	|| (! qi::phrase_parse(iter, end, qi::lit("|||"), standard::space))
	|| (! target.assign(iter, end))
	|| (! qi::phrase_parse(iter, end, qi::lit("|||"), standard::space))
	|| (! alignment.assign(iter, end))
	|| iter != end)
      clear();
  }
  
  friend
  std::istream& operator>>(std::istream& is, Bitext& bitext)
  {
    bitext.clear();
    
    std::string line;
    if (utils::getline(is, line))
      bitext.assign(line);
    
    return is;
  }
  
  friend
  std::ostream& operator<<(std::ostream& os, const Bitext& bitext)
  {
    os << bitext.source
       << " ||| " << bitext.target
       << " ||| " << bitext.alignment;
    return os;
  }
  
};

namespace std
{
  inline
  void swap(Bitext& x, Bitext& y)
  {
    x.swap(y);
  }
};


struct RulePair
{
  typedef std::string phrase_type;
  typedef cicada::Alignment alignment_type;
  typedef double count_type;
  typedef int64_t frequency_type;
  typedef boost::array<frequency_type, 3> freqs_type;
  
  const phrase_type*    source;
  const phrase_type*    target;
  const alignment_type* alignment;
  count_type     count;
  freqs_type     freqs;

  RulePair() : source(0), target(0), alignment(0), count(0), freqs() {}
  RulePair(const phrase_type& __source,
	   const phrase_type& __target,
	   const alignment_type& __alignment,
	   const count_type& __count)
    : source(&__source), target(&__target), alignment(&__alignment), count(__count), freqs() {}

  friend
  size_t hash_value(RulePair const& x)
  {
    typedef utils::hashmurmur3<size_t> hasher_type;
    
    return hasher_type()(x.source, hasher_type()(x.target, (uintptr_t) x.alignment));
  }

  friend
  bool operator==(const RulePair& x, const RulePair& y) 
  {
    return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
  }
  
  friend
  bool operator!=(const RulePair& x, const RulePair& y) 
  {
    return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
  }
  
  friend
  bool operator<(const RulePair& x, const RulePair& y)
  {
    return (RulePair::compare(x.source, y.source)
	    || (!RulePair::compare(y.source, x.source)
		&& (RulePair::compare(x.target, y.target)
		    || (!RulePair::compare(y.target, x.target)
			&& RulePair::compare(x.alignment, y.alignment)))));
  }
  
  friend
  bool operator>(const RulePair& x, const RulePair& y)
  {
    return y < x;
  }
  
  template <typename Tp>
  static inline
  bool compare(const Tp* x, const Tp* y)
  {
    return (x && y && *x < *y) || (!x && y);
  }
};

BOOST_FUSION_ADAPT_STRUCT(RulePair::alignment_type::point_type,
			  (RulePair::alignment_type::index_type, source)
			  (RulePair::alignment_type::index_type, target)
			  )

struct RulePairGenerator
{
  typedef RulePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::frequency_type frequency_type;
  typedef phrase_pair_type::count_type     count_type;
  
  RulePairGenerator() {}
  
  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair) const
  {
    typedef std::ostream_iterator<char> iterator_type;
  
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    iterator_type iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string
			  << " ||| " << -((karma::int_ << '-' << karma::int_) % ' ')
			  << " ||| " << double20 << ' ' << (frequency % ' '),
			  *phrase_pair.source,
			  *phrase_pair.target,
			  *phrase_pair.alignment,
			  phrase_pair.count,
			  phrase_pair.freqs))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };
  
  boost::spirit::karma::real_generator<double, real_precision> double20;
  boost::spirit::karma::int_generator<frequency_type> frequency;
};

//
// forest-based translation rule extraction based on
//
// @InProceedings{mi-huang:2008:EMNLP,
//   author    = {Mi, Haitao  and  Huang, Liang},
//   title     = {Forest-based Translation Rule Extraction},
//   booktitle = {Proceedings of the 2008 Conference on Empirical Methods in Natural Language Processing},
//   month     = {October},
//   year      = {2008},
//   address   = {Honolulu, Hawaii},
//   publisher = {Association for Computational Linguistics},
//   pages     = {206--214},
//   url       = {http://www.aclweb.org/anthology/D08-1022}
// }
//
// and
//
// @InProceedings{galley-EtAl:2006:COLACL,
//    author    = {Galley, Michel  and  Graehl, Jonathan  and  Knight, Kevin  and  Marcu, Daniel  and  DeNeefe, Steve  and  Wang, Wei  and  Thayer, Ignacio},
//    title     = {Scalable Inference and Training of Context-Rich Syntactic Translation Models},
//    booktitle = {Proceedings of the 21st International Conference on Computational Linguistics and 44th Annual Meeting of the Association for Computational Linguistics},
//    month     = {July},
//    year      = {2006},
//    address   = {Sydney, Australia},
//    publisher = {Association for Computational Linguistics},
//    pages     = {961--968},
//    url       = {http://www.aclweb.org/anthology/P06-1121},
//    doi       = {10.3115/1220175.1220296}
// }
//


struct ExtractGHKM
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef RulePair rule_pair_type;
  typedef rule_pair_type::phrase_type phrase_type;
  typedef rule_pair_type::count_type  count_type;


  struct rule_compact_set_type : public utils::hashmurmur3<size_t>
  {
    struct string_hash : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      size_t operator()(const std::string& x) const
      {
	return hasher_type::operator()(x.begin(), x.end(), 0);
      }
    };
    
    typedef utils::indexed_set<phrase_type, string_hash, std::equal_to<phrase_type>, 
			       std::allocator<phrase_type> > phrase_set_type;
    
    typedef phrase_set_type::index_type      index_type;
    typedef phrase_set_type::size_type       size_type;
    typedef phrase_set_type::difference_type difference_type;

    void clear()
    {
      phrases.clear();
    }

    void swap(rule_compact_set_type& x)
    {
      phrases.swap(x.phrases);
    }
    
    index_type insert(const phrase_type& x)
    {
      phrase_set_type::iterator iter = phrases.insert(x).first;
      return iter - phrases.begin();
    }
    
    const phrase_type& operator[](index_type x) const { return phrases[x]; }
    
    size_type size() const { return phrases.size(); }
    bool empty() const { return phrases.empty(); }
    
    phrase_set_type phrases;
  };

  struct alignment_set_type
  {
    typedef rule_pair_type::alignment_type alignment_type;

    typedef utils::indexed_set<alignment_type,
			       boost::hash<alignment_type>,
			       std::equal_to<alignment_type>,
			       std::allocator<alignment_type> > align_set_type;
    
    typedef align_set_type::index_type      index_type;
    typedef align_set_type::size_type       size_type;
    typedef align_set_type::difference_type difference_type;

    void clear()
    {
      aligns.clear();
    }
    
    void swap(alignment_set_type& x)
    {
      aligns.swap(x.aligns);
    }
    
    index_type insert(const alignment_type& x)
    {
      align_set_type::iterator iter = aligns.insert(x).first;
      return iter - aligns.begin();
    }
    
    const alignment_type& operator[](index_type x) const { return aligns[x]; }
    
    size_type size() const { return aligns.size(); }
    bool empty() const { return aligns.empty(); }
    
    align_set_type aligns;
  };
  
  struct RulePairCompact
  {
    rule_compact_set_type::index_type source;
    rule_compact_set_type::index_type target;
    alignment_set_type::index_type    alignment;
    count_type                        count;
    
    RulePairCompact()
      : source(rule_compact_set_type::index_type(-1)),
	target(rule_compact_set_type::index_type(-1)),
	alignment(alignment_set_type::index_type(-1)),
	count(0) {}

    friend
    size_t hash_value(RulePairCompact const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.target, hasher_type()(x.alignment, x.source));
    }
    
    friend
    bool operator==(const RulePairCompact& x, const RulePairCompact& y) 
    {
      return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
    }
    
    friend
    bool operator!=(const RulePairCompact& x, const RulePairCompact& y)
    {
      return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
    }
  };
  
  typedef RulePairCompact rule_pair_compact_type;

  struct rule_pair_compact_unassigned
  {
    rule_pair_compact_type operator()() const
    {
      return rule_pair_compact_type();
    }
  };
  
  typedef utils::compact_set<rule_pair_compact_type,
			     rule_pair_compact_unassigned, rule_pair_compact_unassigned,
			     boost::hash<rule_pair_compact_type>, std::equal_to<rule_pair_compact_type>,
			     std::allocator<rule_pair_compact_type> > rule_pair_compact_set_type;
  
  
  struct rule_pair_set_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef rule_pair_type::alignment_type alignment_type;
    
    struct string_hash : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      size_t operator()(const std::string& x) const
      {
	return hasher_type::operator()(x.begin(), x.end(), 0);
      }
    };
    
    typedef utils::indexed_set<phrase_type, string_hash, std::equal_to<phrase_type>, 
			       std::allocator<phrase_type> > phrase_set_type;

    typedef utils::indexed_set<alignment_type,
			       boost::hash<alignment_type>,
			       std::equal_to<alignment_type>,
			       std::allocator<alignment_type> > alignment_set_type;

    struct rule_pair_unassigned
    {
      rule_pair_type operator()() const
      {
	return rule_pair_type();
      }
    };
    
    typedef utils::compact_set<rule_pair_type,
			       rule_pair_unassigned, rule_pair_unassigned,
			       boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
			       std::allocator<rule_pair_type> > pair_set_type;

    typedef pair_set_type::const_iterator const_iterator;
    typedef pair_set_type::iterator       iterator;
    
    bool empty() const { return rules.empty(); }
    size_type size() const { return rules.size(); }

    const_iterator begin() const { return rules.begin(); }
    const_iterator end() const { return rules.end(); }

    std::pair<iterator, bool> insert(const rule_pair_type& x)
    {
      return rules.insert(rule_pair_type(*sources.insert(*x.source).first,
					 *targets.insert(*x.target).first,
					 *alignments.insert(*x.alignment).first,
					 x.count));
    }
    
    void clear()
    {
      rules.clear();
      sources.clear();
      targets.clear();
      alignments.clear();
    }

    void swap(rule_pair_set_type& x)
    {
      rules.swap(x.rules);
      sources.swap(x.sources);
      targets.swap(x.targets);
      alignments.swap(x.alignments);
    }
    
    pair_set_type      rules;
    phrase_set_type    sources;
    phrase_set_type    targets;
    alignment_set_type alignments;
  };
  
  typedef std::vector<bool, std::allocator<bool> > unique_set_type;
  
  typedef std::pair<rule_compact_set_type::index_type, rule_compact_set_type::index_type> unique_pair_type;
  
  struct unique_pair_unassigned
  {
    unique_pair_type operator()() const
    {
      return unique_pair_type(rule_compact_set_type::index_type(-1),
			      rule_compact_set_type::index_type(-1));
    }
  };
  
  typedef utils::compact_set<unique_pair_type,
			     unique_pair_unassigned, unique_pair_unassigned,
			     utils::hashmurmur3<size_t>, std::equal_to<unique_pair_type>,
			     std::allocator<unique_pair_type> > unique_pair_set_type;
  
  typedef cicada::HyperGraph hypergraph_type;
  typedef cicada::Symbol     word_type;
  typedef cicada::Symbol     symbol_type;
  typedef cicada::Vocab      vocab_type;
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Alignment  alignment_type;
  typedef cicada::TreeRule   tree_rule_type;
  typedef cicada::Rule       rule_type;
  
  typedef alignment_type::point_type point_type;
  
  typedef std::pair<int, int> range_type;
  typedef std::vector<range_type, std::allocator<range_type> > range_set_type;

  typedef hypergraph_type::attribute_set_type attribute_set_type;
  
  typedef attribute_set_type::attribute_type attribute_type;

  struct Span
  {
    typedef std::set<int, std::less<int>, std::allocator<int> > span_type;
    //typedef utils::vector_set<int, std::less<int>, std::allocator<int> > span_type;
    
    typedef span_type::const_iterator const_iterator;
    typedef span_type::iterator       iterator;
    
    Span() : span() {}
    Span(const int first, const int last) : span()
    {
      for (int i = first; i < last; ++ i)
	span.insert(i);
    }
    
    Span& operator|=(const Span& x)
    {
      span.insert(x.span.begin(), x.span.end());
      return *this;
    }

    const_iterator begin() const { return span.begin(); }
    const_iterator end() const { return span.end(); }

    const_iterator lower_bound(int pos) const { return span.lower_bound(pos); }
    
    range_type range() const
    {
      if (span.empty())
	return std::make_pair(0, 0);
      else
	return std::make_pair(*(span.begin()), *(-- span.end()) + 1);
      //return std::make_pair(span.front(), span.back() + 1);
    }
  
    void set(int i)
    {
      span.insert(i);
    }
    
    void clear() { span.clear(); }

    bool intersect(const range_type& range) const
    {
      if (range.first == range.second)
	return false;
      
      if (span.empty())
	return false;
      else
	return span.lower_bound(range.first) != span.lower_bound(range.second);
    }

    bool empty() const { return span.empty(); }
  
    span_type span;
  };

  typedef Span span_type;
  typedef std::vector<span_type, std::allocator<span_type> > span_set_type;

  typedef std::vector<bool, std::allocator<bool> > admissible_set_type;
  
  typedef hypergraph_type::id_type id_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > edge_set_type;
  typedef std::vector<id_type, std::allocator<id_type> > node_set_type;

  typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;

  typedef std::pair<edge_set_type, node_set_type> frontier_type;

  typedef cicada::semiring::Logprob<double> weight_type;
  typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
  
  struct DerivationEdge
  {
    edge_set_type edges;
    node_set_type tails;

    int height;   // height of this rule
    int internal; // # of nodes in the rule including frontiers, but excluding root.
    int compose;  // # of composed rules

    DerivationEdge()
      : edges(), tails(), height(0), internal(0), compose(1) {}
    DerivationEdge(const edge_set_type& __edges, const node_set_type& __tails)
      : edges(__edges), tails(__tails), height(0), internal(0), compose(1) {}
    DerivationEdge(const edge_set_type& __edges, const node_set_type& __tails,
		   const int __height,
		   const int __internal)
      : edges(__edges), tails(__tails), height(__height), internal(__internal), compose(1) {}
    DerivationEdge(const edge_set_type& __edges, const node_set_type& __tails,
		   const int __height,
		   const int __internal,
		   const int __compose)
      : edges(__edges), tails(__tails), height(__height), internal(__internal), compose(__compose) {}
    
    void swap(DerivationEdge& x)
    {
      edges.swap(x.edges);
      tails.swap(x.tails);
      std::swap(height, x.height);
      std::swap(internal, x.internal);
      std::swap(compose, x.compose);
    }
  };
  
  struct DerivationNode
  {
    typedef DerivationEdge derivation_edge_type;
    typedef std::deque<derivation_edge_type, std::allocator<derivation_edge_type> > edge_set_type;

    id_type       node;
    range_type    range;
    edge_set_type edges;
  };
  
  typedef DerivationEdge derivation_edge_type;
  typedef DerivationNode derivation_node_type;

  struct less_derivation_edge_type
  {
    bool operator()(const derivation_edge_type& x, const derivation_edge_type& y) const
    {
      return (x.height < y.height
	      || (!(y.height < x.height)
		  && (x.internal < y.internal
		      || (!(y.internal < x.internal)
			  && x.compose < y.compose))));
    }
  };
  
  typedef utils::chunk_vector<derivation_node_type, 4096 / sizeof(derivation_node_type), std::allocator<derivation_node_type> > derivation_set_type;

  typedef utils::simple_vector<int, std::allocator<int> > index_set_type;
  typedef std::vector<bool, std::allocator<bool> > covered_type;
  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > alignment_map_type;
  typedef std::vector<tree_rule_type, std::allocator<tree_rule_type> > tree_rule_set_type;

  typedef std::vector<symbol_type, std::allocator<symbol_type> > frontier_symbol_set_type;
  
  typedef std::pair<range_type, int> range_pos_type;
  typedef std::vector<range_pos_type, std::allocator<range_pos_type> > range_pos_set_type;

  typedef utils::simple_vector<id_type, std::allocator<id_type> > edge_set_local_type;

  typedef std::vector<rule_pair_compact_type, std::allocator<rule_pair_compact_type> > rule_pair_compact_list_type;

  struct node_set_hash : public utils::hashmurmur3<size_t>
  {
    size_t operator()(const edge_set_local_type& edges) const
    {
      return utils::hashmurmur3<size_t>::operator()(edges.begin(), edges.end(), 0);
    }
  };

  typedef utils::unordered_map<edge_set_local_type, rule_pair_compact_list_type, node_set_hash, std::equal_to<edge_set_local_type>,
			       std::allocator<std::pair<const edge_set_local_type, rule_pair_compact_list_type> > >::type rule_pair_set_span_type;
  
  ExtractGHKM(const symbol_type& __non_terminal,
	      const int __max_sentence_length,
	      const int __max_nodes,
	      const int __max_height,
	      const int __max_compose,
	      const int __max_scope,
	      const double __cutoff,
	      const bool __exhaustive,
	      const bool __constrained,
	      const bool __inverse,
	      const bool __swap_source_target,
	      const bool __project_non_terminal,
	      const bool __collapse_source,
	      const bool __collapse_target)
    : non_terminal(__non_terminal),
      max_sentence_length(__max_sentence_length),
      max_nodes(__max_nodes),
      max_height(__max_height),
      max_compose(__max_compose),
      max_scope(__max_scope),
      cutoff(__cutoff),
      exhaustive(__exhaustive),
      constrained(__constrained),
      inverse(__inverse),
      swap_source_target(__swap_source_target),
      project_non_terminal(__project_non_terminal),
      collapse_source(__collapse_source),
      collapse_target(__collapse_target),
      attr_span_first("span-first"),
      attr_span_last("span-last") {}

  symbol_type non_terminal;
  int max_sentence_length;
  int max_nodes;
  int max_height;
  int max_compose;
  int max_scope;

  double cutoff;
  
  bool exhaustive;
  bool constrained;
  bool inverse;
  bool swap_source_target;
  bool project_non_terminal;
  bool collapse_source;
  bool collapse_target;
  
  attribute_type attr_span_first;
  attribute_type attr_span_last;

  range_set_type span_edges;
  
  range_set_type ranges;
  span_set_type spans;
  span_set_type complements;

  admissible_set_type admissibles;
  
  node_map_type node_map;
  derivation_set_type derivations;

  alignment_map_type alignment_source_target;
  alignment_map_type alignment_target_source;

  weight_set_type weights_inside;
  weight_set_type weights_outside;

  void clear()
  {
    span_edges.clear();
    
    ranges.clear();
    spans.clear();
    complements.clear();
    
    node_map.clear();

    weights_inside.clear();
    weights_outside.clear();
    
    range_set_type(span_edges).swap(span_edges);
    
    range_set_type(ranges).swap(ranges);
    span_set_type(spans).swap(spans);
    span_set_type(complements).swap(complements);
    
    node_map_type(node_map).swap(node_map);
    
    weight_set_type(weights_inside).swap(weights_inside);
    weight_set_type(weights_outside).swap(weights_outside);
  }
  
  template <typename Dumper>
  void operator()(const hypergraph_type& graph,
		  const sentence_type& sentence,
		  const alignment_type& alignment,
		  rule_pair_set_type& rules,
		  const Dumper& dumper)
  {
    if (max_sentence_length > 0 && static_cast<int>(sentence.size()) > max_sentence_length) return;
    
#if 0
    std::cerr << "hypergraph: " << graph << std::endl
	      << "sentence: " << sentence << std::endl
	      << "alignment: " << alignment << std::endl;
#endif
    derivations.clear();
    
    ranges.clear();
    spans.clear();
    complements.clear();
    admissibles.clear();
    node_map.clear();
    
    ranges.reserve(graph.nodes.size());
    spans.reserve(graph.nodes.size());
    complements.reserve(graph.nodes.size());
    admissibles.reserve(graph.nodes.size());
    node_map.reserve(graph.nodes.size());
    
    ranges.resize(graph.nodes.size());
    spans.resize(graph.nodes.size());
    complements.resize(graph.nodes.size());
    admissibles.resize(graph.nodes.size());
    node_map.resize(graph.nodes.size());
    
    //std::cerr << "admissible nodes" << std::endl;
    admissible_nodes(graph, sentence, alignment);
    
    if (max_sentence_length > 0
	&& (static_cast<int>(alignment_source_target.size()) > max_sentence_length
	    || static_cast<int>(alignment_target_source.size()) > max_sentence_length)) return;
    
    //std::cerr << "construct derivations" << std::endl;
    construct_derivations(graph, sentence);
    
    // compute reachable nodes and edges from root...

    //std::cerr << "prune derivations" << std::endl;
    //prune_derivations();
        
    // perform compounds extraction
    //std::cerr << "extraction" << std::endl;
    extract_composed(graph, sentence, rules, dumper);
  }
  
  

  struct Candidate
  {
    const derivation_edge_type* edge;
    derivation_edge_type edge_composed;
    
    index_set_type j;
    
    Candidate(const index_set_type& __j)
      : edge(0), edge_composed(), j(__j) {}
    
    Candidate(const derivation_edge_type& __edge, const index_set_type& __j)
      : edge(&__edge), edge_composed(__edge), j(__j) {}
  };
  typedef Candidate candidate_type;
  typedef utils::chunk_vector<candidate_type, 4096 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;

  struct compare_heap_type : public less_derivation_edge_type
  {
    // we use inverse less, so that when popped from heap, we will grab "less" in back...
    bool operator()(const candidate_type* x, const candidate_type* y) const
    {
      return less_derivation_edge_type::operator()(y->edge_composed, x->edge_composed);
    }
  };
  
  rule_pair_compact_set_type rule_pairs_local;

  template <typename Dumper>
  void extract_composed(const hypergraph_type& graph,
			const sentence_type& sentence,
			rule_pair_set_type& rule_pairs,
			const Dumper& dumper)
  {
    typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
    //typedef utils::b_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type, 512 / sizeof(const candidate_type*)> candidate_heap_type;
    typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;

    rule_pairs_local.clear();
    rules_source.clear();
    rules_target.clear();
    rules_alignment.clear();

    candidate_set_type    candidates;
    candidate_heap_type   cand;
    
    derivation_node_type::edge_set_type derivation_edges_new;
    
    // difficult...
    // we need to keep track of how many nodes / maximum height in the composed rules...
    
    //
    // bottom-up traversal to construct composed rules...
    // actually, we will construct new "composed edges" into new edge set, then, 
    // replace this new edge set with current one, so that we can easily enumerate
    // combination...
    //
    // first, we will implement naive enumeration...
    
    rule_pair_set_span_type rule_pairs_span;
    rule_pair_compact_type rule_pair;

    edge_set_type edges_new;
    node_set_type tails_new;

    //std::cerr << "derivations size: " << derivations.size() << std::endl;

    for (size_t id = 0; id != derivations.size(); ++ id) {
      derivation_node_type& node = derivations[id];
      
      //std::cerr << "node id: " << id << " range: [" << node.range.first << ", " << node.range.second << ")" << std::endl;
           
      derivation_edges_new.clear();
      candidates.clear();
      
      cand.clear();
      cand.reserve(node.edges.size() * 100);

      derivation_node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (derivation_node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const derivation_edge_type& edge = *eiter;
	const index_set_type j(edge.tails.size(), -1);
	
	candidates.push_back(candidate_type(edge, j));
	
	cand.push(&candidates.back());
      }
      
      while (! cand.empty()) {
	const candidate_type* item = cand.top();
	cand.pop();

	const derivation_edge_type& edge_composed = item->edge_composed;
	
	// construct rule pair and insert into rule pair list
	if (construct_rule_pair(graph, sentence, node, edge_composed.edges, edge_composed.tails, rule_pair))
	  rule_pairs_span[edge_set_local_type(edge_composed.edges.begin(), edge_composed.edges.end())].push_back(rule_pair);

	//std::cerr << rule_pair.source << " ||| " << rule_pair.target << std::endl;
	
	// no equality check since, derivation_edges_new are used for "composition"
	// thus, if internal == max_nodes, this cannnot be composed with others!
	// similarly, if compose == max_compose, this cannot be composed with others!
	if ((max_height <= 0 || edge_composed.height < max_height)
	    && (max_nodes <= 0 || edge_composed.internal < max_nodes)
	    && (max_compose <= 0 || edge_composed.compose < max_compose))
	  derivation_edges_new.push_back(edge_composed);
	
	// push-successor...
	
	const derivation_edge_type& edge = *(item->edge);
	
	//candidate_type query(item->j);
	//index_set_type& j = query.j;
	//query.edge = item->edge;
	index_set_type j(item->j);
	
	//std::cerr << "edge: " << edge.tails.size() << std::endl;
	
	// proceed to the next until we find better candidate...
	// @InProceedings{iglesias-EtAl:2009:EACL,
	//  author    = {Iglesias, Gonzalo  and  de Gispert, Adri\`{a}  and  Banga, Eduardo R.  and  Byrne, William},
	//  title     = {Rule Filtering by Pattern for Efficient Hierarchical Translation},
	//  booktitle = {Proceedings of the 12th Conference of the European Chapter of the ACL (EACL 2009)},
	//  month     = {March},
	//  year      = {2009},
	//  address   = {Athens, Greece},
	//  publisher = {Association for Computational Linguistics},
	//  pages     = {380--388},
	//  url       = {http://www.aclweb.org/anthology/E09-1044}
	// }

	for (size_t i = 0; i != j.size(); ++ i) 
	  if (! derivations[edge.tails[i]].edges.empty()) {
	    const int j_i_prev = j[i];
	    ++ j[i];
	    
	    //std::cerr << "i = " << i << " j[i] = " << j[i] << std::endl;

	    int composed_size_prev = edge_composed.compose;
	    if (j_i_prev >= 0)
	      composed_size_prev -= derivations[edge.tails[i]].edges[j_i_prev].compose;

	    int internal_size_prev = edge_composed.internal;
	    if (j_i_prev >= 0)
	      internal_size_prev -= derivations[edge.tails[i]].edges[j_i_prev].internal;
	    
	    for (/**/; j[i] < static_cast<int>(derivations[edge.tails[i]].edges.size()); ++ j[i]) {
	      const int composed_size = composed_size_prev + derivations[edge.tails[i]].edges[j[i]].compose;
	      const int internal_size = internal_size_prev + derivations[edge.tails[i]].edges[j[i]].internal;
	      
	      // early checking!
	      if (max_compose > 0 && composed_size > max_compose) continue;
	      if (max_nodes > 0 && internal_size > max_nodes) continue;
	      
	      index_set_type::const_iterator jiter_begin = j.begin();
	      index_set_type::const_iterator jiter_end   = j.end();
	      node_set_type::const_iterator  titer_begin = edge.tails.begin();
	      edge_set_type::const_iterator  eiter_begin = edge.edges.begin();
	      edge_set_type::const_iterator  eiter_end   = edge.edges.end();
	      
	      //std::cerr << "compose edges" << std::endl;

	      edges_new.clear();
	      tails_new.clear();
	      
	      const std::pair<int, int> rule_stat = compose_edges(graph, jiter_begin, jiter_end, titer_begin, eiter_begin, eiter_end, edges_new, tails_new);
	      
	      // early termination...
	      if (max_height > 0 && rule_stat.first > max_height) break;

#if 0
	      //checking!
	      if (rule_stat.second != internal_size)
		std::cerr << "differ: " << rule_stat.second << ' ' << internal_size << std::endl;

	      std::cerr << "height: " << rule_stat.first << " internal: " << rule_stat.second << std::endl;
	      for (size_t i = 0; i != edges_new.size(); ++ i)
		std::cerr << "rule: " << *(graph.edges[edges_new[i]].rule) << std::endl;
#endif
	      
	      candidates.push_back(candidate_type(edge, j));
	      
	      candidate_type& item_next = candidates.back();
	      
	      item_next.edge_composed.edges = edges_new;
	      item_next.edge_composed.tails = tails_new;
	      item_next.edge_composed.height = rule_stat.first;
	      item_next.edge_composed.internal = rule_stat.second;
	      item_next.edge_composed.compose = composed_size;
	      
	      cand.push(&item_next);
	      
	      break;
	    }
	    
	    if (item->j[i] != -1) break;
	    
	    j[i] = j_i_prev;
	  }
      }
            
      //std::cerr << "finished" << std::endl;
      
      // allocate new derivations!
      node.edges.swap(derivation_edges_new);
      
      // sort for cube-pruning!
      std::sort(node.edges.begin(), node.edges.end(), less_derivation_edge_type());

      // clear candidates, since all the cand is empty
      candidates.clear();
      
      // clear new derivation edges, since we alread swapped!
      derivation_edges_new.clear();
      
      // if we share the same edge set, then, we can easily conclude that it is
      // caused by unaligned word ambiguity...
      
      rule_pair_set_span_type::iterator riter_end = rule_pairs_span.end();
      for (rule_pair_set_span_type::iterator riter = rule_pairs_span.begin(); riter != riter_end; ++ riter) {
	rule_pair_compact_list_type& rule_list = riter->second;
	
	const double factor = 1.0 / rule_list.size();
	
	rule_pair_compact_list_type::iterator liter_end = rule_list.end();
	for (rule_pair_compact_list_type::iterator liter = rule_list.begin(); liter != liter_end; ++ liter) {
	  liter->count *= factor;
	  
	  if (liter->count >= cutoff) {
	    std::pair<rule_pair_compact_set_type::iterator, bool> result = rule_pairs_local.insert(*liter);
	    if (! result.second)
	      const_cast<rule_pair_compact_type&>(*(result.first)).count += liter->count;
	  }
	}
      }
      
      rule_pairs_span.clear();
    }
    
    // try dump when we have spurious allocation...
    dumper(rule_pairs);
    
    uniques_source.clear();
    uniques_target.clear();
    uniques_source.resize(rules_source.size(), true);
    uniques_target.resize(rules_target.size(), true);
    
    uniques_pair.clear();

    size_t iter = 0;
    size_t iter_mask = (1024 * 1024 * 4) - 1;
    
    rule_pair_compact_set_type::const_iterator riter_end = rule_pairs_local.end();
    for (rule_pair_compact_set_type::const_iterator riter = rule_pairs_local.begin(); riter != riter_end; ++ riter, ++ iter) {
      // uncover phrasal representation!
      const bool unique_source = uniques_source[riter->source];
      const bool unique_target = uniques_target[riter->target];

      uniques_source[riter->source] = false;
      uniques_target[riter->target] = false;
      
      std::pair<rule_pair_set_type::iterator, bool> result = rule_pairs.insert(rule_pair_type(rules_source[riter->source],
											      rules_target[riter->target],
											      rules_alignment[riter->alignment],
											      riter->count));
      
      rule_pair_type& rule_pair = const_cast<rule_pair_type&>(*result.first);
      
      if (! result.second)
	rule_pair.count += riter->count;
      
      rule_pair.freqs[0] += uniques_pair.insert(std::make_pair(riter->source, riter->target)).second;
      rule_pair.freqs[1] += unique_source;
      rule_pair.freqs[2] += unique_target;
      
      // dumping...
      if ((iter & iter_mask) == iter_mask)
	dumper(rule_pairs);
      
      //rule_pairs_local.erase(riter ++);
    }
    
    rule_pairs_local.clear();
    rules_source.clear();
    rules_target.clear();
    rules_alignment.clear();
    
    uniques_source.clear();
    uniques_target.clear();
    uniques_pair.clear();
    
    dumper(rule_pairs);
    
    if (rule_pairs.empty()) {
      rule_pair_compact_set_type(rule_pairs_local).swap(rule_pairs_local);
      rule_compact_set_type(rules_source).swap(rules_source);
      rule_compact_set_type(rules_target).swap(rules_target);
      alignment_set_type(rules_alignment).swap(rules_alignment);
      
      unique_set_type(uniques_source).swap(uniques_source);
      unique_set_type(uniques_target).swap(uniques_target);
      unique_pair_set_type(uniques_pair).swap(uniques_pair);
    }
  }

  template <typename IndexIterator, typename TailIterator, typename EdgeIterator, typename Edges, typename Tails>
  std::pair<int, int> compose_edges(const hypergraph_type& graph,
				    IndexIterator& iter, IndexIterator last,
				    TailIterator& tail_iter,
				    EdgeIterator& edge_iter, EdgeIterator edge_last,
				    Edges& edges_new,
				    Tails& tails_new)
  {
    // this should not happen...
    if (edge_iter == edge_last) return std::make_pair(0, 0);
    
    edges_new.push_back(*edge_iter);
    
    const hypergraph_type::edge_type& edge = graph.edges[*edge_iter];
    ++ edge_iter;

    int height = 1;
    int num_tails = edge.tails.size();
    //int num_tails = 0;
    
    hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
      if (edge_iter != edge_last && graph.edges[*edge_iter].head == *titer) {
	const std::pair<int, int> result = compose_edges(graph, iter, last, tail_iter, edge_iter, edge_last, edges_new, tails_new);
	
	height = utils::bithack::max(height, result.first + 1);
	num_tails += result.second;
	
	// early termination...
	if (max_height > 0 && height > max_height)
	  return std::make_pair(height, num_tails);

      } else if (iter != last) {
	if (*iter >= 0) {
	  const derivation_edge_type& edge = derivations[*tail_iter].edges[*iter];
	  
	  height = utils::bithack::max(height, edge.height + 1);
	  num_tails += edge.internal;
	  
	  // early termination...
	  if (max_height > 0 && height > max_height)
	    return std::make_pair(height, num_tails);
	  
	  edges_new.insert(edges_new.end(), edge.edges.begin(), edge.edges.end());
	  tails_new.insert(tails_new.end(), edge.tails.begin(), edge.tails.end());
	} else
	  tails_new.push_back(*tail_iter);
	
	++ iter;
	++ tail_iter;
      } else
	std::cerr << "WARNING: tails size and edge set frontier do not match" << std::endl;
    }
    
    return std::make_pair(height, num_tails);
  }  
  
  // construct rule-pair related data
  frontier_symbol_set_type frontier_symbols;
  point_set_type positions_source;
  point_set_type positions_target;
  tree_rule_set_type trees;
  covered_type covered;
  range_pos_set_type range_pos;
  
  template <typename Tp>
  struct less_first
  {
    bool operator()(const Tp& x, const Tp& y) const
    {
      return x.first < y.first;
    }
  };

  struct ScopeFrontierIterator
  {
    ScopeFrontierIterator(symbol_type& __prev, int& __scope) : prev(__prev), scope(__scope) {}
    
    ScopeFrontierIterator& operator=(const symbol_type& value)
    {
      scope += (prev == vocab_type::EMPTY || prev.is_non_terminal()) && value.is_non_terminal();
      prev = value;
      return *this;
    }
    
    ScopeFrontierIterator& operator*()  { return *this; }
    ScopeFrontierIterator& operator++() { return *this; }
    
    symbol_type& prev;
    int& scope;
  };

  struct CollapseFrontierIterator
  {
    CollapseFrontierIterator(tree_rule_set_type& __trees) : trees(__trees) {}
    
    CollapseFrontierIterator& operator=(const symbol_type& value)
    {
      trees.push_back(value);
      return *this;
    }
    
    CollapseFrontierIterator& operator*()  { return *this; }
    CollapseFrontierIterator& operator++() { return *this; }
    
    tree_rule_set_type& trees;
  };
  
  typedef std::vector<char, std::allocator<char> > buffer_type;
  
  buffer_type buffer_source;
  buffer_type buffer_target;
  
  rule_compact_set_type rules_source;
  rule_compact_set_type rules_target;
  alignment_set_type    rules_alignment;
  unique_set_type       uniques_source;
  unique_set_type       uniques_target;
  unique_pair_set_type  uniques_pair;
  
  bool construct_rule_pair(const hypergraph_type& graph,
			   const sentence_type& sentence,
			   const derivation_node_type& node,
			   const edge_set_type& edges,
			   const node_set_type& tails,
			   rule_pair_compact_type& rule_pair)
  {
    frontier_symbols.clear();
    positions_source.clear();
    covered.clear();
    covered.resize(alignment_source_target.size(), false);
    
    tree_rule_type rule_source;
    weight_type weight = weights_outside[node.node] / weights_inside.back();
    
    //std::cerr << "construct source" << std::endl;
    {
      int frontier_pos = 0;
      int index = 1;
      edge_set_type::const_iterator iter = edges.begin();
      edge_set_type::const_iterator iter_end = edges.end();
      construct_rule(graph, iter, iter_end, rule_source, index, frontier_pos, positions_source, frontier_symbols, covered, weight);
    }

    //std::cerr << "source: " << rule_source << std::endl;
    
    rule_pair.count = weight;
    
    //std::cerr << "construct target" << std::endl;
    
    range_pos.clear();
    trees.clear();
    positions_target.clear();
    positions_target.reserve(sentence.size());
    positions_target.resize(sentence.size());
    
    if (tails.empty()) {
      //std::cerr << "no tails" << std::endl;

      //std::cerr << "sentence: " << sentence << std::endl;
      //std::cerr << "range: " << node.range.first << ", " << node.range.second << std::endl;

      trees.insert(trees.end(), sentence.begin() + node.range.first, sentence.begin() + node.range.second);
      for (int i = node.range.first; i != node.range.second; ++ i)
	positions_target[i] = i - node.range.first;
    } else {
      //std::cerr << "tails: ";
      //std::copy(tails.begin(), tails.end(), std::ostream_iterator<int>(std::cerr, " "));
      //std::cerr << std::endl;
      
      int index = 1;
      range_pos.reserve(tails.size());
      node_set_type::const_iterator titer_end = tails.end();
      for (node_set_type::const_iterator titer = tails.begin(); titer != titer_end; ++ titer, ++ index)
	range_pos.push_back(std::make_pair(derivations[*titer].range, index));
      
      std::sort(range_pos.begin(), range_pos.end(), less_first<range_pos_type>());
      
      int pos_first = node.range.first;
      range_pos_set_type::const_iterator riter_end = range_pos.end();
      for (range_pos_set_type::const_iterator riter = range_pos.begin(); riter != riter_end; ++ riter) {
	const int pos_last = riter->first.first;
	
	int mapped_pos = trees.size();
	for (int i = pos_first; i != pos_last; ++ i, ++ mapped_pos)
	  positions_target[i] = mapped_pos;
	
	trees.insert(trees.end(), sentence.begin() + pos_first, sentence.begin() + pos_last);
	
	if (project_non_terminal)
	  trees.push_back(frontier_symbols[riter->second - 1].non_terminal(riter->second));
	else
	  trees.push_back(non_terminal.non_terminal(riter->second));
	
	pos_first = riter->first.second;
      }
      
      const int pos_last = node.range.second;
      int mapped_pos = trees.size();
      for (int i = pos_first; i != pos_last; ++ i, ++ mapped_pos)
	positions_target[i] = mapped_pos;
      
      trees.insert(trees.end(), sentence.begin() + pos_first, sentence.begin() + pos_last);
    }
    
    // this non-terminal can be derived from the root of rule_target
    tree_rule_type rule_target(project_non_terminal ? rule_source.label : non_terminal, trees.begin(), trees.end());

    //std::cerr << "target: " << rule_target << std::endl;

    //
    // construct word alignment...
    // we can easily reconsruct from positions_{source, target} + covered
    //

    //std::cerr << "construct alignment" << std::endl;
    
    //rule_pair.alignment.clear();
    alignment_type alignment;
    int pos_src = 0;
    
    if (! positions_source.empty()) {
      for (size_t src = 0; src != covered.size(); ++ src)
	if (covered[src]) {
	  point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	  point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	  
	  for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	    alignment.push_back(std::make_pair(positions_source[pos_src], positions_target[*aiter]));
	  
	  ++ pos_src;
	}
    }

    //std::cerr << "alignment: " << alignment << std::endl;

    
    if (swap_source_target) {
      //std::cerr << "swapping" << std::endl;
      
      // inverse alignment
      alignment.inverse();
      std::sort(alignment.begin(), alignment.end());

      rule_source.swap(rule_target);
      
      // re-assign index...
      cicada::sort(rule_source, rule_target);

#if 0
      std::cerr << "swapped source: " << rule_source << std::endl
		<< "swapped target: " << rule_target << std::endl
		<< "swapped alignment: " << alignment << std::endl;
#endif
    }

    //std::cerr << "dumping" << std::endl;
    
    if (max_scope > 0) {
      // compute the maximum scope for the source-side
      
      symbol_type prev(vocab_type::EMPTY);
      int scope = 0;
      
      rule_source.frontier(ScopeFrontierIterator(prev, scope));

      if ((scope + prev.is_non_terminal()) > max_scope)
	return false;
    }

    if (collapse_source) {
      // collapse source-side
      
      trees.clear();
      rule_source.frontier(CollapseFrontierIterator(trees));
      rule_source = tree_rule_type(rule_source.label, trees.begin(), trees.end());
    }
    
    if (collapse_target) {
      // collapse target-side
      
      trees.clear();
      rule_target.frontier(CollapseFrontierIterator(trees));
      rule_target = tree_rule_type(rule_target.label, trees.begin(), trees.end());
    }
    
    buffer_source.clear();
    buffer_target.clear();
    
    boost::iostreams::filtering_ostream os_source;
    boost::iostreams::filtering_ostream os_target;
    os_source.push(boost::iostreams::back_inserter(buffer_source));
    os_target.push(boost::iostreams::back_inserter(buffer_target));
    
    os_source << rule_source;
    os_target << rule_target;
    
    os_source.reset();
    os_target.reset();
    
    rule_pair.source = rules_source.insert(phrase_type(buffer_source.begin(), buffer_source.end()));
    rule_pair.target = rules_target.insert(phrase_type(buffer_target.begin(), buffer_target.end()));
    rule_pair.alignment = rules_alignment.insert(alignment);
    
    return true;
  }
  

  
  template <typename Iterator, typename PosMap, typename Frontiers, typename Covered>
  void construct_rule(const hypergraph_type& graph,
		      Iterator& iter,
		      Iterator last,
		      tree_rule_type& tree_rule,
		      int& index,
		      int& frontier_pos,
		      PosMap& pos_map,
		      Frontiers& frontiers,
		      Covered& covered,
		      weight_type& weight)
  {
    // pre-order traversal...
    
    if (iter == last) {
      tree_rule.label = tree_rule.label.non_terminal();
      return;
    }
    
    const hypergraph_type::edge_type& edge = graph.edges[*iter];

    const int edge_first = span_edges[*iter].first;
    const int edge_last  = span_edges[*iter].second;
    
    for (int pos = edge_first; pos != edge_last; ++ pos)
      covered[pos] = true;
    
    weight *= cicada::operation::weight_function_one<weight_type>()(edge);
    
    tree_rule = tree_rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end());
    
    ++ iter;
    
    
    size_t tail_pos = 0;
    tree_rule_type::iterator titer_end = tree_rule.end();
    for (tree_rule_type::iterator titer = tree_rule.begin(); titer != titer_end; ++ titer) {
      if (titer->label.is_non_terminal()) {
	
	const id_type node_id = edge.tails[tail_pos];
	
	if (iter != last && node_id == graph.edges[*iter].head)
	  construct_rule(graph, iter, last, *titer, index, frontier_pos, pos_map, frontiers, covered, weight);
	else {
	  for (int pos = ranges[node_id].first; pos != ranges[node_id].second; ++ pos)
	    covered[pos] = false;
	  
	  weight *= weights_inside[node_id];
	  
	  titer->label = titer->label.non_terminal(index ++);
	  frontiers.push_back(titer->label.non_terminal());
	  ++ frontier_pos;
	}
	
	++ tail_pos;
      } else {
	pos_map.push_back(frontier_pos);
	++ frontier_pos;
      }
    }
  }
  
  enum color_type {
    white,
    gray,
    black
  };
  
  struct dfs_type
  {
    int node;
    int edge;
    int tail;
    
    dfs_type(const int& _node, const int& _edge, const int& _tail) 
      : node(_node), edge(_edge), tail(_tail) {}
  };

  typedef std::vector<int, std::allocator<int> > reloc_set_type;
  typedef std::vector<color_type, std::allocator<color_type> > color_set_type;
  typedef std::vector<dfs_type, std::allocator<dfs_type> > stack_type;
  
  void prune_derivations()
  {
    // topologically sort...
    // I'm not sure whether we need this... but this is simply to make sure 
    // the extracted GHKM rules are reachable from root of the parse forest
    
    reloc_set_type reloc_node(derivations.size(), -1);
    color_set_type color(derivations.size(), white);
    stack_type stack;
    
    stack.reserve(derivations.size());
    stack.push_back(dfs_type(derivations.size() - 1, 0, 0));
    
    int node_count = 0;
    
    while (! stack.empty()) {
      const dfs_type& dfs = stack.back();
      
      int node_id     = dfs.node;
      size_t pos_edge = dfs.edge;
      size_t pos_tail = dfs.tail;
      
      stack.pop_back();
      
      const derivation_node_type* curr_node = &(derivations[node_id]);
      
      while (pos_edge != curr_node->edges.size()) {
	const derivation_edge_type& curr_edge = curr_node->edges[pos_edge];
	
	if (pos_tail == curr_edge.tails.size()) {
	  // reach end...
	  ++ pos_edge;
	  pos_tail = 0;
	  continue;
	}
	
	const int        tail_node  = curr_edge.tails[pos_tail];
	const color_type tail_color = color[tail_node];
	
	switch (tail_color) {
	case white:
	  ++ pos_tail;
	  stack.push_back(dfs_type(node_id, pos_edge, pos_tail));
	  
	  node_id = tail_node;
	  curr_node = &(derivations[node_id]);
	  
	  color[node_id] = gray;
	  pos_edge = 0;
	  pos_tail = 0;
	  
	  break;
	case black:
	  ++ pos_tail;
	  break;
	case gray:
	  ++ pos_tail;
	  // loop!!!!
	  break;
	}
      }
      
      color[node_id] = black;
      reloc_node[node_id] = node_count ++;
    }
    
    // construct new derivations...
    derivation_set_type derivations_new(node_count);
    for (size_t i = 0; i != derivations.size(); ++ i)
      if (reloc_node[i] >= 0) {
	const derivation_node_type& node_old = derivations[i];
	derivation_node_type& node_new = derivations_new[reloc_node[i]];
	
	node_new.node  = node_old.node;
	node_new.range = node_old.range;
	
	derivation_node_type::edge_set_type::const_iterator eiter_end = node_old.edges.end();
	for (derivation_node_type::edge_set_type::const_iterator eiter = node_old.edges.begin(); eiter != eiter_end; ++ eiter) {
	  
	  node_set_type tails;
	  bool is_valid = true;
	  
	  node_set_type::const_iterator titer_end = eiter->tails.end();
	  for (node_set_type::const_iterator titer = eiter->tails.begin(); is_valid && titer != titer_end; ++ titer) {
	    tails.push_back(reloc_node[*titer]);
	    is_valid = reloc_node[*titer] >= 0;
	  }
	  
	  if (is_valid)
	    node_new.edges.push_back(derivation_edge_type(eiter->edges, tails, eiter->height, eiter->internal));
	}
      }
    
    derivations.swap(derivations_new);
    derivations_new.clear();
  }
  
  struct unassigned_range
  {
    range_type operator()() const { return range_type(0, 0); }
  };

  struct deleted_range
  {
    range_type operator()() const { return range_type(-1, -1); }
  };

  void construct_derivations(const hypergraph_type& graph,
			     const sentence_type& sentence)
  {
    typedef std::deque<frontier_type, std::allocator<frontier_type> > queue_type;
    typedef utils::compact_map<range_type, id_type,
			       unassigned_range, deleted_range,
			       utils::hashmurmur3<size_t>, std::equal_to<range_type>,
			       std::allocator<std::pair<const range_type, id_type> > > range_node_map_type;
    typedef utils::compact_set<range_type,
			       unassigned_range, deleted_range,
			       utils::hashmurmur3<size_t>, std::equal_to<range_type>,
			       std::allocator<range_type> > range_set_type;

    // construc derivations wrt non-aligned words...

    size_t goal_node = size_t(-1);

    queue_type queue;
    range_node_map_type buf;
    
    for (size_t id = 0; id != graph.nodes.size(); ++ id) 
      if (admissibles[id]) {
	const hypergraph_type::node_type& node = graph.nodes[id];
	
	//std::cerr << "admissible node: " << id << std::endl;
	//
	// compute minimal range and maximum outer-range
	//
	
	const range_type range_min(spans[id].range());
	
	span_type::const_iterator riter_first = complements[id].lower_bound(range_min.first);
	span_type::const_iterator riter_last  = complements[id].lower_bound(range_min.second);
	
	const range_type range_max(riter_first == complements[id].begin() ? 0 : *(-- riter_first) + 1,
				   riter_last != complements[id].end() ? *riter_last : static_cast<int>(sentence.size()));
	
	const bool is_goal(id == graph.goal);
	
	queue.clear();
	buf.clear();
	
	if (is_goal) {
	  // insert minimum...
	  buf.insert(std::make_pair(range_type(0, sentence.size()), derivations.size()));
	  node_map[id].push_back(derivations.size());
	  
	  goal_node = derivations.size();
	  
	  derivations.resize(derivations.size() + 1);
	  
	  derivations.back().node = id;
	  derivations.back().range = range_type(0, sentence.size());
	} else {
	  // insert minimum...
	  buf.insert(std::make_pair(range_min, derivations.size()));
	  node_map[id].push_back(derivations.size());
	  
	  derivations.resize(derivations.size() + 1);
	  
	  derivations.back().node = id;
	  derivations.back().range = range_min;
	}

	//std::cerr << "new derivation: " << derivations.back().range.first << ", " << derivations.back().range.second << std::endl;
	
	// construct initial frontiers...
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  queue.resize(queue.size() + 1);
	  frontier_type& frontier = queue.back();
	  
	  frontier.first.push_back(*eiter);
	  
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	    if (! admissibles[*titer])
	      frontier.second.push_back(*titer);
	}

	range_set_type ranges;

	while (! queue.empty()) {
	  // stack operation...
	  const frontier_type frontier = queue.back();
	  queue.pop_back();
	  
	  if (frontier.second.empty()) {
	    // construct rule from "fragment", frontier.first
	    // by enumerating edge and its tails in depth-first manner...
	    // use recursive call for simplicity...
	    
	    // compute "tails" for the fragment
	    
	    node_set_type tails;
	    edge_set_type::const_iterator eiter     = frontier.first.begin();
	    edge_set_type::const_iterator eiter_end = frontier.first.end();
	    const std::pair<int, int> rule_stat     = construct_tails(graph, eiter, eiter_end, tails);
	    
	    // iterate over tails and compute span...
	    // code taken from apply_exact...
	    
	    index_set_type j_ends(tails.size(), 0);
	    index_set_type j(tails.size(), 0);
	    
	    for (size_t i = 0; i != tails.size(); ++ i)
	      j_ends[i] = node_map[tails[i]].size();
	    
	    node_set_type tails_next(tails.size());
	    
	    for (;;) {
	      bool is_valid = true;
	      
	      range_type range(range_min);
	      ranges.clear();
	      
	      for (size_t i = 0; is_valid && i != tails.size(); ++ i) {
		tails_next[i] = node_map[tails[i]][j[i]];
		
		const range_type& range_antecedent = derivations[tails_next[i]].range;
		
		if (! ranges.empty()) {
		  range_set_type::const_iterator riter_end = ranges.end();
		  for (range_set_type::const_iterator riter = ranges.begin(); is_valid && riter != riter_end; ++ riter)
		    is_valid = range_antecedent.second <= riter->first || riter->second <= range_antecedent.first;
		}
		
		if (is_valid) {
		  ranges.insert(range_antecedent);
		  
		  range.first  = utils::bithack::min(range.first,  range_antecedent.first);
		  range.second = utils::bithack::max(range.second, range_antecedent.second);
		}
	      }
	      
	      if (is_valid) {

		if (is_goal) {
		  if (goal_node == size_t(-1)) {
		    goal_node = derivations.size();
		    derivations.resize(goal_node + 1);
		    
		    // range is always [0, sentence.size())
		    derivations.back().node = id;
		    derivations.back().range = range_type(0, sentence.size());
		  }
		  
		  derivations[goal_node].edges.push_back(derivation_edge_type(frontier.first, tails_next, rule_stat.first, rule_stat.second));
		} else {
		  // we will compute all possible ranges...
		  
		  if (exhaustive) {
		    for (int first = range_max.first; first <= range.first; ++ first)
		      for (int last = range.second; last <= range_max.second; ++ last) {
			const range_type range_next(first, last);

			typedef std::pair<range_node_map_type::iterator, bool> result_type;
			
			result_type result = buf.insert(std::make_pair(range_next, derivations.size()));
			
			if (result.second) {
			  node_map[id].push_back(result.first->second);
			  
			  derivations.resize(derivations.size() + 1);
			  
			  derivations.back().node = id;
			  derivations.back().range = range_next;
			}
			
			derivations[result.first->second].edges.push_back(derivation_edge_type(frontier.first, tails_next, rule_stat.first, rule_stat.second));
		      }
		  } else {
		    const range_type& range_next = range;

		    typedef std::pair<range_node_map_type::iterator, bool> result_type;
		    
		    result_type result = buf.insert(std::make_pair(range_next, derivations.size()));
		    
		    if (result.second) {
		      node_map[id].push_back(result.first->second);
		      
		      derivations.resize(derivations.size() + 1);
		      
		      derivations.back().node = id;
		      derivations.back().range = range_next;
		    }
		    
		    derivations[result.first->second].edges.push_back(derivation_edge_type(frontier.first, tails_next, rule_stat.first, rule_stat.second));
		  }
		}
	      }
	      
	      size_t index = 0;
	      for (/**/; index != tails.size(); ++ index) {
		++ j[index];
		if (j[index] < j_ends[index]) break;
		j[index] = 0;
	      }
	      
	      // finished!
	      if (index == tails.size()) break;
	    }
	    
	  } else {
	    // incomplete... futher expand!

	    if (constrained) {
	      edge_set_type::const_iterator eiter     = frontier.first.begin();
	      edge_set_type::const_iterator eiter_end = frontier.first.end();
	      const std::pair<int, int> rule_stat     = rule_statistics(graph, eiter, eiter_end);

	      if ((max_height <= 0 || rule_stat.first <= max_height) && (max_nodes <= 0 || rule_stat.second < max_nodes)) {
		const hypergraph_type::node_type& node = graph.nodes[frontier.second.front()];
		
		hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
		for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
		  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
		
		  queue.resize(queue.size() + 1);
		  frontier_type& frontier_next = queue.back();
		  
		  frontier_next.first.reserve(frontier.first.size() + 1);
		  frontier_next.first.insert(frontier_next.first.end(), frontier.first.begin(), frontier.first.end());
		  frontier_next.first.push_back(*eiter);
		  
		  hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
		  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
		    if (! admissibles[*titer])
		      frontier_next.second.push_back(*titer);
		
		  frontier_next.second.insert(frontier_next.second.end(), frontier.second.begin() + 1, frontier.second.end());
		
		  node_set_type(frontier_next.second).swap(frontier_next.second);
		}
	      }
	    } else {
	      const hypergraph_type::node_type& node = graph.nodes[frontier.second.front()];
	      
	      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
		const hypergraph_type::edge_type& edge = graph.edges[*eiter];
		
		queue.resize(queue.size() + 1);
		frontier_type& frontier_next = queue.back();
		
		frontier_next.first.reserve(frontier.first.size() + 1);
		frontier_next.first.insert(frontier_next.first.end(), frontier.first.begin(), frontier.first.end());
		frontier_next.first.push_back(*eiter);
		
		hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
		for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
		  if (! admissibles[*titer])
		    frontier_next.second.push_back(*titer);
		
		frontier_next.second.insert(frontier_next.second.end(), frontier.second.begin() + 1, frontier.second.end());
		
		node_set_type(frontier_next.second).swap(frontier_next.second);
	      }
	    }
	  }
	}
      }
  }

  template <typename Iterator>
  std::pair<int, int> rule_statistics(const hypergraph_type& graph,
				      Iterator& iter,
				      Iterator last)
  {
    if (iter == last) return std::make_pair(0, 0);
    
    const hypergraph_type::edge_type& edge = graph.edges[*iter];
    ++ iter;
    
    int max_height = 1;
    int num_tails = edge.tails.size();
    //int num_tails = 0;
    
    hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
      if (iter != last && *titer == graph.edges[*iter].head) {
	const std::pair<int, int> result = rule_statistics(graph, iter, last);
	
	max_height = utils::bithack::max(max_height, result.first + 1);
	num_tails += result.second;
      }
    
    return std::make_pair(max_height, num_tails);
  }
    
  template <typename Iterator, typename Tails>
  std::pair<int, int> construct_tails(const hypergraph_type& graph,
				      Iterator& iter,
				      Iterator last,
				      Tails& tails)
  {
    if (iter == last) return std::make_pair(0, 0);
    
    const hypergraph_type::edge_type& edge = graph.edges[*iter];
    ++ iter;
    
    int max_height = 1;
    int num_tails = edge.tails.size();
    //int num_tails = 0;
    
    hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
      if (iter != last && *titer == graph.edges[*iter].head) {
	const std::pair<int, int> result = construct_tails(graph, iter, last, tails);
	
	max_height = utils::bithack::max(max_height, result.first + 1);
	num_tails += result.second;
      } else
	tails.push_back(*titer);
    }
    
    return std::make_pair(max_height, num_tails);
  }
  
  void admissible_nodes(const hypergraph_type& graph,
			const sentence_type& sentence,
			const alignment_type& alignment)
  {
    // span-edge
    span_edges.clear();
    span_edges.reserve(graph.edges.size());
    span_edges.resize(graph.edges.size());
    
    cicada::span_edge(graph, span_edges);
    
    // compute inside/outside score

    //std::cerr << "compute inside outside score" << std::endl;
    weights_inside.clear();
    weights_outside.clear();
    
    weights_inside.reserve(graph.nodes.size());
    weights_outside.reserve(graph.nodes.size());

    weights_inside.resize(graph.nodes.size());
    weights_outside.resize(graph.nodes.size());
    
    cicada::inside(graph, weights_inside, cicada::operation::weight_function_one<weight_type>());
    cicada::outside(graph, weights_inside, weights_outside, cicada::operation::weight_function_one<weight_type>());
    
    alignment_source_target.clear();
    alignment_target_source.clear();
    alignment_target_source.reserve(sentence.size());
    alignment_target_source.resize(sentence.size());
    
    if (inverse) {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {

	const int source = aiter->target;
	const int target = aiter->source;
	
	if (source >= static_cast<int>(alignment_source_target.size()))
	  alignment_source_target.resize(source + 1);
	if (target >= static_cast<int>(alignment_target_source.size()))
	  alignment_target_source.resize(target + 1);
	
	alignment_source_target[source].push_back(target);
	alignment_target_source[target].push_back(source);
      }
    } else {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
	
	if (aiter->source >= static_cast<int>(alignment_source_target.size()))
	  alignment_source_target.resize(aiter->source + 1);
	if (aiter->target >= static_cast<int>(alignment_target_source.size()))
	  alignment_target_source.resize(aiter->target + 1);
      
	alignment_source_target[aiter->source].push_back(aiter->target);
	alignment_target_source[aiter->target].push_back(aiter->source);
      }
    }
    
    // inside-outside!
    
    // first, bottom-up traversal to compute span...
    // we assume that every edge is annoated with its corresponding span information...
    // do we compute by span-forest again...?
    // we also assume that hypergrpah is parse-forest in that the each node share the same span and category...

    //std::cerr << "compute target spans" << std::endl;
    for (size_t id = 0; id != graph.nodes.size(); ++ id) {
      const hypergraph_type::node_type& node = graph.nodes[id];
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const range_type& edge_range = span_edges[*eiter];
	
	// copy range...
	ranges[id] = edge_range;

	if (edge_range.second - 1 >= static_cast<int>(alignment_source_target.size()))
	  alignment_source_target.resize(edge_range.second);

	for (int i = edge_range.first; i != edge_range.second; ++ i) {
	  point_set_type::const_iterator aiter_begin = alignment_source_target[i].begin();
	  point_set_type::const_iterator aiter_end   = alignment_source_target[i].end();
	  
	  for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	    spans[id].set(*aiter);
	}
      }
    }
    
    // second, top-down traversal to compute complement
    //std::cerr << "compute target complements" << std::endl;
    for (int id = graph.nodes.size() - 1; id >= 0; -- id) {

      //std::cerr << "complement: " << id << std::endl;
      const hypergraph_type::node_type& node = graph.nodes[id];
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	hypergraph_type::edge_type::node_set_type::const_iterator titer_begin = edge.tails.begin();
	hypergraph_type::edge_type::node_set_type::const_iterator titer_end   = edge.tails.end();
	
	for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  complements[*titer] |= complements[id];
	  for (hypergraph_type::edge_type::node_set_type::const_iterator iiter = titer_begin; iiter != titer_end; ++ iiter)
	    if (iiter != titer)
	      complements[*titer] |= spans[*iiter];
	}
      }
    }

    // check whether admissible or not...
    // do we also compute non-aligned words...?
    
    //std::cerr << "check admissible" << std::endl;
    for (size_t id = 0; id != graph.nodes.size(); ++ id)
      admissibles[id] = ! spans[id].empty() && ! complements[id].intersect(spans[id].range());
  }
};


struct Task
{
  typedef boost::filesystem::path path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef Bitext bitext_type;
  
  typedef ExtractGHKM::rule_pair_type     rule_pair_type;
  typedef ExtractGHKM::rule_pair_set_type rule_pair_set_type;
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> > queue_type;
  
  Task(queue_type& __queue,
       const path_type& __output,
       const std::string& non_terminal,
       const int max_sentence_length,
       const int max_nodes,
       const int max_height,
       const int max_compose,
       const int max_scope,
       const double cutoff,
       const bool exhaustive,
       const bool constrained,
       const bool inverse,
       const bool swap,
       const bool project,
       const bool collapse_source,
       const bool collapse_target,
       const double __max_malloc)
    : queue(__queue),
      output(__output),
      extractor(non_terminal, max_sentence_length, max_nodes, max_height, max_compose, max_scope, cutoff, exhaustive, constrained, inverse, swap, project, collapse_source, collapse_target),
      max_malloc(__max_malloc) {}
  
  queue_type&   queue;
  path_type     output;
  path_set_type paths;
  
  ExtractGHKM extractor;
  
  double max_malloc;

  struct Dumper
  {
    void operator()(rule_pair_set_type& rule_pairs) const
    {
      if (rule_pairs.size() < 1024 * 1024) return;
      
      const size_t allocated = utils::malloc_stats::used() + rule_pairs.size() * sizeof(void*);
      
      if ((min_counts_size && rule_pairs.size() < min_counts_size && allocated <= (malloc_threshold + (malloc_threshold >> 1)))
	  || allocated <= malloc_threshold) return;
      
      if (! min_counts_size)
	const_cast<size_t&>(min_counts_size) = rule_pairs.size() >> 5;
      
      dump(rule_pairs);
      
      rule_pairs.clear();
      rule_pair_set_type(rule_pairs).swap(rule_pairs);
    }
    
    typedef std::vector<const rule_pair_type*, std::allocator<const rule_pair_type*> > sorted_type;
    
    template <typename Tp>
    struct less_ptr
    {
      bool operator()(const Tp* x, const Tp* y) const
      {
	return *x < *y;
      }
    };
  
    void dump(const rule_pair_set_type& rule_pairs) const
    {
      if (rule_pairs.empty()) return;
    
      // sorting...
      sorted_type sorted(rule_pairs.size());
      {
	sorted_type::iterator siter = sorted.begin();
	rule_pair_set_type::const_iterator citer_end = rule_pairs.end();
	for (rule_pair_set_type::const_iterator citer = rule_pairs.begin(); citer != citer_end; ++ citer, ++ siter)
	  *siter = &(*citer);
      }
      std::sort(sorted.begin(), sorted.end(), less_ptr<rule_pair_type>());
    
      const path_type path_tmp = utils::tempfile::file_name(output / "counts-XXXXXX");
      utils::tempfile::insert(path_tmp);
      const path_type path = path_tmp.string() + ".gz";
      utils::tempfile::insert(path);
    
      utils::compress_ostream os(path, 1024 * 1024);
      sorted_type::const_iterator siter_end = sorted.end();
      for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
	generator(os, *(*siter)) << '\n';
    
      const_cast<path_set_type&>(paths).push_back(path);
    }
    
    Dumper(const path_type& __output,
	   path_set_type& __paths,
	   const size_t __malloc_threshold)
      : output(__output),
	paths(__paths),
	malloc_threshold(__malloc_threshold),
	min_counts_size(0) {}
    
    const path_type& output;
    path_set_type& paths;
    const size_t malloc_threshold;
    size_t min_counts_size;
    
    RulePairGenerator generator;
  };
  
  void operator()()
  {
    bitext_type bitext;
    rule_pair_set_type rule_pairs;
    
    Dumper dumper(output, paths, max_malloc * 1024 * 1024 * 1024);

    const int iter_mask = (1 << 4) - 1;
    
    for (int iter = 0;/**/; ++ iter) {
      bitext.clear();
      queue.pop_swap(bitext);
      
      if (! bitext.source.is_valid()) break;
      
      extractor(bitext.source, bitext.target, bitext.alignment, rule_pairs, dumper);
      
      if (rule_pairs.empty() || (iter & iter_mask) == iter_mask)
	extractor.clear();
    }
    
    dumper.dump(rule_pairs);
  }
};

#endif
