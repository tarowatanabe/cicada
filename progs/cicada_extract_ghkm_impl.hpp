#ifndef __CICADA__EXTRACT_GHKM_IMPL__HPP__
#define __CICADA__EXTRACT_GHKM_IMPL__HPP__ 1


#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <set>

#include <boost/array.hpp>

#include "cicada/sentence.hpp"
#include "cicada/alignment.hpp"
#include "cicada/vocab.hpp"
#include "cicada/tree_rule.hpp"
#include "cicada/hypergraph.hpp"

#include "utils/sgi_hash_map.hpp"
#include "utils/sgi_hash_set.hpp"
#include "utils/chart.hpp"

#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/chunk_vector.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

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
  
  friend
  std::istream& operator>>(std::istream& is, Bitext& bitext)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    using qi::phrase_parse;
    using qi::lit;
    using standard::space;

    bitext.clear();
    std::string line;
    if (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end  = line.end();
      
      if((!bitext.source.assign(iter, end))
	 || (!phrase_parse(iter, end, lit("|||"), space))
	 || (!bitext.target.assign(iter, end))
	 || (!phrase_parse(iter, end, lit("|||"), space))
	 || (!bitext.alignment.assign(iter, end))
	 || iter != end)
	bitext.clear();
    }
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
  
  phrase_type    source;
  phrase_type    target;
  alignment_type alignment;
  count_type     count;

  RulePair() : source(), target(), alignment(), count(0) {}

  friend
  size_t hash_value(RulePair const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    return hasher_type()(x.source.begin(), x.source.end(),
			 hasher_type()(x.target.begin(), x.target.end(),
				       hasher_type()(x.alignment.begin(), x.alignment.end(), 0)));
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
    return (x.source < y.source
	    || (!(y.source < x.source)
		&& (x.target < y.target
		    || (!(y.target < x.target)
			&& x.alignment < y.alignment))));
  }
  
  friend
  bool operator>(const RulePair& x, const RulePair& y)
  {
    return y < x;
  }
};

BOOST_FUSION_ADAPT_STRUCT(RulePair::alignment_type::point_type,
			  (RulePair::alignment_type::index_type, source)
			  (RulePair::alignment_type::index_type, target)
			  )
BOOST_FUSION_ADAPT_STRUCT(
			  RulePair,
			  (RulePair::phrase_type, source)
			  (RulePair::phrase_type, target)
			  (RulePair::alignment_type, alignment)
			  (RulePair::count_type, count)
			  )

struct RulePairGenerator
{
  typedef RulePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::count_type    count_type;
  
  RulePairGenerator() : grammar() {}
  RulePairGenerator(const RulePairGenerator& x) : grammar() {}

  
  template <typename Iterator>
  struct phrase_pair_generator : boost::spirit::karma::grammar<Iterator, phrase_pair_type()>
  {
    phrase_pair_generator() : phrase_pair_generator::base_type(phrase_pair)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::repeat;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      using standard::space;
      
      phrase %= +char_;
      alignment %= -((int_ << '-' << int_) % ' ');
      phrase_pair %= phrase << " ||| " << phrase << " ||| " << alignment << " ||| " << double20;
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, std::string()> phrase;
    boost::spirit::karma::rule<Iterator, alignment_type()> alignment;
    boost::spirit::karma::rule<Iterator, phrase_pair_type()> phrase_pair;
  };

  typedef std::ostream_iterator<char> iterator_type;
  
  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair)
  {
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, phrase_pair))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  phrase_pair_generator<iterator_type> grammar;
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


struct ExtractGHKM
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef RulePair rule_pair_type;
  typedef rule_pair_type::phrase_type phrase_type;
  typedef rule_pair_type::count_type  count_type;
  
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
				  std::allocator<rule_pair_type> > rule_pair_set_type;
#else
 typedef sgi::hash_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
		      std::allocator<rule_pair_type> > rule_pair_set_type;
#endif
  
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

  struct Span
  {
    typedef std::set<int, std::less<int>, std::allocator<int> > span_type;
    
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
  
  struct DerivationEdge
  {
    edge_set_type edges;
    node_set_type tails;

    DerivationEdge() : edges(), tails() {}
    DerivationEdge(const edge_set_type& __edges, const node_set_type& __tails) : edges(__edges), tails(__tails) {}
  };
  
  struct DerivationNode
  {
    typedef DerivationEdge derivation_edge_type;
    typedef std::vector<derivation_edge_type, std::allocator<derivation_edge_type> > edge_set_type;

    id_type       node;
    range_type    range;
    edge_set_type edges;
  };
  
  typedef DerivationEdge derivation_edge_type;
  typedef DerivationNode derivation_node_type;
  
  typedef std::deque<derivation_node_type, std::allocator<derivation_node_type> > derivation_set_type;

  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > alignment_map_type;
  
  
  ExtractGHKM(const symbol_type& __non_terminal,
	      const int __max_nodes,
	      const int __max_height)
    : non_terminal(__non_terminal),
      max_nodes(__max_nodes),
      max_height(__max_height) {}

  symbol_type non_terminal;
  int max_nodes;
  int max_height;
  
  range_set_type ranges;
  span_set_type spans;
  span_set_type complements;
  span_set_type unaligned;

  admissible_set_type admissibles;
  
  node_map_type node_map;
  derivation_set_type derivations;

  alignment_map_type alignment_source_target;
  alignment_map_type alignment_target_source;
  
  void operator()(const hypergraph_type& graph,
		  const sentence_type& sentence,
		  const alignment_type& alignment,
		  rule_pair_set_type& rules)
  {
    ranges.clear();
    spans.clear();
    complements.clear();
    unaligned.clear();
    admissibles.clear();
    derivations.clear();
    node_map.clear();
    
    ranges.resize(graph.nodes.size());
    spans.resize(graph.nodes.size());
    complements.resize(graph.nodes.size());
    unaligned.resize(graph.nodes.size());
    node_map.resize(graph.nodes.size());
    
    admissible_nodes(graph, sentence, alignment);
    
    construct_derivations(graph, sentence);
    
    // compute reachable nodes and edges from root...
    prune_derivations();
    
    // perform extraction
    extract_minimum(graph, sentence, alignment, rules);
    
    // perform compounds extraction
    extract_composed(graph, sentence, alignment, rules);
  }

  template <typename Tp>
  struct less_first
  {
    bool operator()(const Tp& x, const Tp& y) const
    {
      return x.first < y.first;
    }
  };
  
  
  void extract_minimum(const hypergraph_type& graph,
		       const sentence_type& sentence,
		       const alignment_type& alignment,
		       rule_pair_set_type& rule_pairs)
  {
    // easy... but do we "propagate" edges into parent?

    typedef std::pair<range_type, int> range_pos_type;
    typedef std::vector<range_pos_type, std::allocator<range_pos_type> > range_pos_set_type;
    typedef std::vector<tree_rule_type, std::allocator<tree_rule_type> > tree_rule_set_type;
    typedef std::vector<bool, std::allocator<bool> > covered_type;
    
    range_pos_set_type range_pos;
    tree_rule_set_type trees;

    point_set_type positions_source(alignment_source_target.size());
    point_set_type positions_target(alignment_target_source.size());
    covered_type covered(alignment_source_target.size());

    boost::iostreams::filtering_ostream os_source;
    boost::iostreams::filtering_ostream os_target;

    rule_pair_type rule_pair;
    
    for (size_t id = 0; id != derivations.size(); ++ id) {
      const derivation_node_type& node = derivations[id];
      
      derivation_node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (derivation_node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_set_type& edges = eiter->edges;
	const node_set_type& tails = eiter->tails;
	
	// we will re-construct composed rule from edges...
	// edges is ordered by pre-traversal order...
	tree_rule_type rule_source;
	{
	  positions_source.clear();
	  std::fill(covered.begin(), covered.end(), false);
	  int frontier_pos = 0;
	  int index = 1;
	  edge_set_type::const_iterator iter = edges.begin();
	  edge_set_type::const_iterator iter_end = edges.end();
	  construct_rule(graph, iter, iter_end, rule_source, index, frontier_pos, positions_source, covered);
	}
	
	// construct target side rule...
	// tails is ordered by graph... thus, reordering is already replresented by derivations[tails[i]].range...
	
	trees.clear();
	
	if (tails.empty()) {
	  trees.insert(trees.end(), sentence.begin() + node.range.first, sentence.begin() + node.range.second);
	  for (int i = node.range.first; i != node.range.second; ++ i)
	    positions_target[i] = i - node.range.first;
	} else {
	  range_pos.clear();
	  int index = 1;
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
	    trees.push_back(non_terminal.non_terminal(riter->second));
	    
	    pos_first = riter->first.second;
	  }
	  
	  const int pos_last = node.range.second;
	  
	  int mapped_pos = trees.size();
	  for (int i = pos_first; i != pos_last; ++ i, ++ mapped_pos)
	    positions_target[i] = mapped_pos;
	  
	  trees.insert(trees.end(), sentence.begin() + pos_first, sentence.begin() + pos_last);
	}
	
	tree_rule_type rule_target(non_terminal, trees.begin(), trees.end());
	
	//
	// construct word alignment...
	// we can easily reconsruct from positions_{source, target} + covered
	//
	
	rule_pair.alignment.clear();
	int pos_src = 0;

	for (size_t src = 0; src != covered.size(); ++ src)
	  if (covered[src]) {
	    
	    point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	    point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	    
	    for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	      rule_pair.alignment.push_back(std::make_pair(positions_source[pos_src], positions_target[*aiter]));
	    
	    ++ pos_src;
	  }
	
	rule_pair.source.clear();
	rule_pair.target.clear();
	
	os_source.push(boost::iostreams::back_inserter(rule_pair.source));
	os_target.push(boost::iostreams::back_inserter(rule_pair.target));
	
	os_source << rule_source;
	os_target << rule_target;
	
	os_source.pop();
	os_target.pop();
	
	// increment by one... actually, we need to normalize...
	const_cast<rule_pair_type&>(*(rule_pairs.insert(rule_pair).first)).count += 1;
      }
    }
  }
  
  template <typename Iterator, typename PosMap, typename Covered>
  void construct_rule(const hypergraph_type& graph,
		      Iterator& iter,
		      Iterator last,
		      tree_rule_type& tree_rule,
		      int& index,
		      int& frontier_pos,
		      PosMap& pos_map,
		      Covered& covered)
  {
    // pre-order traversal...
    
    if (iter == last) {
      tree_rule.label = tree_rule.label.non_terminal();
      return;
    }
    
    const hypergraph_type::edge_type& edge = graph.edges[*iter];
    
    for (int pos = edge.first; pos != edge.last; ++ pos)
      covered[pos] = true;
    
    tree_rule = tree_rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end());
    
    ++ iter;
    
    if (iter != last) {
      const hypergraph_type::edge_type& edge_next = graph.edges[*iter];
      
      size_t tail_pos = 0;
      tree_rule_type::iterator titer_end = tree_rule.end();
      for (tree_rule_type::iterator titer = tree_rule.begin(); titer != titer_end; ++ titer) {
	if (titer->label.is_non_terminal()) {

	  const id_type node_id = edge.tails[tail_pos];
	  
	  if (node_id == edge_next.head)
	    construct_rule(graph, iter, last, *titer, index, frontier_pos, pos_map, covered);
	  else {
	    for (int pos = ranges[node_id].first; pos != ranges[node_id].second; ++ pos)
	      covered[pos] = false;
	    
	    titer->label = titer->label.non_terminal(index ++);
	    ++ frontier_pos;
	  }
	  
	  ++ tail_pos;
	} else {
	  pos_map.push_back(frontier_pos);
	  ++ frontier_pos;
	}
      }
    } else {
      size_t tail_pos = 0;
      tree_rule_type::iterator titer_end = tree_rule.end();
      for (tree_rule_type::iterator titer = tree_rule.begin(); titer != titer_end; ++ titer, ++ frontier_pos) {
	if (titer->label.is_non_terminal()) {
	  const id_type node_id = edge.tails[tail_pos];
	  
	  for (int pos = ranges[node_id].first; pos != ranges[node_id].second; ++ pos)
	    covered[pos] = false;

	  titer->label = titer->label.non_terminal(index ++);
	  ++ tail_pos;
	} else
	  pos_map.push_back(frontier_pos);
      }
    }
  }
  
  void extract_composed(const hypergraph_type& graph,
			const sentence_type& sentence,
			const alignment_type& alignment,
			rule_pair_set_type& rule_pairs)
  {
    // difficult...
    // we need to keep track of how many nodes / maximum height in the composed ruels.
    
    
    
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
	    is_valid = tails.back() >= 0;
	  }
	  
	  if (is_valid)
	    node_new.edges.push_back(derivation_edge_type(eiter->edges, tails));
	}
      }
    
    derivations.swap(derivations_new);
    derivations_new.clear();
  }
  
  void construct_derivations(const hypergraph_type& graph,
			     const sentence_type& sentence)
  {
    typedef std::vector<int, std::allocator<int> > index_set_type;
    typedef std::deque<frontier_type, std::allocator<frontier_type> > queue_type;
    typedef google::dense_hash_map<range_type, id_type, utils::hashmurmur<size_t>, std::equal_to<range_type> > range_node_map_type;
    typedef google::dense_hash_set<range_type, utils::hashmurmur<size_t>, std::equal_to<range_type> > range_set_type;

    // construc derivations wrt non-aligned words...

    size_t goal_node = size_t(-1);
    
    for (size_t id = 0; id != graph.nodes.size(); ++ id) 
      if (admissibles[id]) {
	const hypergraph_type::node_type& node = graph.nodes[id];

	//
	// compute minimal range and maximum outer-range
	//
	
	const range_type range_min(spans[id].range());
	
	span_type::const_iterator riter_first = complements[id].lower_bound(range_min.first);
	span_type::const_iterator riter_last  = complements[id].lower_bound(range_min.second);
	
	const range_type range_max(riter_first == complements[id].begin() ? 0 : *(-- riter_first) + 1,
				   riter_last != complements[id].end() ? *riter_last : static_cast<int>(sentence.size()));
	
	const bool is_goal(id == graph.goal);
	
	range_node_map_type buf;
	buf.set_empty_key(range_type(0, 0));
	
	queue_type queue;
	
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

	while (! queue.empty()) {
	  const frontier_type& frontier = queue.front();
	  
	  if (frontier.second.empty()) {
	    // construct rule from "fragment", frontier.first
	    // by enumerating edge and its tails in depth-first manner...
	    // use recursive call for simplicity...
	    
	    // compute "tails" for the fragment
	    
	    node_set_type tails;
	    edge_set_type::const_iterator eiter_end = frontier.first.end();
	    for (edge_set_type::const_iterator eiter = frontier.first.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
		if (admissibles[*titer])
		  tails.push_back(*titer);
	    }
	    
	    // iterate over tails and compute span...
	    // code taken from apply_exact...
	    
	    index_set_type j_ends(tails.size(), 0);
	    index_set_type j(tails.size(), 0);
	    
	    for (size_t i = 0; i != tails.size(); ++ i)
	      j_ends[i] = node_map[tails[i]].size();
	    
	    node_set_type tails_next(tails.size());

	    range_set_type ranges;
	    ranges.set_empty_key(range_type(0, 0));
	    
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
		  
		  range.first  = utils::bithack::min(range.first, range_antecedent.first);
		  range.second = utils::bithack::max(range.second, range_antecedent.second);
		}
	      }
	      
	      if (is_valid) {
		// our tail-ranges are valid... and we will consider alternatives from range to range_max...
		
		if (is_goal) {
		  if (goal_node == size_t(-1)) {
		    goal_node = derivations.size();
		    derivations.resize(goal_node + 1);
		    
		    // range is always [0, sentence.size())
		    derivations.back().node = id;
		    derivations.back().range = range_type(0, sentence.size());
		  }
		  
		  derivations[goal_node].edges.push_back(derivation_edge_type(frontier.first, tails_next));
		} else {
		  // we will compute all possible ranges...
		  for (int first = range.first; first <= range_max.first; ++ first)
		    for (int last = range.second; last <= range_max.second; ++ last) {
		      const range_type range_next(first, last);
		      
		      range_node_map_type::iterator biter = buf.find(range_next);
		      if (biter == buf.end()) {
			derivations.resize(derivations.size() + 1);
			
			derivations.back().node = id;
			derivations.back().range = range_next;
			
			node_map[id].push_back(derivations.size() - 1);
			
			biter = buf.insert(std::make_pair(range_next, derivations.size() - 1)).first;
		      }
		      
		      derivations[biter->second].edges.push_back(derivation_edge_type(frontier.first, tails_next));
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
	    const hypergraph_type::node_type& node = graph.nodes[frontier.second.front()];
	    
	    hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	      
	      queue.resize(queue.size() + 1);
	      frontier_type& frontier_next = queue.back();
	      
	      frontier_next.first = frontier.first;
	      frontier_next.first.push_back(*eiter);
	      
	      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
		if (! admissibles[*titer])
		  frontier_next.second.push_back(*titer);
	      
	      frontier_next.second.insert(frontier_next.second.end(), frontier.second.begin() + 1, frontier.second.end());
	    }
	  }
	  
	  queue.pop_front();
	}
      }
  }
  
  void admissible_nodes(const hypergraph_type& graph,
			const sentence_type& sentence,
			const alignment_type& alignment)
  {
    alignment_source_target.clear();
    alignment_target_source.clear();
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      
      if (aiter->source >= static_cast<int>(alignment_source_target.size()))
	alignment_source_target.resize(aiter->source + 1);
      if (aiter->target >= static_cast<int>(alignment_target_source.size()))
	alignment_target_source.resize(aiter->target + 1);
      
      alignment_source_target[aiter->source].push_back(aiter->target);
      alignment_target_source[aiter->target].push_back(aiter->source);
    }
    
    // first, bottom-up traversal to compute span...
    // we assume that every edge is annoated with its corresponding span information...
    // do we compute by span-forest again...?
    // we also assume that hypergrpah is parse-forest in that the each node share the same span and category...
    for (size_t id = 0; id != graph.nodes.size(); ++ id) {
      const hypergraph_type::node_type& node = graph.nodes[id];
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	// copy range...
	ranges[id] = range_type(edge.first, edge.last);

	if (edge.last - 1 >= static_cast<int>(alignment_source_target.size()))
	  alignment_source_target.resize(edge.last);

	for (int i = edge.first; i != edge.last; ++ i) {
	  point_set_type::const_iterator aiter_begin = alignment_source_target[i].begin();
	  point_set_type::const_iterator aiter_end   = alignment_source_target[i].end();
	  
	  for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	    spans[id].set(*aiter);
	}
      }
    }
    
    // second, top-down traversal to compute complement
    for (int id = graph.nodes.size(); id >= 0; -- id) {
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
       const int max_nodes,
       const int max_height,
       const double __max_malloc)
    : queue(__queue),
      output(__output),
      extractor(non_terminal, max_nodes, max_height),
      max_malloc(__max_malloc) {}
  
  queue_type&   queue;
  path_type     output;
  path_set_type paths;
  
  ExtractGHKM extractor;
  RulePairGenerator generator;
  
  double max_malloc;
  
  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };
  
  void dump(const rule_pair_set_type& rule_pairs)
  {
    typedef std::vector<const rule_pair_type*, std::allocator<const rule_pair_type*> > sorted_type;
    
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
    const path_type path = path_tmp.file_string() + ".gz";
    utils::tempfile::insert(path);
    
    utils::compress_ostream os(path, 1024 * 1024);
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      generator(os, *(*siter)) << '\n';
    
    paths.push_back(path);
  }
  
  void operator()()
  {
    bitext_type bitext;
    rule_pair_set_type rule_pairs;
    
    const int iteration_mask = (1 << 10) - 1;
    
    for (int iter = 0;/**/; ++ iter) {
      queue.pop_swap(bitext);
      
      if (! bitext.source.is_valid()) break;
      
      extractor(bitext.source, bitext.target, bitext.alignment, rule_pairs);
      
      if (((iter & iteration_mask) == iteration_mask) && (utils::malloc_stats::used() > size_t(max_malloc * 1024 * 1024 * 1024))) {
	dump(rule_pairs);
	rule_pairs.clear();
      }
    }
    
    if (! rule_pairs.empty()) {
      dump(rule_pairs);
      rule_pairs.clear();
    }
  }
};

#endif
