//
//  Copyright(C) 2011 Graham Neubig <neubig@ar.media.kyoto-u.ac.jp>
//                    Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// Extract hiero rule from pialign derivation
//
// Given ITG, we have only to traverse children and compute holes.
// We will always insert terminal in source-side, but not target-side..
// How to handle counts? Do we use original-hiero-style fractional counts?
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>
#include <memory>
#include <deque>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/remove_epsilon.hpp>
#include <cicada/remove_unary.hpp>
#include <cicada/tree_rule.hpp>

#include <utils/compress_stream.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/std_heap.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


typedef boost::filesystem::path path_type;
typedef std::pair<int, int> span_type;
struct span_pair_type
{
  span_type source;
  span_type target;
  
  span_pair_type()
    : source(), target() {}
  span_pair_type(const span_type& __source, const span_type& __target)
    : source(__source), target(__target) {}
};

typedef cicada::Vocab    vocab_type;
typedef cicada::Symbol   word_type;
typedef cicada::Sentence sentence_type;

struct itg_type
{
  typedef std::vector<itg_type, std::allocator<itg_type> > itg_pair_type;
  
  std::string source;
  std::string target;
  
  itg_pair_type antecedent;
  
  span_pair_type spans;
  
  bool inverse;
  bool block;

  itg_type() : source(), target(), antecedent(), inverse(false), block(false) {}

  void clear()
  {
    source.clear();
    target.clear();
    antecedent.clear();
    
    spans = span_pair_type();

    inverse = false;
    block = false;
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  itg_type,
			  (std::string, source)
			  (std::string, target)
			  (itg_type::itg_pair_type, antecedent)
			  (bool, inverse)
			  (bool, block)
			  )

template <typename Iterator>
struct derivation_parser : boost::spirit::qi::grammar<Iterator, itg_type(), boost::spirit::standard::space_type>
{
  derivation_parser() : derivation_parser::base_type(derivation)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    terminal %= qi::lexeme[+(standard::char_ - standard::space) - "(((" - ")))" - "|||"] | qi::attr("");
    
    derivation_terminal %= "(((" >> terminal >> "|||" >> terminal >> ")))";
    derivation_regular  %= '[' >> qi::attr("") >> qi::attr("") >> derivation_pair >> qi::attr(false) >> ']';
    derivation_inverse  %= '<' >> qi::attr("") >> qi::attr("") >> derivation_pair >> qi::attr(true) >> '>';
    derivation_pair     %= qi::repeat(2)[derivation];
    
    derivation_block %= ('{' >> (qi::hold[derivation_regular] | qi::hold[derivation_inverse] | derivation_terminal) >> '}')[phoenix::at_c<4>(qi::_val) = true];
    derivation_break  %= (qi::hold[derivation_regular] | qi::hold[derivation_inverse] | derivation_terminal)[phoenix::at_c<4>(qi::_val) = false];
    
    derivation %= qi::hold[derivation_block] | derivation_break;
    
    //qi::on_error<qi::fail>(derivation, std::cerr << qi::_4 << phoenix::val(": ") << phoenix::construct<std::string>(qi::_3, qi::_2) << std::endl);
  }
  
  typedef boost::spirit::standard::space_type space_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), space_type>  terminal;

  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_block;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_break;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_terminal;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_regular;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_inverse;
  
  boost::spirit::qi::rule<Iterator, itg_type::itg_pair_type(), space_type> derivation_pair;
};

std::ostream& print_tree(std::ostream& os, const itg_type& itg)
{
  if (itg.block)
    os << '{';

#if 0
  os << ' ' << itg.spans.source.first << ".." << itg.spans.source.second
     << ' ' << itg.spans.target.first << ".." << itg.spans.target.second
     << ' ';
#endif

  if (itg.antecedent.empty())
    os << "((( " << itg.source << " ||| " << itg.target << " )))";
  else {
    os << (itg.inverse ? '<' : '[');
    print_tree(os, itg.antecedent.front());
    print_tree(os, itg.antecedent.back());
    os << (itg.inverse ? '>' : ']');
  }
  if (itg.block)
    os << '}';
  
  return os;
}


void span_derivation_source(itg_type& itg, sentence_type& sentence)
{
  itg.spans.source.first = sentence.size();
  if (itg.antecedent.empty()) {
    if (! itg.source.empty())
      sentence.push_back(itg.source);
  } else {
    span_derivation_source(itg.antecedent.front(), sentence);
    span_derivation_source(itg.antecedent.back(),  sentence);
  }
  itg.spans.source.second = sentence.size();
}

void span_derivation_target(itg_type& itg, sentence_type& sentence)
{
  itg.spans.target.first = sentence.size();
  if (itg.antecedent.empty()) {
    if (! itg.target.empty())
      sentence.push_back(itg.target);
  } else {
    if (! itg.inverse) {
      span_derivation_target(itg.antecedent.front(), sentence);
      span_derivation_target(itg.antecedent.back(),  sentence);
    } else {
      span_derivation_target(itg.antecedent.back(),  sentence);
      span_derivation_target(itg.antecedent.front(), sentence);
    }
  }
  itg.spans.target.second = sentence.size();
}

struct Grammar
{
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  
  typedef std::set<int, std::less<int>, std::allocator<int> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > alignment_type;
  
  typedef std::pair<int, int> point_type;
  typedef std::vector<point_type, std::allocator<point_type> > point_set_type;
};

struct TreeSource : public Grammar
{
  typedef cicada::HyperGraph hypergraph_type;
  
  typedef hypergraph_type::rule_type     rule_type;
  typedef hypergraph_type::rule_ptr_type rule_ptr_type;

  TreeSource(std::ostream& __os,
	     const bool __remove_epsilon,
	     const bool __remove_unary)  
    :  os(__os),
       remove_epsilon(__remove_epsilon),
       remove_unary(__remove_unary)
  {
    rule_type::symbol_set_type rhs(2, vocab_type::X);
    
    rule_binary = rule_type::create(rule_type(vocab_type::X, rhs));
  }
  
  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    alignment.clear();
    alignment.resize(source.size());
    
    hypergraph_type forest;
    
    forest.goal = operator()(itg, source, target, forest, alignment, blocker);
    
    if (forest.is_valid()) {
      forest.topologically_sort();
      if (remove_epsilon)
	cicada::remove_epsilon(forest);
      if (remove_unary)
	cicada::remove_unary(forest);
    }
    
    os << forest << '\n';
  }
  
  template <typename Blocker>
  hypergraph_type::id_type operator()(const itg_type& itg,
				      const sentence_type& source,
				      const sentence_type& target,
				      hypergraph_type& forest,
				      alignment_type& alignment,
				      Blocker blocker)
  {
    if (blocker(itg) || itg.antecedent.empty())  {
      rule_type::symbol_set_type rhs(source.begin() + itg.spans.source.first, source.begin() + itg.spans.source.second);
      
      hypergraph_type::edge_type& edge = forest.add_edge();
      edge.rule = rule_type::create(rhs.empty()
				    ? rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1)
				    : rule_type(vocab_type::X, rhs));
      
      const hypergraph_type::id_type node_id = forest.add_node().id;
      
      forest.connect_edge(edge.id, node_id);
      
      // compute alignment matrix
      if (itg.spans.target.first != itg.spans.target.second)
	for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src)
	  for (int trg = itg.spans.target.first; trg != itg.spans.target.second; ++ trg)
	    alignment[src].insert(trg);
      
      // return parent...
      return node_id;
    } else {
      hypergraph_type::id_type tails[2];
      
      tails[0] = operator()(itg.antecedent.front(), source, target, forest, alignment, blocker);
      tails[1] = operator()(itg.antecedent.back(),  source, target, forest, alignment, blocker);
      
      hypergraph_type::edge_type& edge = forest.add_edge(tails, tails + 2);
      edge.rule = rule_binary;
      
      const hypergraph_type::id_type node_id = forest.add_node().id;
      
      forest.connect_edge(edge.id, node_id);
      
      return  node_id;
    }
  }
  
  std::ostream& os;
  const bool remove_epsilon;
  const bool remove_unary;
  
  rule_ptr_type rule_binary;
};

struct TreeTarget : public Grammar
{
  typedef cicada::HyperGraph hypergraph_type;
  
  typedef hypergraph_type::rule_type     rule_type;
  typedef hypergraph_type::rule_ptr_type rule_ptr_type;

  TreeTarget(std::ostream& __os,
	     const bool __remove_epsilon,
	     const bool __remove_unary)  
    :  os(__os),
       remove_epsilon(__remove_epsilon),
       remove_unary(__remove_unary)
  {
    rule_type::symbol_set_type rhs(2, vocab_type::X);
    
    rule_binary = rule_type::create(rule_type(vocab_type::X, rhs));
  }
  
  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    alignment.clear();
    alignment.resize(source.size());
    
    hypergraph_type forest;
    
    forest.goal = operator()(itg, source, target, forest, alignment, blocker);
    
    if (forest.is_valid()) {
      forest.topologically_sort();
      if (remove_epsilon)
	cicada::remove_epsilon(forest);
      if (remove_unary)
	cicada::remove_unary(forest);
    }
    
    os << forest << '\n';
  }
  
  template <typename Blocker>
  hypergraph_type::id_type operator()(const itg_type& itg,
				      const sentence_type& source,
				      const sentence_type& target,
				      hypergraph_type& forest,
				      alignment_type& alignment,
				      Blocker blocker)
  {
    if (blocker(itg) || itg.antecedent.empty())  {
      rule_type::symbol_set_type rhs(target.begin() + itg.spans.target.first, target.begin() + itg.spans.target.second);
      
      hypergraph_type::edge_type& edge = forest.add_edge();
      edge.rule = rule_type::create(rhs.empty()
				    ? rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1)
				    : rule_type(vocab_type::X, rhs));
      
      const hypergraph_type::id_type node_id = forest.add_node().id;
      
      forest.connect_edge(edge.id, node_id);
      
      // compute alignment matrix
      if (itg.spans.target.first != itg.spans.target.second)
	for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src)
	  for (int trg = itg.spans.target.first; trg != itg.spans.target.second; ++ trg)
	    alignment[src].insert(trg);
      
      // return parent...
      return node_id;
    } else {
      hypergraph_type::id_type tails[2];
      
      tails[  itg.inverse] = operator()(itg.antecedent.front(), source, target, forest, alignment, blocker);
      tails[! itg.inverse] = operator()(itg.antecedent.back(),  source, target, forest, alignment, blocker);
      
      hypergraph_type::edge_type& edge = forest.add_edge(tails, tails + 2);
      edge.rule = rule_binary;
      
      const hypergraph_type::id_type node_id = forest.add_node().id;
      
      forest.connect_edge(edge.id, node_id);
      
      return  node_id;
    }
  }
  
  std::ostream& os;
  const bool remove_epsilon;
  const bool remove_unary;
  
  rule_ptr_type rule_binary;
};

struct GHKMGrammar : public Grammar
{
  //
  // transform itg_type into a paired hypergraph structure by traversing in
  // post-traversal order.
  // (Thus, we need not construct hypergraph explicitly...?)
  //
  
  //
  // we need to map the itg-derivation-post-order-traversal id into spans and rules
  //
  typedef cicada::HyperGraph hypergraph_type;
  
  typedef cicada::TreeRule   tree_rule_type;
  
  typedef hypergraph_type::rule_type     rule_type;
  typedef hypergraph_type::rule_ptr_type rule_ptr_type;

  typedef hypergraph_type::feature_set_type   feature_set_type;
  typedef hypergraph_type::attribute_set_type attribute_set_type;

  typedef attribute_set_type::attribute_type attribute_type;
  typedef feature_set_type::feature_type     feature_type;

  typedef std::pair<hypergraph_type::id_type, hypergraph_type::id_type> id_pair_type;
  
  typedef std::vector<size_t, std::allocator<size_t> > admissible_set_type;
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;
  
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
  
  struct DerivationEdge
  {
    edge_set_type edges;
    tail_set_type tails;
    
    int height;
    int internal;
    int compose;
    
    DerivationEdge() : edges(), tails(), height(0), internal(0), compose(1) {}
    DerivationEdge(const edge_set_type& __edges,
		   const tail_set_type& __tails,
		   const int& __height,
		   const int& __internal)
      : edges(__edges), tails(__tails), height(__height), internal(__internal), compose(1) {}
    DerivationEdge(const edge_set_type& __edges,
		   const tail_set_type& __tails,
		   const int& __height,
		   const int& __internal,
		   const int& __compose)
      : edges(__edges), tails(__tails), height(__height), internal(__internal), compose(__compose) {}

  };
  
  typedef DerivationEdge derivation_edge_type;
  typedef std::deque<derivation_edge_type, std::allocator<derivation_edge_type> > derivation_edge_set_type;
  typedef std::vector<derivation_edge_set_type, std::allocator<derivation_edge_set_type> > derivation_graph_type;
  
  struct DerivationPair
  {
    derivation_edge_type source;
    derivation_edge_type target;

    DerivationPair() : source(), target() {}
    DerivationPair(const derivation_edge_type& __source,
		   const derivation_edge_type& __target)
      : source(__source), target(__target) {}
  };
  typedef DerivationPair derivation_pair_type;
  typedef std::deque<derivation_pair_type, std::allocator<derivation_pair_type> > derivation_pair_set_type;
  typedef std::vector<derivation_pair_set_type, std::allocator<derivation_pair_set_type> > derivation_pair_graph_type;
  
  struct less_derivation_pair_type
  {
    bool operator()(const derivation_pair_type& x, const derivation_pair_type& y) const
    {
      return (x.source.compose < y.source.compose
	      && (!(y.source.compose < x.source.compose)
		  || (x.source.internal < y.source.internal
		      && (!(y.source.internal < x.source.internal)
			  && x.target.internal < y.target.internal))));
    }
  };

  struct rule_pair_type
  {
    tree_rule_type source;
    tree_rule_type target;
    point_set_type alignment;
    
    rule_pair_type() : source(), target(), alignment() {}
    
    void clear()
    {
      source.clear();
      target.clear();
      alignment.clear();
    }

    friend
    std::ostream& operator<<(std::ostream& os, const rule_pair_type& x)
    {
      os << x.source << " ||| " << x.target << " |||";
      point_set_type::const_iterator piter_end = x.alignment.end();
      for (point_set_type::const_iterator piter = x.alignment.begin(); piter != piter_end; ++ piter)
	os << ' ' << piter->first << '-' << piter->second;
      return os;
    }
  };

  
  GHKMGrammar(std::ostream& __os,
	      const int __max_nodes,
	      const int __max_height,
	      const int __max_compose,
	      const int __max_scope,
	      const bool __frontier_source,
	      const bool __frontier_target,
	      const bool __remove_unary)
    : os(__os),
      max_nodes(__max_nodes), max_height(__max_height),
      max_compose(__max_compose), max_scope(__max_scope),
      frontier_source(__frontier_source),
      frontier_target(__frontier_target),
      remove_unary(__remove_unary),
      attr_node_id("node-id")
  {
    rule_type::symbol_set_type rhs(2, vocab_type::X);
    
    rule_binary = rule_type::create(rule_type(vocab_type::X, rhs));
  }
  
  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    alignment.clear();
    alignment.resize(source.size());
    
    // first, compute epsilon-removed tree
    
    spans.clear();
    graph_source.clear();
    graph_target.clear();
    
    const id_pair_type id_pair = operator()(itg, source, target, spans, graph_source, graph_target, alignment, blocker);
    
    graph_source.goal = id_pair.first;
    graph_target.goal = id_pair.second;
    
    admissible_source.clear();
    admissible_target.clear();
    admissible_source.resize(spans.size(), 0);
    admissible_target.resize(spans.size(), 0);
    
    if (graph_source.is_valid()) {
      graph_source.topologically_sort();
      cicada::remove_epsilon(graph_source);
      if (remove_unary)
	cicada::remove_unary(graph_source);
      compute_node_map(graph_source, admissible_source, nodes_map_source);
    }
    if (graph_target.is_valid()) {
      graph_target.topologically_sort();
      cicada::remove_epsilon(graph_target);
      if (remove_unary)
	cicada::remove_unary(graph_target);
      compute_node_map(graph_target, admissible_target, nodes_map_target);
    }
    
    if (! graph_source.is_valid() || ! graph_target.is_valid()) return;

    derivation_source.clear();
    derivation_target.clear();
    derivation_source.resize(spans.size());
    derivation_target.resize(spans.size());
    
    // construct derivation graph
    //std::cerr << "source" << std::endl;
    construct_derivation(graph_source, admissible_target, nodes_map_source, derivation_source);
    
    //std::cerr << "target" << std::endl;
    construct_derivation(graph_target, admissible_source, nodes_map_target, derivation_target);
    
    // extract composed sub-tree pairs from derivation_source and derivation_target
    extract_composed(source, target, alignment);
  }

  // we will extract pairs from derivatins...

  typedef std::vector<int, std::allocator<int> > index_set_type;

  struct Candidate
  {
    derivation_pair_type composed;
    index_set_type j;
    
    Candidate(const derivation_pair_type& __edge, const index_set_type& __j)
      : composed(__edge), j(__j) {}
  };
  typedef Candidate candidate_type;
  typedef utils::chunk_vector<candidate_type, 4096 / sizeof(candidate_type), std::allocator<candidate_type> > candidate_set_type;
  
  struct compare_heap_type
  {
    
    // we use greater, so that when popped from heap, we will grab "less" in back...
    bool operator()(const candidate_type* x, const candidate_type* y) const
    {
      return (x->composed.source.compose > y->composed.source.compose
	      || (! (y->composed.source.compose > x->composed.source.compose)
		  && (x->composed.source.internal > y->composed.source.internal
		      || ( !(y->composed.source.internal > x->composed.source.internal)
			   && x->composed.target.internal > y->composed.target.internal))));
    }
  };
  
  typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > candidate_heap_base_type;
  typedef utils::std_heap<const candidate_type*,  candidate_heap_base_type, compare_heap_type> candidate_heap_type;
  
  void extract_composed(const sentence_type& source,
			const sentence_type& target,
			const alignment_type& alignment)
  {
    derivations.clear();
    derivations.resize(spans.size());
    
    edge_set_type edges_source_new;
    edge_set_type edges_target_new;
    tail_set_type tails_source_new;
    tail_set_type tails_target_new;

    rule_pair_type rule_pair;
    
    // first, compute pairing...
    for (size_t itg_pos = 0; itg_pos != spans.size(); ++ itg_pos)
      if (admissible_source[itg_pos] && admissible_target[itg_pos]) {
	candidates.clear();
	cands.clear();
	
	const derivation_pair_type edge(derivation_source[itg_pos].front(),
					derivation_target[itg_pos].front());

	if (edge.source.tails.size() != edge.target.tails.size())
	  throw std::runtime_error("do not match with tails-size?");
	
	//
	// we need to compute aligns... HOW?
	// we need mapping from target-side into source-side
	//
	index_set_type aligns(edge.target.tails.size());
	for (size_t i = 0; i != aligns.size(); ++ i) {
	  const size_t itg_pos = nodes_map_target[edge.target.tails[i]];
	  
	  size_t pos = 0;
	  for (/**/; pos != edge.source.tails.size() && itg_pos != nodes_map_source[edge.source.tails[pos]]; ++ pos);
	  aligns[i] = pos;
	}
	
	index_set_type j(edge.source.tails.size(), -1);
	candidates.push_back(candidate_type(edge, j));
	cands.push(&candidates.back());
	
	while (! cands.empty()) {
	  const candidate_type* item = cands.top();
	  cands.pop();
	  
	  const derivation_pair_type& edge_composed = item->composed;
	  
	  construct_rule_pair(source, target, alignment, edge_composed, rule_pair);
	  
	  os << rule_pair << " ||| 1" << '\n';

	  // include into derivations[itg_pos]! 
	  if ((max_height <= 0 || (edge_composed.source.height <= max_height && edge_composed.target.height <= max_height))
	      && (max_nodes <= 0 || (edge_composed.source.internal < max_nodes && edge_composed.target.internal < max_nodes)))
	    derivations[itg_pos].push_back(edge_composed);
	  
	  // push-successor...
	  index_set_type j = item->j;
	  
	  for (size_t i = 0; i != j.size(); ++ i)
	    if (! derivations[nodes_map_source[edge.source.tails[i]]].empty()) {
	      ++ j[i];
	      
	      if (j[i] < static_cast<int>(derivations[nodes_map_source[edge.source.tails[i]]].size())) {
		edges_source_new.clear();
		edges_target_new.clear();
		tails_source_new.clear();
		tails_target_new.clear();

		int composed_size_source = edge_composed.source.compose;
		if (j[i] - 1 >= 0)
		  composed_size_source -= derivations[edge.source.tails[i]][j[i] - 1].source.compose;
		composed_size_source += derivations[edge.source.tails[i]][j[i]].source.compose;
		const int& composed_size_target = composed_size_source;

		if (max_compose <= 0 || composed_size_source <= max_compose) {
		  const std::pair<int, int> composed_stat = compose_tails(j, edge.source.tails, edge.target.tails, aligns,
									  edge.source.internal,
									  edge.target.internal,
									  tails_source_new,
									  tails_target_new);
		  
		  if (max_nodes <= 0 || (composed_stat.first <= max_nodes && composed_stat.second <= max_nodes)) {
		    
		    std::pair<int, int> source_stat;
		    std::pair<int, int> target_stat;
		    {
		      edge_set_type::const_iterator  eiter_begin = edge.source.edges.begin();
		      edge_set_type::const_iterator  eiter_end   = edge.source.edges.end();
		      int i = 0;
		      
		      source_stat = compose_edges(j, edge.source.tails, i, eiter_begin, eiter_end, edges_source_new);
		    }

		    {
		      edge_set_type::const_iterator  eiter_begin = edge.target.edges.begin();
		      edge_set_type::const_iterator  eiter_end   = edge.target.edges.end();
		      int i = 0;
		      
		      target_stat = compose_edges(j, edge.target.tails, aligns, i, eiter_begin, eiter_end, edges_target_new);
		    }
		    
		    
		    if (max_height <= 0 || (source_stat.first <= max_height && target_stat.first <= max_height)) {
		      candidates.push_back(candidate_type(derivation_pair_type(derivation_edge_type(edges_source_new,
												    tails_source_new,
												    source_stat.first,
												    source_stat.second,
												    composed_size_source),
									       derivation_edge_type(edges_target_new,
												    tails_target_new,
												    target_stat.first,
												    target_stat.second,
												    composed_size_target)), j));
		      
		      cands.push(&candidates.back());
		    }
		  }
		}
	      }

	      if (item->j[i] != -1) break;
	      
	      -- j[i];
	    }
	}
	
	// finished...
	std::sort(derivations[itg_pos].begin(), derivations[itg_pos].end(), less_derivation_pair_type());
      }
  }
  
  typedef std::vector<bool, std::allocator<bool> > covered_type;
  typedef std::vector<int, std::allocator<int> >   position_set_type;
  typedef std::vector<tree_rule_type, std::allocator<tree_rule_type> > tree_rule_set_type;
  
  struct CollapseFrontierIterator
  {
    CollapseFrontierIterator(tree_rule_set_type& __trees) : trees(__trees) {}
    
    CollapseFrontierIterator& operator=(const word_type& value)
    {
      trees.push_back(value);
      return *this;
    }
    
    CollapseFrontierIterator& operator*()  { return *this; }
    CollapseFrontierIterator& operator++() { return *this; }
    
    tree_rule_set_type& trees;
  };
  
  position_set_type positions_target;
  position_set_type positions_source;

  covered_type   covered;
  position_set_type positions_relative;
  index_set_type aligns;  
  tree_rule_set_type trees;

  void construct_rule_pair(const sentence_type& source,
			   const sentence_type& target,
			   const alignment_type& alignment,
			   const derivation_pair_type& derivation,
			   rule_pair_type& rule_pair)
  {
    rule_pair.clear();

    positions_source.clear();
    positions_target.clear();

    positions_source.resize(source.size(), -1);
    positions_target.resize(target.size(), -1);
    
    {
      positions_relative.clear();
      
      covered.clear();
      covered.resize(source.size(), false);
    
      int index = 0;
      int frontier_pos = 0;
      
      edge_set_type::const_iterator iter     = derivation.source.edges.begin();
      edge_set_type::const_iterator iter_end = derivation.source.edges.end();
      
      construct_rule(iter, iter_end, rule_pair.source, index, frontier_pos, positions_relative, covered);
      
      if (! positions_relative.empty()) {
	position_set_type::const_iterator piter = positions_relative.begin();
	for (size_t i = 0; i != covered.size(); ++ i)
	  if (covered[i]) {
	    positions_source[i] = *piter;
	    ++ piter;
	  }
      }
    }
    
    {
      // conctruct non-terminal alignment...
      aligns.resize(derivation.target.tails.size());
      for (size_t i = 0; i != aligns.size(); ++ i) {
	const size_t itg_pos = nodes_map_target[derivation.target.tails[i]];
	
	size_t pos = 0;
	for (/**/; pos != derivation.source.tails.size() && itg_pos != nodes_map_source[derivation.source.tails[pos]]; ++ pos);
	aligns[i] = pos;
      }

      positions_relative.clear();
      
      covered.clear();
      covered.resize(target.size(), false);
      
      int index = 0;
      int frontier_pos = 0;
      
      edge_set_type::const_iterator iter     = derivation.target.edges.begin();
      edge_set_type::const_iterator iter_end = derivation.target.edges.end();
      
      construct_rule(aligns, iter, iter_end, rule_pair.target, index, frontier_pos, positions_relative, covered);
      
      if (! positions_relative.empty()) {
	position_set_type::const_iterator piter = positions_relative.begin();
	for (size_t i = 0; i != covered.size(); ++ i)
	  if (covered[i]) {
	    positions_target[i] = *piter;
	    ++ piter;
	  }
      }
    }
    
    // construct alignment-vector...
    if (! positions_relative.empty())
      for (size_t src = 0; src != positions_source.size(); ++ src)
	if (positions_source[src] >= 0 && ! alignment[src].empty()) {
	  alignment_type::value_type::const_iterator aiter_begin = alignment[src].begin();
	  alignment_type::value_type::const_iterator aiter_end   = alignment[src].end();
	  
	  for (alignment_type::value_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)  {
	    if (positions_target[*aiter] < 0)
	      throw std::runtime_error("invalid alignment...?");
	    
	    rule_pair.alignment.push_back(std::make_pair(positions_source[src], positions_target[*aiter]));
	  }
	}

    if (frontier_source) {
      trees.clear();
      rule_pair.source.frontier(CollapseFrontierIterator(trees));
      rule_pair.source = tree_rule_type(rule_pair.source.label, trees.begin(), trees.end());
    }

    if (frontier_target) {
      trees.clear();
      rule_pair.target.frontier(CollapseFrontierIterator(trees));
      rule_pair.target = tree_rule_type(rule_pair.target.label, trees.begin(), trees.end());
    }
  }
  
  template <typename Iterator>
  void construct_rule(Iterator& iter,
		      Iterator last,
		      tree_rule_type& tree_rule,
		      int& index,
		      int& frontier_pos,
		      position_set_type& pos_map,
		      covered_type& covered)
  {
    // pre-order traversal...
    
    if (iter == last) {
      // this should not happen...???
      tree_rule.label = tree_rule.label.non_terminal();
      return;
    }
    
    const hypergraph_type::edge_type& edge = graph_source.edges[*iter];
    ++ iter;
    
    const size_t itg_pos = nodes_map_source[edge.head];
    const int edge_first = spans[itg_pos].source.first;
    const int edge_last  = spans[itg_pos].source.second;
    
    for (int pos = edge_first; pos != edge_last; ++ pos)
      covered[pos] = true;
    
    tree_rule = tree_rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end());
    
    size_t tail_pos = 0;
    tree_rule_type::iterator titer_end = tree_rule.end();
    for (tree_rule_type::iterator titer = tree_rule.begin(); titer != titer_end; ++ titer) {
      if (titer->label.is_non_terminal()) {
	const hypergraph_type::id_type node_id = edge.tails[tail_pos];
	++ tail_pos;
	
	if (iter != last && node_id == graph_source.edges[*iter].head)
	  construct_rule(iter, last, *titer, index, frontier_pos, pos_map, covered);
	else {
	  const size_t itg_pos = nodes_map_source[node_id];
	  const int edge_first = spans[itg_pos].source.first;
	  const int edge_last  = spans[itg_pos].source.second;
	  
	  for (int pos = edge_first; pos != edge_last; ++ pos)
	    covered[pos] = false;
	  
	  titer->label = titer->label.non_terminal(index + 1);
	  ++ index;
	  ++ frontier_pos;
	}
      } else {
	pos_map.push_back(frontier_pos);
	++ frontier_pos;
      }
    }
  }

  template <typename Iterator>
  void construct_rule(const index_set_type& aligns,
		      Iterator& iter,
		      Iterator last,
		      tree_rule_type& tree_rule,
		      int& index,
		      int& frontier_pos,
		      position_set_type& pos_map,
		      covered_type& covered)
  {
    // pre-order traversal...
    
    if (iter == last) {
      // this should not happen...???
      tree_rule.label = tree_rule.label.non_terminal();
      return;
    }
    
    const hypergraph_type::edge_type& edge = graph_target.edges[*iter];
    ++ iter;
    
    const size_t itg_pos = nodes_map_target[edge.head];
    const int edge_first = spans[itg_pos].target.first;
    const int edge_last  = spans[itg_pos].target.second;
    
    for (int pos = edge_first; pos != edge_last; ++ pos)
      covered[pos] = true;
    
    tree_rule = tree_rule_type(edge.rule->lhs, edge.rule->rhs.begin(), edge.rule->rhs.end());
    
    size_t tail_pos = 0;
    tree_rule_type::iterator titer_end = tree_rule.end();
    for (tree_rule_type::iterator titer = tree_rule.begin(); titer != titer_end; ++ titer) {
      if (titer->label.is_non_terminal()) {
	const hypergraph_type::id_type node_id = edge.tails[tail_pos];
	++ tail_pos;
	
	if (iter != last && node_id == graph_target.edges[*iter].head)
	  construct_rule(aligns, iter, last, *titer, index, frontier_pos, pos_map, covered);
	else {
	  const size_t itg_pos = nodes_map_target[node_id];
	  const int edge_first = spans[itg_pos].target.first;
	  const int edge_last  = spans[itg_pos].target.second;
	  
	  for (int pos = edge_first; pos != edge_last; ++ pos)
	    covered[pos] = false;
	  
	  titer->label = titer->label.non_terminal(aligns[index] + 1);
	  ++ index;
	  ++ frontier_pos;
	}
      } else {
	pos_map.push_back(frontier_pos);
	++ frontier_pos;
      }
    }
  }

  
  template <typename Iterator>
  std::pair<int, int> compose_edges(const index_set_type& j,
				    const tail_set_type& tails,
				    const index_set_type& aligns,
				    int& i,
				    Iterator& iter,
				    Iterator last,
				    edge_set_type& edges_new)
  {
    // this should not happen...
    if (iter == last) return std::make_pair(0, 0);
    
    edges_new.push_back(*iter);
    const hypergraph_type::edge_type& edge = graph_target.edges[*iter];
    ++ iter;
    
    int height = 1;
    int num_tails = edge.tails.size();
    
    hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
      if (iter != last && graph_target.edges[*iter].head == *titer) {
	const std::pair<int, int> result = compose_edges(j, tails, aligns, i, iter, last, edges_new);
	
	height = utils::bithack::max(height, result.first + 1);
	num_tails += result.second;
      } else if (i != tails.size()) {
	if (j[aligns[i]] >= 0) {
	  const derivation_pair_type& edge = derivations[nodes_map_target[tails[i]]][j[aligns[i]]];
	  edges_new.insert(edges_new.end(), edge.target.edges.begin(), edge.target.edges.end());
	  
	  height = utils::bithack::max(height, edge.target.height + 1);
	  num_tails += edge.target.internal;
	}
	++ i;
      } else 
	throw std::runtime_error("# of tails and # of frontier do not match");
    }
    
  }

  template <typename Iterator>
  std::pair<int, int> compose_edges(const index_set_type& j,
				    const tail_set_type& tails,
				    int& i,
				    Iterator& iter,
				    Iterator last,
				    edge_set_type& edges_new)
  {
    // this should not happen...
    if (iter == last) return std::make_pair(0, 0);
    
    edges_new.push_back(*iter);
    const hypergraph_type::edge_type& edge = graph_source.edges[*iter];
    ++ iter;
    
    int height = 1;
    int num_tails = edge.tails.size();
    
    hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
    for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
      if (iter != last && graph_source.edges[*iter].head == *titer) {
	const std::pair<int, int> result = compose_edges(j, tails, i, iter, last, edges_new);
	
	height = utils::bithack::max(height, result.first + 1);
	num_tails += result.second;
      } else if (i != tails.size()) {
	if (j[i] >= 0) {
	  const derivation_pair_type& edge = derivations[nodes_map_source[tails[i]]][j[i]];
	  edges_new.insert(edges_new.end(), edge.source.edges.begin(), edge.source.edges.end());
	  
	  height = utils::bithack::max(height, edge.source.height + 1);
	  num_tails += edge.source.internal;
	}
	++ i;
      } else
	throw std::runtime_error("# of tails and # of frontier do not match");
    }
    
    return std::make_pair(height, num_tails);
  }

  std::pair<int, int> compose_tails(const index_set_type& j,
				    const tail_set_type& tails_source,
				    const tail_set_type& tails_target,
				    const index_set_type& aligns,
				    int internal_source,
				    int internal_target,
				    tail_set_type& tails_source_new,
				    tail_set_type& tails_target_new)
  {
    for (size_t i = 0; i != j.size(); ++ i) {
      if (j[i] < 0)
	tails_source_new.push_back(tails_source[i]);
      else {
	const derivation_pair_type& edge = derivations[nodes_map_source[tails_source[i]]][j[i]];
	
	tails_source_new.insert(tails_source_new.end(), edge.source.tails.begin(), edge.source.tails.end());
	
	internal_source += edge.source.internal;

	if (max_nodes > 0 && internal_source > max_nodes)
	  return std::make_pair(internal_source, internal_target);
      }
      
      if (j[aligns[i]] < 0)
	tails_target_new.push_back(tails_target[i]);
      else {
	const derivation_pair_type& edge = derivations[nodes_map_target[tails_target[i]]][j[aligns[i]]];
	
	tails_target_new.insert(tails_target_new.end(), edge.target.tails.begin(), edge.target.tails.end());
	
	internal_target += edge.target.internal;
	
	if (max_nodes > 0 && internal_target > max_nodes)
	  return std::make_pair(internal_source, internal_target);
      }
    }
    
    return std::make_pair(internal_source, internal_target);
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

  void construct_derivation(const hypergraph_type& graph,
			    const admissible_set_type& admissible,
			    const node_map_type& nodes_map,
			    derivation_graph_type& derivation)
  {
    for (size_t edge_id = 0; edge_id != graph.edges.size(); ++ edge_id) {
      const hypergraph_type::edge_type& edge = graph.edges[edge_id];
      
      const size_t itg_pos = nodes_map[edge.head];
      
      if (! admissible[itg_pos]) continue;
      
      edge_set_type edges(1, edge_id);
      tail_set_type tails;
      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	if (! admissible[nodes_map[*titer]])
	  tails.push_back(*titer);

      //std::cerr << "edges: " << *edge.rule;
      
      tail_set_type tails_next;
      while (! tails.empty()) {
	const hypergraph_type::node_type& node = graph.nodes[tails.front()];
	
	if (node.edges.size() != 1)
	  throw std::runtime_error("invalid node with multiple edges?");
	
	edges.push_back(node.edges.front());
	
	tails_next.clear();
	
	const hypergraph_type::edge_type& edge = graph.edges[node.edges.front()];
	hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer)
	  if (! admissible[nodes_map[*titer]])
	    tails_next.push_back(*titer);
	
	//std::cerr << " + " << *edge.rule;
	
	tails_next.insert(tails_next.end(), tails.begin() + 1, tails.end());
	
	tails.swap(tails_next);
      }
      
      // we have a set of edges!
      // we will re-construct sub-tree's admissible tails
      
      tails.clear();
      edge_set_type::const_iterator eiter = edges.begin();
      edge_set_type::const_iterator eiter_end = edges.end();
      const std::pair<int, int> rule_stat = construct_tails(graph, eiter, eiter_end, tails);
      
      //std::cerr << " tails: ";
      //std::copy(tails.begin(), tails.end(), std::ostream_iterator<int>(std::cerr, " "));
      //std::cerr << std::endl;
      
      derivation[itg_pos].push_back(derivation_edge_type(edges, tails, rule_stat.first, rule_stat.second));
    }
  }
  
  struct attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
  {
    attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
    attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
    attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
  };

  void compute_node_map(const hypergraph_type& forest, admissible_set_type& admissible, node_map_type& nodes_map)
  {
    nodes_map.resize(forest.nodes.size());
    
    hypergraph_type::edge_set_type::const_iterator eiter_end = forest.edges.end();
    for (hypergraph_type::edge_set_type::const_iterator eiter = forest.edges.begin(); eiter != eiter_end; ++ eiter) {
      const hypergraph_type::edge_type& edge = *eiter;
      
      attribute_set_type::const_iterator aiter = edge.attributes.find(attr_node_id);
      if (aiter == edge.attributes.end())
	throw std::runtime_error("no attributes?");
      
      const int pos = boost::apply_visitor(attribute_integer(), aiter->second);
      if (pos < 0)
	throw std::runtime_error("invalid attribute?");
      
      if (admissible[pos])
	throw std::runtime_error("multiple edges?");

      // we store node-id + 1 so that we can easily investigate asmissibility
      admissible[pos] = edge.head + 1;
      nodes_map[edge.head] = pos;
    }
  }

  template <typename Blocker>
  id_pair_type operator()(const itg_type& itg,
			  const sentence_type& source,
			  const sentence_type& target,
			  span_pair_set_type& spans,
			  hypergraph_type& graph_source,
			  hypergraph_type& graph_target,
			  alignment_type& alignment,
			  Blocker blocker)
  {
    if (blocker(itg) || itg.antecedent.empty())  {
      rule_type::symbol_set_type rhs_source(source.begin() + itg.spans.source.first, source.begin() + itg.spans.source.second);
      rule_type::symbol_set_type rhs_target(target.begin() + itg.spans.target.first, target.begin() + itg.spans.target.second);
      
      hypergraph_type::edge_type& edge_source = graph_source.add_edge();
      hypergraph_type::edge_type& edge_target = graph_target.add_edge();
      
      edge_source.rule = rule_type::create(rhs_source.empty()
					   ? rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1)
					   : rule_type(vocab_type::X, rhs_source));
      edge_target.rule = rule_type::create(rhs_target.empty()
					   ? rule_type(vocab_type::X, &vocab_type::EPSILON, (&vocab_type::EPSILON) + 1)
					   : rule_type(vocab_type::X, rhs_target));

      edge_source.attributes[attr_node_id] = attribute_set_type::int_type(spans.size());
      edge_target.attributes[attr_node_id] = attribute_set_type::int_type(spans.size());
      
      const hypergraph_type::id_type node_id_source = graph_source.add_node().id;
      const hypergraph_type::id_type node_id_target = graph_target.add_node().id;
      
      graph_source.connect_edge(edge_source.id, node_id_source);
      graph_target.connect_edge(edge_target.id, node_id_target);
      
      // compute alignment matrix
      if (itg.spans.target.first != itg.spans.target.second)
	for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src)
	  for (int trg = itg.spans.target.first; trg != itg.spans.target.second; ++ trg)
	    alignment[src].insert(trg);

      // push-back spans...
      spans.push_back(itg.spans);
      
      return std::make_pair(node_id_source, node_id_target);
    } else {
      const id_pair_type id_pair_front = operator()(itg.antecedent.front(), source, target, spans, graph_source, graph_target, alignment, blocker);
      const id_pair_type id_pair_back  = operator()(itg.antecedent.back(),  source, target, spans, graph_source, graph_target, alignment, blocker);
      
      hypergraph_type::id_type tails_source[2];
      hypergraph_type::id_type tails_target[2];

      tails_source[0] = id_pair_front.first;
      tails_source[1] = id_pair_back.first;
      
      tails_target[  itg.inverse] = id_pair_front.second;
      tails_target[! itg.inverse] = id_pair_back.second;
      
      hypergraph_type::edge_type& edge_source = graph_source.add_edge(tails_source, tails_source + 2);
      hypergraph_type::edge_type& edge_target = graph_target.add_edge(tails_target, tails_target + 2);
      
      edge_source.rule = rule_binary;
      edge_target.rule = rule_binary;
      
      edge_source.attributes[attr_node_id] = attribute_set_type::int_type(spans.size());
      edge_target.attributes[attr_node_id] = attribute_set_type::int_type(spans.size());
      
      const hypergraph_type::id_type node_id_source = graph_source.add_node().id;
      const hypergraph_type::id_type node_id_target = graph_target.add_node().id;
      
      graph_source.connect_edge(edge_source.id, node_id_source);
      graph_target.connect_edge(edge_target.id, node_id_target);
      
      // push-back spans...
      spans.push_back(itg.spans);
      
      return std::make_pair(node_id_source, node_id_target);
    }
  }
  
  std::ostream& os;
  
  const int max_nodes;
  const int max_height;
  const int max_compose;
  const int max_scope;
  
  const bool frontier_source;
  const bool frontier_target;
  const bool remove_unary;
  
  span_pair_set_type spans;
  admissible_set_type admissible_source;
  admissible_set_type admissible_target;
  node_map_type nodes_map_source;
  node_map_type nodes_map_target;
  derivation_graph_type derivation_source;
  derivation_graph_type derivation_target;
  
  derivation_pair_graph_type derivations;

  candidate_set_type  candidates;
  candidate_heap_type cands;

  hypergraph_type graph_source;
  hypergraph_type graph_target;

  rule_ptr_type rule_binary;

  const attribute_type attr_node_id;
};

struct HieroGrammar : public Grammar
{

  typedef sentence_type phrase_type;
  
  struct rule_pair_type
  {
    phrase_type source;
    phrase_type target;
    point_set_type alignment;
    
    rule_pair_type() : source(), target(), alignment() {}
    
    void clear()
    {
      source.clear();
      target.clear();
      alignment.clear();
    }

    friend
    std::ostream& operator<<(std::ostream& os, const rule_pair_type& x)
    {
      std::copy(x.source.begin(), x.source.end(), std::ostream_iterator<word_type>(os, " "));
      os << "||| ";
      std::copy(x.target.begin(), x.target.end(), std::ostream_iterator<word_type>(os, " "));
      os << "|||";
      point_set_type::const_iterator piter_end = x.alignment.end();
      for (point_set_type::const_iterator piter = x.alignment.begin(); piter != piter_end; ++ piter)
	os << ' ' << piter->first << '-' << piter->second;
      return os;
    }
  };
  
  typedef std::vector<rule_pair_type, std::allocator<rule_pair_type> > rule_pair_set_type;
  
  HieroGrammar(std::ostream& __os,
	       const int __max_span,
	       const int __max_length)
    : os(__os),
      max_span(__max_span),
      max_length(__max_length) {}
  
  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    span_pair_set_type spans;
    
    alignment.clear();
    alignment.resize(source.size());
    
    operator()(itg, source, target, spans, alignment, blocker);
  }

  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  span_pair_set_type& spans,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    // post-traversal order to collect paired-span
  
    if (itg.spans.source.first == itg.spans.source.second
	|| itg.spans.target.first == itg.spans.target.second)
      return;
    
    const int length_source = itg.spans.source.second - itg.spans.source.first;
    const int length_target = itg.spans.target.second - itg.spans.target.first;
    
    if (blocker(itg) || itg.antecedent.empty()) {
      spans.push_back(itg.spans);
      
      if (itg.spans.target.first != itg.spans.target.second)
	for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src)
	  for (int trg = itg.spans.target.first; trg != itg.spans.target.second; ++ trg)
	    alignment[src].insert(trg);

    } else if (max_span <= 0 || length_source <= max_span) {
      span_pair_set_type spans1;
      span_pair_set_type spans2;
      
      operator()(itg.antecedent.front(), source, target, spans1, alignment, blocker);
      operator()(itg.antecedent.back(),  source, target, spans2, alignment, blocker);
      
      // we will create hiero rule from spans1 and spans2!

      rule_pair_set_type rule_pairs;
      
      // first, single non-terminal
      if (! spans1.empty())
	hiero_rule(source, target, itg.spans, spans1, alignment, rule_pairs);
      
      if (! spans2.empty())
	hiero_rule(source, target, itg.spans, spans2, alignment, rule_pairs);
      
      // second, two non-terminals
      if (! spans1.empty() && ! spans2.empty())
	hiero_rule(source, target, itg.spans, spans1, spans2, alignment, rule_pairs);
    
      spans.swap(spans1);
      spans.insert(spans.end(), spans2.begin(), spans2.end());
      
      spans.push_back(itg.spans);
      
      if (! rule_pairs.empty()) {
	const double weight = 1.0 / rule_pairs.size();
	
	rule_pair_set_type::const_iterator riter_end = rule_pairs.end();
	for (rule_pair_set_type::const_iterator riter = rule_pairs.begin(); riter != riter_end; ++ riter)
	  os << *riter << " ||| " << weight << '\n';
      }
    }
    
    if (max_length <= 0 || (length_source <= max_length && length_target <= max_length)) {
      rule_pair_type rule_pair;
      rule_pair.source.push_back("[x]");
      rule_pair.target.push_back("[x]");
      rule_pair.source.insert(rule_pair.source.end(), source.begin() + itg.spans.source.first, source.begin() + itg.spans.source.second);
      rule_pair.target.insert(rule_pair.target.end(), target.begin() + itg.spans.target.first, target.begin() + itg.spans.target.second);
      
      for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src) {
	index_set_type::const_iterator aiter_begin = alignment[src].begin();
	index_set_type::const_iterator aiter_end   = alignment[src].end();
	
	for (index_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  rule_pair.alignment.push_back(std::make_pair(src - itg.spans.source.first, *aiter - itg.spans.target.first));
      }
      
      os << rule_pair << " ||| " << 1.0 << '\n';
    }
  }
  
  void hiero_rule(const sentence_type& source,
		  const sentence_type& target,
		  const span_pair_type& span,
		  const span_pair_set_type& spans,
		  const alignment_type& alignment,
		  rule_pair_set_type& rule_pairs)
  {
    const int length_source = span.source.second - span.source.first;
    const int length_target = span.target.second - span.target.first;
    
    // we need at least one terminal...
    span_pair_set_type::const_iterator siter_end = spans.end();
    for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)
      if (siter->source != span.source) {
	
	const int length_source1 = siter->source.second - siter->source.first;
	const int length_target1 = siter->target.second - siter->target.first;
	
	if (max_length > 0 && length_source - length_source1 > max_length) continue;
	if (max_length > 0 && length_target - length_target1 > max_length) continue;

	rule_pairs.resize(rule_pairs.size() + 1);
	rule_pair_type& rule_pair = rule_pairs.back();
	
	rule_pair.source.clear();
      
	rule_pair.source.push_back("[x]");
	rule_pair.source.insert(rule_pair.source.end(), source.begin() + span.source.first, source.begin() + siter->source.first);
	rule_pair.source.push_back("[x,1]");
	rule_pair.source.insert(rule_pair.source.end(), source.begin() + siter->source.second, source.begin() + span.source.second);
      
	rule_pair.target.clear();
	rule_pair.target.push_back("[x]");
	rule_pair.target.insert(rule_pair.target.end(), target.begin() + span.target.first, target.begin() + siter->target.first);
	rule_pair.target.push_back("[x,1]");
	rule_pair.target.insert(rule_pair.target.end(), target.begin() + siter->target.second, target.begin() + span.target.second);
	
	for (int src = span.source.first; src != span.source.second; ++ src) 
	  if (is_out_of_span(siter->source, src)) {
	    index_set_type::const_iterator aiter_begin = alignment[src].begin();
	    index_set_type::const_iterator aiter_end   = alignment[src].end();
	    
	    for (index_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	      if (is_out_of_span(siter->target, *aiter)) {
		const int mask_source = - (src    >= siter->source.second);
		const int mask_target = - (*aiter >= siter->target.second);
		
		// -1 for non-terminal label, [x,1]
		const int shift_source = (mask_source & (length_source1 - 1)) + span.source.first;
		const int shift_target = (mask_target & (length_target1 - 1)) + span.target.first;
		
		rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	      }
	  }
      }
  }

  void hiero_rule(const sentence_type& source,
		  const sentence_type& target,
		  const span_pair_type& span,
		  const span_pair_set_type& spans1,
		  const span_pair_set_type& spans2,
		  const alignment_type& alignment,
		  rule_pair_set_type& rule_pairs)
  {
    const int length_source = span.source.second - span.source.first;
    const int length_target = span.target.second - span.target.first;
    
    span_pair_set_type::const_iterator siter1_begin = spans1.begin();
    span_pair_set_type::const_iterator siter1_end = spans1.end();
    span_pair_set_type::const_iterator siter2_begin = spans2.begin();
    span_pair_set_type::const_iterator siter2_end = spans2.end();
    for (span_pair_set_type::const_iterator siter1 = siter1_begin; siter1 != siter1_end; ++ siter1) {
      const span_pair_type& span1 = *siter1;

      const int length_source1 = span1.source.second - span1.source.first;
      const int length_target1 = span1.target.second - span1.target.first;
      
      for (span_pair_set_type::const_iterator siter2 = siter2_begin; siter2 != siter2_end; ++ siter2) 
	if (span1.source.second != siter2->source.first) {
	  const span_pair_type& span2 = *siter2;
	  
	  const int length_source2 = span2.source.second - span2.source.first;
	  const int length_target2 = span2.target.second - span2.target.first;
	  
	  if (max_length > 0 && length_source - length_source1 - length_source2 > max_length) continue;
	  if (max_length > 0 && length_target - length_target1 - length_target2 > max_length) continue;

	  rule_pairs.resize(rule_pairs.size() + 1);
	  rule_pair_type& rule_pair = rule_pairs.back();

	  const bool inverse = ! (span1.target.second <= span2.target.first);
	
	  const span_type& span_source1 = span1.source;
	  const span_type& span_source2 = span2.source;
	
	  const span_type& span_target1 = (inverse ? span2.target : span1.target);
	  const span_type& span_target2 = (inverse ? span1.target : span2.target);
	
	  rule_pair.source.clear();
	  rule_pair.source.push_back("[x]");
	  rule_pair.source.insert(rule_pair.source.end(), source.begin() + span.source.first, source.begin() + span_source1.first);
	  rule_pair.source.push_back("[x,1]");
	  rule_pair.source.insert(rule_pair.source.end(), source.begin() + span_source1.second, source.begin() + span_source2.first);
	  rule_pair.source.push_back("[x,2]");
	  rule_pair.source.insert(rule_pair.source.end(), source.begin() + span_source2.second, source.begin() + span.source.second);
	
	  rule_pair.target.clear();
	  rule_pair.target.push_back("[x]");
	  rule_pair.target.insert(rule_pair.target.end(), target.begin() + span.target.first, target.begin() + span_target1.first);
	  rule_pair.target.push_back(inverse ? "[x,2]" : "[x,1]");
	  rule_pair.target.insert(rule_pair.target.end(), target.begin() + span_target1.second, target.begin() + span_target2.first);
	  rule_pair.target.push_back(inverse ? "[x,1]" : "[x,2]");
	  rule_pair.target.insert(rule_pair.target.end(), target.begin() + span_target2.second, target.begin() + span.target.second);
	  
	  for (int src = span.source.first; src != span.source.second; ++ src) 
	    if (is_out_of_span(span1.source, src) && is_out_of_span(span2.source, src)) {
	      index_set_type::const_iterator aiter_begin = alignment[src].begin();
	      index_set_type::const_iterator aiter_end   = alignment[src].end();
	      
	      for (index_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
		if (is_out_of_span(span1.target, *aiter) && is_out_of_span(span2.target, *aiter)) {
		  const int mask_source1 = - (src >= span1.source.second);
		  const int mask_source2 = - (src >= span2.source.second);
		  
		  const int mask_target1 = - (*aiter >= span1.target.second);
		  const int mask_target2 = - (*aiter >= span2.target.second);
		  
		  // -1 for non-terminal label, [x,1] or [x,2]
		  const int shift_source = (span.source.first
					    + (mask_source1 & (length_source1 - 1))
					    + (mask_source2 & (length_source2 - 1)));
		  const int shift_target = (span.target.first
					    + (mask_target1 & (length_target1 - 1))
					    + (mask_target2 & (length_target2 - 1)));
		  
		  rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
		}
	    }
	}
    }
  }

  bool is_out_of_span(const span_type& span, const int pos) const
  {
    return pos < span.first || span.second <= pos;
  }

  std::ostream& os;
  const int max_span;
  const int max_length;
};


struct BlockerModel
{
  bool operator()(const itg_type& itg) const
  {
    return itg.block;
  }
};

struct BlockerBlock
{
  bool operator()(const itg_type& itg) const
  {
    return ((itg.spans.source.second - itg.spans.source.first) == 1) || ((itg.spans.target.second - itg.spans.target.first) == 1);
  }
};

struct BlockerTerminal
{
  bool operator()(const itg_type& itg) const
  {
    return itg.antecedent.empty();
  }
};

int main(int argc, char** argv)
{
  try {
    namespace po = boost::program_options;
    
    path_type input_file = "-";
    path_type output_file = "-";
    path_type output_source_file;
    path_type output_target_file;
    path_type output_alignment_file;
    int max_length = 7;
    int max_span = 15;
    
    int max_nodes   = 15;
    int max_height  = 4;
    int max_compose = 0;
    int max_scope   = 0;

    bool frontier_source_mode = false;
    bool frontier_target_mode = false;
    
    bool remove_epsilon = false;
    bool remove_unary = false;
    
    bool scfg_mode = false;
    bool ghkm_mode = false;
    bool tree_source_mode = false;
    bool tree_target_mode = false;
    
    bool phrase_mode = false;
    bool block_mode = false;
    bool exhaustive_mode = false;

    int debug = 0;
  
    po::options_description desc("options");
    desc.add_options()
      ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
      ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")
      
      ("source",    po::value<path_type>(&output_source_file),    "output source yield")
      ("target",    po::value<path_type>(&output_target_file),    "output target yield")
      ("alignment", po::value<path_type>(&output_alignment_file), "output word-for-word alignment")
      
      ("max-length", po::value<int>(&max_length)->default_value(max_length), "max terminal length")
      ("max-span",   po::value<int>(&max_span)->default_value(max_span),     "max span")
      
      ("max-nodes",   po::value<int>(&max_nodes)->default_value(max_nodes),     "max nodes")
      ("max-height",  po::value<int>(&max_height)->default_value(max_height),   "max height")
      ("max-compose", po::value<int>(&max_compose)->default_value(max_compose), "max compose")
      ("max-scope",   po::value<int>(&max_scope)->default_value(max_scope),     "max scope")
      ("frontier-source", po::bool_switch(&frontier_source_mode),               "take frontier of source side (string-to-* model)")
      ("frontier-target", po::bool_switch(&frontier_target_mode),               "take frontier of target side (*-to-string model)")
      
      ("remove-epsilon", po::bool_switch(&remove_epsilon), "remove <epsilon> from trees")
      ("remove-unary",   po::bool_switch(&remove_unary),   "remove unary rules from trees")

      ("scfg", po::bool_switch(&scfg_mode), "extract SCFG rules")
      ("ghkm", po::bool_switch(&ghkm_mode), "extract GHKM rules")
      ("tree-source", po::bool_switch(&tree_source_mode), "extract source tree")
      ("tree-target", po::bool_switch(&tree_target_mode), "extract target tree")
    
      ("phrase",     po::bool_switch(&phrase_mode),     "phrase-wise model alignment (many-to-many)")
      ("block",      po::bool_switch(&block_mode),      "block-wise alignment (one-to-many)")
      ("exhaustive", po::bool_switch(&exhaustive_mode), "exhaustive alignment (one-to-one)")
      
      ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
      ("help", "help message");
  
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
    po::notify(vm);
  
    if (vm.count("help")) {
      std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
      return 0;
    }
      
    if (int(phrase_mode) + block_mode + exhaustive_mode > 1)
      throw std::runtime_error("either phrase|block|exhaustive");

    if (int(phrase_mode) + block_mode + exhaustive_mode == 0)
      phrase_mode = true;
    
    if (int(scfg_mode) + ghkm_mode + tree_source_mode + tree_target_mode > 1)
      throw std::runtime_error("either scfg|ghkm|tree-source|tree-target");
    
    if (int(scfg_mode) + ghkm_mode + tree_source_mode + tree_target_mode == 0)
      scfg_mode = true;
    
    typedef boost::spirit::istream_iterator iter_type;
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    std::auto_ptr<std::ostream> os_src(! output_source_file.empty() ? new utils::compress_ostream(output_source_file, 1024 * 1024) : 0);
    std::auto_ptr<std::ostream> os_trg(! output_target_file.empty() ? new utils::compress_ostream(output_target_file, 1024 * 1024) : 0);
    std::auto_ptr<std::ostream> os_align(! output_alignment_file.empty() ? new utils::compress_ostream(output_alignment_file, 1024 * 1024) : 0);

    is.unsetf(std::ios::skipws);
    os.precision(20);
    
    derivation_parser<iter_type> parser;
    
    iter_type iter(is);
    iter_type end;
  
    itg_type itg;
    sentence_type source;
    sentence_type target;

    Grammar::alignment_type alignment;

    HieroGrammar scfg_grammar(os, max_span, max_length);
    GHKMGrammar  ghkm_grammar(os, max_nodes, max_height, max_compose, max_scope, frontier_source_mode, frontier_target_mode, remove_unary);
    TreeSource   tree_source(os, remove_epsilon, remove_unary);
    TreeTarget   tree_target(os, remove_epsilon, remove_unary);
    
    while (iter != end) {
      itg.clear();
    
      if (! boost::spirit::qi::phrase_parse(iter, end, parser, boost::spirit::standard::space, itg))
	throw std::runtime_error("parsing failed");
      
      source.clear();
      target.clear();

      span_derivation_source(itg, source);
      span_derivation_target(itg, target);
      
      if (os_src.get())
	*os_src << source << '\n';
      if (os_trg.get())
	*os_trg << target << '\n';
      
      //print_tree(std::cout, itg) << std::endl;
      
      if (scfg_mode) {
	if (phrase_mode)
	  scfg_grammar(itg, source, target, alignment, BlockerModel());
	else if (block_mode)
	  scfg_grammar(itg, source, target, alignment, BlockerBlock());
	else
	  scfg_grammar(itg, source, target, alignment, BlockerTerminal());
      } else if (tree_source_mode) {
	if (phrase_mode)
	  tree_source(itg, source, target, alignment, BlockerModel());
	else if (block_mode)
	  tree_source(itg, source, target, alignment, BlockerBlock());
	else
	  tree_source(itg, source, target, alignment, BlockerTerminal());
      } else if (tree_target_mode) {
	if (phrase_mode)
	  tree_target(itg, source, target, alignment, BlockerModel());
	else if (block_mode)
	  tree_target(itg, source, target, alignment, BlockerBlock());
	else
	  tree_target(itg, source, target, alignment, BlockerTerminal());	
      } else {
	if (phrase_mode)
	  ghkm_grammar(itg, source, target, alignment, BlockerModel());
	else if (block_mode)
	  ghkm_grammar(itg, source, target, alignment, BlockerBlock());
	else
	  ghkm_grammar(itg, source, target, alignment, BlockerTerminal());
      }
      
      // dump alignment...
      if (os_align.get()) {
	bool initial = true;
	for (size_t src = 0; src != alignment.size(); ++ src) 
	  if (! alignment[src].empty()) {
	    const Grammar::alignment_type::value_type& align = alignment[src];
	    
	    Grammar::alignment_type::value_type::const_iterator aiter_end = align.end();
	    for (Grammar::alignment_type::value_type::const_iterator aiter = align.begin(); aiter != aiter_end; ++ aiter) {
	      if (! initial)
		*os_align << ' ';
	      *os_align << src << '-' << *aiter;
	      
	      initial = false;
	    }
	  }
	*os_align << '\n';
      }
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
