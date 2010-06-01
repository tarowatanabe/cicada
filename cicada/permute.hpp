// -*- mode: c++ -*-

#ifndef __CICADA__PERMUTE__HPP__
#define __CICADA__PERMUTE__HPP__ 1

#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>

// we will add more permuted edges to the existing hypergraph
// with features, rule=1 score...

namespace cicada
{
  struct Permute
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef feature_set_type::feature_type    feature_type;
    
    typedef std::vector<int, std::allocator<int> > permutation_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > phrase_type;
    
            
    void operator()(const hypergraph_type& source, hypergraph_type& target, const int permute_size)
    {
      target.clear();
      
      permutation_type permutation;

      tails_type tails;

      phrase_type source_phrase;
      phrase_type source_non_terminals;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type& node = target.add_node();
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	
	  hypergraph_type::edge_type& edge = target.add_edge(edge_source.tails.begin(), edge_source.tails.end());
	  
	  edge.rule = edge_source.rule;
	  
	  edge.features = edge_source.features;
	  edge.features[rule_feature(*edge.rule)] = 1.0;
	  edge.features["rule:original"] = 1.0;
	  
	  target.connect_edge(edge.id, node.id);
	  
	  if (edge_source.tails.size() > 1) {
	    
	    permutation.clear();
	    for (int i = 0; i < edge_source.tails.size(); ++ i)
	      permutation.push_back(i);

	    source_non_terminals.clear();
	    rule_type::symbol_set_type::const_iterator siter_end = edge_source.rule->source.end();
	    for (rule_type::symbol_set_type::const_iterator siter = edge_source.rule->source.begin(); siter != siter_end; ++ siter)
	      if (siter->is_non_terminal())
		source_non_terminals.push_back(*siter);
	    
	    source_phrase.clear();
	    source_phrase.insert(source_phrase.end(), edge_source.rule->source.begin(), edge_source.rule->source.end());
	    
	    tails.clear();
	    tails.insert(tails.end(), edge_source.tails.begin(), edge_source.tails.end());
	    
	    // perform permutation...
	    // minus for all permutation
	    // zero for no-permutation
	    if (permute_size < 0 || permute_size >= 1)
	      while (std::next_permutation(permutation.begin(), permutation.end())) {
		
		if (! is_valid_permutation(permutation, permute_size)) continue;
		
		// permute nodes...
		for (int i = 0; i < tails.size(); ++ i)
		  tails[i] = edge_source.tails[permutation[i]];
		
		// permute source-phrase
		int non_terminal_pos = 0;
		phrase_type::iterator siter_end = source_phrase.end();
		for (phrase_type::iterator siter = source_phrase.begin(); siter != siter_end; ++ siter)
		  if (siter->is_non_terminal()) {
		    *siter = source_non_terminals[permutation[non_terminal_pos]];
		    
		    ++ non_terminal_pos;
		  }
		
		hypergraph_type::edge_type& edge = target.add_edge(tails.begin(), tails.end());
		
		edge.rule.reset(new rule_type(*edge_source.rule));
		edge.rule->source = rule_type::symbol_set_type(source_phrase.begin(), source_phrase.end());
		
		edge.features = edge_source.features;
		edge.features[rule_feature(*edge.rule)] = 1.0;
		
		target.connect_edge(edge.id, node.id);
	      }
	  }
	}
      }
      
      target.goal = source.goal;
    }
    

  private:
    bool is_valid_permutation(const permutation_type& permutation, const int permute_size)
    {
      if (permute_size < 0) return true;
      
      for (int index = 0; index < permutation.size(); ++ index)
	if (utils::bithack::abs(index - permutation[index]) > permute_size)
	  return false;
      return true;
    }

    feature_type rule_feature(const rule_type& rule)
    {
      std::string rule_string("rule:");
      rule_string += static_cast<const std::string&>(rule.lhs) + "->";
      
      rule_type::symbol_set_type::const_iterator siter_end = rule.source.end();
      for (rule_type::symbol_set_type::const_iterator siter = rule.source.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal())
	  rule_string += static_cast<const std::string&>(*siter);
	else
	  rule_string += '<' + static_cast<const std::string&>(*siter) + '>';
      
      return rule_string;
    }

  };
  
  inline
  void permute(const HyperGraph& source, HyperGraph& target, const int permute_size)
  {
    Permute()(source, target, permute_size);
  }
  
  inline
  void permute(HyperGraph& graph, const int permute_size)
  {
    HyperGraph target;
    
    permute(graph, target, permute_size);
    
    graph.swap(target);
  }
  
};

#endif
