// -*- mode: c++ -*-

#ifndef __CICADA__PERMUTE__HPP__
#define __CICADA__PERMUTE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>
#include <utils/bit_vector.hpp>
#include <utils/sgi_hash_set.hpp>

#include <boost/lexical_cast.hpp>

// we will add more permuted edges to the existing hypergraph
// with features, rule=1 score...

namespace cicada
{

  struct PermuteNoFeature
  {
    template <typename Features, typename Rule, typename Permutation>
    void operator()(Features& features, const Rule& rule, const Permutation& permutation)
    {
      
    }
  };

  struct PermuteFeature
  {
    
    std::vector<std::string, std::allocator<std::string> > non_terminal_symbols;
    
    template <typename Features, typename Rule, typename Permutation>
    void operator()(Features& features, const Rule& rule, const Permutation& permutation)
    {
      typedef Rule rule_type;
      
      non_terminal_symbols.clear();
      
      int non_terminal_pos = 0;
      typename rule_type::symbol_set_type::const_iterator siter_end = rule.source.end();
      for (typename rule_type::symbol_set_type::const_iterator siter = rule.source.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  non_terminal_symbols.push_back('[' + static_cast<const std::string&>(*siter).substr(1, siter->size() - 2) + '_' + boost::lexical_cast<std::string>(non_terminal_pos) + ']');
	  ++ non_terminal_pos;
	}
      
      std::string rule_string("permute:");
      rule_string += static_cast<const std::string&>(rule.lhs) + "->";
      
      int permutation_pos = 0;
      for (typename rule_type::symbol_set_type::const_iterator siter = rule.source.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  rule_string += non_terminal_symbols[permutation[permutation_pos]];
	  ++ permutation_pos;
	} else
	  rule_string += '<' + static_cast<const std::string&>(*siter) + '>';
      
      features[rule_string] = 1.0;
    }
  };

  template <typename FeatureFunction>
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
    
    typedef utils::bit_vector<1024> coverage_type;
    
    Permute(FeatureFunction __function,
	    const int __permute_size)
      : function(__function),
	permute_size(__permute_size) {}
    
    
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      target.clear();
      
      permutation_type permutation;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type& node = target.add_node();
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	
	  if (! edge_source.rule->target.empty())
	    throw std::runtime_error("we do not suppor synchronous-permutation");
	  
	  hypergraph_type::edge_type& edge = target.add_edge(edge_source.tails.begin(), edge_source.tails.end());
	  
	  permutation.clear();
	  for (int i = 0; i < edge_source.tails.size(); ++ i)
	    permutation.push_back(i);
	  
	  edge.rule = edge_source.rule;
	  
	  edge.features = edge_source.features;

	  function(edge.features, *edge_source.rule, permutation);
	  
	  target.connect_edge(edge.id, node.id);
	  
	  if (edge_source.tails.size() > 1) {
	    
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
	    if (permute_size < 0 || permute_size >= 1) {
	      
	      if (permute_size < 0 || permutation.size() <= permute_size + 2) {
		while (std::next_permutation(permutation.begin(), permutation.end())) {
		  
		  if (! is_valid_permutation(permutation)) continue;
		  
		  make_permuted_edge(node, edge_source, target, permutation);
		}
	      } else {
		// we will permute by traversing via transducer... actually, it is almost the same as word-based permutation with "skip"
		
		coverage_type coverage;
		permutation.clear();
		
		permute_recursive(node, edge_source, target, permutation, coverage);
	      }
	    }
	  }
	}
      }
      
      target.goal = source.goal;
    }
    

  private:

    void permute_recursive(const hypergraph_type::node_type& node,
			   const hypergraph_type::edge_type& edge_source,
			   hypergraph_type& graph,
			   const permutation_type& __permutation,
			   const coverage_type& __coverage)
    {
      const int arity = edge_source.tails.size();

      if (__permutation.size() < arity) {
	permutation_type permutation(__permutation);
	coverage_type    coverage(__coverage);

	const int index = permutation.size();
	
	permutation.push_back(0);
	for (int pos = utils::bithack::max(index - permute_size, 0); pos != utils::bithack::min(index + permute_size + 1, arity); ++ pos)
	  if (! coverage[pos]) {
	    permutation.back() = pos;
	    coverage.set(pos, true);
	    
	    permute_recursive(node, edge_source, graph, permutation, coverage);
	    
	    coverage.set(pos, false);
	  }
	
      } else {
	const permutation_type& permutation = __permutation;

	if (! is_valid_permutation(permutation)) return;
	
	make_permuted_edge(node, edge_source, graph, permutation);
      }
    }
      
    void make_permuted_edge(const hypergraph_type::node_type& node,
			    const hypergraph_type::edge_type& edge_source,
			    hypergraph_type& graph,
			    const permutation_type& permutation)
    {
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
      
      hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
      
      edge.rule.reset(new rule_type(*edge_source.rule));
      edge.rule->source = rule_type::symbol_set_type(source_phrase.begin(), source_phrase.end());
      
      if (*edge.rule == *edge_source.rule)
	edge.rule = edge_source.rule;
      
      edge.features = edge_source.features;
      
      function(edge.features, *edge_source.rule, permutation);
      
      graph.connect_edge(edge.id, node.id);
    }

    bool is_valid_permutation(const permutation_type& permutation)
    {
      if (permute_size < 0) return true;
      
      int differences = 0;
      for (int index = 0; index < permutation.size(); ++ index) {
	const int abs = utils::bithack::abs(index - permutation[index]);
	if (abs > permute_size)
	  return false;
	differences += abs;
      }
      return differences > 0;
    }


    FeatureFunction function;

    int permute_size;

    phrase_type source_phrase;
    phrase_type source_non_terminals;
    
    tails_type tails;
  };

  template <typename Function>
  inline
  void permute(const HyperGraph& source, HyperGraph& target, Function function, const int permute_size)
  {
    Permute<Function> permutation(function, permute_size);
    
    permutation(source, target);
  }
  
  inline
  void permute(const HyperGraph& source, HyperGraph& target, const int permute_size)
  {
    Permute<PermuteFeature> permutation(PermuteFeature(), permute_size);
    
    permutation(source, target);
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
