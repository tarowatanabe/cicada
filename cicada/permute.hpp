// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PERMUTE__HPP__
#define __CICADA__PERMUTE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/bithack.hpp>
#include <utils/bit_vector.hpp>
#include <utils/hashmurmur.hpp>

#include <boost/lexical_cast.hpp>

#include <google/dense_hash_set>

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
      typename rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
      for (typename rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  non_terminal_symbols.push_back('[' + static_cast<const std::string&>(*siter).substr(1, siter->size() - 2) + '_' + boost::lexical_cast<std::string>(non_terminal_pos) + ']');
	  ++ non_terminal_pos;
	}
      
      std::string rule_string("permute:");
      rule_string += static_cast<const std::string&>(rule.lhs) + "->";
      
      int permutation_pos = 0;
      for (typename rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  rule_string += non_terminal_symbols[permutation[permutation_pos]];
	  ++ permutation_pos;
	} else
	  rule_string += '<' + static_cast<const std::string&>(*siter) + '>';
      
      features[rule_string] = 1.0;
    }
  };

  template <typename Weights>
  struct PermuteFeatureCollapsed
  {
    std::vector<std::string, std::allocator<std::string> > non_terminal_symbols;
    
    const Weights& weights;
    const typename Weights::feature_type feature_name;

    PermuteFeatureCollapsed(const Weights& __weights)
      : weights(__weights), feature_name("permute:collapsed") {}
    
    template <typename Features, typename Rule, typename Permutation>
    void operator()(Features& features, const Rule& rule, const Permutation& permutation)
    {
      typedef Rule rule_type;
      

      non_terminal_symbols.clear();
      
      int non_terminal_pos = 0;
      typename rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
      for (typename rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  non_terminal_symbols.push_back('[' + static_cast<const std::string&>(*siter).substr(1, siter->size() - 2) + '_' + boost::lexical_cast<std::string>(non_terminal_pos) + ']');
	  ++ non_terminal_pos;
	}
      
      std::string rule_string("permute:");
      rule_string += static_cast<const std::string&>(rule.lhs) + "->";
      
      int permutation_pos = 0;
      for (typename rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  rule_string += non_terminal_symbols[permutation[permutation_pos]];
	  ++ permutation_pos;
	} else
	  rule_string += '<' + static_cast<const std::string&>(*siter) + '>';

      if (Weights::feature_type::exists(rule_string))
	features[feature_name] += weights[rule_string];
    }
  };

  struct PermuteFilter
  {
    template <typename Category>
    bool operator()(const Category& x) const
    {
      return false;
    };
  };

  template <typename FeatureFunction, typename Filter>
  struct Permute
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    typedef hypergraph_type::feature_set_type feature_set_type;
    
    typedef feature_set_type::feature_type    feature_type;
    
    typedef std::vector<int, std::allocator<int> > permutation_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > phrase_type;
    
    typedef utils::bit_vector<1024> coverage_type;
    typedef utils::bit_vector<1024> mask_type;
    

    struct rule_ptr_hash
    {
      size_t operator()(const rule_ptr_type& x) const
      {
	return (x ? hash_value(*x) : size_t(0));
      }
    };

    struct rule_ptr_equal
    {
      bool operator()(const rule_ptr_type& x, const rule_ptr_type& y) const
      {
	return (x == y || (x && y && *x == *y));
      }
    };

    typedef google::dense_hash_set<rule_ptr_type, rule_ptr_hash, rule_ptr_equal> rule_set_type;
    
    
    
    Permute(FeatureFunction __function,
	    Filter __filter,
	    const int __permute_size)
      : function(__function),
	filter(__filter),
	permute_size(__permute_size),
	rules() { rules.set_empty_key(rule_ptr_type()); }
    
    
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      rules.clear();

      target = source;
      
      permutation_type permutation;
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type& node = target.nodes[node_source.id];
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  hypergraph_type::edge_type& edge = target.edges[edge_source.id];

	  permutation.clear();
	  for (size_t i = 0; i != edge_source.tails.size(); ++ i)
	    permutation.push_back(i);
	  
	  function(edge.features, *edge.rule, permutation);
	  
	  if (edge_source.tails.size() > 1) {
	    mask_type mask;
	    
	    bool has_constraint = false;
	    int pos = 0;
	    permuted_non_terminals.clear();
	    rule_type::symbol_set_type::const_iterator siter_end = edge_source.rule->rhs.end();
	    for (rule_type::symbol_set_type::const_iterator siter = edge_source.rule->rhs.begin(); siter != siter_end; ++ siter, ++ pos)
	      if (siter->is_non_terminal()) {
		permuted_non_terminals.push_back(*siter);
		
		const bool constraint = filter(siter->non_terminal());
		
		mask.set(pos, constraint);
		has_constraint |= constraint;
	      }
	    
	    permuted_phrase.clear();
	    permuted_phrase.insert(permuted_phrase.end(), edge_source.rule->rhs.begin(), edge_source.rule->rhs.end());
	    
	    tails.clear();
	    tails.insert(tails.end(), edge_source.tails.begin(), edge_source.tails.end());
	    
	    // perform permutation...
	    // minus for all permutation
	    // zero for no-permutation
	    if (permute_size < 0 || permute_size >= 1) {
	      
	      if (! has_constraint && (permute_size < 0 || static_cast<int>(permutation.size()) <= permute_size + 2)) {
		while (std::next_permutation(permutation.begin(), permutation.end())) {
		  
		  if (! is_valid_permutation(permutation)) continue;
		  
		  make_permuted_edge(node, edge_source, target, permutation);
		}
	      } else {
		// we will permute by traversing via transducer... actually, it is almost the same as word-based permutation with "skip"
		
		coverage_type coverage;
		permutation.clear();
		
		permute_recursive(node, edge_source, target, permutation, coverage, mask);
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
			   const coverage_type& __coverage,
			   const mask_type& mask)
    {
      const int arity = edge_source.tails.size();

      if (static_cast<int>(__permutation.size()) < arity) {
	permutation_type permutation(__permutation);
	coverage_type    coverage(__coverage);
	
	const int index = permutation.size();

	permutation.push_back(0);
	
	if (! coverage[index]) {
	  permutation.back() = index;
	  coverage.set(index, true);
	  
	  permute_recursive(node, edge_source, graph, permutation, coverage, mask);
	  
	  coverage.set(index, false);
	}
	
	if (! mask[index]) {
	  const int pos_first = (permute_size <= 0 ? 0     : utils::bithack::max(index - permute_size, 0));
	  const int pos_last  = (permute_size <= 0 ? arity : utils::bithack::min(index + permute_size + 1, arity));
	  
	  for (int pos = index - 1; pos >= pos_first && ! mask[pos]; -- pos)
	    if (! coverage[pos]) {
	      permutation.back() = pos;
	      coverage.set(pos, true);
	      
	      permute_recursive(node, edge_source, graph, permutation, coverage, mask);
	      
	      coverage.set(pos, false);
	    }
	  
	  for (int pos = index + 1; pos < pos_last && ! mask[pos]; ++ pos)
	    if (! coverage[pos]) {
	      permutation.back() = pos;
	      coverage.set(pos, true);
	      
	      permute_recursive(node, edge_source, graph, permutation, coverage, mask);
	      
	      coverage.set(pos, false);
	    }
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
      for (size_t i = 0; i != tails.size(); ++ i)
	tails[i] = edge_source.tails[permutation[i]];
      
      // permute source-phrase
      int non_terminal_pos = 0;
      phrase_type::iterator siter_end = permuted_phrase.end();
      for (phrase_type::iterator siter = permuted_phrase.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal()) {
	  *siter = permuted_non_terminals[permutation[non_terminal_pos]];
	  
	  ++ non_terminal_pos;
	}
      
      hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());

      if (permuted_phrase.size() == edge_source.rule->rhs.size()
	  && std::equal(permuted_phrase.begin(), permuted_phrase.end(), edge_source.rule->rhs.begin()))
	edge.rule = edge_source.rule;
      else {
	edge.rule.reset(new rule_type(*edge_source.rule));
	edge.rule->rhs = rule_type::symbol_set_type(permuted_phrase.begin(), permuted_phrase.end());
	edge.rule = *(rules.insert(edge.rule).first);
      }
      
      edge.features = edge_source.features;
      edge.attributes = edge_source.attributes;
      
      function(edge.features, *edge_source.rule, permutation);
      
      graph.connect_edge(edge.id, node.id);
    }

    bool is_valid_permutation(const permutation_type& permutation)
    {
      if (permute_size < 0) return true;
      
      int differences = 0;
      for (int index = 0; index < static_cast<int>(permutation.size()); ++ index) {
	const int abs = utils::bithack::abs(index - permutation[index]);
	if (abs > permute_size)
	  return false;
	differences += abs;
      }
      return differences > 0;
    }
    
    
    FeatureFunction function;
    Filter          filter;

    int permute_size;

    phrase_type permuted_phrase;
    phrase_type permuted_non_terminals;
    
    tails_type tails;

    rule_set_type rules;
  };

  template <typename Function, typename Filter>
  inline
  void permute(const HyperGraph& source, HyperGraph& target, Function function, Filter filter, const int permute_size)
  {
    Permute<Function, Filter> permutation(function, filter, permute_size);
    
    permutation(source, target);
  }

  template <typename Function>
  inline
  void permute(const HyperGraph& source, HyperGraph& target, Function function, const int permute_size)
  {
    Permute<Function, PermuteFilter> permutation(function, PermuteFilter(), permute_size);
    
    permutation(source, target);
  }
  
  inline
  void permute(const HyperGraph& source, HyperGraph& target, const int permute_size)
  {
    Permute<PermuteFeature, PermuteFilter> permutation(PermuteFeature(), PermuteFilter(), permute_size);
    
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
