// -*- mode: c++ -*-

#ifndef __CICADA__BINARIZE__HPP__
#define __CICADA__BINARIZE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <utils/bithack.hpp>

#include <boost/lexical_cast.hpp>

namespace cicada
{
  
  struct BinarizeNoFeature
  {
    template <typename Features, typename Rule>
    void operator()(Features& features, const Rule& rule, const Rule& binarized)
    {
      
    }
  };

  struct BinarizeFeature
  {
    template <typename Features, typename Rule>
    void operator()(Features& features, const Rule& rule, const Rule& binarized)
    {
      // we know the parent of binarized constituent...
      features["binarize:" + static_cast<const std::string&>(rule.lhs) + static_cast<const std::string&>(binarized.lhs)] = 1.0;
    }
  };
  
  template <typename Weights>
  struct BinarizeFeatureCollapsed
  {
    const Weights& weights;
    const typename Weights::feature_type feature_name;
    
    BinarizeFeatureCollapsed(const Weights& __weights)
      : weights(__weights), feature_name("binarize:collapsed") {}

    template <typename Features, typename Rule>
    void operator()(Features& features, const Rule& rule, const Rule& binarized)
    {
      // we know the parent of binarized constituent...
      const std::string feature = "binarize:" + static_cast<const std::string&>(rule.lhs) + static_cast<const std::string&>(binarized.lhs);
      
      if (Weights::feature_type::exists(feature))
	features[feature_name] = weights[feature];
    }
  };

  

  template <typename FeatureFunction>
  struct Binarize
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    typedef hypergraph_type::feature_set_type feature_set_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > phrase_type;

    typedef std::vector<bool, std::allocator<bool> > removed_type;
    
    struct filter
    {
      const removed_type& removed;
      
      filter(const removed_type& __removed) : removed(__removed) {}
      
      template <typename Edge>
      bool operator()(const Edge& edge) const
      {
	return removed[edge.id];
      }
    };
    
    Binarize(FeatureFunction __function, int __binarize_size)
      : function(__function), binarize_size(__binarize_size) {}
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      target = source;

      phrase_type source_rule(2);
      hypergraph_type::edge_type::node_set_type tails(2);
      
      // we will traverse source-side in order to avoid confusion with newly created nodes...

      removed_type removed(source.edges.size(), false);
      
      hypergraph_type::node_set_type::const_iterator niter_end = source.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = source.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node_source = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node_source.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node_source.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge_source = source.edges[*eiter];
	  
	  if (! edge_source.rule->target.empty())
	    throw std::runtime_error("we do not suppor synchronous-binarization");
	  
	  if (edge_source.tails.size() >= 3 && edge_source.tails.size() >= binarize_size) {
	    
	    if (edge_source.tails.size() != edge_source.rule->arity)
	      throw std::runtime_error("we do not support terminal-mixed rules (aka Hiero rules)");
	    
	    removed[edge_source.id] = true;
	    
	    hypergraph_type::id_type head = edge_source.head;
	    symbol_type non_terminal_head = edge_source.rule->lhs;
	    
	    const int arity = edge_source.tails.size();
	    
	    for (int i = 0; i < arity - 2; ++ i) {
	      
	      const symbol_type non_terminal_new = merge_non_terminals(edge_source.rule->source.begin(), edge_source.rule->source.end() - i - 1);
	      
	      hypergraph_type::node_type& node_new = target.add_node();
	      
	      tails.front() = node_new.id;
	      tails.back() = edge_source.tails[arity - i - 1];

	      source_rule.front() = non_terminal_new;
	      source_rule.back()  = edge_source.rule->source[arity - i - 1];
	      
	      hypergraph_type::edge_type& edge_new = target.add_edge(tails.begin(), tails.end());
	      edge_new.rule.reset(new rule_type(non_terminal_head,
						rule_type::symbol_set_type(source_rule.begin(), source_rule.end()), 
						rule_type::symbol_set_type(),
						2));
	      
	      target.connect_edge(edge_new.id, head);

	      function(edge_new.features, *edge_source.rule, *edge_new.rule);
	      
	      head = node_new.id;
	      non_terminal_head = non_terminal_new;
	    }
	    
	    hypergraph_type::edge_type& edge_new = target.add_edge(edge_source.tails.begin(), edge_source.tails.begin() + 2);
	    
	    source_rule.front() = edge_source.rule->source[0];
	    source_rule.back()  = edge_source.rule->source[1];
	    
	    edge_new.rule.reset(new rule_type(non_terminal_head,
					      rule_type::symbol_set_type(source_rule.begin(), source_rule.end()), 
					      rule_type::symbol_set_type(),
					      2));
	    // assign features here...
	    edge_new.features = edge_source.features;
	    
	    function(edge_new.features, *edge_source.rule, *edge_new.rule);
	    
	    target.connect_edge(edge_new.id, head);
	  } 
	}
      }
      
      // further resize...
      removed.resize(target.edges.size(), false);
      
      hypergraph_type graph_removed;
      
      topologically_sort(target, graph_removed, filter(removed));
      
      target.swap(graph_removed);
    }

  private:
    template <typename Iterator>
    symbol_type merge_non_terminals(Iterator first, Iterator last)
    {
      // create new binary rule by merging non-terminals at pos to the end
      
      std::string non_terminal;
      
      bool initial = true;
      for (/**/; first != last; ++ first) {
	if (first->is_terminal())
	  throw std::runtime_error("we do not support terminal mixed rules");
	
	const std::string& raw = static_cast<const std::string&>(*first);
	
	if (! initial)
	  non_terminal += '|';
	
	non_terminal += raw.substr(1, raw.size() - 2);
	
	initial = false;
      }
      
      return '[' + non_terminal + ']';
    }

  private:
    FeatureFunction function;
    
    int binarize_size;
  };
  
  template <typename Function>
  inline
  void binarize(const HyperGraph& source, HyperGraph& target, Function function, const int binarize_size=0)
  {
    Binarize<Function> binarizer(function, binarize_size);
    
    binarizer(source, target);
  }
  
  inline
  void binarize(const HyperGraph& source, HyperGraph& target, const int binarize_size=0)
  {
    Binarize<BinarizeFeature> binarizer(BinarizeFeature(), binarize_size);
    
    binarizer(source, target);
  }

  inline
  void binarize(HyperGraph& source, const int binarize_size=0)
  {
    HyperGraph target;

    Binarize<BinarizeFeature> binarizer(BinarizeFeature(), binarize_size);
    
    binarizer(source, target);
    
    source.swap(target);
  }

};

#endif
