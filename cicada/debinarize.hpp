// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DEBINARIZE__HPP__
#define __CICADA__DEBINARIZE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/chart.hpp>

namespace cicada
{
  struct Debinarize
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      // debinarization by stripping off the ^ from syntactic categories...
      
      // bottom-up topological order to find binarised antecedents...
      if (! source.is_valid()) {
	target.clear();
	return;
      }
      
      target = source;
      
      hypergraph_type::node_set_type::const_iterator niter_end = target.nodes.end();
      for (hypergraph_type::node_set_type::conts_iterator niter = target.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = target.edges[*eiter];
	  
	  // search for antecedent nodes, and seek the binarized label..
	  // if found, try merge! 
	  
	  
	}
      }
    }
  };
};

#endif
