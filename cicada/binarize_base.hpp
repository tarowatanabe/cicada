// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BINARIZE_BASE__HPP__
#define __CICADA__BINARIZE_BASE__HPP__ 1

#include <stdexcept>
#include <algorithm>
#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/sort.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/chart.hpp>

namespace cicada
{
  
  struct BinarizeBase
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > phrase_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > context_type;

    typedef utils::chart<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_chart_type;
    typedef utils::chart<symbol_type, std::allocator<symbol_type> > label_chart_type;

    typedef std::vector<int, std::allocator<int> > position_set_type;
    
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
  };
};

#endif
