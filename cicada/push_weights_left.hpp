// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PUSH_WEIGHTS_LEFT__HPP__
#define __CICADA__PUSH_WEIGHTS_LEFT__HPP__ 1

#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <utility>

#include <cicada/hypergraph.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sort_topologically.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>
#include <utils/bithack.hpp>


namespace cicada
{
  
  struct PushWeightsLeft
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    void operator()(const hypergraph_type& source, hypergraph_type& target)
    {
      target.clear();

      if (! source.is_valid())
	return;
      
      
    }
  };
  
  
  inline
  void push_weights_left(const HyperGraph& source, HyperGraph& target)
  {
    PushWeightsLeft()(source, target);
  }
  
  inline
  void push_weights_left(HyperGraph& graph)
  {
    HyperGraph x;
    push_weights_left(graph, x);
    graph.swap(x);
  }
};

#endif
