// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __DEPENDENCY__PERMUTATION__HPP__
#define __DEPENDENCY__PERMUTATION__HPP__ 1

#include <vector>
#include <utility>

#include "kuhn_munkres.hpp"

struct DependencyPermutation
{
  typedef std::pair<int, int> dep_type;
  typedef std::vector<dep_type, std::allocator<dep_type> > dep_set_type;

  template <typename Dependency>
  struct insert_dependency
  {
    Dependency& dependency;
    
    insert_dependency(Dependency& __dependency) : dependency(__dependency) {}

    template <typename Edge>
    insert_dependency& operator=(const Edge& edge)
    {
      // first is head, second is dependent
      
      dependency[edge.second] = edge.first + 1;
      
      return *this;
    }
    
    insert_dependency& operator*() { return *this; }
    insert_dependency& operator++() { return *this; }
    insert_dependency operator++(int) { return *this; }
  };
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    kuhn_munkres_assignment(scores, insert_dependency<Dependency>(dependency));
  }
};

#endif
