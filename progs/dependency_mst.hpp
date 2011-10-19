// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __DEPENDENCY__MST__HPP__
#define __DEPENDENCY__MST__HPP__ 1

#include <vector>
#include <utility>

#include "cicada/optimize/mst.hpp"

// minimum-spanning-tree algorithm

struct DependencyMST
{
  typedef std::pair<int, int> dep_type;
  typedef std::vector<dep_type, std::allocator<dep_type> > dep_set_type;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    
    
  }
};

struct DependencyMSTSingleRoot
{
  typedef std::pair<int, int> dep_type;
  typedef std::vector<dep_type, std::allocator<dep_type> > dep_set_type;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    
    
  }
};

#endif
