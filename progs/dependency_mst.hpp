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
  typedef std::vector<int, std::allocator<int> > dep_type;

  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    
    mst.viterbi_forest(scores, deps);
    
    for (size_t i = 0; i != sentence_size; ++ i)
      dependency[i] = (deps[i + 1] ? deps[i + 1] - 1 : 0);
  }

  cicada::optimize::MST mst;
  dep_type deps;
};

struct DependencyMSTSingleRoot
{
  typedef std::vector<int, std::allocator<int> > dep_type;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    
    mst.viterbi_tree(scores, deps);
    
    for (size_t i = 0; i != sentence_size; ++ i)
      dependency[i] = (deps[i + 1] ? deps[i + 1] - 1 : 0);
  }

  cicada::optimize::MST mst;
  dep_type deps;
};

#endif
