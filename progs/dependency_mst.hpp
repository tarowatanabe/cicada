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
  typedef boost::numeric::ublas::matrix<double> matrix_type;

  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    // we will transform scores into MST's score format:
    // diagonal contains probability into root...
    
    const size_t sentence_size = dependency.size();
    const size_t matrix_size   = scores.size1();

    if (sentence_size + 1 != matrix_size)
      throw std::runtime_error("scores and dependency size do not match");
    
    matrix.resize(matrix_size - 1, matrix_size - 1);
    matrix.clear();
    
    for (size_t head = 1; head != matrix_size; ++ head)
      for (size_t dep = 1; dep != matrix_size; ++ dep)
	matrix(head - 1, dep - 1) = scores(head, dep);
    
    // assign root score at diagonal
    for (size_t dep = 1; dep != matrix_size; ++ dep)
      matrix(dep - 1, dep - 1) = scores(0, dep);
    
    mst.viterbi_forest(matrix, dependency);
  }

  matrix_type matrix;
  cicada::optimize::MST mst;
};

struct DependencyMSTSingleRoot
{
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    const size_t matrix_size   = scores.size1();

    if (sentence_size + 1 != matrix_size)
      throw std::runtime_error("scores and dependency size do not match");
    
    matrix.resize(matrix_size - 1, matrix_size - 1);
    matrix.clear();
    
    for (size_t head = 1; head != matrix_size; ++ head)
      for (size_t dep = 1; dep != matrix_size; ++ dep)
	matrix(head - 1, dep - 1) = scores(head, dep);
    
    // assign root score at diagonal
    for (size_t dep = 1; dep != matrix_size; ++ dep)
      matrix(dep - 1, dep - 1) = scores(0, dep);
    
    mst.viterbi_tree(matrix, dependency);
  }
  
  matrix_type matrix;
  cicada::optimize::MST mst;
};

#endif
