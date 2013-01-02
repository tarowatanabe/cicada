// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// A Sparse feature vector representation using the eigen's sparse vector
//

#ifndef __CICADA__FEATURE_VECTOR_SPARSE__HPP__
#define __CICADA__FEATURE_VECTOR_SPARSE__HPP__ 1

#include <cicada/feature.hpp>

#include <Eigen/SparseCore>

namespace cicada
{
  // forward declaration...
  template <typename Tp, typename Alloc >
  class WeightVector;

  class FeatureVectorCompact;
  
  template <typename Tp, typename Alloc >
  class FeatureVectorLinear;
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class FeatureVectorSparse
  {
  public:
    typedef cicada::Feature feature_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  private:
    typedef std::pair<const feature_type, Tp> __value_type;
    typedef Eigen::SparseVector<Tp> vector_type;
    
  public:
    
    
    
    
  };
};

#include <cicada/weight_vector.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/feature_vector_compact.hpp>
#include <cicada/feature_vector_linear.hpp>


#endif
