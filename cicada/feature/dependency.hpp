// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__DEPENDENCY__HPP__
#define __CICADA__FEATURE__DEPENDENCY__HPP__ 1

#include <cicada/cluster_stemmer.hpp>
#include <cicada/feature_function.hpp>

namespace cicada
{
  namespace feature
  {
    class Dependency : public FeatureFunction
    {
    public:
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
    public:
      Dependency();
      Dependency(const Dependency& x);
      Dependency& operator=(const Dependency& x);
      
    private:
      const attribute_type attr_dependency_pos;
      const attribute_type attr_dependency_head;
      const attribute_type attr_dependency_dependent;
    };
  };
};

#endif
