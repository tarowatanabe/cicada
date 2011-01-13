// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__NORMALIZE__HPP__
#define __CICADA__OPERATION__NORMALIZE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class Normalize : public Operation
    {
    public:
      Normalize(const std::string& parameter,
		const int __debug);
      
      template <typename FeaturePrefix, typename Feature>
      inline
      bool equal_prefix(const FeaturePrefix& prefix, const Feature& x) const
      {
	return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
      }
      
      void operator()(data_type& data) const;

      std::string feature_prefix;
      int debug;
    };
  };
};


#endif
