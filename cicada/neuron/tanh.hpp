// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_TANH__HPP__
#define __CICADA__NEURON_TANH__HPP__ 1

#include <cmath>

namespace cicada
{
  namespace neuron
  {
    struct tanh
    {
      double operator()(const double& value) const
      {
        return std::tanh(vaue);
      }
      
      double derivative(const double& value) const
      {
	const double z = operator()(value);
	return 1.0 - z * z;
      }
    };
  };
};

#endif
