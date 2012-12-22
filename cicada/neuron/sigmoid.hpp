// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_SIGMOID__HPP__
#define __CICADA__NEURON_SIGMOID__HPP__ 1

#include <cmath>

namespace cicada
{
  namespace neuron
  {
    struct sigmoid
    {
      double operator()(const double& value) const
      {
	return 1.0 / (1.0 + std::exp(- value));
      }
      
      double derivative(const double& value) const
      {
	const double z = operator()(value);
	return (1.0 - z) * z;
      }
    };
  };
};

#endif
