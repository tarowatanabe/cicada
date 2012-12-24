// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_MIN__HPP__
#define __CICADA__NEURON_MIN__HPP__ 1

#include <cmath>

#include <cicada/neuron/layer.hpp>

#include <Eigen/Core>

namespace cicada
{
  namespace neuron
  {
    class Min : public Layer
    {
    public:
      Min(const bool __dimension=false) : dimension(__dimension) {}
      
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);

    private:
      typedef Eigen::Matrix<int, Eigen::Dynamic, 1> index_set_type;
      
    private:
      index_set_type indices;
      bool dimension;
    };
  };
};

#endif
