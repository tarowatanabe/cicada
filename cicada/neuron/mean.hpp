// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_MEAN__HPP__
#define __CICADA__NEURON_MEAN__HPP__ 1

#include <cmath>

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Mean : public Layer
    {
    public:
      Mean(const bool __dimension=false) : dimension(__dimension) {}

    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output) {}
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Mean(*this)); }
      virtual std::ostream& write(std::ostream& os) const;
    private:
      bool dimension;
    };
  };
};

#endif
