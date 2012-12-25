// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_MIN__HPP__
#define __CICADA__NEURON_MIN__HPP__ 1

#include <cmath>

#include <vector>

#include <cicada/neuron/layer.hpp>

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
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output) {}
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Min(*this)); }

    private:
      typedef std::vector<int, std::allocator<int> > index_set_type;
      
    private:
      index_set_type indices;
      bool dimension;
    };
  };
};

#endif
