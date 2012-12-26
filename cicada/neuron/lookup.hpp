// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_LOOKUP__HPP__
#define __CICADA__NEURON_LOOKUP__HPP__ 1

//
// lookup-table layer which map input, a sequence of id, into a sequence of tensor
//

#include <cicada/neuron/layer.hpp>

#include <vector>

namespace cicada
{
  namespace neuron
  {
    class Lookup : public Layer
    {
    public:
      Lookup(size_type __size=1) : size(__size) {}
      Lookup(const tensor_type& __weight) : size(__weight.rows()), weight(__weight) {}
      
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Lookup(*this)); }
      
    public:
      size_type     size;
      tensor_type   weight;
      tensor_type   gradient_weight;
    };
  };
};

#endif
