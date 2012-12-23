// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_LAYER__HPP__
#define __CICADA__NEURON_LAYER__HPP__ 1

//
// Base class for NN layer
// We assume that input/outputs are simply one-dimensional (thus, we do not handle image, but "sentence")
//
// We will hold all the data in single precision, not double precision for compatibility
// with GPU (via tensor_type).
// For GPU implementation, we will use openCL...
//

#include <boost/shared_ptr.hpp>

namespace cicada
{
  namespace neuron
  {
    class Layer
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef Layer                         layer_type;
      typedef boost::shared_ptr<layer_type> layer_ptr_type;
      
    public:
      Layer() {}
      virtual ~Layer() {}
      
    public:
      // we want to allow sparse input.... currently, we support only dense input....
      virtual void forward(const tensor_type& data_input) = 0;
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output) = 0;
      
    public:
      // we will explicitly comment out unused data structure...!
      
      //tensor_type data_input;
      tensor_type data_output;
      
      tensor_type gradient_input;
      //tensor_type gradient_output;
    };
  };
};

#endif
