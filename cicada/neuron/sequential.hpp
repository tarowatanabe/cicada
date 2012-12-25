// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_SEQUENTIAL__HPP__
#define __CICADA__NEURON_SEQUENTIAL__HPP__ 1

#include <vector>

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Sequential : public Layer
    {
    public:
      void clear() { layers.clear(); }
      
      void push_back(const layer_ptr_type& layer) { layers.push_back(layer); }
      
      inline const layer_ptr_type& operator[](size_type pos) const { return layers[pos]; }
      inline       layer_ptr_type& operator[](size_type pos)       { return layers[pos]; }
      
      size_type size() const { return layers.size(); }
      bool empty() const { return layers.empty(); }
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual layer_ptr_type clone() const;
    private:
      typedef std::vector<layer_ptr_type, std::allocator<layer_ptr_type> > layer_set_type;
      
    private:
      layer_set_type layers;
    };
  };
};

#endif
