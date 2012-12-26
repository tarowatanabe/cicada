// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_PARALLEL__HPP__
#define __CICADA__NEURON_PARALLEL__HPP__ 1

#include <vector>

#include <cicada/neuron/layer.hpp>

namespace cicada
{
  namespace neuron
  {
    class Parallel : public Layer
    {
    public:
      Parallel(const bool __dimension_input=false, const bool __dimension_output=false)
	: dimension_input(__dimension_input), dimension_output(__dimension_output) {}

      template <typename Iterator>
      Parallel(Iterator first,
	       Iterator last,
	       const bool __dimension_input,
	       const bool __dimension_output)
	: layers(first, last),
	  dimension_input(__dimension_input),
	  dimension_output(__dimension_output) {}
      
    public:
      void clear() { layers.clear(); }
      
      void push_back(const layer_ptr_type& layer) { layers.push_back(layer); }
      void push_back(Layer* layer) { layers.push_back(layer_ptr_type(layer)); }
      
      inline const layer_ptr_type& operator[](size_type pos) const { return layers[pos]; }
      inline       layer_ptr_type& operator[](size_type pos)       { return layers[pos]; }
      
      inline const layer_ptr_type& front() const { return layers.front(); }
      inline       layer_ptr_type& front()       { return layers.front(); }
      
      inline const layer_ptr_type& back() const { return layers.back(); }
      inline       layer_ptr_type& back()       { return layers.back(); }
      
      size_type size() const { return layers.size(); }
      bool empty() const { return layers.empty(); }
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual layer_ptr_type clone() const;
      virtual std::ostream& write(std::ostream& os) const;
    private:
      typedef std::vector<layer_ptr_type, std::allocator<layer_ptr_type> > layer_set_type;
      typedef std::vector<size_type, std::allocator<size_type> > size_set_type;
      
    private:
      layer_set_type layers;
      size_set_type  sizes;
      bool dimension_input;
      bool dimension_output;
    };
  };
};

#endif
