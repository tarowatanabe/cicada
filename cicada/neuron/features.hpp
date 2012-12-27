// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NEURON_FEATURES__HPP__
#define __CICADA__NEURON_FEATURES__HPP__ 1

//
// features-table layer which map feature-fector input as a table (perform reordering)
//

#include <vector>

#include <boost/functional/hash/hash.hpp>

#include <cicada/neuron/layer.hpp>
#include <cicada/feature.hpp>

#include <utils/indexed_set.hpp>

namespace cicada
{
  namespace neuron
  {
    class Features : public Layer
    {
    public:
      typedef cicada::Feature feature_type;
      
    private:
      typedef utils::indexed_set<feature_type,
				 boost::hash<feature_type>,
				 std::equal_to<feature_type>,
				 std::allocator<feature_type > > feature_map_type;

    public:
      typedef feature_map_type::iterator       iterator;
      typedef feature_map_type::const_iterator const_iterator;
      typedef feature_map_type::index_type     index_type;
      
    public:
      Features() {}
      
      template <typename Iterator>
      Features(Iterator first, Iterator last) : features(first, last) {}
      
    public:
      const feature_type& operator[](index_type x) { return features[x]; }
      
      template <typename Iterator>
      void insert(Iterator first, Iterator last) { features.insert(first, last); }
      
      std::pair<iterator, bool> insert(const feature_type& x) { return features.insert(x); }
      const_iterator find(const feature_type& x) const { return features.find(x); }
      
      const_iterator begin() const { return features.begin(); }
      const_iterator end() const { return features.end(); }
      
      void clear() { features.clear(); }
      
      size_type size() const { return features.size(); }
      bool empty() const { return features.empty(); }
      
    public:
      virtual void forward(const tensor_type& data_input);
      virtual void backward(const tensor_type& data_input, const tensor_type& gradient_output);
      virtual void accumulate(const tensor_type& data_input, const tensor_type& gradient_output) {}
      virtual layer_ptr_type clone() const { return layer_ptr_type(new Features(*this)); }
      virtual std::ostream& write(std::ostream& os) const;

      
    private:
      feature_map_type features;
    };
  };
};

#endif
