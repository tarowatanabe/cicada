// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__OUTPUT__HPP__
#define __CICADA__OPERATION__OUTPUT__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/semiring.hpp>
#include <cicada/kbest.hpp>
#include <cicada/graphviz.hpp>
#include <cicada/inside_outside.hpp>

#include <cicada/operation/functional.hpp>
#include <cicada/operation/traversal.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/sgi_hash_map.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  namespace operation
  {


    class Output : public Operation
    {
    public:
      Output(const std::string& parameter, output_data_type& __output_data, const int __debug);

      void assign(const weight_set_type& __weights);
  
      void clear();
  
      void operator()(data_type& data) const;
      
      output_data_type& output_data;
  
      path_type file;
      path_type directory;
  
      const weight_set_type* weights;
      bool weights_one;
      int  kbest_size;
      bool kbest_unique;
      
      std::string insertion_prefix;

      bool yield_string;
      bool yield_tree;
      bool yield_alignment;

      bool graphviz;
      bool statistics;
  
      int debug;
    };

  };
};


#endif
