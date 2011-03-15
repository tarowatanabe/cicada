// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__OUTPUT__HPP__
#define __CICADA__OPERATION__OUTPUT__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

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
  
      const weights_path_type* weights;
      const weight_set_type*   weights_assigned;
      bool weights_one;
      int  kbest_size;
      bool kbest_unique;
      
      std::string insertion_prefix;

      bool yield_string;
      bool yield_terminal_pos;
      bool yield_tree;
      bool yield_graphviz;
      bool yield_treebank;
      bool yield_alignment;
      bool yield_span;

      bool debinarize;
      bool graphviz;
      bool statistics;
      bool lattice_mode;
      bool forest_mode;
      bool no_id;
  
      int debug;
    };

  };
};


#endif
