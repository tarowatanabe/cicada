// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__EXPECTED_NGRAM__HPP__
#define __CICADA__OPERATION__EXPECTED_NGRAM__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class ExpectedNGram : public Operation
    {
    public:
      ExpectedNGram(const std::string& parameter, const int __debug);

      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);
  
      int order;
      bool bos_eos;
      
      const weights_path_type* weights;
      const weight_set_type*   weights_assigned;
      bool weights_one;
      
      double scale;

      int debug;
    };

  };
};


#endif
