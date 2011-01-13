// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__APPLY__HPP__
#define __CICADA__OPERATION__APPLY__HPP__ 1

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {

    class Apply : public Operation
    {
    public:
      Apply(const std::string& parameter,
	    const model_type& __model,
	    const int __debug);
      
      void operator()(data_type& data) const;

      void assign(const weight_set_type& __weights);

      model_type model_local;

      const model_type& model;
      const weight_set_type* weights;
      int size;
      bool weights_one;
  
      bool exact;
      bool prune;
      bool grow;
      bool incremental;
  
      bool forced;
  
      int debug;
    };

  };
};


#endif
