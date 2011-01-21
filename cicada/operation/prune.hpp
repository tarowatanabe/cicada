// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PRUNE__HPP__
#define __CICADA__OPERATION__PRUNE__HPP__ 1


#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class Prune : public Operation
    {
    public:
      Prune(const std::string& parameter, const int __debug);

      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);

      const weight_set_type* weights;
  
      size_t kbest;
      double beam;
      double density;
      double scale;
  
      bool weights_one;

      bool semiring_tropical;
      bool semiring_logprob;
      bool semiring_log;
  
      int debug;
    };

  };
};


#endif
