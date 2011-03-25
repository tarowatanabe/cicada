// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__GENERATE__HPP__
#define __CICADA__OPERATION__GENERATE__HPP__ 1

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class GenerateEarley : public Operation
    {
    public:
      GenerateEarley(const std::string& parameter,
		     const grammar_type& __grammar,
		     const std::string& __goal,
		     const int __debug);
  
      void operator()(data_type& data) const;
      
      int depth;
      int width;
      
      int debug;
    };

  };
};


#endif
