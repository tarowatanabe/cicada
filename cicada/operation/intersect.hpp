// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__INTERSECT__HPP__
#define __CICADA__OPERATION__INTERSECT__HPP__ 1

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class Intersect : public Operation
    {
    public:
      Intersect(const std::string& parameter, const int __debug);

      void operator()(data_type& data) const;

      bool lattice_mode;
      bool target_mode;
  
      int debug;
    };

  };
};


#endif
