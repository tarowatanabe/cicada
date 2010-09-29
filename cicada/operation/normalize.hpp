// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__NORMALIZE__HPP__
#define __CICADA__OPERATION__NORMALIZE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/generate.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>

namespace cicada
{
  namespace operation
  {
    class Normalize : public Operation
    {
    public:
      Normalize(const std::string& parameter,
		const int __debug)
	: debug(__debug)
      {
	
      }
      

      void operator()(data_type& data) const
      {
	
      }
    };
  };
};


#endif
