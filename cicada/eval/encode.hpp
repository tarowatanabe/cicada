// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__ENCODE__HPP__
#define __CICADA__EVAL__ENCODE__HPP__ 1

#include <sstream>
#include <iterator>

#include <utils/base64.hpp>

namespace cicada
{
  namespace eval
  {
    
    struct escape_base64
    {
      escape_base64(const double& __value) : value(__value) {}
      const double& value;
      
      friend
      std::ostream& operator<<(std::ostream& os, const escape_base64& x)
      {
	const std::string base64 = utils::encode_base64(x.value);
	
	os << '\"';
	std::string::const_iterator iter_end = base64.end();
	for (std::string::const_iterator iter = base64.begin(); iter != iter_end; ++ iter)
	  if (*iter == '/')
	    os << "\\/";
	  else
	    os << *iter;
	os << '\"';
	
	return os;
      }
    };
    
    inline
    escape_base64 escaper(const double& x)
    {
      return escape_base64(x);
    }
    
  };
};

#endif
