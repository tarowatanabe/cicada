// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__LEXICAL_CAST__HPP__
#define __UTILS__LEXICAL_CAST__HPP__ 1

#include <boost/lexical_cast.hpp>

namespace utils
{
  template <typename Target, typename Source>
  inline
  Target lexical_cast(const Source& arg)
  {
    return boost::lexical_cast<Target>(arg);
  }

  template <>
  inline
  std::string lexical_cast<std::string, bool>(const bool& arg)
  {
    return (arg ? "true" : "false");
  }
  
  template <>
  inline
  bool lexical_cast<bool, std::string>(const std::string& arg)
  {
    if (strcasecmp(arg.c_str(), "true") == 0)
      return true;
    else if (strcasecmp(arg.c_str(), "false") == 0)
      return false;
    else if (strcasecmp(arg.c_str(), "yes") == 0)
      return true;
    else if (strcasecmp(arg.c_str(), "no") == 0)
      return false;
    else if (strcasecmp(arg.c_str(), "nil") == 0)
      return false;
    else if (atoi(arg.c_str()) > 0)
      return true;
    else
      return false;
  }
};

#endif
