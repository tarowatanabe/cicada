// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__LEXICAL_CAST__HPP__
#define __UTILS__LEXICAL_CAST__HPP__ 1

#include <boost/lexical_cast.hpp>

#include <utils/piece.hpp>

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
  bool lexical_cast<bool, utils::ipiece>(const utils::ipiece& arg)
  {
    if (arg == "true")
      return true;
    else if (arg == "false")
      return false;
    else if (arg == "yes")
      return true;
    else if (arg == "no")
      return false;
    else if (arg == "nil")
      return false;
    else if (atol(arg.c_str()) > 0)
      return true;
    else
      return false;
  }
  
  template <>
  inline
  utils::piece lexical_cast<utils::piece, bool>(const bool& arg)
  {
    return (arg ? "true" : "false");
  }

  template <>
  inline
  bool lexical_cast<bool, utils::piece>(const utils::piece& arg)
  {
    return lexical_cast<bool>(utils::ipiece(arg));
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
    return lexical_cast<bool>(utils::ipiece(arg));
  }

};

#endif
