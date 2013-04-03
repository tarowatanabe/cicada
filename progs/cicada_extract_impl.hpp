// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXTRACT_IMPL__HPP__
#define __CICADA__EXTRACT_IMPL__HPP__ 1

#include <cstddef>
#include <stdint.h>

#include <iostream>
#include <stdexcept>

struct Statistic
{
  typedef int64_t count_type;
  
  Statistic(const count_type __bitext=0) : bitext(__bitext) {}
  
  count_type bitext;

  Statistic& operator+=(const Statistic& x)
  {
    bitext += x.bitext;
    return *this;
  }

  Statistic& operator-=(const Statistic& x)
  {
    bitext -= x.bitext;
    return *this;
  }
  
  friend
  std::ostream& operator<<(std::ostream& os, const Statistic& x)
  {
    os << x.bitext << '\n';
    return os;
  }

  friend
  std::istream& operator>>(std::istream& is, Statistic& x)
  {
    is >> x.bitext;

    if (x.bitext <= 0)
      throw std::runtime_error("invalid statistic");

    return is;
  }
};

#endif
