// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BASE64__HPP__
#define __UTILS__BASE64__HPP__ 1

#include <string>
#include <algorithm>

#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>

namespace utils
{
  template <typename Tp>
  inline
  std::string encode_base64(const Tp& x)
  {
    using namespace boost::archive::iterators;

    typedef base64_from_binary<transform_width<const char*, 6, 8> > encoder_type;

    std::string encoded;
    std::copy(encoder_type((const char*) &x), encoder_type(((const char*) &x) + sizeof(x)), std::back_inserter(encoded));
  
    return encoded;
  }

  template <typename Tp>
  inline
  Tp decode_base64(const std::string& x)
  {
    using namespace boost::archive::iterators;
  
    typedef transform_width<binary_from_base64<std::string::const_iterator>, 8, 6> decoder_type;
  
    Tp value;
  
    char* iter = (char*) &value;
    char* iter_end = iter + sizeof(Tp);

    decoder_type decoder(x.begin());
    for (/**/; iter != iter_end; ++ iter, ++ decoder)
      *iter = *decoder;
  
    return value;
  }
};

#endif
