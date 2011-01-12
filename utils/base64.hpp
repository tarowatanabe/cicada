// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BASE64__HPP__
#define __UTILS__BASE64__HPP__ 1

#include <string>
#include <algorithm>

#include <utils/bithack.hpp>

#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>

namespace utils
{
  struct __base64_impl
  {
    
    
    
  };

  template <typename Tp>
  inline
  std::string encode_base64(const Tp& x)
  {
    static const char* enc64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    
    const unsigned char* buf = reinterpret_cast<const unsigned char*>(&x);
    const size_t size = sizeof(Tp);

    std::string encoded;
    
    size_t curr = 0;
    while (curr < size) {
      const size_t len = utils::bithack::min(size, size_t(3));
      
      encoded += enc64[buf[curr + 0] >> 2];
      encoded += enc64[((buf[curr + 0] & 0x03) << 4) | (len > 1 ? ((buf[curr + 1] & 0xf0) >> 4) : 0)];
      encoded += (len > 1 ? enc64[((buf[curr + 1] & 0x0f) << 2) | ((buf[curr + 2] & 0xc0) >> 6) ] : '=');
      encoded += (len > 2 ? enc64[ buf[curr + 2] & 0x3f ] : '=');
      curr += len;
    }

    return encoded;
    
#if 0
    using namespace boost::archive::iterators;

    typedef base64_from_binary<transform_width<const char*, 6, 8> > encoder_type;

    std::string encoded;
    std::copy(encoder_type((const char*) &x), encoder_type(((const char*) &x) + sizeof(Tp)), std::back_inserter(encoded));
  
    return encoded;
#endif
  }

  template <typename Tp, typename Iterator>
  inline
  Iterator encode_base64(const Tp& x, Iterator iter)
  {
    static const char* enc64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    
    const unsigned char* buf = reinterpret_cast<const unsigned char*>(&x);
    const size_t size = sizeof(Tp);
    
    size_t curr = 0;
    while (curr < size) {
      const size_t len = utils::bithack::min(size, size_t(3));
      
      *iter = enc64[buf[curr + 0] >> 2];
      ++ iter;
      
      *iter = enc64[((buf[curr + 0] & 0x03) << 4) | ((buf[curr + 1] & 0xf0) >> 4)];
      ++ iter;
      
      *iter = (len > 1 ? enc64[((buf[curr + 1] & 0x0f) << 2) | ((buf[curr + 2] & 0xc0) >> 6) ] : '=');
      ++ iter;
      
      *iter = (len > 2 ? enc64[ buf[curr + 2] & 0x3f ] : '=');
      ++ iter;
      
      curr += len;
    }

    return iter;
    
#if 0
    using namespace boost::archive::iterators;
    
    typedef base64_from_binary<transform_width<const char*, 6, 8> > encoder_type;
    
    std::copy(encoder_type((const char*) &x), encoder_type(((const char*) &x) + sizeof(Tp)), iter);
    return iter;
#endif
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
