// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BASE64__HPP__
#define __UTILS__BASE64__HPP__ 1

#include <string>
#include <algorithm>
#include <iterator>
#include <stdexcept>

#include <utils/bithack.hpp>

#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>

namespace utils
{
  struct struct_base64_decoder_base
  {
    static inline
    unsigned char table(size_t pos)
    {
      static const char __table[] = {
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,62,-1,-1,-1,63,
	52,53,54,55,56,57,58,59,60,61,-1,-1,-1,-1,-1,-1,
	-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,
	15,16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,
	-1,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
	41,42,43,44,45,46,47,48,49,50,51,-1,-1,-1,-1,-1};
      
      if (pos > 127)
	throw std::runtime_error("invalid base64 character");

      const char value = __table[pos];
      
      if (value < 0)
	throw std::runtime_error("invalid base64 character");
      
      return value;
    }
  };

  template <size_t Size>
  struct struct_base64_decoder : public struct_base64_decoder_base
  {
    static const size_t i_max    = Size / 3;
    static const size_t i_remain = Size % 3;
    static const size_t input_last = (i_max * 4) + (i_remain == 0 ? size_t(0) : i_remain + 1);
 
    template <typename Iterator>
    static inline
    Iterator decode(unsigned char* buf, Iterator iter, Iterator last)
    {
      if (iter + input_last > last)
	throw std::runtime_error("base64 decoder: insufficient buffer");

      for (size_t i = 0; i != i_max; ++ i, buf += 3, iter += 4) {
	buf[0] = (table(*(iter + 0)) << 2) | (table(*(iter + 1)) >> 4);
	buf[1] = (table(*(iter + 1)) << 4) | (table(*(iter + 2)) >> 2);
	buf[2] = ((table(*(iter + 2)) << 6) & 0xc0) | table(*(iter + 3));
      }
      
      return struct_base64_decoder<i_remain>::decode(buf, iter, last);
    }
  };
  
  template <>
  struct struct_base64_decoder<2> : public struct_base64_decoder_base
  {
    template <typename Iterator>
    static inline
    Iterator decode(unsigned char* buf, Iterator iter, Iterator last)
    {
      if (iter + 3 > last)
	throw std::runtime_error("base64 decoder: insufficient buffer");

      buf[0] = (table(*(iter + 0)) << 2) | (table(*(iter + 1)) >> 4);
      buf[1] = (table(*(iter + 1)) << 4) | (table(*(iter + 2)) >> 2);
      
      return iter + 3;
    }
  };

  template <>
  struct struct_base64_decoder<1> : public struct_base64_decoder_base
  {
    template <typename Iterator>
    static inline
    Iterator decode(unsigned char* buf, Iterator iter, Iterator last)
    {
      if (iter + 2 > last)
	throw std::runtime_error("base64 decoder: insufficient buffer");
      
      buf[0] = (table(*(iter + 0)) << 2) | (table(*(iter + 1)) >> 4);
      return iter + 2;
    }
  };
  
  template <>
  struct struct_base64_decoder<0> : public struct_base64_decoder_base
  {
    template <typename Iterator>
    static inline
    Iterator decode(unsigned char* buf, Iterator iter, Iterator last)
    {
      return iter;
    }
  };


  struct struct_base64_encoder_base
  {
    static inline
    char table(size_t pos)
    {
      static const char* __table = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
      return __table[pos];
    }
  };
  
  template <size_t Size>
  struct struct_base64_encoder : public struct_base64_encoder_base
  {
    static const size_t i_max    = Size / 3;
    static const size_t i_remain = Size % 3;
 
    template <typename Iterator>
    static inline
    Iterator encode(const unsigned char* buf, Iterator iter)
    {
      for (size_t i = 0; i != i_max; ++ i, buf += 3) {
	*iter = table(buf[0] >> 2);                                     ++ iter;
	*iter = table(((buf[0] & 0x03) << 4) | ((buf[1] & 0xf0) >> 4)); ++ iter;
	*iter = table(((buf[1] & 0x0f) << 2) | ((buf[2] & 0xc0) >> 6)); ++ iter;
	*iter = table(buf[2] & 0x3f);                                   ++ iter;
      }
      
      return struct_base64_encoder<i_remain>::encode(buf, iter);
    }
  };

  template <>
  struct struct_base64_encoder<2> : public struct_base64_encoder_base
  {
    template <typename Iterator>
    static inline
    Iterator encode(const unsigned char* buf, Iterator iter)
    {
      *iter = table(buf[0] >> 2);                                     ++ iter;
      *iter = table(((buf[0] & 0x03) << 4) | ((buf[1] & 0xf0) >> 4)); ++ iter;
      *iter = table((buf[1] & 0x0f) << 2);                            ++ iter;

      return iter;
    }
  };

  template <>
  struct struct_base64_encoder<1> : public struct_base64_encoder_base
  {
    template <typename Iterator>
    static inline
    Iterator encode(const unsigned char* buf, Iterator iter)
    {
      *iter = table(buf[0] >> 2);          ++ iter;
      *iter = table((buf[0] & 0x03) << 4); ++ iter;

      return iter;
    }
  };

  template <>
  struct struct_base64_encoder<0> : public struct_base64_encoder_base
  {
    template <typename Iterator>
    static inline
    Iterator encode(const unsigned char* buf, Iterator iter)
    {
      return iter;
    }
  };
  
  
  template <typename Tp, typename Iterator>
  inline
  Iterator encode_base64(const Tp& x, Iterator iter)
  {
    
    return struct_base64_encoder<sizeof(Tp)>::encode(reinterpret_cast<const unsigned char*>(&x), iter);

#if 0
    static const char* enc64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    
    const unsigned char* buf = reinterpret_cast<const unsigned char*>(&x);
    const size_t size = sizeof(Tp);
    const size_t i_max = size / 3;
    
    for (size_t i = 0; i != i_max; ++ i, buf += 3) {
      *iter = enc64[buf[0] >> 2];                                     ++ iter;
      *iter = enc64[((buf[0] & 0x03) << 4) | ((buf[1] & 0xf0) >> 4)]; ++ iter;
      *iter = enc64[((buf[1] & 0x0f) << 2) | ((buf[2] & 0xc0) >> 6)]; ++ iter;
      *iter = enc64[buf[2] & 0x3f];                                   ++ iter;
    }
    
    switch (size % 3) {
    case 2:
      *iter = enc64[buf[0] >> 2];                                     ++ iter;
      *iter = enc64[((buf[0] & 0x03) << 4) | ((buf[1] & 0xf0) >> 4)]; ++ iter;
      *iter = enc64[(buf[1] & 0x0f) << 2];                            ++ iter;
      break;
    case 1:
      *iter = enc64[buf[0] >> 2];          ++ iter;
      *iter = enc64[(buf[0] & 0x03) << 4]; ++ iter;
      break;
    }
    
    return iter;
#endif
    
    
#if 0
    using namespace boost::archive::iterators;
    
    typedef base64_from_binary<transform_width<const char*, 6, 8> > encoder_type;
    
    std::copy(encoder_type((const char*) &x), encoder_type(((const char*) &x) + sizeof(Tp)), iter);
    return iter;
#endif
  }

  template <typename Tp>
  inline
  std::string encode_base64(const Tp& x)
  {
    std::string encoded;
    encode_base64(x, std::back_inserter(encoded));
    return encoded;
  }
  
  
  template <typename Tp>
  inline
  Tp decode_base64(const std::string& x)
  {
    Tp value;
    
    struct_base64_decoder<sizeof(Tp)>::decode(reinterpret_cast<unsigned char*>(&value), x.begin(), x.end());

    return value;
    
#if 0
    static const unsigned char table[] = {
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,62,-1,-1,-1,63,
      52,53,54,55,56,57,58,59,60,61,-1,-1,-1,-1,-1,-1,
      -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,
      15,16,17,18,19,20,21,22,23,24,25,-1,-1,-1,-1,-1,
      -1,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
      41,42,43,44,45,46,47,48,49,50,51,-1,-1,-1,-1,-1};
    
    Tp value;
    
    std::string::const_iterator iter = x.begin();
    std::string::const_iterator iter_end = x.end();
    char* buf = reinterpret_cast<char*>(&value);
    
    const size_t size = sizeof(Tp);
    const size_t i_max = size / 3;
    
    for (size_t i = 0; i != i_max; ++ i, buf += 3, iter += 4) {
      buf[0] = (table[*(iter + 0)] << 2) | (table[*(iter + 1)] >> 4);
      buf[1] = (table[*(iter + 1)] << 4) | (table[*(iter + 2)] >> 2);
      buf[2] = ((table[*(iter + 2)] << 6) & 0xc0) | table[*(iter + 3)];
    }
    
    switch (size % 3) {
    case 2:
      buf[0] = (table[*(iter + 0)] << 2) | (table[*(iter + 1)] >> 4);
      buf[1] = (table[*(iter + 1)] << 4) | (table[*(iter + 2)] >> 2);
      break;
    case 1:
      buf[0] = (table[*(iter + 0)] << 2) | (table[*(iter + 1)] >> 4);
      break;
    }
    
    return value;
#endif
    
#if 0
    using namespace boost::archive::iterators;
  
    typedef transform_width<binary_from_base64<std::string::const_iterator>, 8, 6> decoder_type;
  
    Tp value;
  
    char* iter = (char*) &value;
    char* iter_end = iter + sizeof(Tp);
    
    decoder_type decoder(x.begin());
    for (/**/; iter != iter_end; ++ iter, ++ decoder)
      *iter = *decoder;
    
    return value;
#endif
  }
};

#endif
