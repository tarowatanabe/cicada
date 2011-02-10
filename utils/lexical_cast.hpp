// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__LEXICAL_CAST__HPP__
#define __UTILS__LEXICAL_CAST__HPP__ 1

#include <stdexcept>

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/lexical_cast.hpp>

#include <utils/piece.hpp>

namespace utils
{

  template <typename Target>
  Target __lexical_cast_unsigned(const utils::piece& arg)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    utils::piece::const_iterator iter = arg.begin();
    utils::piece::const_iterator iter_end = arg.end();
    
    Target parsed;
    
    qi::uint_parser<Target> parser;
    
    if (! qi::phrase_parse(iter, iter_end, parser, standard::space, parsed) || iter != iter_end)
      throw std::bad_cast();
    
    return parsed;
  }

  template <typename Target>
  Target __lexical_cast_signed(const utils::piece& arg)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    utils::piece::const_iterator iter = arg.begin();
    utils::piece::const_iterator iter_end = arg.end();
    
    Target parsed;
    
    qi::int_parser<Target> parser;
    
    if (! qi::phrase_parse(iter, iter_end, parser, standard::space, parsed) || iter != iter_end)
      throw std::bad_cast();
    
    return parsed;
  }

  template <typename Target>
  Target __lexical_cast_real(const utils::piece& arg)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    utils::piece::const_iterator iter = arg.begin();
    utils::piece::const_iterator iter_end = arg.end();
    
    Target parsed;
    
    qi::real_parser<Target> parser;
    
    if (! qi::phrase_parse(iter, iter_end, parser, standard::space, parsed) || iter != iter_end)
      throw std::bad_cast();
    
    return parsed;
  }

  template<class T>
  struct __lexical_cast_array_to_pointer_decay
  {
    typedef T type;
  };
  
  template<class T, std::size_t N>
  struct __lexical_cast_array_to_pointer_decay<T[N]>
  {
    typedef const T * type;
  };
  
  
  template <typename Target, typename Source>
  inline
  Target __lexical_cast(const Source& arg)
  {
    return boost::lexical_cast<Target>(arg);
  }
  
  template <typename Target, typename Source>
  inline
  Target lexical_cast(const Source& arg)
  {
    typedef typename __lexical_cast_array_to_pointer_decay<Source>::type src;
    
    return __lexical_cast<Target, src>(arg);
  }
  
  template <>
  inline
  int64_t __lexical_cast<int64_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_signed<int64_t>(arg);
  }
  
  template <>
  inline
  int32_t __lexical_cast<int32_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_signed<int32_t>(arg);
  }

  template <>
  inline
  int16_t __lexical_cast<int16_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_signed<int16_t>(arg);
  }

  template <>
  inline
  int8_t __lexical_cast<int8_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_signed<int8_t>(arg);
  }

  template <>
  inline
  uint64_t __lexical_cast<uint64_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_unsigned<uint64_t>(arg);
  }
  
  template <>
  inline
  uint32_t __lexical_cast<uint32_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_unsigned<uint32_t>(arg);
  }

  template <>
  inline
  uint16_t __lexical_cast<uint16_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_unsigned<uint16_t>(arg);
  }

  template <>
  inline
  uint8_t __lexical_cast<uint8_t, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_unsigned<uint8_t>(arg);
  }
  

  template <>
  inline
  float __lexical_cast<float, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_real<float>(arg);
  }

  template <>
  inline
  double __lexical_cast<double, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_real<double>(arg);
  }

  template <>
  inline
  long double __lexical_cast<long double, utils::piece>(const utils::piece& arg)
  {
    return __lexical_cast_real<long double>(arg);
  }
  
  template <>
  inline
  bool __lexical_cast<bool, utils::piece>(const utils::piece& arg)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    utils::piece::const_iterator iter = arg.begin();
    utils::piece::const_iterator iter_end = arg.end();
    
    int parsed = 0;
    const bool result = qi::phrase_parse(iter, iter_end,
					 qi::no_case["true"] [phoenix::ref(parsed) = 1] 
					 || qi::no_case["yes"] [phoenix::ref(parsed) = 1] 
					 || qi::no_case["no"] [phoenix::ref(parsed) = 0] 
					 || qi::no_case["nil"] [phoenix::ref(parsed) = 0] 
					 || qi::int_ [phoenix::ref(parsed) = qi::_1],
					 standard::space);
    
    return result && iter == iter_end && parsed > 0;
  }
  
  template <>
  inline
  utils::piece __lexical_cast<utils::piece, bool>(const bool& arg)
  {
    return (arg ? "true" : "false");
  }
  
  template <>
  inline
  std::string __lexical_cast<std::string, bool>(const bool& arg)
  {
    return (arg ? "true" : "false");
  }
  
#define __LEXICAL_CAST_DEFINE__(trg, src) 	   \
  template <>                                      \
  inline					   \
  trg __lexical_cast<trg, src>(const src& arg)	   \
  {						   \
    return __lexical_cast<trg, utils::piece>(arg);	\
  } 
  
  
  typedef const char* __lexical_cast_char_pointer;

  __LEXICAL_CAST_DEFINE__(bool, std::string)
  __LEXICAL_CAST_DEFINE__(bool, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(int64_t, std::string)
  __LEXICAL_CAST_DEFINE__(int64_t, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(int32_t, std::string)
  __LEXICAL_CAST_DEFINE__(int32_t, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(int16_t, std::string)
  __LEXICAL_CAST_DEFINE__(int16_t, __lexical_cast_char_pointer)
  
  __LEXICAL_CAST_DEFINE__(int8_t, std::string)
  __LEXICAL_CAST_DEFINE__(int8_t, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(uint64_t, std::string)
  __LEXICAL_CAST_DEFINE__(uint64_t, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(uint32_t, std::string)
  __LEXICAL_CAST_DEFINE__(uint32_t, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(uint16_t, std::string)
  __LEXICAL_CAST_DEFINE__(uint16_t, __lexical_cast_char_pointer)
  
  __LEXICAL_CAST_DEFINE__(uint8_t, std::string)
  __LEXICAL_CAST_DEFINE__(uint8_t, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(float, std::string)
  __LEXICAL_CAST_DEFINE__(float, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(double, std::string)
  __LEXICAL_CAST_DEFINE__(double, __lexical_cast_char_pointer)

  __LEXICAL_CAST_DEFINE__(long double, std::string)
  __LEXICAL_CAST_DEFINE__(long double, __lexical_cast_char_pointer)

};

#endif
