// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// TODO: better implementation by differentiating via type-traits...
//

#ifndef __UTILS__LEXICAL_CAST__HPP__
#define __UTILS__LEXICAL_CAST__HPP__ 1

#include <stdexcept>
#include <iterator>

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits.hpp>

#include <utils/piece.hpp>

namespace utils
{

  namespace impl
  {
    template <typename Tp, bool isFloat, bool isSigned>
    struct __lexical_cast_parser {};

    template <typename Tp>
    struct __lexical_cast_parser<Tp, true, false>
    {
      typedef boost::spirit::qi::real_parser<Tp> parser_type;
    };
    
    template <typename Tp>
    struct __lexical_cast_parser<Tp, false, true>
    {
      typedef boost::spirit::qi::int_parser<Tp> parser_type;
    };
    
    template <typename Tp>
    struct __lexical_cast_parser<Tp, false, false>
    {
      typedef boost::spirit::qi::uint_parser<Tp> parser_type;
    };
    
    
    template <typename Target>
    Target __lexical_cast_parse(const utils::piece& arg)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
    
      utils::piece::const_iterator iter = arg.begin();
      utils::piece::const_iterator iter_end = arg.end();
      
      Target parsed;
      typename __lexical_cast_parser<Target, boost::is_float<Target>::value, boost::is_signed<Target>::value>::parser_type parser;
    
      if (! qi::phrase_parse(iter, iter_end, parser, standard::space, parsed) || iter != iter_end)
	throw std::bad_cast();
      
      return parsed;
    }

    template <typename Tp, bool isFloat, bool isSigned>
    struct __lexical_cast_generator {};
    
    template <typename Tp>
    struct __lexical_cast_generator<Tp, true, false>
    {
      struct __policy : boost::spirit::karma::real_policies<Tp>
      {
	static unsigned int precision(Tp)
	{
	  return std::numeric_limits<Tp>::digits10 + 1;
	}
      };

      typedef boost::spirit::karma::real_generator<Tp, __policy> generator_type;
    };
    
    template <typename Tp>
    struct __lexical_cast_generator<Tp, false, true>
    {
      typedef boost::spirit::karma::int_generator<Tp> generator_type;
    };
    
    template <typename Tp>
    struct __lexical_cast_generator<Tp, false, false>
    {
      typedef boost::spirit::karma::uint_generator<Tp> generator_type;
    };
    
    template <typename Source>
    std::string __lexical_cast_generate(const Source& arg)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      std::string generated;
      typename __lexical_cast_generator<Source, boost::is_float<Source>::value, boost::is_signed<Source>::value>::generator_type generator;
      
      std::back_insert_iterator<std::string> iter(generated);

      if (! karma::generate(iter, generator, arg))
	throw std::bad_cast();
      
      return generated;
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
    
    template <>
    inline
    int64_t __lexical_cast<int64_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<int64_t>(arg);
    }
  
    template <>
    inline
    int32_t __lexical_cast<int32_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<int32_t>(arg);
    }

    template <>
    inline
    int16_t __lexical_cast<int16_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<int16_t>(arg);
    }

    template <>
    inline
    int8_t __lexical_cast<int8_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<int8_t>(arg);
    }

    template <>
    inline
    uint64_t __lexical_cast<uint64_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<uint64_t>(arg);
    }
  
    template <>
    inline
    uint32_t __lexical_cast<uint32_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<uint32_t>(arg);
    }

    template <>
    inline
    uint16_t __lexical_cast<uint16_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<uint16_t>(arg);
    }

    template <>
    inline
    uint8_t __lexical_cast<uint8_t, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<uint8_t>(arg);
    }
  

    template <>
    inline
    float __lexical_cast<float, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<float>(arg);
    }

    template <>
    inline
    double __lexical_cast<double, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<double>(arg);
    }

    template <>
    inline
    long double __lexical_cast<long double, utils::piece>(const utils::piece& arg)
    {
      return __lexical_cast_parse<long double>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, int64_t>(const int64_t& arg)
    {
      return __lexical_cast_generate<int64_t>(arg);
    }
  
    template <>
    inline
    std::string __lexical_cast<std::string, int32_t>(const int32_t& arg)
    {
      return __lexical_cast_generate<int32_t>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, int16_t>(const int16_t& arg)
    {
      return __lexical_cast_generate<int16_t>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, int8_t>(const int8_t& arg)
    {
      return __lexical_cast_generate<int8_t>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, uint64_t>(const uint64_t& arg)
    {
      return __lexical_cast_generate<uint64_t>(arg);
    }
  
    template <>
    inline
    std::string __lexical_cast<std::string, uint32_t>(const uint32_t& arg)
    {
      return __lexical_cast_generate<uint32_t>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, uint16_t>(const uint16_t& arg)
    {
      return __lexical_cast_generate<uint16_t>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, uint8_t>(const uint8_t& arg)
    {
      return __lexical_cast_generate<uint8_t>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, float>(const float& arg)
    {
      return __lexical_cast_generate<float>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, double>(const double& arg)
    {
      return __lexical_cast_generate<double>(arg);
    }

    template <>
    inline
    std::string __lexical_cast<std::string, long double>(const long double& arg)
    {
      return __lexical_cast_generate<long double>(arg);
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
    std::string __lexical_cast<std::string, bool>(const bool& arg)
    {
      return (arg ? "true" : "false");
    }
  
#define __LEXICAL_CAST_PARSE__(trg, src)		\
    template <>						\
    inline						\
    trg __lexical_cast<trg, src>(const src& arg)	\
    {							\
      return __lexical_cast<trg, utils::piece>(arg);	\
    } 

  
    typedef const char* __lexical_cast_char_pointer;

    __LEXICAL_CAST_PARSE__(bool, std::string)
    __LEXICAL_CAST_PARSE__(bool, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(int64_t, std::string)
    __LEXICAL_CAST_PARSE__(int64_t, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(int32_t, std::string)
    __LEXICAL_CAST_PARSE__(int32_t, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(int16_t, std::string)
    __LEXICAL_CAST_PARSE__(int16_t, __lexical_cast_char_pointer)
  
    __LEXICAL_CAST_PARSE__(int8_t, std::string)
    __LEXICAL_CAST_PARSE__(int8_t, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(uint64_t, std::string)
    __LEXICAL_CAST_PARSE__(uint64_t, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(uint32_t, std::string)
    __LEXICAL_CAST_PARSE__(uint32_t, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(uint16_t, std::string)
    __LEXICAL_CAST_PARSE__(uint16_t, __lexical_cast_char_pointer)
  
    __LEXICAL_CAST_PARSE__(uint8_t, std::string)
    __LEXICAL_CAST_PARSE__(uint8_t, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(float, std::string)
    __LEXICAL_CAST_PARSE__(float, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(double, std::string)
    __LEXICAL_CAST_PARSE__(double, __lexical_cast_char_pointer)

    __LEXICAL_CAST_PARSE__(long double, std::string)
    __LEXICAL_CAST_PARSE__(long double, __lexical_cast_char_pointer)
    
  };
  
  template <typename Target, typename Source>
  inline
  Target lexical_cast(const Source& arg)
  {
    typedef typename impl::__lexical_cast_array_to_pointer_decay<Source>::type src;
    
    return impl::__lexical_cast<Target, src >(arg);
  }

};

#endif
