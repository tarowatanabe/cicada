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

    template <typename Tp, size_t Size, typename Alloc=std::allocator<Tp> >
    class __lexical_cast_buffer
    {
    private:
      typedef std::vector<Tp, Alloc> buffer_type;
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef Tp  value_type;
      typedef Tp* pointer;
      
      typedef Tp*       iterator;
      typedef const Tp* const_iterator;

      typedef Tp&       reference;
      typedef const Tp& const_reference;
      
    public:
      __lexical_cast_buffer() : buffer_static(), buffer_dynamic(), first(buffer_static), last(buffer_static) {}
      
    private:
      __lexical_cast_buffer(const __lexical_cast_buffer&) {}
      __lexical_cast_buffer& operator=(const __lexical_cast_buffer& x) { return *this; }

    public:
      iterator begin() { return first; }
      const_iterator begin() const { return first; }
      
      iterator end() { return last; }
      const_iterator end() const { return last; }

      bool empty() const { return first == last; }
      size_type size() const { return last - first; }

      void clear()
      {
	buffer_dynamic.clear();
	first = buffer_static;
	last = buffer_static;
      }

      void push_back(const value_type& x)
      {
	const size_type __size = size();
	
	if (__size < Size) {
	  *last = x;
	  ++ last;
	} else {
	  if (__size == Size)
	    buffer_dynamic = buffer_type(first, last);
	  buffer_dynamic.push_back(x);
	  
	  first = &(*buffer_dynamic.begin());
	  last  = &(*buffer_dynamic.end());
	}
      }
      
    private:
      Tp          buffer_static[Size];
      buffer_type buffer_dynamic;

      iterator first;
      iterator last;
    };
    
    
    template <typename Source>
    std::string __lexical_cast_generate(const Source& arg)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      typedef __lexical_cast_buffer<char, 64> buffer_type;
      
      buffer_type buffer;
      
      typename __lexical_cast_generator<Source, boost::is_float<Source>::value, boost::is_signed<Source>::value>::generator_type generator;
      
      std::back_insert_iterator<buffer_type > iter(buffer);
      if (! karma::generate(iter, generator, arg))
	throw std::bad_cast();
      
      return std::string(buffer.begin(), buffer.end());

#if 0
      // a little bit dangerous, though...
      char buffer[256];
      typename __lexical_cast_generator<Source, boost::is_float<Source>::value, boost::is_signed<Source>::value>::generator_type generator;
      
      char* iter = buffer;
      if (! karma::generate(iter, generator, arg))
	throw std::bad_cast();
      *iter = '\0';
      
      return buffer;
#endif
     
#if 0 
      std::string generated;
      typename __lexical_cast_generator<Source, boost::is_float<Source>::value, boost::is_signed<Source>::value>::generator_type generator;
      
      std::back_insert_iterator<std::string> iter(generated);
      if (! karma::generate(iter, generator, arg))
	throw std::bad_cast();
      
      return generated;
#endif
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

    template <typename Target, typename Source, bool TargetArithmetic, bool SourceArithmetic>
    struct __lexical_cast
    {
      static inline
      Target cast(const Source& arg)
      {
	return boost::lexical_cast<Target>(arg);
      }
    };

    template <typename Target>
    struct __lexical_cast<Target, utils::piece, true, false>
    {
      static inline
      Target cast(const utils::piece& arg)
      {
	return __lexical_cast_parse<Target>(arg);
      }
    };
    
    template <>
    struct __lexical_cast<bool, utils::piece, true, false>
    {
      static inline
      bool cast(const utils::piece& arg)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	namespace phoenix = boost::phoenix;
    
	utils::piece::const_iterator iter = arg.begin();
	utils::piece::const_iterator iter_end = arg.end();
	
	bool parsed = false;
	const bool result = qi::phrase_parse(iter, iter_end,
					     qi::no_case["true"] [phoenix::ref(parsed) = true] 
					     || qi::no_case["yes"] [phoenix::ref(parsed) = true] 
					     || qi::no_case["no"] [phoenix::ref(parsed) = false] 
					     || qi::no_case["nil"] [phoenix::ref(parsed) = false] 
					     || qi::int_ [phoenix::ref(parsed) = (qi::_1 > 0)],
					     standard::space);
	
	return result && iter == iter_end && parsed;
      }
    };
    
    typedef const char* __lexical_cast_char_pointer;

    template <typename Target>
    struct __lexical_cast<Target, __lexical_cast_char_pointer, true, false>
    {
      static inline
      Target cast(const __lexical_cast_char_pointer& arg)
      {
	return __lexical_cast<Target, utils::piece, true, false>::cast(arg);
      }
    };

    template <typename Target>
    struct __lexical_cast<Target, std::string, true, false>
    {
      static inline
      Target cast(const std::string& arg)
      {
	return __lexical_cast<Target, utils::piece, true, false>::cast(arg);
      }
    };
    
    template <typename Source>
    struct __lexical_cast<std::string, Source, false, true>
    {
      static inline
      std::string cast(const Source& arg)
      {
	return __lexical_cast_generate<Source>(arg);
      }
    };

    template <>
    struct __lexical_cast<std::string, bool, false, true>
    {
      static inline
      std::string cast(const bool& arg)
      {
	return (arg ? "true" : "false");
      }
    };
    
  };
  
  template <typename Target, typename Source>
  inline
  Target lexical_cast(const Source& arg)
  {
    typedef typename impl::__lexical_cast_array_to_pointer_decay<Source>::type src;
    
    return impl::__lexical_cast<Target, src, boost::is_arithmetic<Target>::value, boost::is_arithmetic<src>::value>::cast(arg);
  }

};

#endif
