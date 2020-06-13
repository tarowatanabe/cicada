//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__PYTHON_STRING_PARSER__HPP__
#define __UTILS__PYTHON_STRING_PARSER__HPP__ 1

#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace utils
{
  template <typename Iterator>
  struct python_string_parser : boost::spirit::qi::grammar<Iterator, std::string()>
  {
    typedef uint32_t uchar_type;

    struct push_raw_func
    {
#if BOOST_PHOENIX_VERSION >= 0x3000
      template<class>
#else
      template<class, class>
#endif
      struct result {
	typedef void type;
      };
      
      void operator()(std::string& result, const uchar_type c) const
      {
	result += c;
      }
    };
    
    struct push_utf8_func
    {
#if BOOST_PHOENIX_VERSION >= 0x3000
      template<class>
#else
      template<class, class>
#endif
      struct result {
	typedef void type;
      };

      void operator()(std::string& result, const uchar_type c) const
      {
	if (c <= 0x7f)
	  result += c;
	else if (c <= 0x7ff) {
	  const char buffer[2] = {static_cast<char>(((c >> 6) & 0x1f) | 0xc0),
				  static_cast<char>((c & 0x3f) | 0x80)};
	  result.append(buffer, 2);
	} else if (c <= 0xffff) {
	  const char buffer[3] = {static_cast<char>(((c >> 12) & 0x0f) | 0xe0),
				  static_cast<char>(((c >> 6)  & 0x3f) | 0x80),
				  static_cast<char>((c & 0x3f) | 0x80)};
	  result.append(buffer, 3);
	} else if (c <= 0x1fffff) {
	  const char buffer[4] = {static_cast<char>(((c >> 18) & 0x07) | 0xf0),
				  static_cast<char>(((c >> 12) & 0x3f) | 0x80),
				  static_cast<char>(((c >> 6)  & 0x3f) | 0x80),
				  static_cast<char>((c & 0x3f) | 0x80)};
	  result.append(buffer, 4);
	} else if (c <= 0x3ffffff) {
	  const char buffer[5] = {static_cast<char>(((c >> 24) & 0x03) | 0xf8),
				  static_cast<char>(((c >> 18) & 0x3f) | 0x80),
				  static_cast<char>(((c >> 12) & 0x3f) | 0x80),
				  static_cast<char>(((c >> 6)  & 0x3f) | 0x80),
				  static_cast<char>((c & 0x3f) | 0x80)};
	  result.append(buffer, 5);
	} else {
	  const char buffer[6] = {static_cast<char>(((c >> 30) & 0x01) | 0xfc),
				  static_cast<char>(((c >> 24) & 0x3f) | 0x80),
				  static_cast<char>(((c >> 18) & 0x3f) | 0x80),
				  static_cast<char>(((c >> 12) & 0x3f) | 0x80),
				  static_cast<char>(((c >> 6)  & 0x3f) | 0x80),
				  static_cast<char>((c & 0x3f) | 0x80)};
	  result.append(buffer, 6);
	}
      }
    };
    
    struct push_escaped_func
    {
#if BOOST_PHOENIX_VERSION >= 0x3000
      template<class>
#else
      template<class, class>
#endif
       struct result {
	 typedef void type;
       };
      
      void operator()(std::string& result, const uchar_type c) const
      {
	switch (c) {
	case '\\': result += '\\'; break;
	case '\'': result += '\''; break;
	case '\"': result += '\"'; break;
	case 'a':  result += '\a'; break;
	case 'b':  result += '\b'; break;
	case 'f':  result += '\f'; break;
	case 'n':  result += '\n'; break;
	case 'r':  result += '\r'; break;
	case 't':  result += '\t'; break;
	case 'v':  result += '\v'; break;
	}
      }
    };
    
    python_string_parser() : python_string_parser::base_type(string)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      escaped = '\\' >> (('u' >> hex4)   [push_utf8(qi::_r1, qi::_1)]
			 | ('U' >> hex8) [push_utf8(qi::_r1, qi::_1)]
			 | ('x' >> hex2) [push_raw(qi::_r1, qi::_1)]
			 | standard::char_("abfnrtv\\\"'") [push_escaped(qi::_r1, qi::_1)]
			 | oct [push_raw(qi::_r1, qi::_1)]);
      
      string = (('\"' >> *(escaped(qi::_val) | (~standard::char_('\"'))[qi::_val += qi::_1]) >> '\"')
		| ('\'' >> *(escaped(qi::_val) | (~standard::char_('\''))[qi::_val += qi::_1]) >> '\''));
    }
    
    boost::phoenix::function<push_raw_func> const     push_raw;
    boost::phoenix::function<push_utf8_func> const    push_utf8;
    boost::phoenix::function<push_escaped_func> const push_escaped;
    
    boost::spirit::qi::uint_parser<uchar_type, 8, 3, 3>  oct;
    boost::spirit::qi::uint_parser<uchar_type, 16, 2, 2> hex2;
    boost::spirit::qi::uint_parser<uchar_type, 16, 4, 4> hex4;
    boost::spirit::qi::uint_parser<uchar_type, 16, 8, 8> hex8;
    
    boost::spirit::qi::rule<Iterator, void(std::string&)> escaped;
    boost::spirit::qi::rule<Iterator, std::string()>      string;
  };
};

#endif
