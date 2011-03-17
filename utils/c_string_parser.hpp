//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__C_STRING_PARSER__HPP__
#define __UTILS__C_STRING_PARSER__HPP__ 1

#include <string>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace utils
{
  template <typename Iterator>
  struct c_string_parser : boost::spirit::qi::grammar<Iterator, std::string()>
  {
    typedef uint32_t uchar_type;
    
    struct push_raw_func
    {
      template<class, class>
      struct result {
	typedef void type;
      };
      
      void operator()(std::string& result, const uchar_type c) const
      {
	result += c;
      }
    };
    
    struct push_escaped_func
    {
       template<class, class>
       struct result {
	 typedef void type;
       };
      
      void operator()(std::string& result, const uchar_type c) const
      {
	switch (c) {
	case '\\': result += '\\'; break;
	case '\'': result += '\''; break;
	case '\"': result += '\"'; break;
	case '\?': result += '\?'; break;
	case 'a':  result += '\a'; break;
	case 'b':  result += '\b'; break;
	case 'f':  result += '\f'; break;
	case 'n':  result += '\n'; break;
	case 'r':  result += '\r'; break;
	case 't':  result += '\t'; break;
	case 'v':  result += '\v'; break;
	case '0':  result += '\0'; break;
	}
      }
    };
    
    c_string_parser() : c_string_parser::base_type(string)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      escaped = '\\' >> (('x' >> hex2) [push_raw(qi::_r1, qi::_1)]
			 | standard::char_("abfnrtv0\\\"'?") [push_escaped(qi::_r1, qi::_1)]
			 | oct [push_raw(qi::_r1, qi::_1)]);
      
      string = '\"' >> *(escaped(qi::_val) | (~standard::char_('\"'))[qi::_val += qi::_1]) >> '\"';
    }

    boost::phoenix::function<push_raw_func> const     push_raw;
    boost::phoenix::function<push_escaped_func> const push_escaped;
    
    boost::spirit::qi::uint_parser<uchar_type, 8, 3, 3>  oct;
    boost::spirit::qi::uint_parser<uchar_type, 16, 2, 2> hex2;
    
    boost::spirit::qi::rule<Iterator, void(std::string&)> escaped;
    boost::spirit::qi::rule<Iterator, std::string()>      string;
  };
};

#endif
