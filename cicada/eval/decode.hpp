// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__DECODE__HPP__
#define __CICADA__EVAL__DECODE__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <string>

#include <boost/spirit/include/qi.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <utils/base64.hpp>

namespace cicada
{
  namespace eval
  {    
    template <typename Iterator>
    struct double_base64_parser : boost::spirit::qi::grammar<Iterator, double(), boost::spirit::standard::space_type>
    {
      class double_base64_type : public std::string
      {
      public:
	operator double() const { return utils::decode_base64<double>(static_cast<const std::string&>(*this)); }
      };
      
      double_base64_parser() : double_base64_parser::base_type(decoded)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	escape_char.add
	  ("\\\"", '\"')
	  ("\\\\", '\\')
	  ("\\/", '/')
	  ("\\b", '\b')
	  ("\\f", '\f')
	  ("\\n", '\n')
	  ("\\r", '\r')
	  ("\\t", '\t');
	
	encoded %= '\"' >> qi::lexeme[+(escape_char | ~standard::char_('\"'))] >> '\"';
	decoded %= encoded;
      }
      
      typedef boost::spirit::standard::space_type space_type;

      boost::spirit::qi::symbols<char, char> escape_char;
      boost::spirit::qi::rule<Iterator, double_base64_type(), space_type> encoded;
      boost::spirit::qi::rule<Iterator, double(), space_type> decoded;
    };
    
  };
};

#endif
