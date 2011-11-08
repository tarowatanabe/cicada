//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DOUBLE_BASE64_PARSER__HPP__
#define __UTILS__DOUBLE_BASE64_PARSER__HPP__ 1

#include <string>

#include <boost/spirit/include/qi.hpp>

#include <utils/base64.hpp>

namespace utils
{
  template <typename Iterator>
  struct double_base64_parser : boost::spirit::qi::grammar<Iterator, double()>
  {
    class double_base64_type : public std::string
    {
    public:
      operator double() const { return utils::decode_base64<double>(static_cast<const std::string&>(*this)); }
    };
    
    double_base64_parser() : double_base64_parser::base_type(double_base64)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      double_token %= qi::repeat(11)[standar::char_ - standard::space];
      double_base64 %= double_token;
    }
    
    boost::spirit::qi::rule<Iterator, double_base64_type()> double_token;
    boost::spirit::qi::rule<Iterator, double()> double_base64;
  };
  
};

#endif
