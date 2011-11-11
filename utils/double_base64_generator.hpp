//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DOUBLE_BASE64_GENERATOR__HPP__
#define __UTILS__DOUBLE_BASE64_GENERATOR__HPP__ 1

#include <string>
#include <iterator>

#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <utils/base64.hpp>

namespace utils
{
  template <typename Iterator>
  struct double_base64_generator : boost::spirit::karma::grammar<Iterator, double()>
  {
    struct double_base64_func
    {
      template<class, class>
      struct result {
	typedef void type;
      };

      void operator()(std::string& x, const double& value) const
      {
	utils::encode_base64(value, std::back_inserter(x));
      }
    };
    
    double_base64_generator() : double_base64_generator::base_type(double_token)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      double_token = standard::string[double_base64(karma::_1, karma::_val)];
    }
    
    boost::phoenix::function<double_base64_func> const double_base64;
    boost::spirit::karma::rule<Iterator, double()> double_token;
  };
    
};

#endif
