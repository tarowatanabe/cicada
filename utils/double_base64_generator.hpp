//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__DOUBLE_BASE64_GENERATOR__HPP__
#define __UTILS__DOUBLE_BASE64_GENERATOR__HPP__ 1

#include <string>

#include <boost/spirit/include/karma.hpp>

#include <utils/base64.hpp>

namespace utils
{
  template <typename Iterator>
  struct double_base64_geneartor : boost::spirit::karma::grammar<Iterator, double()>
  {
    class double_base64_type : public std::string
    {
    public:
      double_base64_type(const double& x) : std::string(utils::encode_base64(x)) {}
    };
    
    
  };
    
};

#endif
