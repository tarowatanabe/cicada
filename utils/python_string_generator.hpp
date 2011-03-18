//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__PYTHON_STRING_GENERATOR__HPP__
#define __UTILS__PYTHON_STRING_GENERATOR__HPP__ 1

#include <string>

#include <boost/spirit/include/karma.hpp>

namespace utils
{
  template <typename Iterator>
  struct python_string_generator : boost::spirit::karma::grammar<Iterator, std::string()>
  {
    typedef uint32_t uchar_type;
    
    python_string_generator() : python_string_generator::base_type(string)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      string = ('\"'
		<< *(&standard::char_('\a') << "\\a"
		     | &standard::char_('\b') << "\\b"
		     | &standard::char_('\f') << "\\f"
		     | &standard::char_('\n') << "\\n"
		     | &standard::char_('\r') << "\\r"
		     | &standard::char_('\t') << "\\t"
		     | &standard::char_('\v') << "\\v"
		     | &standard::char_('\"') << "\\\""
		     | &standard::char_('\'') << "\\'"
		     | &standard::char_('\\') << "\\\\"
		     | standard::char_)
		<< '\"');
    }
    
    boost::spirit::karma::rule<Iterator, std::string()> string;
  };
};

#endif
