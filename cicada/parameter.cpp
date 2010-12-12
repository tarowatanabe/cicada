//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include "parameter.hpp"

#include <boost/spirit/include/qi.hpp>

#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

namespace cicada
{
  
  typedef std::pair<std::string, std::string> value_parsed_type;
  typedef std::vector<value_parsed_type, std::allocator<value_parsed_type> > value_parsed_set_type;
  
  typedef std::pair<std::string, value_parsed_set_type> parameter_parsed_type;

  template <typename Iterator>
  struct parameter_parser : boost::spirit::qi::grammar<Iterator, parameter_parsed_type(), boost::spirit::standard::space_type>
  {
    parameter_parser() : parameter_parser::base_type(parameter)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      escape_char.add
	("\\\"", '\"')
	("\\\\", '\\')
	("\\a", '\a')
	("\\b", '\b')
	("\\f", '\f')
	("\\n", '\n')
	("\\r", '\r')
	("\\t", '\t')
	("\\v", '\v');
      
      escaped %= '"' >> qi::lexeme[*(escape_char | ~standard::char_('"'))] >> '"';
      
      param %= qi::lexeme[+(standard::char_ - standard::space - ':')];
      key   %= qi::lexeme[+(standard::char_ - standard::space - '=')];
      value %= qi::lexeme[+(standard::char_ - standard::space - ',')];
      key_values %= ((qi::hold[escaped] | key) >> '=' >> (qi::hold[escaped] | value)) % ',';
      parameter %= (qi::hold[escaped] | param) >> -(':' >> key_values);
    }
    
    typedef boost::spirit::standard::space_type space_type;
    
    boost::spirit::qi::symbols<char, char> escape_char;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           escaped;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           param;
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           key;
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           value;
    boost::spirit::qi::rule<Iterator, value_parsed_set_type(), space_type> key_values;
    boost::spirit::qi::rule<Iterator, parameter_parsed_type(), space_type> parameter;
  };
  
  void Parameter::parse(const std::string& parameter)
  {
    typedef std::string::const_iterator iter_type;
    typedef parameter_parser<iter_type> parser_type;
    
    __attr.clear();
    __values.clear();
    
    parser_type parser;

    parameter_parsed_type parsed;

    iter_type iter     = parameter.begin();
    iter_type iter_end = parameter.end();
    
    const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::space, parsed);
    
    if (! result || iter != iter_end)
      throw std::runtime_error(std::string("parameter parsing failed: ") + parameter);

    __attr   = parsed.first;
    __values.insert(__values.end(), parsed.second.begin(), parsed.second.end());
  }
};
