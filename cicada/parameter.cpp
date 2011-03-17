//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iterator>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include "parameter.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include "utils/c_string_parser.hpp"
#include "utils/c_string_generator.hpp"

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
      
      param %= qi::lexeme[+(standard::char_ - standard::space - ':')];
      key   %= qi::lexeme[+(standard::char_ - standard::space - '=')];
      value %= qi::lexeme[+(standard::char_ - standard::space - ',')];
      key_values %= ((qi::hold[escaped] | key) >> '=' >> (qi::hold[escaped] | value)) % ',';
      parameter %= (qi::hold[escaped] | param) >> -(':' >> key_values);
    }
    
    typedef boost::spirit::standard::space_type space_type;
    
    utils::c_string_parser<Iterator> escaped;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           param;
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           key;
    boost::spirit::qi::rule<Iterator, std::string(), space_type>           value;
    boost::spirit::qi::rule<Iterator, value_parsed_set_type(), space_type> key_values;
    boost::spirit::qi::rule<Iterator, parameter_parsed_type(), space_type> parameter;
  };
  
  void Parameter::parse(const utils::piece& parameter)
  {
    typedef utils::piece::const_iterator iter_type;
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
  
  std::ostream& operator<<(std::ostream& os, const Parameter& x)
  {
    typedef std::ostream_iterator<char> iterator_type;

    utils::c_string_generator<iterator_type> escaped;
    boost::spirit::karma::rule<iterator_type, value_parsed_type()> value;
    
    value %= escaped << '=' << escaped;
    
    iterator_type iter(os);
    boost::spirit::karma::generate(iter, escaped, x.__attr);
    if (! x.__values.empty()) {
      os << ':';
      boost::spirit::karma::generate(iter, value % ',', x.__values);
    }
    
    return os;
  }
};
