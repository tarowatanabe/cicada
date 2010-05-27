
#include "parameter.hpp"

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>


namespace cicada
{
  void Parameter::parse(const std::string& parameter)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
	
    using qi::phrase_parse;
    using qi::lexeme;
    using standard::char_;
    using qi::_1;
    using standard::space;
    
    using phoenix::push_back;
    using phoenix::ref;
    
    typedef std::string::const_iterator iter_type;
    
    qi::rule<iter_type, attribute_type(), standard::space_type> rule_param     = (lexeme[+(char_ - space - ',')]);
    qi::rule<iter_type, key_type(), standard::space_type>       rule_key       = (lexeme[+(char_ - space - ',' - '=')]);
    qi::rule<iter_type, value_type(), standard::space_type>     rule_key_value = (rule_key >> '=' >> rule_key);
    
    __attr.clear();
    __values.clear();

    std::string::const_iterator iter = parameter.begin();
    std::string::const_iterator iter_end = parameter.end();
    
    const bool result = phrase_parse(iter, iter_end,
				     (
				      rule_param[ref(__attr) = _1]
				      >> *(',' >> rule_key_value[push_back(ref(__values), _1)])
				      ),
				     space);
    if (! result || iter != iter_end)
      throw std::runtime_error(std::string("parameter parsing failed: ") + parameter);
  }
};
