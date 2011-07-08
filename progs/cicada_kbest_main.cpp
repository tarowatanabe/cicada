#include <iostream>
#include <iterator>

#include "cicada_kbest_impl.hpp"

int main(int argc, char** argv)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  parser_type parser;
  kbest_feature_type kbest;
  
  std::cin.unsetf(std::ios::skipws);
  
  iter_type iter(std::cin);
  iter_type iter_end;
  
  while (iter != iter_end) {
    boost::fusion::get<1>(kbest).clear();
    boost::fusion::get<2>(kbest).clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
      if (iter != iter_end)
	throw std::runtime_error("kbest parsing failed");

    std::copy(boost::fusion::get<1>(kbest).begin(), boost::fusion::get<1>(kbest).end(), std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << "|||";
    
    feature_parsed_set_type::const_iterator fiter_end = boost::fusion::get<2>(kbest).end();
    for (feature_parsed_set_type::const_iterator fiter = boost::fusion::get<2>(kbest).begin(); fiter != fiter_end; ++ fiter)
      std::cout << ' ' << fiter->first << "=" << fiter->second;
    std::cout << std::endl;
  }
}
