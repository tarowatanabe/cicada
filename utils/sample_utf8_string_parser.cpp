#include "utils/utf8_string_parser.hpp"

#include <iostream>

int main(int argc, char** argv)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
    
  utils::utf8_string_parser<std::string::const_iterator> parser;
  
  std::string line;
  std::string parsed;
  while (std::getline(std::cin, line)) {
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    parsed.clear();
    if (qi::phrase_parse(iter, end, parser, standard::space, parsed))
      std::cout << "parsed: " << parsed << std::endl;
    else
      std::cout << "failed: " << std::string(iter, end) << std::endl;
  }
  
}
