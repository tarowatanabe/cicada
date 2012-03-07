#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>

#include "utils/double_base64_parser.hpp"
#include "utils/double_base64_generator.hpp"

int main(int argc, char** argv)
{
  namespace qi = boost::spirit::qi;
  namespace karma = boost::spirit::karma;
  namespace standard = boost::spirit::standard;

  typedef std::ostream_iterator<char> oiterator_type;
  
  utils::double_base64_parser<std::string::const_iterator> parser;
  utils::double_base64_generator<oiterator_type> generator;

  srandom(time(0));
  
  for (int i = 0; i != 1024 * 16; ++ i) {
    const double value = double(random()) / double(random());

    std::ostringstream encoded_stream;
    
    karma::generate(oiterator_type(encoded_stream), generator, value);
    
    const std::string encoded = encoded_stream.str();
    
    std::string::const_iterator iter = encoded.begin();
    std::string::const_iterator end = encoded.end();
    
    double decoded;
    qi::parse(iter, end, parser, decoded);
    
    if (decoded != value)
      std::cerr << "value differ" << std::endl;
  }
}
