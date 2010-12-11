#include <iostream>

#include "cicada_text_impl.hpp"

int main(int argc, char** argv)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  parser_type parser;
  id_sentence_type id_sentence;
  
  std::cin.unsetf(std::ios::skipws);
  
  iter_type iter(std::cin);
  iter_type iter_end;
  
  while (iter != iter_end) {
    id_sentence.second.clear();
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
      if (iter != iter_end)
	throw std::runtime_error("refset parsing failed");
    
    std::cout << id_sentence.first << " ||| " << id_sentence.second << std::endl;
  }
}
