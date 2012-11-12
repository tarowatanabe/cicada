

#include "cicada_extract_score_impl.hpp"

void dump(std::ostream& os, const RootCountParser::root_count_type& root_count)
{
  os << "label: " << root_count.label << std::endl;
  os << "counts: ";
  std::copy(root_count.counts.begin(), root_count.counts.end(), std::ostream_iterator<double>(os, " "));
  os << std::endl;
  os << "observed: " << root_count.observed << std::endl;
}

int main(int argc, char** argv)
{
  
  RootCountParser root_parser;
  
  RootCountParser::root_count_type root_count;
  
  if (! root_parser("Good ||| 5 ||| 5", root_count))
    std::cout << "parsing failed" << std::endl;
  dump(std::cout, root_count);
  
  if (! root_parser("|||| ||| 5 ||| 5", root_count))
    std::cout << "parsing failed" << std::endl;
  dump(std::cout, root_count);

  if (! root_parser("||| ||| 5 ||| 5", root_count))
    std::cout << "parsing failed" << std::endl;
  dump(std::cout, root_count);

  if (! root_parser("||| ||| 5 ||| 5", root_count))
    std::cout << "parsing failed" << std::endl;
  dump(std::cout, root_count);

  PhrasePairParser phrase_pair_parser;
  PhrasePairParser::phrase_pair_type phrase_pair;
  
  if (! phrase_pair_parser("Good morning |||| bad ||| good ||| 0-0 ||| 5", phrase_pair))
    std::cout << "parsing failed" << std::endl;
  
  PhrasePairGenerator()(std::cout, phrase_pair) << std::endl;
}
