
#include "head_finder.hpp"
#include "head/collins.hpp"

namespace cicada
{
  
  const char* HeadFinder::lists()
  {
    static const char* desc = "\
lower: matching by lower-case\n\
stemmer: matching by stemming algorithm\n\
\talgorithm=[stemmer spec]\n\
wordnet: matching by wordnet synsets\n\
\tpath=[path to wordnet database]\n   \
";
    return desc;
  }

  HeadFinder& HeadFinder::create(const std::string& parameter)
  {
    
  }
  
};
