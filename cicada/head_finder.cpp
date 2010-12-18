
#include "head_finder.hpp"
#include "head/chinese.hpp"
#include "head/collins.hpp"

namespace cicada
{
  
  const char* HeadFinder::lists()
  {
    static const char* desc = "\
collins: Collins head finder\n\
chinese: Chinese head finder\n\
";
    return desc;
  }

  HeadFinder& HeadFinder::create(const std::string& parameter)
  {
    
  }
  
};
