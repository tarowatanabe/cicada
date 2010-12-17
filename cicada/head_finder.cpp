
#include "head_finder.hpp"
#include "head/collins.hpp"

namespace cicada
{
  
  const char* HeadFinder::lists()
  {
    static const char* desc = "\
collins: collins head finder\n\
";
    return desc;
  }

  HeadFinder& HeadFinder::create(const std::string& parameter)
  {
    
  }
  
};
