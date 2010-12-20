
#include <sstream>
#include <iterator>

#include "f.hpp"

namespace cicada
{
  namespace eval
  {
    std::string F::description() const
    {
      std::ostringstream stream;
      stream << __description() << ": " << score()
	     << ' ' << (norm_hyp != 0.0 ? match_hyp / norm_hyp : 0.0)
	     << '|' << (norm_ref != 0.0 ? match_ref / norm_ref : 0.0);
      
      return stream.str();
    }
    
  };
  
};
