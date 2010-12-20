#include <sstream>
#include <iterator>

#include "combined.hpp"

namespace cicada
{
  namespace eval
  {
    std::string Combined::description() const
    {
      std::ostringstream stream;
      stream << "combined: " << score() << ' ';
      
      stream << '{';
      if (! scores.empty()) {
	score_ptr_set_type::const_iterator siter_end = scores.end() - 1;
	for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter)
	  stream << (*siter)->description() << ", ";
	stream << scores.back()->description();
      }
      stream << '}';
      
      return stream.str();
    }
    
  };
  
};
