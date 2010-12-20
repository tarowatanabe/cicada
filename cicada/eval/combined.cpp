#include <sstream>
#include <iterator>

#include "combined.hpp"

#include "utils/base64.hpp"

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

    inline
    std::string escape_base64(const std::string& x)
    {
      std::string result;
      
      std::string::const_iterator iter_end = x.end();
      for (std::string::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	if (*iter == '/')
	  result += "\\/";
	else
	  result += *iter;
      
      return result;
    }
    
    
    std::string Combined::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"combined\",";
      stream << "\"score\":[";
      if (! scores.empty()) {
	for (size_t i = 0; i != scores.size() - 1; ++ i)
	  stream << (*scores[i]) << ',';
	stream << (*scores.back());
      }
      stream << "],";
      stream << "\"weight\":[";
      if (! weights.empty()) {
	for (size_t i = 0; i != weights.size() - 1; ++ i)
	  stream << "\"" << escape_base64(utils::encode_base64(weights[i])) << "\",";
	stream << "\"" << escape_base64(utils::encode_base64(weights.back())) << "\"";
      }
      stream << "]";
      stream << '}';
      
      return stream.str();
    }    
  };
  
};
