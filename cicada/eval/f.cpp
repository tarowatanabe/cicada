
#include <sstream>
#include <iterator>

#include "f.hpp"

#include "utils/base64.hpp"

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

    std::string F::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"" << __description() << "\",";
      stream << "\"reference\":[";
      stream << "\"" << escape_base64(utils::encode_base64(match_ref)) << "\",";
      stream << "\"" << escape_base64(utils::encode_base64(norm_ref)) << "\"";
      stream << "],";
      stream << "\"hypothesis\":[";
      stream << "\"" << escape_base64(utils::encode_base64(match_hyp)) << "\",";
      stream << "\"" << escape_base64(utils::encode_base64(norm_hyp)) << "\"";
      stream << "]";
      stream << '}';
      return stream.str();
    }
    
    
  };
  
};
