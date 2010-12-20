
#include <sstream>
#include <iterator>

#include "bleu.hpp"

#include "utils/base64.hpp"

namespace cicada
{
  namespace eval
  {
    std::string Bleu::description() const
    {
      std::vector<double, std::allocator<double> > precisions;
      
      const double penalty = std::min(1.0 - length_reference / length_hypothesis, 0.0);
      
      double smooth = 0.5;
      double score = 0.0;
      
      for (size_t n = 0; n < ngrams_hypothesis.size(); ++ n) {
	const double p = (ngrams_reference[n] > 0
			  ? (ngrams_hypothesis[n] > 0 ? ngrams_hypothesis[n] : smooth) / ngrams_reference[n]
			  : 0.0);
	
	precisions.push_back(p);
	
	score += p > 0.0 ? std::log(p) : 0.0;
	smooth *= 0.5;
      }
      
      score /= ngrams_hypothesis.size();
      score += penalty;
      
      std::ostringstream stream;
      stream << "bleu: " << std::exp(score) << ' ';
      if (! precisions.empty()) {
	std::copy(precisions.begin(), precisions.end() - 1, std::ostream_iterator<double>(stream, "|"));
	stream << precisions.back();
      }
      stream << " penalty: " << std::exp(penalty);
      
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
    
    std::string Bleu::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"bleu\",";
      stream << "\"reference\":[";
      stream << "\"" << escape_base64(utils::encode_base64(length_reference)) << "\"";
      for (size_t i = 0; i != ngrams_reference.size(); ++ i)
	stream << ",\"" << escape_base64(utils::encode_base64(ngrams_reference[i])) << "\"";
      stream << "],";
      stream << "\"hypothesis\":[";
      stream << "\"" << escape_base64(utils::encode_base64(length_hypothesis)) << "\"";
      for (size_t i = 0; i != ngrams_hypothesis.size(); ++ i)
	stream << ",\"" << escape_base64(utils::encode_base64(ngrams_hypothesis[i])) << "\"";
      stream << "]";
      stream << '}';
      
      return stream.str();
    }
  };
};
