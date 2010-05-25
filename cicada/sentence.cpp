#include <iterator>

#include "sentence.hpp"

#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>

namespace cicada
{
  
  void Sentence::assign(const std::string& x)
  {
    typedef boost::tokenizer<utils::space_separator> tokenizer_type;
    
    tokenizer_type tokenizer(x);
    __sent.assign(tokenizer.begin(), tokenizer.end());
  }

  std::ostream& operator<<(std::ostream& os, const Sentence& x)
  {
    if (! x.empty()) {
      std::copy(x.__sent.begin(), x.__sent.end() - 1, std::ostream_iterator<Sentence::word_type>(os, " "));
      os << x.__sent.back();
    }
    return os;
  }

  std::istream& operator>>(std::istream& is, Sentence& x)
  {
    typedef boost::tokenizer<utils::space_separator> tokenizer_type;
    
    std::string line;
    
    x.clear();
    if (std::getline(is, line)) {
      tokenizer_type tokenizer(line);
      x.__sent.assign(tokenizer.begin(), tokenizer.end());
    }
    
    return is;
  }
};
