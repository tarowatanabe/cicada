//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "symbol_vector.hpp"

#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>


namespace cicada
{
  
  std::ostream& operator<<(std::ostream& os, const SymbolVector& x)
  {
    if (! x.empty()) {
      std::copy(x.begin(), x.end() - 1, std::ostream_iterator<SymbolVector::symbol_type>(os, " "));
      os << x.back();
    }
    return os;
  }
  
  std::istream& operator>>(std::istream& is, SymbolVector& x)
  {
    typedef boost::tokenizer<utils::space_separator> tokenizer_type;
    
    std::string line;
    x.clear();
    if (std::getline(is, line)) {
      tokenizer_type tokenizer(line);
      x.assign(tokenizer.begin(), tokenizer.end());
    }
    
    return is;
  }
  
};

