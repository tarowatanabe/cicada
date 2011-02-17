//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include "symbol_vector.hpp"

#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>


namespace cicada
{
  
  std::ostream& operator<<(std::ostream& os, const SymbolVector& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    iterator_type iter(os);
    
    if (! karma::generate(iter, -(standard::string % ' '), x.__impl))
      throw std::runtime_error("sentence generation failed...?");
    
    return os;
  }
  
  std::istream& operator>>(std::istream& is, SymbolVector& x)
  {
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
    
    std::string line;
    x.clear();
    if (std::getline(is, line)) {
      utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
      x.assign(tokenizer.begin(), tokenizer.end());
    }
    
    return is;
  }
  
};

