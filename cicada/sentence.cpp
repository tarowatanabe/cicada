//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <iterator>

#include "sentence.hpp"

#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>

namespace cicada
{

  bool Sentence::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    typedef std::string::const_iterator iter_type;

    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    qi::rule<iter_type, std::string(), standard::space_type> word = qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    
    clear();
    
    return qi::phrase_parse(iter, end, *(word[phoenix::push_back(phoenix::ref(__sent), qi::_1)]), standard::space);
  }
  
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
