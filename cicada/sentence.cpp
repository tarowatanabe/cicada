//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <iterator>

#include "sentence.hpp"

namespace cicada
{
  bool Sentence::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    typedef std::string::const_iterator iter_type;
    
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    qi::rule<iter_type, std::string(), standard::space_type> word = qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    
    return qi::phrase_parse(iter, end, *(word), standard::space, __sent);
  }
  
  void Sentence::assign(const utils::piece& x)
  {
    std::string::const_iterator iter(x.begin());
    std::string::const_iterator end(x.end());
    
    const bool result = assign(iter, end);
    
    if (! result || iter != end)
      throw std::runtime_error("sentence parsing failed");
  }

  std::ostream& operator<<(std::ostream& os, const Sentence& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    iterator_type iter(os);
    
    if (! karma::generate(iter, -(standard::string % ' '), x.__sent))
      throw std::runtime_error("sentence generation failed...?");
    
    return os;
  }

  std::istream& operator>>(std::istream& is, Sentence& x)
  {
    std::string line;
    
    x.clear();
    if (std::getline(is, line))
      x.assign(line);
    
    return is;
  }
};
