//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iterator>
#include <algorithm>
#include <string>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include "dependency.hpp"

namespace cicada
{
  
  bool Dependency::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    return qi::phrase_parse(iter, end, *(qi::lexeme[qi::int_]), standard::space, __dep);
  }

  void Dependency::assign(const utils::piece& line)
  {
    std::string::const_iterator iter(line.begin());
    std::string::const_iterator end(line.end());
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("dependency format error");
  }
  
  std::ostream& operator<<(std::ostream& os, const Dependency& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! karma::generate(iter, -(karma::int_ % ' '), x.__dep))
      throw std::runtime_error("dependency generation failed...?");
    
    return os;
  }
  
  std::istream& operator>>(std::istream& is, Dependency& x)
  {
    std::string line;
    x.clear();
    if (std::getline(is, line))
      x.assign(line);
    
    return is;
  }
  
};
