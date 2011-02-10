//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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

#include "alignment.hpp"

#include "utils/lexical_cast.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Alignment::point_type,
			  (cicada::Alignment::index_type, source)
			  (cicada::Alignment::index_type, target)
			  )


namespace cicada
{
  
  bool Alignment::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    return qi::phrase_parse(iter, end, *(qi::lexeme[qi::int_ >> '-' >> qi::int_]), standard::space, __align);
  }

  void Alignment::assign(const utils::piece& line)
  {
    std::string::const_iterator iter(line.begin());
    std::string::const_iterator end(line.end());
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("alignment format error");
  }
  
  void Alignment::point_type::assign(const utils::piece& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string::const_iterator iter(x.begin());
    std::string::const_iterator end(x.end());
    
    const bool result = qi::phrase_parse(iter, end, qi::lexeme[qi::int_ >> '-' >> qi::int_], standard::space, *this);
    if (! result || iter != end)
      throw std::runtime_error("invalid point format? " + x);
  }

  std::ostream& operator<<(std::ostream& os, const Alignment::point_type& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! karma::generate(iter, karma::int_ << '-' << karma::int_, x))
      throw std::runtime_error("point generation failed...?");
    
    return os;
  }
  
  std::istream& operator>>(std::istream& is, Alignment::point_type& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string point;
    
    if (is >> point)
      x.assign(point);
    else {
      x.source = 0;
      x.target = 0;
    }
    return is;
  }
  
  
  std::ostream& operator<<(std::ostream& os, const Alignment& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! karma::generate(iter, -((karma::int_ << '-' << karma::int_) % ' '), x.__align))
      throw std::runtime_error("alignment generation failed...?");
    
    return os;
  }
  
  std::istream& operator>>(std::istream& is, Alignment& x)
  {
    std::string line;
    x.clear();
    if (std::getline(is, line))
      x.assign(line);
    
    return is;
  }
  
};
