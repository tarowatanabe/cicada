//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
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
    
    return qi::phrase_parse(iter, end, *(qi::int_ >> '-' >> qi::int_), standard::space, __align);
  }

  void Alignment::assign(const std::string& line)
  {
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("alignment format error");
  }
  

  std::ostream& operator<<(std::ostream& os, const Alignment::point_type& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, karma::int_ << '-' << karma::int_, x))
      throw std::runtime_error("point generation failed...?");
    
    return os;
  }
  
  std::istream& operator>>(std::istream& is, Alignment::point_type& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    std::string point;

    if (is >> point) {
      std::string::const_iterator iter = point.begin();
      std::string::const_iterator end = point.end();
      
      const bool result = qi::phrase_parse(iter, end, qi::int_ >> '-' >> qi::int_, standard::space, x);
      if (! result || iter != end)
	throw std::runtime_error("invalid point format? " + point);
    } else {
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
    
    if (! boost::spirit::karma::generate(iter, -((karma::int_ << '-' << karma::int_) % ' '), x.__align))
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
