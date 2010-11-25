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
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include "alignment.hpp"

#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>

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
    namespace phoenix = boost::phoenix;
    
    using qi::int_;
    using standard::space;
    
    clear();
    
    return qi::phrase_parse(iter, end,
			    *(int_ >> '-' >> int_),
			    space,
			    __align);
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
    os << x.source << '-' << x.target;
    return os;
  }
  
  std::istream& operator>>(std::istream& is, Alignment::point_type& x)
  {
    std::string point;

    if (is >> point) {
      std::string::size_type pos = point.find('-');
      if (pos == std::string::npos)
	throw std::runtime_error(std::string("invalid format? ") + point);
      
      x.source = atoi(point.substr(0, pos).c_str());
      x.target = atoi(point.substr(pos + 1).c_str());
    } else {
      x.source = 0;
      x.target = 0;
    }
    return is;
  }
  
  
  std::ostream& operator<<(std::ostream& os, const Alignment& x)
  {
    if (! x.empty()) {
      std::copy(x.begin(), x.end() - 1, std::ostream_iterator<Alignment::point_type>(os, " "));
      os << x.back();
    }
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
