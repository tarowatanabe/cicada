//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//
 
#include <iterator>
#include <algorithm>
#include <string>
#include <stdexcept>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/thread_specific_ptr.hpp"

#include "span_vector.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::SpanVector::Span,
			  (cicada::SpanVector::index_type, first)
			  (cicada::SpanVector::index_type, last)
			  (cicada::SpanVector::label_type, label)
			  )

namespace cicada
{

  template <typename Iterator, typename Container>
  struct span_vector_parser : boost::spirit::qi::grammar<Iterator, Container()>
  {
    span_vector_parser() : span_vector_parser::base_type(spans)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      label %= qi::lexeme[+(standard::char_ - standard::space)];
      span  %= qi::int_ >> '-' >> qi::int_ >> -(':' >> label);
      spans %= qi::omit[*standard::space] >> -(span % (+standard::space)) >> qi::omit[*standard::space];
    }
    
    boost::spirit::qi::rule<Iterator, std::string()> label;
    boost::spirit::qi::rule<Iterator, typename Container::value_type()> span;
    boost::spirit::qi::rule<Iterator, Container()> spans;
  };

  namespace span_vector_impl
  {
    typedef std::vector<SpanVector::span_type, std::allocator<SpanVector::span_type> > spans_type;
    typedef span_vector_parser<std::string::const_iterator, spans_type> grammar_type;
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif

    static grammar_type& instance()
    {
#ifdef HAVE_TLS
      if (! __grammar_tls) {
	__grammar.reset(new grammar_type());
	__grammar_tls = __grammar.get();
      }
      
      return *__grammar_tls;
#else
      if (! __grammar.get())
	__grammar.reset(new grammar_type());
      
      return *__grammar;
#endif
    }
  };
  
  bool SpanVector::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    return qi::parse(iter, end, span_vector_impl::instance(), __spans);
  }
  
  void SpanVector::assign(const utils::piece& line)
  {
    clear();
    
    std::string::const_iterator iter(line.begin());
    std::string::const_iterator end(line.end());
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("span vector parsing failed");
  }

  std::ostream& operator<<(std::ostream& os, const SpanVector::span_type& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! karma::generate(iter, karma::int_ << '-' << karma::int_ << -karma::buffer[':' << standard::string], x))
      throw std::runtime_error("span generation failed...?");
    
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const SpanVector& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! karma::generate(iter, -((karma::int_ << '-' << karma::int_ << -karma::buffer[':' << standard::string]) % ' '), x))
      throw std::runtime_error("span vector generation failed...?");
    
    return os;
  }

  void SpanVector::span_type::assign(const utils::piece& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string::const_iterator iter(x.begin());
    std::string::const_iterator end(x.end());
    
    qi::rule<std::string::const_iterator, std::string()> label;
    
    label %= qi::lexeme[+(standard::char_ - standard::space)];
    
    const bool result = qi::parse(iter, end,
				  qi::omit[*standard::space] >> qi::int_ >> '-' >> qi::int_ >> -(':' >> label) >> qi::omit[*standard::space],
				  *this);
    
    if (! result || iter != end)
      throw std::runtime_error("span format error: " + x);
  }
  
  std::istream& operator>>(std::istream& is, SpanVector::span_type& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string token;
    
    x.clear();
    if (is >> token)
      x.assign(token);
    return is;
  }
  
  std::istream& operator>>(std::istream& is, SpanVector& x)
  {
    std::string line;
    
    x.clear();
    if (std::getline(is, line)) {
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! x.assign(iter, end) || iter != end)
	throw std::runtime_error("span format error");
    }
    return is;
  }
  
};
