 
#include <iterator>
#include <algorithm>
#include <string>
#include <stdexcept>

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

#include <boost/thread.hpp>

#include "utils/config.hpp"

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
  struct span_vector_parser : boost::spirit::qi::grammar<Iterator, Container(), boost::spirit::standard::space_type>
  {
    span_vector_parser() : span_vector_parser::base_type(spans)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::phrase_parse;
      using qi::int_;
      using qi::lexeme;
      using standard::space;
      using standard::char_;
      
      label %= lexeme[+(char_ - space)];
      spans %= *(int_ >> '-' >> int_ >> -(':' >> label));
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> label;
    boost::spirit::qi::rule<Iterator, Container(), boost::spirit::standard::space_type> spans;
  };
  

  bool SpanVector::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    typedef span_vector_parser<std::string::const_iterator, spans_type> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
    
    if (! __grammar_tls) {
      __grammar.reset(new grammar_type());
      __grammar_tls = __grammar.get();
    }
    
    grammar_type& grammar = *__grammar_tls;
#else
    static boost::thread_specific_ptr<grammar_type > __grammar;
    if (! __grammar.get())
      __grammar.reset(new grammar_type());
    
    grammar_type& grammar = *__grammar;
#endif
    
    return phrase_parse(iter, end,
			grammar,
			boost::spirit::standard::space,
			__spans);
  }


  std::ostream& operator<<(std::ostream& os, const SpanVector::span_type& x)
  {
    if (x.label.empty())
      os << x.first << '-' << x.last;
    else
      os << x.first << '-' << x.last << ':' << x.label;
    return os;
  }

  std::ostream& operator<<(std::ostream& os, const SpanVector& x)
  {
    if (! x.empty()) {
      std::copy(x.begin(), x.end() - 1, std::ostream_iterator<SpanVector::span_type>(os, " "));
      os << x.back();
    }
    return os;
  }

  std::istream& operator>>(std::istream& is, SpanVector::span_type& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    using qi::phrase_parse;
    using qi::int_;
    using qi::lexeme;
    using standard::space;
    using standard::char_;

    std::string token;
    
    x.clear();
    if (is >> token) {
      std::string::const_iterator iter = token.begin();
      std::string::const_iterator end = token.end();
      
      qi::rule<std::string::const_iterator, std::string(), standard::space_type> label;

      label %= lexeme[+(char_ - space)];
      
      const bool result = phrase_parse(iter, end,
				       int_ >> '-' >> int_ >> -(':' >> label),
				       space,
				       x);
      
      if (! result || iter != end)
	throw std::runtime_error("span format error");
    }
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
