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
#include <stdexcept>

#include "sentence_vector.hpp"

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/thread_specific_ptr.hpp"

namespace cicada
{
  
  typedef std::vector<cicada::Sentence, std::allocator<cicada::Sentence> > sentence_set_type;

  template <typename Iterator>
  struct sentence_vector_parser : boost::spirit::qi::grammar<Iterator, sentence_set_type(), boost::spirit::standard::space_type>
  {
    sentence_vector_parser() : sentence_vector_parser::base_type(sentence_vector)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      word = qi::lexeme[+(standard::char_ - standard::space) - "|||"];
      sentence = *word;
      sentence_vector = -(sentence % "|||");
    }
    typedef boost::spirit::standard::space_type space_type;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type> word;
    boost::spirit::qi::rule<Iterator, cicada::Sentence(), space_type> sentence;
    boost::spirit::qi::rule<Iterator, sentence_set_type(), space_type> sentence_vector;
  };

  void SentenceVector::assign(const std::string& x)
  {
    std::string::const_iterator iter = x.begin();
    std::string::const_iterator end = x.end();
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("sentence vector format parsing failed...");
  }
  
  namespace sentence_vector_impl {
    typedef sentence_vector_parser<std::string::const_iterator> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif
  };
  
  bool SentenceVector::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    using namespace sentence_vector_impl;

#ifdef HAVE_TLS
    if (! __grammar_tls) {
      __grammar.reset(new grammar_type());
      __grammar_tls = __grammar.get();
    }
    
    grammar_type& grammar = *__grammar_tls;
#else
    if (! __grammar.get())
      __grammar.reset(new grammar_type());
    
    grammar_type& grammar = *__grammar;
#endif
    
    clear();
    if (iter == end) return true;
    
    return boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, __sents);
  }
  
  std::ostream& operator<<(std::ostream& os, const SentenceVector& x)
  {
    if (! x.empty()) {
      std::copy(x.begin(), x.end() - 1, std::ostream_iterator<Sentence>(os, " ||| "));
      os << x.back();
    }
    return os;
  }
  
  std::istream& operator>>(std::istream& is, SentenceVector& x)
  {
    std::string line;
    
    x.clear();
    if (std::getline(is, line))
      x.assign(line);
    return is;
  }
};
