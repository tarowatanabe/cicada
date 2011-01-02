//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TEXT__IMPL__HPP__
#define __CICADA__TEXT__IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

#include <boost/fusion/adapted.hpp>

#include <cicada/sentence.hpp>

typedef std::pair<int, cicada::Sentence> id_sentence_type;

template <typename Iterator>
struct cicada_sentence_parser : boost::spirit::qi::grammar<Iterator, id_sentence_type(), boost::spirit::standard::blank_type>
{
  cicada_sentence_parser() : cicada_sentence_parser::base_type(id_sentence)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    word        %= qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    sentence    %= *word;
    id_sentence %= qi::int_ >> "|||" >> sentence >> -qi::omit["|||" >> *(standard::char_ - qi::eol)] >> (qi::eol | qi::eoi);
  };
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>      word;
  boost::spirit::qi::rule<Iterator, cicada::Sentence(), blank_type> sentence;
  boost::spirit::qi::rule<Iterator, id_sentence_type(), blank_type> id_sentence;
};

#endif

