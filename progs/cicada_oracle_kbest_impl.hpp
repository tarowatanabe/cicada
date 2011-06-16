//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__ORACLE_KBEST_IMPL__HPP__
#define __CICADA__ORACLE_KBEST_IMPL__HPP__ 1

#include <vector>
#include <deque>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <utils/simple_vector.hpp>

#include "cicada/sentence.hpp"
#include "cicada/eval.hpp"
#include "cicada/feature.hpp"
#include "cicada/symbol.hpp"

typedef cicada::Sentence sentence_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type  score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;

typedef size_t size_type;
typedef std::vector<std::string, std::allocator<std::string> > tokens_type;
typedef std::pair<std::string, double> feature_parsed_type;
typedef std::vector<feature_parsed_type, std::allocator<feature_parsed_type> > feature_parsed_set_type;
typedef boost::fusion::tuple<size_type, tokens_type, feature_parsed_set_type> kbest_feature_type;

struct hypothesis_type
{
  typedef cicada::Symbol  word_type;
  typedef cicada::Feature feature_type;
  typedef std::pair<feature_type, double> feature_value_type;
  
  typedef utils::simple_vector<word_type, std::allocator<word_type> >                   sentence_type;
  typedef utils::simple_vector<feature_value_type, std::allocator<feature_value_type> > feature_set_type;
  
  hypothesis_type() : sentence(), features(), score() {}
  hypothesis_type(const kbest_feature_type& x)
    : sentence(boost::fusion::get<1>(x).begin(), boost::fusion::get<1>(x).end()),
      features(boost::fusion::get<2>(x).begin(), boost::fusion::get<2>(x).end()),
      score()
  {
    std::sort(features.begin(), features.end());
  }
  
  sentence_type    sentence;
  feature_set_type features;
  score_ptr_type   score;
};

typedef std::vector<hypothesis_type, std::allocator<hypothesis_type> > hypothesis_set_type;
typedef std::vector<hypothesis_set_type, std::allocator<hypothesis_set_type> > hypothesis_map_type;

typedef std::vector<const hypothesis_type*, std::allocator<const hypothesis_type*> > oracle_set_type;
typedef std::vector<oracle_set_type, std::allocator<oracle_set_type> > oracle_map_type;

template <typename Iterator>
struct kbest_feature_parser : boost::spirit::qi::grammar<Iterator, kbest_feature_type(), boost::spirit::standard::blank_type>
{
  kbest_feature_parser() : kbest_feature_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    remains %= *qi::lexeme[+(standard::char_ - standard::space)];
    
    feature %= qi::lexeme[+(standard::char_ - standard::space - '=')] >> '=' >> qi::double_;
    features %= *feature;
    
    kbest %= size >> "|||" >> tokens >> "|||" >> features >> -("|||" >> remains) >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  
  boost::spirit::qi::rule<Iterator, feature_parsed_type(), blank_type>  feature;
  boost::spirit::qi::rule<Iterator, feature_parsed_set_type(), blank_type> features;

  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  
  boost::spirit::qi::rule<Iterator, kbest_feature_type(), blank_type>  kbest;
};


#endif
