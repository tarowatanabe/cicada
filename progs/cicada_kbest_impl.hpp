//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__KBEST_IMPL__HPP__
#define __CICADA__KBEST_IMPL__HPP__ 1

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
#include <utils/base64.hpp>
#include <utils/double_base64_parser.hpp>
#include <utils/double_base64_generator.hpp>

#include "cicada/sentence.hpp"
#include "cicada/eval.hpp"
#include "cicada/feature.hpp"
#include "cicada/symbol.hpp"
#include "cicada/feature_vector_compact.hpp"

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
  //typedef utils::simple_vector<feature_value_type, std::allocator<feature_value_type> > feature_set_type;
  typedef cicada::FeatureVectorCompact feature_set_type;
  
  hypothesis_type() : sentence(), features(), score(), loss(0) {}
  hypothesis_type(const kbest_feature_type& x)
    : sentence(boost::fusion::get<1>(x).begin(), boost::fusion::get<1>(x).end()),
      features(boost::fusion::get<2>(x).begin(), boost::fusion::get<2>(x).end()),
      score(),
      loss(0)
  {
    
  }
  template <typename IteratorSentence, typename IteratorFeature>
  hypothesis_type(IteratorSentence sfirst, IteratorSentence slast,
		  IteratorFeature  ffirst, IteratorFeature  flast)
    : sentence(sfirst, slast),
      features(ffirst, flast),
      score(),
      loss(0)
  {
    
  }
  
  sentence_type    sentence;
  feature_set_type features;
  score_ptr_type   score;
  double           loss;
};


inline
size_t hash_value(hypothesis_type const& x)
{
  typedef utils::hashmurmur<size_t> hasher_type;
  
  return hasher_type()(x.sentence.begin(), x.sentence.end(), hasher_type()(x.features.begin(), x.features.end(), 0));
}

inline
bool operator==(const hypothesis_type& x, const hypothesis_type& y)
{
  return x.sentence == y.sentence && x.features == y.features;
}

inline
bool operator<(const hypothesis_type& x, const hypothesis_type& y)
{
  return x.sentence < y.sentence || (!(y.sentence < x.sentence) && x.features < y.features);
}

typedef std::vector<hypothesis_type, std::allocator<hypothesis_type> > hypothesis_set_type;
typedef std::vector<hypothesis_set_type, std::allocator<hypothesis_set_type> > hypothesis_map_type;

template <typename Iterator>
struct kbest_feature_parser : boost::spirit::qi::grammar<Iterator, kbest_feature_type(), boost::spirit::standard::blank_type>
{
  kbest_feature_parser() : kbest_feature_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
        
    // TODO: we want to handle longest character sequences... HOW?
    feature %= qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi)) >> (standard::char_ - standard::space))] >> '=' >> qi::double_;
    features %= -(feature % (+standard::space));
    
    kbest %= size >> "|||" >> tokens >> -("|||" >> features) >> -("|||" >> qi::double_) >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
    
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  
  boost::spirit::qi::rule<Iterator, feature_parsed_type()>  feature;
  boost::spirit::qi::rule<Iterator, feature_parsed_set_type()> features;
  
  boost::spirit::qi::rule<Iterator, kbest_feature_type(), blank_type>  kbest;
};

template <typename Iterator>
struct kbest_feature_generator : boost::spirit::karma::grammar<Iterator, kbest_feature_type()>
{
  kbest_feature_generator() : kbest_feature_generator::base_type(kbest) 
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    tokens %= standard::string % ' ';
    features %= (standard::string << '=' << double20) % ' ';
    kbest %= size << " ||| " << tokens << -karma::buffer[" ||| " << features];
  }
  
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };
  
  boost::spirit::karma::uint_generator<size_type>              size;
  boost::spirit::karma::real_generator<double, real_precision> double20;
  
  boost::spirit::karma::rule<Iterator, tokens_type()>             tokens;
  boost::spirit::karma::rule<Iterator, feature_parsed_set_type()> features;
  boost::spirit::karma::rule<Iterator, kbest_feature_type()>      kbest;
    

};

#endif
