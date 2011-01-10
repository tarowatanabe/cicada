//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// from conll dependency to hypergraph conversion...
//
// we may use POS or semantic-role as our label...
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

// this is important!
#define FUSION_MAX_VECTOR_SIZE 15

#include <boost/variant.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>

#include <cicada/hypergraph.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

struct conll_type
{
  typedef size_t size_type;
  typedef boost::variant<size_type, std::string> phead_type;

  size_type   id;
  std::string form;
  std::string lemma;
  std::string cpostag;
  std::string postag;
  std::string feats;
  size_type   head;
  std::string deprel;
  phead_type  phead;
  std::string pdeprel;
};

template <typename Iterator>
struct conll_parser : boost::spirit::qi::grammar<Iterator, conll_set_type(), boost::spirit::standard::blank_type>
{
  conll_parser() : conll_parser::base_type(conlls)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    token %= qi::lexeme[+(standard::char_ - standard::space)];
    
    conll  %= size >> token >> token >> token >> token >> token >> size >> token >> (size | token) >> token >> qi::eol;
    conlls %= *conll >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>            size;
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>    token;
  boost::spirit::qi::rule<Iterator, conll_type(), blank_type>     conll;
  boost::spirit::qi::rule<Iterator, conll_set_type(), blank_type> conlls;
  
};
