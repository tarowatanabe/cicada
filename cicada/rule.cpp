
#include "rule.hpp"

#define BOOST_SPIRIT_THREADSAFE

#include <boost/numeric/conversion/bounds.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"

namespace cicada
{
  
  typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
  typedef std::pair<std::string, double> score_parsed_type;
  typedef std::vector<score_parsed_type, std::allocator<score_parsed_type> > scores_parsed_type;
  typedef boost::fusion::tuple<std::string, phrase_parsed_type, phrase_parsed_type, scores_parsed_type > rule_parsed_type;

  template <typename Iterator>
  struct rule_grammar_parser : boost::spirit::qi::grammar<Iterator, rule_parsed_type(), boost::spirit::standard::space_type>
  {
    
    rule_grammar_parser() : rule_grammar_parser::base_type(rule_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
    
      using qi::phrase_parse;
      using qi::lexeme;
      using qi::no_skip;
      using qi::repeat;
      using qi::hold;
      using qi::lit;
      using qi::inf;
      using qi::attr;
      using standard::char_;
      using qi::double_;
      using qi::_1;
      using standard::space;
      
      lhs    %= (lexeme[char_('[') >> +(char_ - space - ']') >> char_(']')]);
      phrase %= *(lexeme[+(char_ - space) - "|||"]);
      score  %= lexeme[+(char_ - space - '=')] >> '=' >> double_;
      scores %= +score;
      
      rule_grammar %= (hold[lhs >> "|||"] | attr("")) >> phrase >> "|||" >> phrase >> -("|||" >> scores);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> lhs;
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, score_parsed_type(), boost::spirit::standard::space_type>  score;
    boost::spirit::qi::rule<Iterator, scores_parsed_type(), boost::spirit::standard::space_type> scores;
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), boost::spirit::standard::space_type> rule_grammar;
  };

  void Rule::assign(const std::string& x)
  {
    clear();
    
    if (x.empty()) return;
    
    typedef rule_grammar_parser<std::string::const_iterator> grammar_type;
    
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
    
    hypergraph_parser<grammar_type>& grammar = *__grammar;
#endif
    
      
    rule_parsed_type rule_parsed;
    
    std::string::const_iterator iter = x.begin();
    std::string::const_iterator iter_end = x.end();
    
    const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, grammar, boost::spirit::standard::space, rule_parsed);
    if (! result || iter != iter_end)
      throw std::runtime_error("rule parsing failed...");
    
    lhs = boost::fusion::get<0>(rule_parsed);
    if (lhs.empty())
      lhs = vocab_type::X;
    
    source = Rule::symbol_set_type(boost::fusion::get<1>(rule_parsed).begin(), boost::fusion::get<1>(rule_parsed).end());
    target = Rule::symbol_set_type(boost::fusion::get<2>(rule_parsed).begin(), boost::fusion::get<2>(rule_parsed).end());
    features = Rule::feature_set_type(boost::fusion::get<3>(rule_parsed).begin(), boost::fusion::get<3>(rule_parsed).end());
    arity = source.arity();
    
    if (! target.empty() && arity != target.arity())
      throw std::runtime_error("rule parsing failed because of different arity...");
  }
  
  std::istream& operator>>(std::istream& is, Rule& x)
  {
    x.clear();
    
    std::string line;
    if (std::getline(is, line) && ! line.empty()) {
      
      typedef rule_grammar_parser<std::string::const_iterator> grammar_type;
      
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
      
      hypergraph_parser<grammar_type>& grammar = *__grammar;
#endif
      
      rule_parsed_type rule_parsed;
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator iter_end = line.end();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, grammar, boost::spirit::standard::space, rule_parsed);
      if (! result || iter != iter_end)
	throw std::runtime_error("rule parsing failed...");
      
      x.lhs = boost::fusion::get<0>(rule_parsed);
      if (x.lhs.empty())
	x.lhs = Rule::vocab_type::X;
      
      x.source = Rule::symbol_set_type(boost::fusion::get<1>(rule_parsed).begin(), boost::fusion::get<1>(rule_parsed).end());
      x.target = Rule::symbol_set_type(boost::fusion::get<2>(rule_parsed).begin(), boost::fusion::get<2>(rule_parsed).end());
      x.features = Rule::feature_set_type(boost::fusion::get<3>(rule_parsed).begin(), boost::fusion::get<3>(rule_parsed).end());
      x.arity = x.source.arity();
      
      if (! x.target.empty() && x.arity != x.target.arity())
	throw std::runtime_error("rule parsing failed because of different arity...");
    }
    return is;
  }
  
  
  std::ostream& operator<<(std::ostream& os, const Rule& x)
  {
    os << x.lhs << " ||| " << x.source << " ||| " << x.target;
    if (! x.features.empty()) {
      os << " |||";
      Rule::feature_set_type::const_iterator fiter_end = x.features.end();
      for (Rule::feature_set_type::const_iterator fiter = x.features.begin(); fiter != fiter_end; ++ fiter)
	os << ' ' << fiter->first << '=' << fiter->second;
    }
    
    return os;
  }
  
};
