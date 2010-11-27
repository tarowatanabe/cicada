//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include "rule.hpp"

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
#include "utils/array_power2.hpp"

namespace cicada
{
  
  Rule::rule_ptr_type Rule::create(const Rule& x)
  {
    typedef utils::array_power2<rule_ptr_type, 1024 * 128, std::allocator<rule_ptr_type> > cache_type;

#ifdef HAVE_TLS
    static __thread cache_type* __cache_tls = 0;
    static boost::thread_specific_ptr<cache_type> __cache;
    
    if (! __cache_tls) {
      __cache.reset(new cache_type());
      __cache_tls = __cache.get();
    }
    cache_type& cache = *__cache_tls;
#else
    static boost::thread_specific_ptr<cache_type> __cache;
    if (! __cache.get())
      __cache.reset(new cache_type());
    cache_type& cache = *__cache;
#endif
    
    const size_t cache_pos = hash_value(x) & (cache.size() - 1);
    if (*cache[cache_pos] != x)
      cache[cache_pos].reset(new Rule(x));
    return cache[cache_pos];
  }
  
  typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
  
  typedef boost::fusion::tuple<std::string, phrase_parsed_type > rule_parsed_type;

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
      using qi::hold;
      using qi::attr;
      using standard::char_;
      using qi::double_;
      using standard::space;
      
      lhs    %= (lexeme[char_('[') >> +(char_ - space - ']') >> char_(']')]);
      phrase %= *(lexeme[+(char_ - space) - "|||"]);
      
      rule_grammar %= (hold[lhs >> "|||"] | attr("")) >> phrase;
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> lhs;
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), boost::spirit::standard::space_type> rule_grammar;
  };

  bool Rule::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
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
    
    grammar_type& grammar = *__grammar;
#endif

    clear();
      
    rule_parsed_type rule_parsed;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, rule_parsed);
    if (result) {
      lhs = boost::fusion::get<0>(rule_parsed);
      if (lhs.empty())
	lhs = vocab_type::X;
      
      rhs = Rule::symbol_set_type(boost::fusion::get<1>(rule_parsed).begin(), boost::fusion::get<1>(rule_parsed).end());
    }

    return result;
  }

  void Rule::assign(const std::string& x)
  {
    clear();
    
    if (x.empty()) return;

    std::string::const_iterator iter = x.begin();
    std::string::const_iterator end = x.end();
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("rule format parsing failed...\"" + x + "\"");
  }
  
  std::istream& operator>>(std::istream& is, Rule& x)
  {
    x.clear();
    
    std::string line;
    if (std::getline(is, line) && ! line.empty())
      x.assign(line);
    
    return is;
  }
  
  
  std::ostream& operator<<(std::ostream& os, const Rule& x)
  {
    os << x.lhs << " ||| " << x.rhs;
    return os;
  }
  
};
