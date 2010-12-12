//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "grammar_mutable.hpp"
#include "parameter.hpp"
#include "attribute_vector.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/compress_stream.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"

namespace std
{
  std::ostream& operator<<(std::ostream& os, const std::pair<std::string, double>& x)
  {
    os << x.first << "=" << x.second;
    return os;
  }
};

namespace cicada
{
  
  class GrammarMutableImpl
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol    symbol_type;
    typedef cicada::Symbol    word_type;
    typedef cicada::Feature   feature_type;
    typedef cicada::Attribute attribute_type;
    typedef cicada::Vocab     vocab_type;

    typedef Transducer::rule_type          rule_type;
    typedef Transducer::rule_ptr_type      rule_ptr_type;
    typedef Transducer::rule_pair_type     rule_pair_type;
    typedef Transducer::rule_pair_set_type rule_pair_set_type;
    
    typedef utils::compact_trie_dense<symbol_type, rule_pair_set_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				      std::allocator<std::pair<const symbol_type, rule_pair_set_type> > > trie_type;
    typedef trie_type::id_type id_type;

    typedef std::vector<feature_type, std::allocator<feature_type> >     feature_name_set_type;
    typedef std::vector<attribute_type, std::allocator<attribute_type> > attribute_name_set_type;

    GrammarMutableImpl() : trie(symbol_type()) {}
    
    void read(const std::string& parameter);
    
    id_type root() const { return trie.root(); }
    id_type next(id_type node, const symbol_type& symbol) const { return trie.find(node, symbol); }
    bool has_next(id_type node) const { return ! trie.empty(node); }
    const rule_pair_set_type& rules(id_type node) { return trie[node]; }

    void insert(const std::string& pattern);
    void insert(const rule_pair_type& rule)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
      
      if (rule.source->rhs.empty()) return;
      
      // no checking for non-terminal index!
      
      sequence_type source_index(rule.source->rhs.begin(), rule.source->rhs.end());
      sequence_type::iterator siter_end = source_index.end();
      for (sequence_type::iterator siter = source_index.begin(); siter != siter_end; ++ siter)
	*siter = siter->non_terminal();
      
      const id_type id = trie.insert(source_index.begin(), source_index.end());
      
      trie[id].push_back(rule);
    }
    
    
    void clear() { trie.clear(); }

    void read();

  private:
    trie_type trie;

    feature_name_set_type   feature_names_default;
    attribute_name_set_type attribute_names_default;
  };
  
  typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
  typedef std::pair<std::string, double> score_parsed_type;
  typedef std::vector<score_parsed_type, std::allocator<score_parsed_type> > scores_parsed_type;
  
  typedef std::pair<std::string, AttributeVector::data_type> attr_parsed_type;
  typedef std::vector<attr_parsed_type, std::allocator<attr_parsed_type> > attrs_parsed_type;
  
  typedef boost::fusion::tuple<std::string, phrase_parsed_type, phrase_parsed_type, scores_parsed_type, attrs_parsed_type> rule_parsed_type;
  
  template <typename Iterator>
  struct rule_grammar_parser_mutable : boost::spirit::qi::grammar<Iterator, rule_parsed_type(), boost::spirit::standard::space_type>
  {
    
    rule_grammar_parser_mutable() : rule_grammar_parser_mutable::base_type(rule_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
    
      escape_char.add
	("\\\"", '\"')
	("\\\\", '\\')
	("\\/", '/')
	("\\b", '\b')
	("\\f", '\f')
	("\\n", '\n')
	("\\r", '\r')
	("\\t", '\t')
	("\\u0020", ' ');
      
      lhs %= (qi::lexeme[standard::char_('[') >> +(standard::char_ - standard::space - ']') >> standard::char_(']')]);
      phrase %= *(qi::lexeme[+(standard::char_ - standard::space) - "|||"]);

      score %= (qi::hold[qi::lexeme[+(standard::char_ - standard::space - '=')] >> '='] | qi::attr("")) >> qi::double_;
      scores %= *score;
      
      data_value %= ('\"' >> qi::lexeme[*(escape_char | (standard::char_ - '\"'))] >> '\"');
      data %= data_value | double_dot | int64_;
      
      attribute %= (qi::hold[qi::lexeme[+(standard::char_ - standard::space - '=')] >> '='] | qi::attr("")) >> data;
      attributes %= *attribute;
      
      rule_grammar %= (qi::hold[lhs >> "|||"] | qi::attr("")) >> phrase >> "|||" >> phrase >> -("|||" >> scores) >> -("|||" >> attributes);
    }
  
    typedef boost::spirit::standard::space_type space_type;

    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1> int64_;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > double_dot;
    
    boost::spirit::qi::symbols<char, char> escape_char;

    boost::spirit::qi::rule<Iterator, std::string(), space_type> lhs;
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), space_type> phrase;
    
    boost::spirit::qi::rule<Iterator, score_parsed_type(), space_type>  score;
    boost::spirit::qi::rule<Iterator, scores_parsed_type(), space_type> scores;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type>                data_value;
    boost::spirit::qi::rule<Iterator, AttributeVector::data_type(), space_type> data;
    boost::spirit::qi::rule<Iterator, attr_parsed_type(), space_type>           attribute;
    boost::spirit::qi::rule<Iterator, attrs_parsed_type(), space_type>          attributes;
    
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), space_type> rule_grammar;
  };

  
  void GrammarMutableImpl::insert(const std::string& pattern)
  {
    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
    
    typedef rule_grammar_parser_mutable<std::string::const_iterator> rule_parser_type;
#ifdef HAVE_TLS
    static __thread rule_parser_type* __rule_parser_tls = 0;
    static boost::thread_specific_ptr<rule_parser_type > __rule_parser;
    
    if (! __rule_parser_tls) {
      __rule_parser.reset(new rule_parser_type());
      __rule_parser_tls = __rule_parser.get();
    }
    
    rule_parser_type& rule_parser = *__rule_parser_tls;
#else
    static boost::thread_specific_ptr<rule_parser_type > __rule_parser;
    if (! __rule_parser.get())
      __rule_parser.reset(new rule_parser_type());
    
    rule_parser_type& rule_parser = *__rule_parser;
#endif
    
    rule_parsed_type rule_parsed;

    std::string::const_iterator iter_begin = pattern.begin();
    std::string::const_iterator iter_end = pattern.end();
    std::string::const_iterator iter = iter_begin;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule_parsed);
    
    if (! result || iter != iter_end)
      throw std::runtime_error(std::string("rule parsing failed: ") + pattern);    
    
    const symbol_type lhs(boost::fusion::get<0>(rule_parsed).empty() ? vocab_type::X : symbol_type(boost::fusion::get<0>(rule_parsed)));

    rule_type rule_source(lhs, boost::fusion::get<1>(rule_parsed).begin(), boost::fusion::get<1>(rule_parsed).end());
    rule_type rule_target(lhs, boost::fusion::get<2>(rule_parsed).begin(), boost::fusion::get<2>(rule_parsed).end());
    
    if (! rule_target.rhs.empty() && rule_source.rhs.arity() != rule_target.rhs.arity())
      throw std::runtime_error("arity do not match");
    
    cicada::sort(rule_source, rule_target);
    
    // we will transform into "plain" non-terminal, so that we can query!
    sequence_type source_index(rule_source.rhs.size());
    sequence_type::iterator iiter = source_index.begin();
    rule_type::symbol_set_type::const_iterator siter_end = rule_source.rhs.end();
    for (rule_type::symbol_set_type::const_iterator siter = rule_source.rhs.begin(); siter != siter_end; ++ siter, ++ iiter)
      *iiter = siter->non_terminal();
    
    // we do not check duplicates!
    rule_pair_set_type& rules = trie[trie.insert(source_index.begin(), source_index.end())];
    rules.push_back(rule_pair_type(rule_type::create(rule_source), rule_type::create(rule_target)));
    
    // add features...
    const scores_parsed_type& scores = boost::fusion::get<3>(rule_parsed);
    
    int feature = 0;
    scores_parsed_type::const_iterator fiter_end = scores.end();
    for (scores_parsed_type::const_iterator fiter = scores.begin(); fiter != fiter_end; ++ fiter)
      if (fiter->first.empty()) {
	// default name!
	
	if (feature >= static_cast<int>(feature_names_default.size()))
	  feature_names_default.resize(feature + 1);
	if (feature_names_default[feature].empty())
	  feature_names_default[feature] = "rule-table-" + boost::lexical_cast<std::string>(feature);
	
	rules.back().features[feature_names_default[feature]] = fiter->second;
	
	++ feature;
      } else
	rules.back().features[fiter->first] = fiter->second;
    
    int attribute = 0;
    attrs_parsed_type::const_iterator aiter_end = boost::fusion::get<4>(rule_parsed).end();
    for (attrs_parsed_type::const_iterator aiter = boost::fusion::get<4>(rule_parsed).begin(); aiter != aiter_end; ++ aiter) 
      if (aiter->first.empty()) {
	
	if (attribute >= static_cast<int>(attribute_names_default.size()))
	  attribute_names_default.resize(attribute + 1);
	if (attribute_names_default[attribute].empty())
	  attribute_names_default[attribute] = "rule-table-" + boost::lexical_cast<std::string>(attribute);
	
	rules.back().attributes[attribute_names_default[attribute]] = aiter->second;
	
	++ attribute;
      } else
	rules.back().attributes[aiter->first] = aiter->second;
  }
  
  void GrammarMutableImpl::read(const std::string& parameter)
  {
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;

    typedef cicada::Parameter parameter_type;
    
    
    const parameter_type param(parameter);
    
    if (param.name() != "-" && ! boost::filesystem::exists(param.name()))
      throw std::runtime_error("no grammar file: " + param.name());
    
    feature_name_set_type   feature_names;
    attribute_name_set_type attribute_names;
    parameter_type::iterator piter_end = param.end();
    for (parameter_type::iterator piter = param.begin(); piter != piter_end; ++ piter) {
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;

      {
	std::string::const_iterator iter = piter->first.begin();
	std::string::const_iterator iter_end = piter->first.end();
	
	int feature_id = -1;
	const bool result = qi::parse(iter, iter_end, "feature" >> qi::int_[phoenix::ref(feature_id) = qi::_1]);
	if (result && iter == iter_end && feature_id >= 0) {
	  if (feature_id >= int(feature_names.size()))
	    feature_names.resize(feature_id + 1);
	  
	  feature_names[feature_id] = piter->second;
	  continue;
	}
      }
      
      {
	std::string::const_iterator iter = piter->first.begin();
	std::string::const_iterator iter_end = piter->first.end();
	
	int attribute_id = -1;
	const bool result = qi::parse(iter, iter_end, "attribute" >> qi::int_[phoenix::ref(attribute_id) = qi::_1]);
	if (result && iter == iter_end && attribute_id >= 0) {
	  if (attribute_id >= int(attribute_names.size()))
	    attribute_names.resize(attribute_id + 1);
	  
	  attribute_names[attribute_id] = piter->second;
	  continue;
	}
      }
      
      throw std::runtime_error("unsupported key: " + piter->first);
    }

    typedef rule_grammar_parser_mutable<std::string::const_iterator> rule_parser_type;
    
#ifdef HAVE_TLS
    static __thread rule_parser_type* __rule_parser_tls = 0;
    static boost::thread_specific_ptr<rule_parser_type > __rule_parser;
    
    if (! __rule_parser_tls) {
      __rule_parser.reset(new rule_parser_type());
      __rule_parser_tls = __rule_parser.get();
    }
    
    rule_parser_type& rule_parser = *__rule_parser_tls;
#else
    static boost::thread_specific_ptr<rule_parser_type > __rule_parser;
    if (! __rule_parser.get())
      __rule_parser.reset(new rule_parser_type());
    
    rule_parser_type& rule_parser = *__rule_parser;
#endif

    
    utils::compress_istream is(param.name(), 1024 * 1024);
    std::string line;
    rule_parsed_type rule_parsed;
    
    sequence_type source_index;

    while (std::getline(is, line)) {
      if (line.empty()) continue;
      
      boost::fusion::get<0>(rule_parsed).clear();
      boost::fusion::get<1>(rule_parsed).clear();
      boost::fusion::get<2>(rule_parsed).clear();
      boost::fusion::get<3>(rule_parsed).clear();
      boost::fusion::get<4>(rule_parsed).clear();
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator iter_end = line.end();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule_parsed);
      
      if (! result || iter != iter_end)
	throw std::runtime_error("rule parsing failed: " + line);
      
      const symbol_type lhs(boost::fusion::get<0>(rule_parsed).empty() ? vocab_type::X : symbol_type(boost::fusion::get<0>(rule_parsed)));

      rule_type rule_source(lhs, boost::fusion::get<1>(rule_parsed).begin(), boost::fusion::get<1>(rule_parsed).end());
      rule_type rule_target(lhs, boost::fusion::get<2>(rule_parsed).begin(), boost::fusion::get<2>(rule_parsed).end());
      
      if (! rule_target.rhs.empty() && rule_source.rhs.arity() != rule_target.rhs.arity())
      throw std::runtime_error("arity do not match");
      
      cicada::sort(rule_source, rule_target);
      
      source_index.resize(rule_source.rhs.size());
      sequence_type::iterator iiter = source_index.begin();
      rule_type::symbol_set_type::const_iterator siter_end = rule_source.rhs.end();
      for (rule_type::symbol_set_type::const_iterator siter = rule_source.rhs.begin(); siter != siter_end; ++ siter, ++ iiter)
	*iiter = siter->non_terminal();
      
      // we do not check duplicates!
      rule_pair_set_type& rules = trie[trie.insert(source_index.begin(), source_index.end())];
      rules.push_back(rule_pair_type(rule_type::create(rule_source), rule_type::create(rule_target)));
      
      // add features...
      const scores_parsed_type& scores = boost::fusion::get<3>(rule_parsed);
      
      int feature = 0;
      scores_parsed_type::const_iterator fiter_end = scores.end();
      for (scores_parsed_type::const_iterator fiter = scores.begin(); fiter != fiter_end; ++ fiter) 
	if (fiter->first.empty()) {
	  
	  if (feature < int(feature_names.size()) && ! feature_names[feature].empty())
	    rules.back().features[feature_names[feature]] = fiter->second;
	  else {
	    // default name!
	    if (feature >= static_cast<int>(feature_names_default.size()))
	      feature_names_default.resize(feature + 1);
	    if (feature_names_default[feature].empty())
	      feature_names_default[feature] = "rule-table-" + boost::lexical_cast<std::string>(feature);
	    
	    rules.back().features[feature_names_default[feature]] = fiter->second;
	  }
	  
	  ++ feature;
	} else
	  rules.back().features[fiter->first] = fiter->second;
      
      int attribute = 0;
      attrs_parsed_type::const_iterator aiter_end = boost::fusion::get<4>(rule_parsed).end();
      for (attrs_parsed_type::const_iterator aiter = boost::fusion::get<4>(rule_parsed).begin(); aiter != aiter_end; ++ aiter) 
	if (aiter->first.empty()) {
	  
	  if (attribute < int(attribute_names.size()) && ! attribute_names[attribute].empty())
	    rules.back().attributes[attribute_names[attribute]] = aiter->second;
	  else {
	    // default name!
	    if (attribute >= static_cast<int>(attribute_names_default.size()))
	      attribute_names_default.resize(attribute + 1);
	    if (attribute_names_default[attribute].empty())
	      attribute_names_default[attribute] = "rule-table-" + boost::lexical_cast<std::string>(attribute);
	    
	    rules.back().attributes[attribute_names_default[attribute]] = aiter->second;
	  }
	  
	  ++ attribute;
	} else
	  rules.back().attributes[aiter->first] = aiter->second;
    }
  }
  
  void GrammarMutable::read(const std::string& parameter)
  {
    pimpl->read(parameter);
  }

  void GrammarMutable::clear()
  {
    pimpl->clear();
  }

  void GrammarMutable::insert(const std::string& pattern)
  {
    pimpl->insert(pattern);
  }

  void GrammarMutable::insert(const rule_pair_type& rule)
  {
    pimpl->insert(rule);
  }
  
  GrammarMutable::GrammarMutable()
    : pimpl(new impl_type()) {}
  GrammarMutable::GrammarMutable(const std::string& parameter)
    : pimpl(new impl_type())
  {
    pimpl->read(parameter);
  }
  GrammarMutable::~GrammarMutable() { std::auto_ptr<impl_type>(tmp); }
  
  GrammarMutable::GrammarMutable(const GrammarMutable& x)
    : pimpl(new impl_type(*x.pimpl)) {}
  
  GrammarMutable& GrammarMutable::operator=(const GrammarMutable& x)
  {
    *pimpl = *x.pimpl;
    return *this;
  }
  
  GrammarMutable::transducer_ptr_type GrammarMutable::clone() const
  {
    return transducer_ptr_type(new GrammarMutable(*this));
  }
  
  
  GrammarMutable::id_type GrammarMutable::root() const
  {
    return pimpl->root();
  }
  
  GrammarMutable::id_type GrammarMutable::next(const id_type& node, const symbol_type& symbol) const
  {
    return pimpl->next(node, symbol.non_terminal());
  }
  
  bool GrammarMutable::has_next(const id_type& node) const
  {
    return pimpl->has_next(node);
  }
  
  const GrammarMutable::rule_pair_set_type& GrammarMutable::rules(const id_type& node) const
  {
    static const rule_pair_set_type __empty;
    return (node == pimpl->root() ? __empty : pimpl->rules(node));
  }
  
};
