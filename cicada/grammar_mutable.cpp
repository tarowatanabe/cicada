#include <memory>

#include "grammar_mutable.hpp"
#include "parameter.hpp"

#include "utils/compact_trie.hpp"
#include "utils/compress_stream.hpp"

#include <boost/lexical_cast.hpp>

#define BOOST_SPIRIT_THREADSAFE

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
    
    typedef cicada::Symbol  symbol_type;
    typedef cicada::Symbol  word_type;
    typedef cicada::Feature feature_type;
    typedef cicada::Vocab   vocab_type;

    typedef Transducer::rule_type     rule_type;
    typedef Transducer::rule_ptr_type rule_ptr_type;
    typedef Transducer::rule_set_type rule_set_type;
    

    typedef utils::compact_trie<symbol_type, rule_set_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				std::allocator<std::pair<const symbol_type, rule_set_type> > > trie_type;
    typedef trie_type::id_type id_type;
    
    void read(const std::string& parameter);
    
    id_type root() const { return trie.root(); }
    id_type next(id_type node, const symbol_type& symbol) const { return trie.find(node, symbol); }
    bool has_next(id_type node) const { return ! trie.empty(node); }
    const rule_set_type& rules(id_type node) { return trie[node]; }

    void insert(const std::string& pattern);
    void insert(const rule_type& rule)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
      
      if (rule.source.empty()) return;
      
      sequence_type source_index(rule.source.begin(), rule.source.end());
      sequence_type::iterator siter_end = source_index.end();
      for (sequence_type::iterator siter = source_index.begin(); siter != siter_end; ++ siter)
	*siter = siter->non_terminal();      
      
      const id_type id = trie.insert(source_index.begin(), source_index.end());
      
      trie[id].push_back(boost::shared_ptr<rule_type>(new rule_type(rule)));
      trie[id].back()->sort_source_index();
    }

    void insert(const rule_ptr_type& rule)
    {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
      
      if (rule->source.empty()) return;
      
      sequence_type source_index(rule->source.begin(), rule->source.end());
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
  };
  
  
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
      using qi::attr;
      using qi::hold;
      using standard::char_;
      using qi::double_;
      using qi::_1;
      using standard::space;
      
      lhs %= (lexeme[char_('[') >> +(char_ - space - ']') >> char_(']')]);
      phrase %= *(lexeme[+(char_ - space) - "|||"]);

      score %= (hold[lexeme[+(char_ - space - '=')] >> '='] | attr("")) >> double_;
      scores %= +score;
      
      rule_grammar %= (hold[lhs >> "|||"] | attr("")) >> phrase >> "|||" >> phrase >> -("|||" >> scores);
    }
  
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> lhs;
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, score_parsed_type(), boost::spirit::standard::space_type>  score;
    boost::spirit::qi::rule<Iterator, scores_parsed_type(), boost::spirit::standard::space_type> scores;
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), boost::spirit::standard::space_type> rule_grammar;
  };

  
  void GrammarMutableImpl::insert(const std::string& pattern)
  {
    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
    

    rule_parsed_type rule;
    
    typedef rule_grammar_parser<std::string::const_iterator> rule_parser_type;
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


    std::string::const_iterator iter_begin = pattern.begin();
    std::string::const_iterator iter_end = pattern.end();
    std::string::const_iterator iter = iter_begin;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);
    
    if (! result || iter != iter_end)
      throw std::runtime_error(std::string("rule parsing failed: ") + pattern);    
    
    const symbol_type lhs(boost::fusion::get<0>(rule).empty() ? vocab_type::X : symbol_type(boost::fusion::get<0>(rule)));
    const rule_type::symbol_set_type source(boost::fusion::get<1>(rule).begin(), boost::fusion::get<1>(rule).end());
    const rule_type::symbol_set_type target(boost::fusion::get<2>(rule).begin(), boost::fusion::get<2>(rule).end());
    
    const size_type arity = source.arity();
    if (! target.empty() && arity != target.arity())
      throw std::runtime_error("arity do not match");

    rule_ptr_type rule_ptr(new rule_type(lhs, source, target, arity));
    rule_ptr->sort_source_index();
    
    // we will transform into "plain" non-terminal, so that we can query!
    sequence_type source_index(source.size());
    sequence_type::iterator iiter = source_index.begin();
    rule_type::symbol_set_type::const_iterator siter_end = source.end();
    for (rule_type::symbol_set_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter, ++ iiter)
      *iiter = siter->non_terminal();
    
    // we do not check duplicates!
    trie[trie.insert(source_index.begin(), source_index.end())].push_back(rule_ptr);
    
    // add features...
    const scores_parsed_type& scores = boost::fusion::get<3>(rule);
      
    int feature = 0;
    scores_parsed_type::const_iterator fiter_end = scores.end();
    for (scores_parsed_type::const_iterator fiter = scores.begin(); fiter != fiter_end; ++ fiter)
      if (fiter->first.empty()) {
	// default name!
	const std::string name = std::string("rule-table-") + boost::lexical_cast<std::string>(feature);
	
	rule_ptr->features[name] = fiter->second;
	
	++ feature;
      } else
	rule_ptr->features[fiter->first] = fiter->second;
  }
  
  void GrammarMutableImpl::read(const std::string& parameter)
  {
    typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;

    typedef cicada::Parameter parameter_type;
    
    
    const parameter_type param(parameter);
    
    if (param.name() != "-" && ! boost::filesystem::exists(param.name()))
      throw std::runtime_error("no grammar file: " + param.name());
    
    feature_name_set_type feature_names;
    parameter_type::iterator piter_end = param.end();
    for (parameter_type::iterator piter = param.begin(); piter != piter_end; ++ piter) {

      std::string::const_iterator iter = piter->first.begin();
      std::string::const_iterator iter_end = piter->first.end();

      int feature_id = -1;
      
      const bool result = boost::spirit::qi::parse(iter, iter_end,
						   "feature" >> boost::spirit::qi::int_[boost::phoenix::ref(feature_id) = boost::spirit::qi::_1]);
      if (result && iter == iter_end && feature_id >= 0) {
	if (feature_id >= feature_names.size())
	  feature_names.resize(feature_id + 1);
	feature_names[feature_id] = piter->second;
      } else
	throw std::runtime_error("unsupported key: " + piter->first);
    }

    typedef rule_grammar_parser<std::string::const_iterator> rule_parser_type;
    
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
    rule_parsed_type rule;
    
    sequence_type source_index;

    while (std::getline(is, line)) {
      if (line.empty()) continue;
      
      boost::fusion::get<0>(rule).clear();
      boost::fusion::get<1>(rule).clear();
      boost::fusion::get<2>(rule).clear();
      boost::fusion::get<3>(rule).clear();
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator iter_end = line.end();

      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);
      
      if (! result || iter != iter_end)
	throw std::runtime_error("rule parsing failed: " + line);
      
      const symbol_type lhs(boost::fusion::get<0>(rule).empty() ? vocab_type::X : symbol_type(boost::fusion::get<0>(rule)));
      const rule_type::symbol_set_type source(boost::fusion::get<1>(rule).begin(), boost::fusion::get<1>(rule).end());
      const rule_type::symbol_set_type target(boost::fusion::get<2>(rule).begin(), boost::fusion::get<2>(rule).end());

#if 0
      // debug...
      std::cerr << "lhs: " << lhs;
      std::cerr << " source: ";
      std::copy(source.begin(), source.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << "target: ";
      std::copy(target.begin(), target.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << "features: ";
      std::copy(boost::fusion::get<3>(rule).begin(), boost::fusion::get<3>(rule).end(), std::ostream_iterator<std::pair<std::string, double> >(std::cerr, " "));
      std::cerr << std::endl;
#endif
      
      const size_type arity = source.arity();
      if (! target.empty() && arity != target.arity())
	throw std::runtime_error("arity do not match");
      
      rule_ptr_type rule_ptr(new rule_type(lhs, source, target, arity));
      rule_ptr->sort_source_index();
      
      source_index.resize(source.size());
      sequence_type::iterator iiter = source_index.begin();
      rule_type::symbol_set_type::const_iterator siter_end = source.end();
      for (rule_type::symbol_set_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter, ++ iiter)
	*iiter = siter->non_terminal();

      // we do not check duplicates!
      trie[trie.insert(source_index.begin(), source_index.end())].push_back(rule_ptr);
      
      // add features...
      const scores_parsed_type& scores = boost::fusion::get<3>(rule);
      
      int feature = 0;
      scores_parsed_type::const_iterator fiter_end = scores.end();
      for (scores_parsed_type::const_iterator fiter = scores.begin(); fiter != fiter_end; ++ fiter) 
	if (fiter->first.empty()) {
	  
	  if (feature < feature_names.size() && ! feature_names[feature].empty())
	    rule_ptr->features[feature_names[feature]] = fiter->second;
	  else {
	    // default name!
	    const std::string name = std::string("rule-table-") + boost::lexical_cast<std::string>(feature);
	    
	    rule_ptr->features[name] = fiter->second;
	  }
	  
	  ++ feature;
	} else
	  rule_ptr->features[fiter->first] = fiter->second;
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

  void GrammarMutable::insert(const rule_type& rule)
  {
    pimpl->insert(rule);
  }

  void GrammarMutable::insert(const rule_ptr_type& rule)
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
  
  const GrammarMutable::rule_set_type& GrammarMutable::rules(const id_type& node) const
  {
    static const rule_set_type __empty;
    return (node == pimpl->root() ? __empty : pimpl->rules(node));
  }
  
};
