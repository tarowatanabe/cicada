//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/numeric/conversion/bounds.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/array_power2.hpp"
#include "utils/thread_specific_ptr.hpp"

#include "tree_rule.hpp"

struct cicada_treebank_type
{
  typedef std::vector<cicada_treebank_type, std::allocator<cicada_treebank_type> > antecedents_type;
  
  std::string label;
  antecedents_type antecedents;
  
  cicada_treebank_type() {}
  cicada_treebank_type(const std::string& __label) : label(__label) {}
  
  void clear()
  {
    label.clear();
    antecedents.clear();
  }
  
  void transform(cicada::TreeRule& x) const
  {
    // pre-order traversal...
    x.label = label;
    x.antecedents = cicada::TreeRule::antecedent_set_type(antecedents.size());
    
    for (size_t i = 0; i != x.antecedents.size(); ++ i)
      antecedents[i].transform(x.antecedents[i]);
  }
  
};


BOOST_FUSION_ADAPT_STRUCT(
			  cicada_treebank_type,
			  (std::string, label)
			  (cicada_treebank_type::antecedents_type, antecedents)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::TreeRule,
			  (cicada::TreeRule::label_type, label)
			  (cicada::TreeRule::antecedent_set_type, antecedents)
			  )

namespace cicada
{

  typedef utils::array_power2<TreeRule::rule_ptr_type, 1024 * 32, std::allocator<TreeRule::rule_ptr_type> > cache_type;
  
#ifdef HAVE_TLS
  static __thread cache_type* __tree_rule_cache_tls = 0;
  static boost::thread_specific_ptr<cache_type> __tree_rule_cache;
#else
  static utils::thread_specific_ptr<cache_type> __tree_rule_cache;
#endif
  
  TreeRule::rule_ptr_type TreeRule::create(const TreeRule& x)
  {
#ifdef HAVE_TLS
    if (! __tree_rule_cache_tls) {
      __tree_rule_cache.reset(new cache_type());
      __tree_rule_cache_tls = __tree_rule_cache.get();
    }
    cache_type& cache = *__tree_rule_cache_tls;
#else
    if (! __tree_rule_cache.get())
      __tree_rule_cache.reset(new cache_type());
    cache_type& cache = *__tree_rule_cache;
#endif
    
    const size_t cache_pos = hash_value(x) & (cache.size() - 1);
    if (! cache[cache_pos] || *cache[cache_pos] != x)
      cache[cache_pos].reset(new TreeRule(x));
    
    return cache[cache_pos];
  }

  
  template <typename Iterator>
  struct tree_rule_parser_grammar : boost::spirit::qi::grammar<Iterator, cicada_treebank_type(), boost::spirit::standard::space_type>
  {
    tree_rule_parser_grammar() : tree_rule_parser_grammar::base_type(root)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
    
      escaped_char.add
	("\\\\", '\\')
	("\\(", '(')
	("\\)", ')');
      
      label %= qi::lexeme[+(escaped_char | (standard::char_ - standard::space - '(' - ')' - '\\')) - "|||"];
      tree_rule %= label >> (qi::hold['(' >> +tree_rule >> ')'] | qi::eps);
      root %= -tree_rule;
    }
    
    boost::spirit::qi::symbols<char, char>        escaped_char;
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type>   label;
    boost::spirit::qi::rule<Iterator, cicada_treebank_type(), boost::spirit::standard::space_type> tree_rule;
    boost::spirit::qi::rule<Iterator, cicada_treebank_type(), boost::spirit::standard::space_type> root;
  };

  template <typename Iterator>
  struct tree_rule_generator_grammar : boost::spirit::karma::grammar<Iterator, TreeRule()>
  {
    tree_rule_generator_grammar() : tree_rule_generator_grammar::base_type(tree_rule)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
    
      escaped_char.add
	('\\', "\\\\")
	('(',  "\\(")
	(')',  "\\)");
      
      label %= *(escaped_char | standard::char_);
      antecedents %= tree_rule % ' ';
      tree_rule %= label << (karma::buffer['(' << antecedents << ')'] | karma::eps);
    }
    
    boost::spirit::karma::symbols<char, const char*>                       escaped_char;
    boost::spirit::karma::rule<Iterator, TreeRule::label_type()>           label;
    boost::spirit::karma::rule<Iterator, TreeRule::antecedent_set_type()>  antecedents;
    boost::spirit::karma::rule<Iterator, TreeRule()>                       tree_rule;
  };
  
  namespace tree_rule_parser_impl
  {
    typedef tree_rule_parser_grammar<std::string::const_iterator> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif
    
    static grammar_type& instance()
    {
#ifdef HAVE_TLS
      if (! __grammar_tls) {
	__grammar.reset(new grammar_type());
	__grammar_tls = __grammar.get();
      }
      
      return *__grammar_tls;
#else
      if (! __grammar.get())
	__grammar.reset(new grammar_type());
      
      return *__grammar;
#endif
    }
  };


  bool TreeRule::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    cicada_treebank_type treebank;
    
    if (qi::phrase_parse(iter, end, tree_rule_parser_impl::instance(), standard::space, treebank)) {
      treebank.transform(*this);
      return true;
    } else
      return false;
  }
 
  void TreeRule::assign(const std::string& x)
  {
    clear();
    
    if (x.empty()) return;

    std::string::const_iterator iter = x.begin();
    std::string::const_iterator end = x.end();
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("tree rule format parsing failed..." + x);
  }
  
  std::istream& operator>>(std::istream& is, TreeRule& x)
  {
    x.clear();
    
    std::string line;
    if (std::getline(is, line) && ! line.empty())
      x.assign(line);
    
    return is;
  }

  namespace tree_rule_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef tree_rule_generator_grammar<iterator_type> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif

    static grammar_type& instance()
    {
#ifdef HAVE_TLS
      if (! __grammar_tls) {
	__grammar.reset(new grammar_type());
	__grammar_tls = __grammar.get();
      }
      
      return *__grammar_tls;
#else
      if (! __grammar.get())
	__grammar.reset(new grammar_type());
      
      return *__grammar;
#endif
    }
  };

  std::ostream& operator<<(std::ostream& os, const TreeRule& x)
  {
    namespace karma = boost::spirit::karma;
    
    tree_rule_generator_impl::iterator_type iter(os);
    
    if (! karma::generate(iter, tree_rule_generator_impl::instance(), x))
      throw std::runtime_error("failed tree-rule generation!");
    
    return os;
  }
  
  template <typename Index>
  inline
  void sort_source(TreeRule& rule, Index& index, int& pos)
  {
    if (rule.label.is_non_terminal() && rule.antecedents.empty()) {
      int non_terminal_pos = rule.label.non_terminal_index();
      if (non_terminal_pos == 0)
	non_terminal_pos = pos;
      
      if (non_terminal_pos >= static_cast<int>(index.size()))
	index.resize(non_terminal_pos + 1);
      index[non_terminal_pos] = pos;
      
      rule.label = rule.label.non_terminal(pos);
      
      ++ pos;
    }
    
    TreeRule::antecedent_set_type::iterator iter_end = rule.antecedents.end();
    for (TreeRule::antecedent_set_type::iterator iter = rule.antecedents.begin(); iter != iter_end; ++ iter)
      sort_source(*iter, index, pos);
  }
  
  template <typename Index>
  inline
  void sort_target(TreeRule& rule, const Index& index, int& pos)
  {
    if (rule.label.is_non_terminal() && rule.antecedents.empty()) {
      if (pos >= static_cast<int>(index.size()))
	throw std::runtime_error("sort tree failed: output index");
      
      rule.label = rule.label.non_terminal(index[pos]);
      ++ pos;
    }
    
    TreeRule::antecedent_set_type::iterator iter_end = rule.antecedents.end();
    for (TreeRule::antecedent_set_type::iterator iter = rule.antecedents.begin(); iter != iter_end; ++ iter)
      sort_target(*iter, index, pos);
  }

  void sort(TreeRule& x, TreeRule& y)
  {
    typedef std::vector<int, std::allocator<int> > index_type;
    
    // we will sort in pre-order traversal...
    index_type index;
    
    int pos = 1;
    sort_source(x, index, pos);
    
    pos = 1;
    sort_target(y, index, pos);
  }
};
