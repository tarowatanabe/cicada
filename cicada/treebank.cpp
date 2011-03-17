//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>


#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iterator>

#include "treebank.hpp"

#include "utils/thread_specific_ptr.hpp"


namespace cicada
{
  // TODO: escaping + de-normalizing...
  
  template <typename Iterator>
  struct treebank_label_generator : boost::spirit::karma::grammar<Iterator, std::string()>
  {
    treebank_label_generator() : treebank_label_generator::base_type(label)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      label %= +(&standard::char_('(') << "-LRB-"
		 | &standard::char_(')') << "-RRB-"
		 | &standard::char_('[') << "-LSB-"
		 | &standard::char_(']') << "-RSB-"
		 | &standard::char_('{') << "-LCB-"
		 | &standard::char_('}') << "-RCB-"
		 | &standard::char_('/') << "\\/"
		 | &standard::char_('*') << "\\*"
		 | standard::char_);
    };
    
    boost::spirit::karma::rule<Iterator, std::string()>     label;
  };

  namespace treebank_label_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef treebank_label_generator<iterator_type> grammar_type;
    
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
  

  struct TreeBankDump
  {
    typedef HyperGraph hypergraph_type;
    typedef hypergraph_type::rule_type rule_type;

    TreeBankDump(std::ostream& __os) : os(__os), label_generator(treebank_label_generator_impl::instance()) {}
    
    void operator()(const hypergraph_type& graph) const
    {
      if (graph.is_valid())
	traverse(graph, graph.goal);
      else
	os << "(())";
    }
    
    void traverse(const hypergraph_type& graph, const hypergraph_type::id_type& id) const
    {
      namespace karma = boost::spirit::karma;
      
      const hypergraph_type::node_type& node = graph.nodes[id];
      
      if (node.edges.empty()) return;
      
      const hypergraph_type::edge_type& edge = graph.edges[node.edges.front()];
      const rule_type& rule = *edge.rule;

      const utils::piece lhs_piece(static_cast<const std::string&>(rule.lhs));

      os << '(';
      if (lhs_piece == "[PERIOD]")
	os << '.';
      else if (lhs_piece == "[COMMA]")
	os << ',';
      else if (lhs_piece == "[COLON]")
	os << ':';
      else if (lhs_piece == "[SEMICOLON]")
	os << ';';
      else
	os << rule.lhs.non_terminal_strip();

      int non_terminal_pos = 0;
      rule_type::symbol_set_type::const_iterator riter_end = rule.rhs.end();
      for (rule_type::symbol_set_type::const_iterator riter = rule.rhs.begin(); riter != riter_end; ++ riter) {
	
	os << ' ';
	
	if (riter->is_non_terminal()) {
	  const int __non_terminal_index = riter->non_terminal_index();
	  const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	  ++ non_terminal_pos;
	  
	  traverse(graph, edge.tails[antecedent_index]);
	} else {
	  if (! karma::generate(treebank_label_generator_impl::iterator_type(os), label_generator, *riter))
	    throw std::runtime_error("invalid label");
	}
      }
      os << ')';
    }
    
    std::ostream& os;
    treebank_label_generator_impl::grammar_type& label_generator;
  };

  std::ostream& treebank(std::ostream& os, const HyperGraph& hypergraph)
  {
    if (! hypergraph.is_valid()) return os;
    
    TreeBankDump dumper(os);
    
    dumper(hypergraph);
    
    return os;
  }
};
