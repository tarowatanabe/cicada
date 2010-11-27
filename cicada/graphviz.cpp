//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>


#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iterator>

#include "graphviz.hpp"


BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Rule,
			  (cicada::Rule::symbol_type,     lhs)
			  (cicada::Rule::symbol_set_type, rhs)
			  )

namespace cicada
{

  template <typename Iterator>
  struct graphviz_label_generator : boost::spirit::karma::grammar<Iterator, std::string()>
  {
    graphviz_label_generator() : graphviz_label_generator::base_type(label)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
      using karma::buffer;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('{', "\\{")
	('}', "\\}")
	('<', "\\<")
	('>', "\\>")
	('|', "\\|")
	(' ',  "\\ ")
	('/', "\\/")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t");
      
      label %= *(escape_char | ~char_('\"'));
    }
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    
    boost::spirit::karma::rule<Iterator, std::string()> label;
  };
  
  template <typename Iterator>
  struct graphviz_rule_generator : boost::spirit::karma::grammar<Iterator, cicada::Rule()>
  {
    typedef cicada::Rule                 rule_type;
    typedef rule_type::symbol_type       symbol_type;
    typedef rule_type::symbol_set_type   symbol_set_type;
    
    graphviz_rule_generator() : graphviz_rule_generator::base_type(rule)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
      using karma::buffer;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('{', "\\{")
	('}', "\\}")
	('<', "\\<")
	('>', "\\>")
	('|', "\\|")
	(' ',  "\\ ")
	('/', "\\/")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t");
      
      lhs %= *(escape_char | ~char_('\"'));
      phrase %= -(lhs % "\\ ");
      
      rule %= lhs << " | " << phrase;
    }
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    
    boost::spirit::karma::rule<Iterator, symbol_type()>      lhs;
    boost::spirit::karma::rule<Iterator, symbol_set_type()>  phrase;
    boost::spirit::karma::rule<Iterator, rule_type()>        rule;
  };

  template <typename Iterator>
  struct graphviz_tail_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::edge_type::node_set_type()>
  {
    
    typedef cicada::HyperGraph::edge_type::node_set_type node_set_type;

    graphviz_tail_generator() : graphviz_tail_generator::base_type(tail)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
      using karma::buffer;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      tail %= -(int_ % "\\ ");
    }
    
    boost::spirit::karma::rule<Iterator, node_set_type()> tail;
  };
  
  
  template <typename Iterator>
  struct graphviz_feature_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::feature_set_type()>
  {
    typedef cicada::Rule                 rule_type;
    typedef rule_type::symbol_type       symbol_type;
    typedef rule_type::symbol_set_type   symbol_set_type;

    typedef cicada::HyperGraph::feature_set_type  feature_set_type;
    typedef feature_set_type::value_type value_type;
    
    graphviz_feature_generator() : graphviz_feature_generator::base_type(features)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
      using karma::buffer;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('{', "\\{")
	('}', "\\}")
	('<', "\\<")
	('>', "\\>")
	('|', "\\|")
	(' ',  "\\ ")
	('/', "\\/")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t");
      
      // left adjusted newlines
      features %= -(((+(escape_char | ~char_('\"')) << ":\\ " << double_) % "\\l") << "\\l");
    }
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    boost::spirit::karma::rule<Iterator, feature_set_type()> features;
    
  };


  
  std::ostream& graphviz(std::ostream& os, const HyperGraph& hypergraph)
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef std::ostream_iterator<char> iterator_type;
    
    typedef graphviz_rule_generator<iterator_type>    rule_grammar_type;
    typedef graphviz_tail_generator<iterator_type>    tail_grammar_type;
    typedef graphviz_feature_generator<iterator_type> feature_grammar_type;
    
    rule_grammar_type    rule_grammar;
    tail_grammar_type    tail_grammar;
    feature_grammar_type feature_grammar;
    
    os << "digraph { rankdir=BT;" << '\n';
    
    hypergraph_type::node_set_type::const_iterator niter_end = hypergraph.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = hypergraph.nodes.begin(); niter != niter_end; ++ niter) {
      const node_type& node = *niter;
      
      os << " node_" << node.id << " [label=\"" << node.id << "\", shape=circle, height=0.1, width=0.1];" << '\n';
      
      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = hypergraph.edges[*eiter];

	if (edge.rule) {
	  os << "  edge_" << edge.id << " [label=\"{";
	  
	  iterator_type iter_rule(os);
	  boost::spirit::karma::generate(iter_rule, rule_grammar, *edge.rule);
	  
	  os << " | ";
	  
	  iterator_type iter_tail(os);
	  boost::spirit::karma::generate(iter_tail, tail_grammar, edge.tails);
	  
	  os << "}";
	  
	  if (! edge.features.empty()) {
	    os << " | ";
	    iterator_type iter_feature(os);
	    boost::spirit::karma::generate(iter_feature, feature_grammar, edge.features);
	  }
	  
	  os << "\", shape=record];" << '\n';
	} else
	  os << "  edge_" << edge.id << " [label=\"\", shape=rect];" << '\n';
	
	os << "    edge_" << edge.id << " -> node_" << node.id << ';' << '\n';
	edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	  os << "    node_" << *niter << " -> edge_" << edge.id << ';' << '\n';
	
      }
    }
    
    os << '}' << '\n';
    
    return os;
  }


  std::ostream& graphviz(std::ostream& os, const Lattice& lattice)
  {
    typedef Lattice lattice_type;
    
    typedef std::ostream_iterator<char> iterator_type;

    typedef graphviz_label_generator<iterator_type>   label_grammar_type;
    typedef graphviz_feature_generator<iterator_type> feature_grammar_type;

    label_grammar_type   label_grammar;
    feature_grammar_type feature_grammar;

    os << "digraph { " << '\n';

    int id_edge = 0;
    for (size_t id = 0; id != lattice.size(); ++ id) {
      os << " node_" << id << " [label=\"" << id << "\", shape=circle, height=0.1, width=0.1];" << '\n';
      
      lattice_type::arc_set_type::const_iterator aiter_end = lattice[id].end();
      for (lattice_type::arc_set_type::const_iterator aiter = lattice[id].begin(); aiter != aiter_end; ++ aiter) {
	const lattice_type::arc_type& arc = *aiter;
	
	os << "   edge_" << id_edge << " [label=\"";

	iterator_type iter_label(os);
	boost::spirit::karma::generate(iter_label, label_grammar, arc.label);
	
	if (! arc.features.empty()) {
	  os << " | ";
	  
	  iterator_type iter_feature(os);
	  boost::spirit::karma::generate(iter_feature, feature_grammar, arc.features);
	}
	
	os << "\", shape=record];" << '\n';
	
	os << "   node_" << id << " -> edge_" << id_edge << ';' << '\n';
	os << "   edge_" << id_edge << " -> node_" << (id + arc.distance) << ';' << '\n';
	++ id_edge;
      }
    }
    
    os << " node_" << lattice.size() << " [label=\"" << lattice.size() << "\", shape=circle, height=0.1, width=0.1];" << '\n';
    
    os << '}' << '\n';
    
    return os;
  }
};
