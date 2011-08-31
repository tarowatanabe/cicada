//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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

#include "graphviz.hpp"

#include "utils/thread_specific_ptr.hpp"

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
    graphviz_label_generator() : graphviz_label_generator::base_type(string)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      string = *(&standard::char_('\t') << "\\t"
		 | &standard::char_('\n') << "\\n"
		 | &standard::char_('\r') << "\\r"
		 | &standard::char_('\"') << "\\\""
		 | &standard::char_('\\') << "\\\\"
		 | &standard::char_('{') << "\\{"
		 | &standard::char_('}') << "\\}"
		 | &standard::char_('<') << "\\<"
		 | &standard::char_('>') << "\\>"
		 | &standard::char_('|') << "\\|"
		 | &standard::char_(' ') << "\\ "
		 | &standard::char_('/') << "\\/"
		 | standard::char_);
    }    
    boost::spirit::karma::rule<Iterator, std::string()> string;
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
      
      phrase %= -(lhs % "\\ ");
      
      rule %= lhs << " | " << phrase;
    }
    
    graphviz_label_generator<Iterator> lhs;
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
      
      tail %= -(karma::int_ % "\\ ");
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
      
      // left adjusted newlines
      features %= -(((string << ":\\ " << karma::double_) % "\\l") << "\\l");
    }
    
    graphviz_label_generator<Iterator> string;
    boost::spirit::karma::rule<Iterator, feature_set_type()> features;
  };

  template <typename Iterator>
  struct graphviz_attribute_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::attribute_set_type()>
  {
    typedef cicada::Rule                 rule_type;
    typedef rule_type::symbol_type       symbol_type;
    typedef rule_type::symbol_set_type   symbol_set_type;

    typedef cicada::HyperGraph::attribute_set_type  attribute_set_type;
    typedef attribute_set_type::value_type value_type;
    
    graphviz_attribute_generator() : graphviz_attribute_generator::base_type(attributes)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      // left adjusted newlines
      attributes %= -(((string << ":\\ " << (int64_ | double10 | "\\\"" << string << "\\\"")) % "\\l") << "\\l");
    }
    
    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 10;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double10;
    boost::spirit::karma::int_generator<AttributeVector::int_type, 10, false> int64_;

    graphviz_label_generator<Iterator> string;
    boost::spirit::karma::rule<Iterator, attribute_set_type()> attributes;
  };
  
  
  namespace graphviz_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    typedef graphviz_label_generator<iterator_type>     grammar_label_type;
    typedef graphviz_rule_generator<iterator_type>      grammar_rule_type;
    typedef graphviz_tail_generator<iterator_type>      grammar_tail_type;
    typedef graphviz_feature_generator<iterator_type>   grammar_feature_type;
    typedef graphviz_attribute_generator<iterator_type> grammar_attribute_type;
    
#ifdef HAVE_TLS
    static __thread grammar_label_type*     __grammar_label_tls = 0;
    static __thread grammar_rule_type*      __grammar_rule_tls = 0;
    static __thread grammar_tail_type*      __grammar_tail_tls = 0;
    static __thread grammar_feature_type*   __grammar_feature_tls = 0;
    static __thread grammar_attribute_type* __grammar_attribute_tls = 0;
    
    static boost::thread_specific_ptr<grammar_label_type >     __grammar_label;
    static boost::thread_specific_ptr<grammar_rule_type >      __grammar_rule;
    static boost::thread_specific_ptr<grammar_tail_type >      __grammar_tail;
    static boost::thread_specific_ptr<grammar_feature_type >   __grammar_feature;
    static boost::thread_specific_ptr<grammar_attribute_type > __grammar_attribute;
#else
    static utils::thread_specific_ptr<grammar_label_type >     __grammar_label;
    static utils::thread_specific_ptr<grammar_rule_type >      __grammar_rule;
    static utils::thread_specific_ptr<grammar_tail_type >      __grammar_tail;
    static utils::thread_specific_ptr<grammar_feature_type >   __grammar_feature;
    static utils::thread_specific_ptr<grammar_attribute_type > __grammar_attribute;
#endif
    
    static grammar_label_type& instance_label()
    {
#ifdef HAVE_TLS
      if (! __grammar_label_tls) {
	__grammar_label.reset(new grammar_label_type());
	__grammar_label_tls = __grammar_label.get();
      }
      return * __grammar_label_tls;
#else
      if (! __grammar_label.get())
	__grammar_label.reset(new grammar_label_type());
      
      return *__grammar_label;
#endif
    }

    static grammar_rule_type& instance_rule()
    {
#ifdef HAVE_TLS
      if (! __grammar_rule_tls) {
	__grammar_rule.reset(new grammar_rule_type());
	__grammar_rule_tls = __grammar_rule.get();
      }
      return * __grammar_rule_tls;
#else
      if (! __grammar_rule.get())
	__grammar_rule.reset(new grammar_rule_type());
      
      return *__grammar_rule;
#endif
    }
    
    static grammar_tail_type& instance_tail()
    {
#ifdef HAVE_TLS
      if (! __grammar_tail_tls) {
	__grammar_tail.reset(new grammar_tail_type());
	__grammar_tail_tls = __grammar_tail.get();
      }
      return * __grammar_tail_tls;
#else
      if (! __grammar_tail.get())
	__grammar_tail.reset(new grammar_tail_type());
      
      return *__grammar_tail;
#endif
    }

    static grammar_feature_type& instance_feature()
    {
#ifdef HAVE_TLS
      if (! __grammar_feature_tls) {
	__grammar_feature.reset(new grammar_feature_type());
	__grammar_feature_tls = __grammar_feature.get();
      }
      return * __grammar_feature_tls;
#else
      if (! __grammar_feature.get())
	__grammar_feature.reset(new grammar_feature_type());
      
      return *__grammar_feature;
#endif
    }

    static grammar_attribute_type& instance_attribute()
    {
#ifdef HAVE_TLS
      if (! __grammar_attribute_tls) {
	__grammar_attribute.reset(new grammar_attribute_type());
	__grammar_attribute_tls = __grammar_attribute.get();
      }
      return * __grammar_attribute_tls;
#else
      if (! __grammar_attribute.get())
	__grammar_attribute.reset(new grammar_attribute_type());
      
      return *__grammar_attribute;
#endif
    }

  };

  
  std::ostream& graphviz(std::ostream& os, const HyperGraph& hypergraph)
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;

    typedef graphviz_impl::iterator_type iterator_type;
    
    graphviz_impl::grammar_rule_type&    grammar_rule        = graphviz_impl::instance_rule();
    graphviz_impl::grammar_tail_type&    grammar_tail        = graphviz_impl::instance_tail();
    graphviz_impl::grammar_feature_type& grammar_feature     = graphviz_impl::instance_feature();
    graphviz_impl::grammar_attribute_type& grammar_attribute = graphviz_impl::instance_attribute();
    
    os << "digraph { rankdir=BT; ordering=in;";
    
    hypergraph_type::node_set_type::const_iterator niter_end = hypergraph.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = hypergraph.nodes.begin(); niter != niter_end; ++ niter) {
      const node_type& node = *niter;
      
      os << " node_" << node.id << " [label=\"" << node.id << "\", shape=circle, height=0.1, width=0.1];";
      
      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const edge_type& edge = hypergraph.edges[*eiter];

	if (edge.rule) {
	  os << " edge_" << edge.id << " [label=\"{";
	  
	  boost::spirit::karma::generate(iterator_type(os), grammar_rule, *edge.rule);
	  if (! edge.tails.empty()) {
	    os << " | ";
	    boost::spirit::karma::generate(iterator_type(os), grammar_tail, edge.tails);
	  }
	  
	  os << "}";
	  
	  if (! edge.features.empty()) {
	    os << " | ";
	    boost::spirit::karma::generate(iterator_type(os), grammar_feature, edge.features);
	  }

	  if (! edge.attributes.empty()) {
	    os << " | ";
	    boost::spirit::karma::generate(iterator_type(os), grammar_attribute, edge.attributes);
	  }
	  
	  os << "\", shape=record];";
	} else
	  os << " edge_" << edge.id << " [label=\"\", shape=rect];";
	
	os << " edge_" << edge.id << " -> node_" << node.id << ';';
	edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	  os << " node_" << *niter << " -> edge_" << edge.id << ';';
	
      }
    }
    
    os << '}';
    
    return os;
  }


  std::ostream& graphviz(std::ostream& os, const Lattice& lattice)
  {
    typedef Lattice lattice_type;
    
    typedef graphviz_impl::iterator_type iterator_type;
    
    graphviz_impl::grammar_label_type&   grammar_label   = graphviz_impl::instance_label();
    graphviz_impl::grammar_feature_type& grammar_feature = graphviz_impl::instance_feature();

    os << "digraph { ordering=out;";
    
    int id_edge = 0;
    for (size_t id = 0; id != lattice.size(); ++ id) {
      os << " node_" << id << " [label=\"" << id << "\", shape=circle, height=0.1, width=0.1];";
      
      lattice_type::arc_set_type::const_iterator aiter_end = lattice[id].end();
      for (lattice_type::arc_set_type::const_iterator aiter = lattice[id].begin(); aiter != aiter_end; ++ aiter) {
	const lattice_type::arc_type& arc = *aiter;
	
	os << " edge_" << id_edge << " [label=\"";
	
	boost::spirit::karma::generate(iterator_type(os), grammar_label, arc.label);
	if (! arc.features.empty()) {
	  os << " | ";
	  boost::spirit::karma::generate(iterator_type(os), grammar_feature, arc.features);
	}
	
	os << "\", shape=record];";
	
	os << " node_" << id << " -> edge_" << id_edge << ';';
	os << " edge_" << id_edge << " -> node_" << (id + arc.distance) << ';';
	++ id_edge;
      }
    }
    
    os << " node_" << lattice.size() << " [label=\"" << lattice.size() << "\", shape=circle, height=0.1, width=0.1];";
    
    os << '}';
    
    return os;
  }
};
