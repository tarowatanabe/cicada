//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include "hypergraph.hpp"
#include "sort_topologically.hpp"
#include "unite.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/thread_specific_ptr.hpp"
#include "utils/json_string_parser.hpp"
#include "utils/json_string_generator.hpp"
#include "utils/dense_hash_map.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Rule,
			  (cicada::Rule::symbol_type,     lhs)
			  (cicada::Rule::symbol_set_type, rhs)
			  )

namespace cicada
{
  
  // statics...
  const HyperGraph::id_type HyperGraph::invalid = HyperGraph::id_type(-1);
  
  void HyperGraph::topologically_sort()
  {
    cicada::topologically_sort(*this);
  }
  
  void HyperGraph::unite(const HyperGraph& x)
  {
    cicada::unite(*this, x);
  }
  
  typedef std::string rule_parsed_type;
  typedef std::vector<rule_parsed_type> rule_parsed_set_type;
  
  typedef std::vector<int> tail_node_set_type;
  
  typedef std::pair<std::string, double> feature_parsed_type;
  typedef std::vector<feature_parsed_type> feature_parsed_set_type;
  
  typedef AttributeVector::data_type attribute_data_type;
  typedef std::pair<std::string, attribute_data_type> attribute_parsed_type;
  typedef std::vector<attribute_parsed_type> attribute_parsed_set_type;

  typedef boost::tuple<tail_node_set_type, feature_parsed_set_type, attribute_parsed_set_type, unsigned int> edge_parsed_type;
  
  typedef std::vector<edge_parsed_type> node_parsed_type;
  typedef std::vector<node_parsed_type> node_parsed_set_type;
  
  typedef boost::tuple<rule_parsed_set_type, node_parsed_set_type> hypergraph_parsed_type;
  
  template <typename Iterator>
  struct hypergraph_parser : boost::spirit::qi::grammar<Iterator, hypergraph_parsed_type(), boost::spirit::standard::space_type>
  {
    typedef HyperGraph hypergraph_type;
    typedef hypergraph_type::rule_type rule_type; 
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;
    typedef std::vector<rule_ptr_type, std::allocator<rule_ptr_type> > rule_ptr_set_type;
    
    hypergraph_parser() : hypergraph_parser::base_type(hypergraph)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      rule_string_action = rule_string [add_rule(rules)];
      rule_string_set = '[' >> -(rule_string_action % ',')  >> ']';
      
      tail_node_set %= '[' >> -(qi::uint_ % ',') >> ']';
      
      feature       %= feature_name >> qi::lit(':') >> qi::double_;
      feature_set   %= '{' >> -(feature % ',' )  >> '}';
      
      // attributes...
      data %= data_string | double_dot | int64_;
      attribute %= key >> ':' >> data;
      attribute_set %= '{' >> -(attribute % ',') >> '}';
      
      edge %= ('{'
	       >> -(qi::lit("\"tail\"")      >> ':' >> tail_node_set >> ',')
	       >> -(qi::lit("\"feature\"")   >> ':' >> feature_set >> ',')
	       >> -(qi::lit("\"attribute\"") >> ':' >> attribute_set >> ',')
	       >> qi::lit("\"rule\"")        >> ':' >> qi::uint_
	       >> -(',' >> qi::lit("\"first\"") >> ':' >> qi::omit[qi::int_]) // backward compatibility
	       >> -(',' >> qi::lit("\"last\"")  >> ':' >> qi::omit[qi::int_]) // backward compatibility
	       >> '}');
      
      edge_action = edge [add_edge(graph, rules)];
      
      // here, we will call add_node...
      node = (qi::lit('[')[add_node(graph)] >> -(edge_action % ',')  >> ']');
      
      node_set = ('[' >> -(node % ',') >> ']');
      
      hypergraph = ('{'
		    >> qi::lit("\"rules\"") >> ':' >> rule_string_set
		    >> ',' >> qi::lit("\"nodes\"") >> ':' >> node_set
		    >> -(',' >> qi::lit("\"goal\"")  >> ':' >> qi::uint_ [finish_graph(graph)])
		    >> '}');
    }
    
    struct finish_graph
    {
      finish_graph(hypergraph_type& __graph)
	: _graph(&__graph) {}

      void operator()(const int& goal, boost::spirit::qi::unused_type, boost::spirit::qi::unused_type) const
      {
	hypergraph_type& graph   = const_cast<hypergraph_type&>(*_graph);
	
	graph.goal = goal;
      }
      
      hypergraph_type* _graph;
    };

    
    struct add_node
    {
      add_node(hypergraph_type& __graph)
	: _graph(&__graph) {}
      
      void operator()(boost::spirit::qi::unused_type, boost::spirit::qi::unused_type, boost::spirit::qi::unused_type) const
      {
	hypergraph_type& graph   = const_cast<hypergraph_type&>(*_graph);
	graph.add_node();

	//std::cerr << "node: " << graph.nodes.back().id << std::endl;
      }
      
      hypergraph_type* _graph;
    };


    struct add_edge
    {
      add_edge(hypergraph_type& __graph,
	       rule_ptr_set_type& __rules)
	: _graph(&__graph), _rules(&__rules) {}

      void operator()(const edge_parsed_type& edge_parsed, boost::spirit::qi::unused_type, boost::spirit::qi::unused_type) const
      {
	hypergraph_type& graph   = const_cast<hypergraph_type&>(*_graph);
	rule_ptr_set_type& rules = const_cast<rule_ptr_set_type&>(*_rules);
	
	// here, we have not set-up in-edges (yet)
	hypergraph_type::edge_type& edge = graph.add_edge();
	
	//std::cerr << "edge: " << edge.id << std::endl;
	
	edge.tails = hypergraph_type::edge_type::node_set_type(boost::fusion::get<0>(edge_parsed).begin(), boost::fusion::get<0>(edge_parsed).end());
	edge.features.insert(boost::fusion::get<1>(edge_parsed).begin(), boost::fusion::get<1>(edge_parsed).end());
	edge.attributes.insert(boost::fusion::get<2>(edge_parsed).begin(), boost::fusion::get<2>(edge_parsed).end());
	edge.rule = rules[boost::fusion::get<3>(edge_parsed)];
	
	graph.connect_edge(edge.id, graph.nodes.back().id);
      }
      
      hypergraph_type* _graph;
      rule_ptr_set_type* _rules;
    };

    struct add_rule
    {
      add_rule(rule_ptr_set_type& __rules)
	: rules(&__rules) {}
      
      void operator()(const std::string& pattern, boost::spirit::qi::unused_type, boost::spirit::qi::unused_type) const
      {
	const_cast<rule_ptr_set_type&>(*rules).push_back(rule_type::create(rule_type(pattern)));

	//std::cerr << "add rule: " << *(*rules).back() << std::endl;
      }
      
      rule_ptr_set_type* rules;
    };
    
    hypergraph_type   graph;
    rule_ptr_set_type rules;

    typedef boost::spirit::standard::space_type space_type;
    
    utils::json_string_parser<Iterator> rule_string;
    
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), space_type>     rule_string_action;
    boost::spirit::qi::rule<Iterator, rule_parsed_set_type(), space_type> rule_string_set;
    
    boost::spirit::qi::rule<Iterator, tail_node_set_type(), space_type>      tail_node_set;
    
    utils::json_string_parser<Iterator> feature_name;
    boost::spirit::qi::rule<Iterator, feature_parsed_type(), space_type>     feature;
    boost::spirit::qi::rule<Iterator, feature_parsed_set_type(), space_type> feature_set;
    
    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1> int64_;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > double_dot;
    
    utils::json_string_parser<Iterator> key;
    utils::json_string_parser<Iterator> data_string;
    
    boost::spirit::qi::rule<Iterator, AttributeVector::data_type(), space_type> data;
    boost::spirit::qi::rule<Iterator, attribute_parsed_type(), space_type>      attribute;
    boost::spirit::qi::rule<Iterator, attribute_parsed_set_type(), space_type>  attribute_set;
    
    
    boost::spirit::qi::rule<Iterator, edge_parsed_type(), space_type> edge;
    boost::spirit::qi::rule<Iterator, edge_parsed_type(), space_type> edge_action;
    
    boost::spirit::qi::rule<Iterator, node_parsed_type(), space_type>        node;
    
    boost::spirit::qi::rule<Iterator, node_parsed_set_type(), space_type>    node_set;
    
    boost::spirit::qi::rule<Iterator, hypergraph_parsed_type(), space_type> hypergraph;
  };

  void HyperGraph::assign(const utils::piece& x)
  {
    std::string::const_iterator iter(x.begin());
    std::string::const_iterator end(x.end());
    
    const bool result = assign(iter, end);
    
    if (! result || iter != end)
      throw std::runtime_error("hypergraph format error");
  }

  namespace hypergraph_parser_impl
  {
    typedef hypergraph_parser<std::string::const_iterator> grammar_type;

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


  bool HyperGraph::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();

    hypergraph_parser_impl::grammar_type& grammar = hypergraph_parser_impl::instance();
    
    grammar.rules.clear();
    grammar.rules.push_back(HyperGraph::rule_ptr_type());
    grammar.graph.clear();
    
    const bool result = qi::phrase_parse(iter, end, grammar, standard::space);
    
    if (result)
      grammar.graph.swap(*this);

    grammar.rules.clear();
    grammar.graph.clear();
    
    return result;
  }
  
  std::istream& operator>>(std::istream& is, HyperGraph& x)
  {
    std::string line;
    
    x.clear();
    if (std::getline(is, line))
      x.assign(line);
    
    return is;
  }

  // do we use karma...?

  template <typename Iterator>
  struct rule_generator : boost::spirit::karma::grammar<Iterator, cicada::Rule()>
  {
    typedef cicada::Rule                 rule_type;
    typedef rule_type::symbol_type       symbol_type;
    typedef rule_type::symbol_set_type   symbol_set_type;
        
    rule_generator() : rule_generator::base_type(rule)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      lhs %= *(&standard::char_('\b') << "\\b"
	       | &standard::char_('\t') << "\\t"
	       | &standard::char_('\n') << "\\n"
	       | &standard::char_('\f') << "\\f"
	       | &standard::char_('\r') << "\\r"
	       | &standard::char_('\"') << "\\\""
	       | &standard::char_('\\') << "\\\\"
	       | &standard::char_('/') << "\\/"
	       | standard::char_);
      
      phrase %= -(lhs % ' ');
      
      rule %= lhs << " ||| " << phrase;
    }
    
    
    boost::spirit::karma::rule<Iterator, symbol_type()>      lhs;
    boost::spirit::karma::rule<Iterator, symbol_set_type()>  phrase;
    boost::spirit::karma::rule<Iterator, rule_type()>        rule;
  };

  typedef std::pair<Feature, double> feature_gen_type;
  typedef std::vector<feature_gen_type, std::allocator<feature_gen_type> > feature_generated_type;

  template <typename Iterator>
  struct features_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::feature_set_type()>
  {
    typedef cicada::HyperGraph::feature_set_type feature_set_type;
    //typedef feature_generated_type feature_set_type;
    
    features_generator() : features_generator::base_type(features)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      features %= (name << ':' << double10) % ',';
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 10;
      }
    };

    boost::spirit::karma::real_generator<double, real_precision> double10;
    utils::json_string_generator<Iterator>                       name;
    
    boost::spirit::karma::rule<Iterator, feature_set_type()> features;
  };

  namespace hypergraph_rule_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef rule_generator<iterator_type> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif
      
    static grammar_type& instance()
    {
#ifdef HAVE_TLD
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

  namespace hypergraph_feature_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    typedef features_generator<iterator_type> grammar_type;

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
  
  std::ostream& operator<<(std::ostream& os, const HyperGraph& graph)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
      
    typedef HyperGraph hypergraph_type;
    typedef hypergraph_type::rule_type rule_type;

    typedef google::dense_hash_map<const rule_type*, int, utils::hashmurmur<size_t>, std::equal_to<const rule_type*> > rule_unique_map_type;
    
    os << '{';
    
    rule_unique_map_type rules_unique;
    rules_unique.set_empty_key(0);
    
    {
      os << "\"rules\"" << ": " << '[';

      // dump rule part...

      hypergraph_rule_generator_impl::grammar_type& grammar = hypergraph_rule_generator_impl::instance();
      
      bool initial_rule = true;
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {

	const hypergraph_type::node_type& node = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  if (edge.rule) {
	    const rule_type& rule = *(edge.rule);

	    rule_unique_map_type::iterator riter = rules_unique.find(&rule);
	    
	    if (riter == rules_unique.end()) {
	      // + 1 for none-rule which will be zero-rule-id
	      const int rule_id = rules_unique.size() + 1;
	      
	      rules_unique.insert(std::make_pair(&rule, rule_id));
	      
	      if (! initial_rule)
		os << ", ";
	      os << '\"';
	      karma::generate(hypergraph_rule_generator_impl::iterator_type(os), grammar, rule);
	      os << '\"';
	      
	      initial_rule = false;
	    }
	  }
	}
      }
      
      os << ']';
    }
    
    os << ", ";
    
    {
      os << "\"nodes\"" << ": " << '[';
      
      hypergraph_feature_generator_impl::grammar_type& grammar = hypergraph_feature_generator_impl::instance();

      //feature_generated_type features;
      
      // dump nodes...
      bool initial_node = true;
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	if (! initial_node)
	  os << ", ";
	initial_node = false;
	
	os << '[';
	
	bool initial_edge = true;
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = niter->edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = niter->edges.begin(); eiter != eiter_end; ++ eiter) {
	  if (! initial_edge)
	    os << ", ";
	  initial_edge = false;
	  
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  
	  os << '{';

	  if (! edge.tails.empty()) {
	    typedef std::ostream_iterator<char> iterator_type;
	    
	    os << "\"tail\":[";
	    karma::generate(iterator_type(os), karma::uint_ % ',', edge.tails);
	    os << "],";
	  }
	  
#if 0
	  features.clear();
	  hypergraph_type::feature_set_type::const_iterator fiter_end = edge.features.end();
	  for (hypergraph_type::feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	    if (fiter->second != 0.0 && ! fiter->first.empty())
	      features.push_back(*fiter);
#endif
	  
	  if (! edge.features.empty()) {
	    os << "\"feature\":{";
	    karma::generate(hypergraph_feature_generator_impl::iterator_type(os), grammar, edge.features);
	    os << "},";
	  }
	  
	  if (! edge.attributes.empty())
	    os << "\"attribute\":" << edge.attributes << ',';
	  
	  os << "\"rule\":" << (! edge.rule ? 0 : rules_unique.find(&(*edge.rule))->second);
	  
	  os << '}';
	}
	
	os << ']';
      }
      
      os << ']';
    }
    
    // dump goal...
    if (graph.is_valid())
      os << ", \"goal\":" << graph.goal;
    
    os << '}';
    
    return os;
  }
  
};
