//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include "hypergraph.hpp"
#include "sort_topologically.hpp"
#include "unite.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <boost/functional/hash/hash.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/array_power2.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/thread_specific_ptr.hpp"
#include "utils/json_string_parser.hpp"
#include "utils/json_string_generator.hpp"
#include "utils/compact_map.hpp"
#include "utils/getline.hpp"

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
      namespace phoenix = boost::phoenix;
      
      rule_string_action = rule_string [add_rule(rules)];
      rule_string_set = '[' >> -(rule_string_action % ',')  >> ']';
      
      tail_node_set %= '[' >> -(qi::uint_ % ',') >> ']';
      
      feature       %= feature_name >> qi::lit(':') >> qi::double_;
      feature_set   %= '{' >> -(feature % ',' )  >> '}';
      
      // attributes...
      data = (data_string [qi::_val = data_string_cache(qi::_1)]
              | data_double [qi::_val = qi::_1]
	      | data_int [qi::_val = qi::_1]);
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
	
	edge.features.rehash(boost::fusion::get<1>(edge_parsed).size());
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


    struct attribute_string_cache_type : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      typedef std::string string_type;
      
      typedef utils::array_power2<string_type, 1024 * 4, std::allocator<string_type> > string_set_type;
      
      template <typename >
      struct result { typedef const string_type& type; };
      
      const string_type& operator()(const string_type& x) const
      {
	const size_t pos = hasher_type::operator()(x.begin(), x.end(), 0) & (caches.size() - 1);
	
	string_type& ret = const_cast<string_type&>(caches[pos]);
	if (ret != x)
	  ret = x;
	
	return ret;
      }
      
      string_set_type caches;
    };

    
    utils::json_string_parser<Iterator> rule_string;
    
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), space_type>     rule_string_action;
    boost::spirit::qi::rule<Iterator, rule_parsed_set_type(), space_type> rule_string_set;
    
    boost::spirit::qi::rule<Iterator, tail_node_set_type(), space_type> tail_node_set;
    
    utils::json_string_parser<Iterator>                                      feature_name;
    boost::spirit::qi::rule<Iterator, feature_parsed_type(), space_type>     feature;
    boost::spirit::qi::rule<Iterator, feature_parsed_set_type(), space_type> feature_set;
    
    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1>                      data_int;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > data_double;
    utils::json_string_parser<Iterator>                                                      data_string;
    boost::phoenix::function<attribute_string_cache_type>                                    data_string_cache;
    
    utils::json_string_parser<Iterator>                                         key;
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
    static utils::thread_specific_ptr<grammar_type > __grammar;
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
    if (utils::getline(is, line))
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

  template <typename Iterator>
  struct features_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::feature_set_type()>
  {
    typedef cicada::HyperGraph::feature_set_type feature_set_type;
    
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

  template <typename Iterator>
  struct attributes_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::attribute_set_type()>
  {
    typedef cicada::HyperGraph::attribute_set_type attribute_set_type;
    
    attributes_generator() : attributes_generator::base_type(attributes)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      attribute_data %= data_int | data_double | data_string;
      attribute %= attribute_key << ':' << attribute_data;
      
      attributes %= attribute % ',';
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 10;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision>                 data_double;
    boost::spirit::karma::int_generator<attribute_set_type::int_type, 10, false> data_int;
    utils::json_string_generator<Iterator>                                       data_string;
    
    utils::json_string_generator<Iterator>                                 attribute_key;
    boost::spirit::karma::rule<Iterator, attribute_set_type::data_type()>  attribute_data;
    boost::spirit::karma::rule<Iterator, attribute_set_type::value_type()> attribute;
    boost::spirit::karma::rule<Iterator, attribute_set_type()>             attributes;
  };
  

  namespace hypergraph_rule_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef rule_generator<iterator_type> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static utils::thread_specific_ptr<grammar_type > __grammar;
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
    static utils::thread_specific_ptr<grammar_type > __grammar;
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

  namespace hypergraph_attribute_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    
    typedef attributes_generator<iterator_type> grammar_type;

#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static utils::thread_specific_ptr<grammar_type > __grammar;
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
  
  template <typename Tp>
  struct ptr_unassigned
  {
    const Tp* operator()() const { return 0; }
  };

  template <typename Tp>
  struct ptr_hash
  {
    size_t operator()(const Tp* x) const
    {
      return (x ? hash_value(*x) : size_t(0));
    }
  };
  
  template <typename Tp>
  struct ptr_equal
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return (x == y) || (x && y && *x == *y);
    }
  };

  
  std::ostream& operator<<(std::ostream& os, const HyperGraph& graph)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    typedef std::ostream_iterator<char> iterator_type;
    
    typedef HyperGraph hypergraph_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef utils::compact_map<const rule_type*, int,
			       ptr_unassigned<rule_type>, ptr_unassigned<rule_type>,
			       ptr_hash<rule_type>, ptr_equal<rule_type> > rule_unique_map_type;
    
    karma::generate(iterator_type(os), '{');
    
    rule_unique_map_type rules_unique(graph.edges.size());
    
    {
      karma::generate(iterator_type(os), karma::lit("\"rules\": ["));
      
      // dump rule part...

      hypergraph_rule_generator_impl::grammar_type& grammar = hypergraph_rule_generator_impl::instance();
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {

	const hypergraph_type::node_type& node = *niter;
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  if (edge.rule) {
	    const rule_type& rule = *(edge.rule);
	    
	    std::pair<rule_unique_map_type::iterator, bool> result = rules_unique.insert(std::make_pair(&rule,
													rules_unique.size() + 1));
	    
	    if (result.second) {
	      if (rules_unique.size() != 1)
		karma::generate(iterator_type(os), karma::lit(", "));
	      
	      karma::generate(hypergraph_rule_generator_impl::iterator_type(os),
			      '\"' << grammar << '\"',
			      rule);
	    }
	  }
	}
      }
      
      karma::generate(iterator_type(os), ']');
    }
    
    karma::generate(iterator_type(os), karma::lit(", "));
    
    {
      karma::generate(iterator_type(os), karma::lit("\"nodes\": ["));
      
      hypergraph_feature_generator_impl::grammar_type&   grammar_feature   = hypergraph_feature_generator_impl::instance();
      hypergraph_attribute_generator_impl::grammar_type& grammar_attribute = hypergraph_attribute_generator_impl::instance();
      
      // dump nodes...
      bool initial_node = true;
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	if (! initial_node)
	  karma::generate(iterator_type(os), karma::lit(", "));
	initial_node = false;
	
	karma::generate(iterator_type(os), '[');
	
	bool initial_edge = true;
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = niter->edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = niter->edges.begin(); eiter != eiter_end; ++ eiter) {
	  if (! initial_edge)
	    karma::generate(iterator_type(os), karma::lit(", "));
	  initial_edge = false;
	  
	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  karma::generate(iterator_type(os), '{');
	  
	  if (! edge.tails.empty())
	    karma::generate(iterator_type(os),
			    "\"tail\":[" << (karma::uint_ % ',') << "],",
			    edge.tails);
	  
	  if (! edge.features.empty())
	    karma::generate(hypergraph_feature_generator_impl::iterator_type(os),
			    "\"feature\":{" << grammar_feature << "},",
			    edge.features);
	  
	  if (! edge.attributes.empty())
	    karma::generate(hypergraph_attribute_generator_impl::iterator_type(os),
			    "\"attribute\":{" << grammar_attribute << "},",
			    edge.attributes);
	  
	  karma::generate(iterator_type(os), "\"rule\":" << karma::uint_, ! edge.rule ? 0 : rules_unique.find(&(*edge.rule))->second);
	  
	  karma::generate(iterator_type(os), '}');
	}
	
	karma::generate(iterator_type(os), ']');
      }
      
      karma::generate(iterator_type(os), ']');
    }
    
    // dump goal...
    if (graph.is_valid())
      karma::generate(iterator_type(os), karma::lit(", \"goal\": ") << karma::uint_generator<HyperGraph::id_type>(), graph.goal);
    
    karma::generate(iterator_type(os), '}');
    
    return os;
  }
  
};
