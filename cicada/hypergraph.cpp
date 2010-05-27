
#include "hypergraph.hpp"



#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/sgi_hash_map.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Rule,
			  (cicada::Rule::symbol_type,     lhs)
			  (cicada::Rule::symbol_set_type, source)
			  (cicada::Rule::symbol_set_type, target)
			  (cicada::Rule::feature_set_type, features)
			  )

namespace cicada
{

  // topologically sort...
  
  struct TopologicallySortImpl
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
 
    enum color_type {
      white,
      gray,
      black
    };

    struct dfs_type
    {
      id_type node;
      int edge;
      int tail;
      
      dfs_type(const id_type& _node, const int& _edge, const int& _tail) 
	: node(_node), edge(_edge), tail(_tail) {}
    };
    
    typedef std::vector<int, std::allocator<int> > reloc_set_type;
    typedef std::vector<color_type, std::allocator<color_type> > color_set_type;
    typedef std::vector<dfs_type, std::allocator<dfs_type> > stack_type;

    struct no_filter_edge
    {
      bool operator()(const edge_type& edge) const
      {
	return false;
      }
    };
    
    template <typename Filter>
    void operator()(hypergraph_type& x, Filter filter)
    {
      reloc_set_type reloc_node(x.nodes.size(), -1);
      reloc_set_type reloc_edge(x.edges.size(), -1);
      color_set_type color(x.nodes.size(), white);
      stack_type stack;
      
      stack.reserve(x.nodes.size());
      stack.push_back(dfs_type(x.goal, 0, 0));
      int node_count = 0;
      int edge_count = 0;
      
      while (! stack.empty()) {
	const dfs_type& dfs = stack.back();
	id_type node_id = dfs.node;
	int pos_edge = dfs.edge;
	int pos_tail = dfs.tail;
	
	stack.pop_back();
	
	const node_type* curr_node = &(x.nodes[node_id]);
	while (pos_edge != curr_node->in_edges.size()) {
	  const edge_type& curr_edge = x.edges[curr_node->in_edges[pos_edge]];
	  
	  if (pos_tail == curr_edge.tail_nodes.size() || filter(curr_edge)) {
	    // reach end: proceed to the next edge with pos_tail initialized
	    ++ pos_edge;
	    pos_tail = 0;
	    continue;
	  }
	  
	  const id_type tail_node = curr_edge.tail_nodes[pos_tail];
	  const color_type tail_color = color[tail_node];
	  
	  switch (tail_color) {
	  case white:
	    ++ pos_tail;
	    stack.push_back(dfs_type(node_id, pos_edge, pos_tail));
	    
	    node_id = tail_node;
	    curr_node = &(x.nodes[node_id]);
	    
	    color[node_id] = gray;
	    pos_edge = 0;
	    pos_tail = 0;
	    
	    break;
	  case black:
	    ++ pos_tail;
	    break;
	  case gray:
	    throw std::runtime_error("detected cycle!");
	    break;
	  }
	}
	
	color[node_id] = black;
	reloc_node[node_id] = node_count ++;
	for (int i = 0; i < curr_node->in_edges.size(); ++ i)
	  if (! filter(x.edges[curr_node->in_edges[i]]))
	    reloc_edge[curr_node->in_edges[i]] = edge_count ++;
      }
      
      // sorted graph!
      HyperGraph sorted;
      
      // construct edges...
      for (int i = 0; i < reloc_edge.size(); ++ i)
	if (reloc_edge[i] >= 0) {
	  const edge_type& edge_old = x.edges[i];

	  const id_type edge_id = sorted.edges.size();
	  
	  sorted.edges.push_back(edge_old);
	  
	  edge_type& edge_new = sorted.edges.back();
	  
	  edge_new.id = edge_id;
	  
	  edge_new.head_node = reloc_node[edge_new.head_node];
	  edge_type::node_set_type::iterator niter_end = edge_new.tail_nodes.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tail_nodes.begin(); niter != niter_end; ++ niter)
	    *niter = reloc_node[*niter];
	  
	  reloc_edge[i] = edge_id;
	}
      
      
      // construct reverse node-map ...
      reloc_set_type reloc_map_node(node_count, -1);
      for (int i = 0; i < x.nodes.size(); ++ i)
	if (reloc_node[i] >= 0)
	  reloc_map_node[reloc_node[i]] = i;
      
      // construct nodes
      for (int i = 0; i < reloc_map_node.size(); ++ i) {
	const node_type& node_old = x.nodes[reloc_map_node[i]];
	node_type& node_new = sorted.add_node();
	
	node_type::edge_set_type::const_iterator eiter_end = node_old.in_edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node_old.in_edges.begin(); eiter != eiter_end; ++ eiter)
	  if (reloc_edge[*eiter] >= 0)
	    node_new.in_edges.push_back(reloc_edge[*eiter]);
	
	eiter_end = node_old.out_edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node_old.out_edges.begin(); eiter != eiter_end; ++ eiter)
	  if (reloc_edge[*eiter] >= 0)
	    node_new.in_edges.push_back(reloc_edge[*eiter]);
      }
      
      sorted.goal = sorted.nodes.size() - 1;
      
      sorted.swap(x);
    }
  };
  
  void HyperGraph::topologically_sort()
  {
    TopologicallySortImpl()(*this, TopologicallySortImpl::no_filter_edge());
  }
  
  // comput union
  struct UnionImpl
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    typedef Rule   rule_type;
    
    void operator()(hypergraph_type& x, const hypergraph_type& y)
    {
      if (&x == &y || y.nodes.empty() || y.goal == hypergraph_type::invalid) return;
      if (x.nodes.empty() || x.goal == hypergraph_type::invalid) {
	x = y;
	return;
      }
      
      // check if we share the same goal... otherwise, we will create new edge and fill...
      //
      // 
      const symbol_type& goal_x = x.edges[x.nodes[x.goal].in_edges.front()].rule->lhs;
      const symbol_type& goal_y = y.edges[y.nodes[y.goal].in_edges.front()].rule->lhs;

      if (goal_x == goal_y) {
	const int y_node_offset = x.nodes.size();
	const int y_edge_offset = x.edges.size();
	
	// -1 to adjust goal
	x.nodes.resize(x.nodes.size() + y.nodes.size() - 1);
	x.edges.resize(x.edges.size() + y.edges.size());
	
	// merge nodes
	for (id_type id = 0; id < y.nodes.size(); ++ id)
	  if (id != y.goal) {
	    const id_type id_new = id + y_node_offset - (id > y.goal);
	    
	    const node_type& node_old = y.nodes[id];
	    node_type& node_new = x.nodes[id_new];
	    
	    node_new = node_old;
	    
	    node_new.id =  id_new;
	    node_type::edge_set_type::iterator eiter_end = node_new.in_edges.end();
	    for (node_type::edge_set_type::iterator eiter = node_new.in_edges.begin(); eiter != eiter_end; ++ eiter)
	      *eiter += y_edge_offset;
	    
	    eiter_end = node_new.out_edges.end();
	    for (node_type::edge_set_type::iterator eiter = node_new.out_edges.begin(); eiter != eiter_end; ++ eiter)
	      *eiter += y_edge_offset;
	  }
	
	// merge edges
	for (id_type id = 0; id < y.edges.size(); ++ id) {
	  const edge_type& edge_old = y.edges[id];
	  edge_type& edge_new = x.edges[id + y_edge_offset];
	  
	  edge_new = edge_old;
	  edge_new.id = id + y_edge_offset;
	  
	  if (edge_new.head_node == y.goal) {
	    edge_new.head_node = x.goal;
	    x.nodes[x.goal].in_edges.push_back(edge_new.id);
	  } else
	    edge_new.head_node += y_node_offset - (edge_new.head_node > y.goal);
	  
	  edge_type::node_set_type::iterator niter_end = edge_new.tail_nodes.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tail_nodes.begin(); niter != niter_end; ++ niter)
	    *niter += y_node_offset - (*niter > y.goal);
	}
      } else {
	const int y_node_offset = x.nodes.size();
	const int y_edge_offset = x.edges.size();
	
	// +1 to adjust goal
	x.nodes.resize(x.nodes.size() + y.nodes.size() + 1);
	// +2 to adjust edges toward goal
	x.edges.resize(x.edges.size() + y.edges.size() + 2);
	
	// merge nodes
	for (id_type id = 0; id < y.nodes.size(); ++ id) {
	  const id_type id_new = id + y_node_offset;
	  
	  const node_type& node_old = y.nodes[id];
	  node_type& node_new = x.nodes[id_new];
	  
	  node_new = node_old;
	    
	  node_new.id =  id_new;
	  node_type::edge_set_type::iterator eiter_end = node_new.in_edges.end();
	  for (node_type::edge_set_type::iterator eiter = node_new.in_edges.begin(); eiter != eiter_end; ++ eiter)
	    *eiter += y_edge_offset;
	  
	  eiter_end = node_new.out_edges.end();
	  for (node_type::edge_set_type::iterator eiter = node_new.out_edges.begin(); eiter != eiter_end; ++ eiter)
	    *eiter += y_edge_offset;
	}
	
	// merge edges
	for (id_type id = 0; id < y.edges.size(); ++ id) {
	  const edge_type& edge_old = y.edges[id];
	  edge_type& edge_new = x.edges[id + y_edge_offset];
	  
	  edge_new = edge_old;
	  edge_new.id = id + y_edge_offset;
	  
	  edge_new.head_node += y_node_offset;
	  edge_type::node_set_type::iterator niter_end = edge_new.tail_nodes.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tail_nodes.begin(); niter != niter_end; ++ niter)
	    *niter += y_node_offset;
	}
	
	// create new node...
	node_type& goal_node = x.nodes.back();
	edge_type& goal_edge_x = x.edges[x.edges.size() - 2];
	edge_type& goal_edge_y = x.edges[x.edges.size() - 1];
	
	goal_node.id = x.nodes.size() - 1;
	goal_edge_x.id = x.edges.size() - 2;
	goal_edge_y.id = x.edges.size() - 1;

	goal_edge_x.rule.reset(new rule_type(vocab_type::GOAL,
					     rule_type::symbol_set_type(1, goal_x.non_terminal(1)),
					     rule_type::symbol_set_type(1, goal_x.non_terminal(1)),
					     1));
	goal_edge_y.rule.reset(new rule_type(vocab_type::GOAL,
					     rule_type::symbol_set_type(1, goal_y.non_terminal(1)),
					     rule_type::symbol_set_type(1, goal_y.non_terminal(1)),
					     1));
	
	goal_edge_x.head_node = goal_node.id;
	goal_edge_y.head_node = goal_node.id;
	
	goal_edge_x.tail_nodes = edge_type::node_set_type(1, x.goal);
	goal_edge_y.tail_nodes = edge_type::node_set_type(1, y.goal + y_node_offset);
	
	x.nodes[x.goal].out_edges.push_back(goal_edge_x.id);
	x.nodes[y.goal + y_node_offset].out_edges.push_back(goal_edge_y.id);
	
	goal_node.in_edges.push_back(goal_edge_x.id);
	goal_node.in_edges.push_back(goal_edge_y.id);
	
      }
    }
  };
  
  void HyperGraph::unite(const HyperGraph& x)
  {
    UnionImpl()(*this, x);
  }
  
  typedef std::string rule_parsed_type;
  typedef std::vector<rule_parsed_type> rule_parsed_set_type;
  
  typedef std::vector<int> tail_node_set_type;
  typedef std::pair<std::string, double> feature_parsed_type;
  typedef std::vector<feature_parsed_type> feature_parsed_set_type;
  typedef boost::tuple<tail_node_set_type, feature_parsed_set_type, int> edge_parsed_type;
  
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
      
      using qi::phrase_parse;
      using qi::lexeme;
      using qi::omit;
      using qi::repeat;
      using qi::lit;
      using qi::eps;
      using qi::inf;
      using qi::attr;
      using standard::char_;
      using qi::double_;
      using qi::int_;
      
      using namespace qi::labels;
      
      using standard::space;
      
      escape_char.add
	("\\\"", '\"')
	("\\\\", '\\')
	("\\/", '/')
	("\\b", '\b')
	("\\f", '\f')
	("\\n", '\n')
	("\\r", '\r')
	("\\t", '\t');
      
      rule_string %= ('\"' >> lexeme[*(escape_char | ~char_('\"'))] >> '\"') ;
      rule_string_action = rule_string [add_rule(rules)];
      rule_string_set = '[' >> -(rule_string_action % ',')  >> ']';
      
      tail_node_set %= '[' >> -(int_ % ',') >> ']';
      
      feature       %= '\"' >> lexeme[*(escape_char | ~char_('\"'))] >> lit('\"') >> lit(':') >> double_;
      feature_set   %= '{' >> -(feature % ',' )  >> '}';
      
      edge %= ('{' >> lit("\"tail\"")    >> ':' >> tail_node_set >> ','
	       >> lit("\"feature\"") >> ':' >> feature_set >> ','
	       >> lit("\"rule\"")    >> ':' >> int_ >> '}');
      
      edge_action = edge [add_edge(*graph, rules)];
      
      // here, we will call add_node...
      node = (lit('[')[add_node(*graph)] >> -(edge_action % ',')  >> ']');
      
      node_set = ('[' >> -(node % ',') >> ']');
      
      hypergraph = ('{'
		    >> lit("\"rules\"") >> ':' >> rule_string_set >> ','
		    >> lit("\"nodes\"") >> ':' >> node_set >> ','
		    >> lit("\"goal\"") >> ':' >> int_ [finish_graph(*graph)]
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
	
	// fix-up out-edges...
	hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	  
	  hypergraph_type::edge_type::node_set_type::const_iterator niter_end = eiter->tail_nodes.end();
	  for (hypergraph_type::edge_type::node_set_type::const_iterator niter = eiter->tail_nodes.begin(); niter != niter_end; ++ niter)
	    graph.nodes[*niter].out_edges.push_back(eiter->id);
	}
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
	
	edge.tail_nodes = hypergraph_type::edge_type::node_set_type(boost::fusion::get<0>(edge_parsed).begin(), boost::fusion::get<0>(edge_parsed).end());
	edge.features.insert(boost::fusion::get<1>(edge_parsed).begin(), boost::fusion::get<1>(edge_parsed).end());
	edge.rule = rules[boost::fusion::get<2>(edge_parsed)];
	
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
	const_cast<rule_ptr_set_type&>(*rules).push_back(rule_ptr_type(new rule_type(pattern)));
      }
      
      rule_ptr_set_type* rules;
    };
    
    hypergraph_type* graph;
    rule_ptr_set_type rules;

    typedef boost::spirit::standard::space_type space_type;
    
    
    boost::spirit::qi::symbols<char, char> escape_char;
    
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), space_type>     rule_string;
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), space_type>     rule_string_action;
    
    boost::spirit::qi::rule<Iterator, rule_parsed_set_type(), space_type> rule_string_set;
    
    boost::spirit::qi::rule<Iterator, tail_node_set_type(), space_type>      tail_node_set;
    boost::spirit::qi::rule<Iterator, feature_parsed_type(), space_type>     feature;
    boost::spirit::qi::rule<Iterator, feature_parsed_set_type(), space_type> feature_set;
    
    boost::spirit::qi::rule<Iterator, edge_parsed_type(), space_type> edge;
    boost::spirit::qi::rule<Iterator, edge_parsed_type(), space_type> edge_action;
    
    boost::spirit::qi::rule<Iterator, node_parsed_type(), space_type>        node;
    
    boost::spirit::qi::rule<Iterator, node_parsed_set_type(), space_type>    node_set;
    
    boost::spirit::qi::rule<Iterator, hypergraph_parsed_type(), space_type> hypergraph;
  };

  
  std::istream& operator>>(std::istream& is, HyperGraph& x)
  {
    typedef hypergraph_parser<std::string::const_iterator> grammar_type;
    
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
    
    std::string line;
    
    
    x.clear();
    if (std::getline(is, line)) {
      grammar.rules.clear();
      grammar.rules.push_back(HyperGraph::rule_ptr_type());
      grammar.graph = &x;
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space);
      if (! result || iter != end) {
	x.clear();
	grammar.rules.clear();
      }
    }
    
    return is;
  }

  // do we use karma...?
  
  template <typename Iterator>
  struct rule_generator : boost::spirit::karma::grammar<Iterator, cicada::Rule()>
  {
    typedef cicada::Rule                rule_type;
    typedef rule_type::symbol_type      symbol_type;
    typedef rule_type::symbol_set_type  symbol_set_type;
    typedef rule_type::feature_set_type feature_set_type;
    
    rule_generator() : rule_generator::base_type(rule)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('/', "\\/")
	('\b', "\\b")
	('\f', "\\f")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t");
      
      lhs %= *(escape_char | ~char_('\"'));
      phrase %= -(lhs % ' ');
      
      features %= (+(escape_char | ~char_('\"')) << '=' << double_) % ' ';
      
      rule %= lhs << " ||| " << phrase << " ||| " << phrase << -( " ||| " << features);
    }
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    
    boost::spirit::karma::rule<Iterator, symbol_type()>      lhs;
    boost::spirit::karma::rule<Iterator, symbol_set_type()>  phrase;
    boost::spirit::karma::rule<Iterator, feature_set_type()> features;
    boost::spirit::karma::rule<Iterator, rule_type()>        rule;
  };

  template <typename Iterator>
  struct features_generator : boost::spirit::karma::grammar<Iterator, cicada::HyperGraph::feature_set_type()>
  {
    typedef cicada::Rule                rule_type;
    typedef rule_type::feature_set_type feature_set_type;
    
    features_generator() : features_generator::base_type(features)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('/', "\\/")
	('\b', "\\b")
	('\f', "\\f")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t");
      
      features %= ('\"' << (+(escape_char | ~char_('\"')) << '\"' << ':' << double_) % ", ");
    }
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    boost::spirit::karma::rule<Iterator, feature_set_type()> features;
  };
  
  std::ostream& operator<<(std::ostream& os, const HyperGraph& graph)
  {
    typedef HyperGraph hypergraph_type;
    typedef hypergraph_type::rule_type rule_type;

#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<const rule_type*, int, utils::hashmurmur<size_t>, std::equal_to<const rule_type*>,
      std::allocator<std::pair<const rule_type*, int> > > rule_unique_map_type;
#else
    typedef sfi::hash_map<const rule_type*, int, utils::hashmurmur<size_t>, std::equal_to<const rule_type*>,
			  std::allocator<std::pair<const rule_type*, int> > > rule_unique_map_type;
#endif
    
    typedef std::back_insert_iterator<std::string> iterator_type;
    
    os << '{';

    rule_unique_map_type rules_unique;
    
    {
      os << "\"rules\"" << ": " << '[';

      // dump rule part...
      
      rule_generator<iterator_type>     rule_grammar;
      
      
      bool initial_rule = true;
      std::string output_rule;
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) 
	if (eiter->rule) {
	  const rule_type& rule = *(eiter->rule);
	  
	  rule_unique_map_type::iterator riter = rules_unique.find(&rule);
	  if (riter == rules_unique.end()) {
	    
	    // + 1 for none-rule which will be zero-rule-id
	    const int rule_id = rules_unique.size() + 1;
	    
	    rules_unique.insert(std::make_pair(&rule, rule_id));
	    
	    output_rule.clear();
	    iterator_type iter(output_rule);
	    
	    boost::spirit::karma::generate(iter, rule_grammar, rule);
	    
	    if (! initial_rule)
	      os << ", ";
	    os << '\"' << output_rule << '\"';
	    
	    initial_rule = false;
	  }
	}
      os << ']';
    }
    
    os << ", ";
    
    {
      os << "\"nodes\"" << ": " << '[';

      features_generator<iterator_type> features_grammar;
      std::string output_features;
      
      // dump nodes...
      bool initial_node = true;
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	if (! initial_node)
	  os << ", ";
	initial_node = false;

	os << '[';
	
	bool initial_edge = true;
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = niter->in_edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = niter->in_edges.begin(); eiter != eiter_end; ++ eiter) {
	  if (! initial_edge)
	    os << ", ";
	  initial_edge = false;

	  const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	  
	  os << "{\"tail\": [";
	  if (! edge.tail_nodes.empty()) {
	    std::copy(edge.tail_nodes.begin(), edge.tail_nodes.end() - 1, std::ostream_iterator<hypergraph_type::id_type>(os, ","));
	    os << edge.tail_nodes.back();
	  }
	  os << "], ";
	  os << "\"feature\": {";
	  // dump features!
	  
	  output_features.clear();
	  iterator_type iter(output_features);
	  boost::spirit::karma::generate(iter, features_grammar, edge.features);
	  
	  os << output_features;
	  
	  os << "}, ";
	  os << "\"rule\": ";
	  os << (! edge.rule ? 0 : rules_unique.find(&(*edge.rule))->second);
	  os << '}';
	}
	
	os << ']';
      }
      
      os << ']';
    }
    
    os << ", ";
    
    {
      // dump goal...
      os << "\"goal\": " << graph.goal;
    }
    os << '}';
    
    return os;
  }
  
};
