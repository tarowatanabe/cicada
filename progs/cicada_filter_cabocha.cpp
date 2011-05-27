//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "cicada/hypergraph.hpp"
#include "cicada/sentence.hpp"

#include "utils/program_options.hpp"
#include "utils/space_separator.hpp"
#include "utils/lexical_cast.hpp"

typedef boost::filesystem::path path_type;
typedef boost::tokenizer<utils::space_separator> tokenizer_type;
typedef std::vector<std::string, std::allocator<std::string> > tokens_type;

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Rule       rule_type;
typedef cicada::Symbol     symbol_type;

typedef std::pair<std::string, std::string> terminal_type;
typedef std::vector<terminal_type, std::allocator<terminal_type> > terminal_set_type;

struct node_type
{
  node_type(int __pos, int __head) : id(hypergraph_type::invalid), pos(__pos), head(__head) {}

  hypergraph_type::id_type id;
  
  int pos;
  int head;
  
  symbol_type cat;
  terminal_set_type terminals;
};
typedef std::vector<node_type, std::allocator<node_type> > node_set_type;

struct equal_pos
{
  equal_pos(int __pos) : pos(__pos) {}
  
  bool operator()(const node_type& x) const
  {
    return x.pos == pos;
  }

  int pos;
};


path_type input_file = "-";
path_type output_file = "-";

std::string goal = "[s]";
std::string non_terminal = "[x]";

bool head_mode = false;
bool pos_mode = false;
bool leaf_mode = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
							
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file);

    node_set_type nodes;
 
    hypergraph_type graph;
    
    std::string line;
    tokens_type tokens;

    symbol_type __goal(symbol_type(goal).non_terminal());
    symbol_type __non_terminal(symbol_type(non_terminal).non_terminal());
    
    if (! __goal.is_non_terminal())
      throw std::runtime_error(goal + " is not a non-terminal");
    
    if (! __non_terminal.is_non_terminal())
      throw std::runtime_error(non_terminal + " is not a non-terminal");
    
    size_t lineno = -1;
    while (std::getline(is, line)) {
      ++ lineno;
      tokenizer_type tokenizer(line);
      
      tokens.clear();
      tokens.insert(tokens.end(), tokenizer.begin(), tokenizer.end());
      
      if (tokens.empty()) continue;
      
      if (tokens.size() == 1) {
	if (tokens.front() != "EOS")
	  throw std::runtime_error("invalid cabocha F1 format: no EOS: " + utils::lexical_cast<std::string>(lineno));
	
	
	if (leaf_mode) {
	  bool initial = true;
	  node_set_type::const_iterator niter_end = nodes.end();
	  for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	    
	    terminal_set_type::const_iterator titer_end = niter->terminals.end();
	    for (terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer) {
	      if (!initial)
		os << ' ';
	      
	      os << titer->first;
	      initial = false;
	    }
	  }
	  os << '\n';
	} else {
	  graph.clear();

	  std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tails;
	  std::vector<symbol_type, std::allocator<symbol_type> > symbols;
	
	  // handle terminals...
	  node_set_type::iterator niter_end = nodes.end();
	  for (node_set_type::iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	    
	    hypergraph_type::node_type& head = graph.add_node();
	    
	    niter->id = head.id;
	    tails.clear();
	    symbols.clear();

	    if (head_mode) {
	      int pos = 0;
	      terminal_set_type::const_iterator titer_end = niter->terminals.end();
	      for (terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer, ++ pos) {
		
		hypergraph_type::node_type& node = graph.add_node();
		
		tails.push_back(node.id);
		symbols.push_back('[' + titer->second + "*]");
		
		hypergraph_type::edge_type& edge = graph.add_edge();
		edge.rule = rule_type::create(rule_type(symbols.back(), rule_type::symbol_set_type(1, titer->first)));
		graph.connect_edge(edge.id, node.id);
		
		if (pos != niter->head) {
		  hypergraph_type::node_type& node = graph.add_node();

		  const symbol_type symbol_new = '[' + titer->second + ']';
		  
		  hypergraph_type::edge_type& edge = graph.add_edge(&tails.back(), (&tails.back()) + 1);
		  edge.rule = rule_type::create(rule_type(symbol_new, rule_type::symbol_set_type(1, symbols.back())));
		  graph.connect_edge(edge.id, node.id);
		  
		  tails.back() = node.id;
		  symbols.back() = symbol_new;
		}
	      }
	    } else {
	      terminal_set_type::const_iterator titer_end = niter->terminals.end();
	      for (terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer) {
		
		hypergraph_type::node_type& node = graph.add_node();
		
		tails.push_back(node.id);
		symbols.push_back('[' + titer->second + ']');
		
		hypergraph_type::edge_type& edge = graph.add_edge();
		edge.rule = rule_type::create(rule_type(symbols.back(), rule_type::symbol_set_type(1, titer->first)));
		graph.connect_edge(edge.id, node.id);
	      }
	    }
	    
	    niter->cat = '[' + niter->terminals[niter->head].second + ']';
	    
	    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
	    edge.rule = rule_type::create(rule_type(niter->cat, symbols.begin(), symbols.end()));
	    graph.connect_edge(edge.id, head.id);
	  }
	  
	  // handle dependency...
	  // start from root...
	  
	  node_set_type::const_iterator niter = std::find_if(nodes.begin(), nodes.end(), equal_pos(-1));
	  if (niter == nodes.end())
	    throw std::runtime_error("invalid nodes!");
	  
	  // breadth first search...
	  typedef std::deque<std::pair<int, int>, std::allocator<std::pair<int, int> > > queue_type;

	  hypergraph_type::node_type& root = graph.add_node();
	  graph.goal = root.id;

	  queue_type queue;
	  queue.push_back(std::make_pair(root.id, niter - nodes.begin()));
	  
	  while (! queue.empty()) {
	    const int parent_id = queue.front().first;
	    const int node_id = queue.front().second;
	    queue.pop_front();

	    tails.clear();
	    symbols.clear();
	    for (int id = 0; id < static_cast<int>(nodes.size()); ++ id) {
	      if (nodes[id].pos == node_id) {
		tails.push_back(graph.add_node().id);
		symbols.push_back(pos_mode ? nodes[id].cat : __non_terminal);
		
		queue.push_back(std::make_pair(tails.back(), id));
	      } else if (id == node_id) {
		tails.push_back(nodes[id].id);
		symbols.push_back(pos_mode ? nodes[id].cat : __non_terminal);
	      }
	    }
	    
	    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
	    edge.rule = rule_type::create(rule_type(static_cast<hypergraph_type::id_type>(parent_id) == graph.goal
						    ? __goal
						    : (pos_mode ? nodes[node_id].cat : __non_terminal),
						    symbols.begin(), symbols.end()));
	    graph.connect_edge(edge.id, parent_id);
	  }
	  
	  if (! graph.nodes.empty() && graph.goal != hypergraph_type::invalid)
	    graph.topologically_sort();

	  os << graph << '\n';
	}
	
	nodes.clear();
	// end of processing...
	
	// we will dump!
	
      } else if (tokens.size() == 5) {
	if (tokens.front() != "*")
	  throw std::runtime_error("invalid cabocha F1 format: no star: " + utils::lexical_cast<std::string>(lineno));
	
	const int index = utils::lexical_cast<int>(tokens[1]);
	if (index != static_cast<int>(nodes.size()))
	  throw std::runtime_error("invalid cabocha F1 format: node size do not match: " + utils::lexical_cast<std::string>(lineno));
	
	nodes.push_back(node_type(atoi(tokens[2].c_str()), atoi(tokens[3].c_str())));
	
      } else if (tokens.size() == 3) {
	boost::tokenizer<boost::char_separator<char> > tokenizer(tokens[1], boost::char_separator<char>(","));
	tokens_type poss(tokenizer.begin(), tokenizer.end());
	
	if (poss.size() < 2)
	  throw std::runtime_error("invalid cabocha F1 format: invalid POS: " + utils::lexical_cast<std::string>(lineno));
	
	if (poss[1] == "*")
	  nodes.back().terminals.push_back(std::make_pair(tokens[0], poss[0]));
	else
	  nodes.back().terminals.push_back(std::make_pair(tokens[0], poss[0] + '-' + poss[1]));
      } else
	throw std::runtime_error("invalid cabocha F1 format: # of columns do not match: " + utils::lexical_cast<std::string>(lineno));
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
  

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")
    
    ("goal",         po::value<std::string>(&goal)->default_value(goal),                 "goal symbol")
    ("non-terminal", po::value<std::string>(&non_terminal)->default_value(non_terminal), "non-terminal symbol")
    
    ("head",    po::bool_switch(&head_mode), "use non-terminal for head word")
    ("pos",     po::bool_switch(&pos_mode),  "use pos as non-terminal")
    ("leaf",    po::bool_switch(&leaf_mode), "collect leaf nodes only")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

