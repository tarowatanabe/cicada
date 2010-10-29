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

bool leaf = false;

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
    
    while (std::getline(is, line)) {
      tokenizer_type tokenizer(line);
      
      tokens.clear();
      tokens.insert(tokens.end(), tokenizer.begin(), tokenizer.end());
      
      if (tokens.empty()) continue;
      
      if (tokens.size() == 1) {
	if (tokens.front() != "EOS")
	  throw std::runtime_error("invalid cabocha F1 format: no EOS");
	
	
	if (leaf) {
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
	    
	    terminal_set_type::const_iterator titer_end = niter->terminals.end();
	    for (terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer) {
	      
	      hypergraph_type::node_type& node = graph.add_node();
	      
	      tails.push_back(node.id);
	      symbols.push_back('[' + titer->second + ']');
	      
	      hypergraph_type::edge_type& edge = graph.add_edge();
	      edge.rule.reset(new rule_type(symbols.back(), rule_type::symbol_set_type(1, titer->first)));
	      graph.connect_edge(edge.id, node.id);
	    }

	    niter->cat = symbols[niter->head];
	    
	    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
	    edge.rule.reset(new rule_type(niter->cat, rule_type::symbol_set_type(symbols.begin(), symbols.end())));
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
		symbols.push_back(nodes[id].cat);
		
		queue.push_back(std::make_pair(tails.back(), id));
	      } else if (id == node_id) {
		tails.push_back(nodes[id].id);
		symbols.push_back(nodes[id].cat);
	      }
	    }
	    
	    hypergraph_type::edge_type& edge = graph.add_edge(tails.begin(), tails.end());
	    edge.rule.reset(new rule_type(static_cast<hypergraph_type::id_type>(parent_id) == graph.goal ? symbol_type("[root]") : nodes[node_id].cat,
					  rule_type::symbol_set_type(symbols.begin(), symbols.end())));
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
	  throw std::runtime_error("invalid cabocha F1 format: no star");
	
	const int index = boost::lexical_cast<int>(tokens[1]);
	if (index != static_cast<int>(nodes.size()))
	  throw std::runtime_error("invalid cabocha F1 format: node size do not match");
	
	nodes.push_back(node_type(atoi(tokens[2].c_str()), atoi(tokens[3].c_str())));
	
      } else if (tokens.size() == 3) {
	boost::tokenizer<boost::char_separator<char> > tokenizer(tokens[1], boost::char_separator<char>(","));
	tokens_type poss(tokenizer.begin(), tokenizer.end());
	
	if (poss.size() < 2)
	  throw std::runtime_error("invalid cabocha F1 format: invalid POS");
	
	nodes.back().terminals.push_back(std::make_pair(tokens[0], poss[0] + '-' + poss[1]));
      } else
	throw std::runtime_error("invalid cabocha F1 format: # of columns do not match");
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
    ("leaf",      po::bool_switch(&leaf),    "collect leaf nodes only")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

