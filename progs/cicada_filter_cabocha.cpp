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

typedef std::vector<int, std::allocator<int> > dependency_type;
typedef std::vector<int, std::allocator<int> > offset_set_type;

typedef std::vector<int, std::allocator<int> > index_set_type;
typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_map_type;

typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;
typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
typedef sentence_type phrase_type;

path_type input_file = "-";
path_type output_file = "-";

std::string goal = "[s]";
std::string non_terminal = "[x]";

bool head_mode = false;
bool pos_mode = false;
bool leaf_mode = false;
bool dependency_mode = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
							
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file);

    node_set_type     nodes;
    dependency_type   dependency;
    offset_set_type   offsets;
    terminal_set_type terminals;

    hypergraph_type     hypergraph;
    dependency_map_type dependency_map;
    node_map_type       node_map;
    tail_set_type       tails;
    phrase_type         non_terminals;
    phrase_type         phrase;
    
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
	      
	      if (pos_mode)
		os << titer->first << "|[" << titer->second << "]";
	      else
		os << titer->first;
	      
	      initial = false;
	    }
	  }
	  os << '\n';
	} else {
	  // we will convert bunsetsu dependency into word-dependency...
	  dependency.clear();
	  offsets.clear();
	  terminals.clear();
	  
	  node_set_type::const_iterator niter_end = nodes.end();
	  for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	    const size_t head_pos = dependency.size() + niter->head;
	    
	    int pos = 0;
	    terminal_set_type::const_iterator titer_end = niter->terminals.end();
	    for (terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer, ++ pos) {
	      if (pos == niter->head) {
		offsets.push_back(dependency.size());
		dependency.push_back(-1);
	      } else 
		dependency.push_back(head_pos + 1);
	    }
	    
	    terminals.insert(terminals.end(), niter->terminals.begin(), niter->terminals.end());
	  }
	  
	  //
	  // second iteration to perform bunsets-wise dependency, but shifted by the bunsets length
	  //
	  int index = 0;
	  for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	    int pos = 0;
	    terminal_set_type::const_iterator titer_end = niter->terminals.end();
	    for (terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer, ++ pos) {
	      if (pos == niter->head) {
		if (niter->pos < 0)
		  dependency[index] = 0;
		else
		  dependency[index] = offsets[niter->pos] + 1;
	      }
	      
	      ++ index;
	    }
	  }

	  if (dependency_mode) {
	    if (! dependency.empty()) {
	      std::copy(dependency.begin(), dependency.end() - 1, std::ostream_iterator<int>(os, " "));
	      os << dependency.back();
	    }
	    os << '\n';
	  } else {
	    hypergraph.clear();
	    
	    if (dependency.empty())
	      os << hypergraph << '\n';
	    else {
	      dependency_map.clear();
	      dependency_map.resize(dependency.size() + 1);
	      
	      node_map.clear();
	      node_map.resize(dependency.size() + 1, hypergraph_type::invalid); 
	      
	      non_terminals.clear();
	      non_terminals.push_back(__goal);
	      
	      for (size_t i = 0; i != dependency.size(); ++ i) {
		dependency_map[dependency[i]].push_back(i + 1);
		
		if (pos_mode)
		  non_terminals.push_back('[' + terminals[i].second + ']');
		else
		  non_terminals.push_back(__non_terminal);
	      }
	      
	      if (! dependency_map[0].empty()) {
		tails.clear();
		phrase.clear();
		
		for (size_t i = 0; i != dependency_map[0].size(); ++ i) {
		  const int antecedent = dependency_map[0][i];
		  
		  if (node_map[antecedent] == hypergraph_type::invalid)
		    node_map[antecedent] = hypergraph.add_node().id;
		  
		  tails.push_back(node_map[antecedent]);
		  phrase.push_back(non_terminals[antecedent]);
		}
		
		if (node_map[0] == hypergraph_type::invalid)
		  node_map[0] = hypergraph.add_node().id;
		
		hypergraph_type::edge_type& edge = hypergraph.add_edge(tails.begin(), tails.end());
		edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(non_terminals[0], phrase.begin(), phrase.end()));
		
		hypergraph.connect_edge(edge.id, node_map[0]);
		hypergraph.goal = node_map[0];
	      }
	      
	      for (size_t id = 1; id != dependency_map.size(); ++ id) {
		tails.clear();
		phrase.clear();
		
		index_set_type::const_iterator iiter_begin = dependency_map[id].begin();
		index_set_type::const_iterator iiter_end   = dependency_map[id].end();
		index_set_type::const_iterator iiter_lex   = std::lower_bound(iiter_begin, iiter_end, id);
		
		for (index_set_type::const_iterator iiter = iiter_begin; iiter != iiter_lex; ++ iiter) {
		  const size_t antecedent = *iiter;
		  
		  if (node_map[antecedent] == hypergraph_type::invalid)
		    node_map[antecedent] = hypergraph.add_node().id;
		  
		  tails.push_back(node_map[antecedent]);
		  phrase.push_back(non_terminals[antecedent]);
		}
		
		if (head_mode) {
		  // create a new node and edge to terminal(s)!
		  
		  const symbol_type lhs = '[' + non_terminals[id].non_terminal_strip() + "*]";
		  tails.push_back(hypergraph.add_node().id);
		  phrase.push_back(lhs);
		  
		  hypergraph_type::edge_type& edge = hypergraph.add_edge();
		  edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(lhs, hypergraph_type::rule_type::symbol_set_type(1, terminals[id - 1].first)));
		  
		  hypergraph.connect_edge(edge.id, tails.back());
		} else
		  phrase.push_back(terminals[id - 1].first);
		
		for (index_set_type::const_iterator iiter = iiter_lex; iiter != iiter_end; ++ iiter) {
		  const size_t antecedent = *iiter;
		  
		  if (node_map[antecedent] == hypergraph_type::invalid)
		    node_map[antecedent] = hypergraph.add_node().id;
		  
		  tails.push_back(node_map[antecedent]);
		  phrase.push_back(non_terminals[antecedent]);
		}
		
		if (node_map[id] == hypergraph_type::invalid)
		  node_map[id] = hypergraph.add_node().id;
		
		const symbol_type& lhs = non_terminals[id];
			       
		hypergraph_type::edge_type& edge = hypergraph.add_edge(tails.begin(), tails.end());
		edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(lhs, phrase.begin(), phrase.end()));
		
		hypergraph.connect_edge(edge.id, node_map[id]);
	      }
	      
	      if (! hypergraph.nodes.empty() && hypergraph.is_valid())
		hypergraph.topologically_sort();
	      
	      os << hypergraph << '\n';
	    }
	  }
	}
	
	nodes.clear();
	// end of processing...
	
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
	
	if (poss.size() < 2) {
	  poss.resize(2);
	  poss[0] = "UNK";
	  poss[1] = "*";
	}
	
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
    
    ("head",       po::bool_switch(&head_mode),       "use non-terminal for head word")
    ("pos",        po::bool_switch(&pos_mode),        "use pos as non-terminal")
    ("leaf",       po::bool_switch(&leaf_mode),       "collect leaf nodes only")
    ("dependency", po::bool_switch(&dependency_mode), "collect dependency only")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

