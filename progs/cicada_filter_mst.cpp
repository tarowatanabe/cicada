//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// from MST dependency to hypergraph conversion...
//
// we may use POS or semantic-role as our label...
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

// this is important!
#define FUSION_MAX_VECTOR_SIZE 15

#include <boost/variant.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>

#include <cicada/hypergraph.hpp>
#include <cicada/sentence.hpp>
#include <cicada/vocab.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

struct mst_type
{
  typedef size_t size_type;
  typedef std::string label_type;
  

  typedef std::vector<label_type, std::allocator<label_type> > label_set_type;
  typedef std::vector<size_type, std::allocator<size_type> >   position_set_type;
  
  
  label_set_type words;
  label_set_type poss;
  label_set_type labels;
  position_set_type positions;

  bool verify() const
  {
    return words.size() == poss.size() && words.size() == labels.size() && words.size() == positions.size();
  }
  
  void clear()
  {
    words.clear();
    poss.clear();
    labels.clear();
    positions.clear();
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  mst_type,
			  (mst_type::label_set_type,    words)
			  (mst_type::label_set_type,    poss)
			  (mst_type::label_set_type,    labels)
			  (mst_type::position_set_type, positions)
			  )

template <typename Iterator>
struct mst_parser : boost::spirit::qi::grammar<Iterator, mst_type(), boost::spirit::standard::blank_type>
{
  mst_parser() : mst_parser::base_type(mst)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    token %= qi::lexeme[+(standard::char_ - standard::space)];
    labels %= (+token);
    
    positions %= (+qi::int_);
    
    mst %= (labels >> qi::eol
	    >> labels >> qi::eol
	    >> labels >> qi::eol
	    >> positions >> qi::eol
	    >> (qi::eol || qi::eoi));
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type> token;
  boost::spirit::qi::rule<Iterator, mst_type::label_set_type(), blank_type> labels;
  boost::spirit::qi::rule<Iterator, mst_type::position_set_type(), blank_type> positions;
  
  boost::spirit::qi::rule<Iterator, mst_type(), blank_type> mst;
};


typedef boost::filesystem::path path_type;
typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Vocab      vocab_type;
typedef cicada::Symbol     word_type;
typedef cicada::Symbol     symbol_type;

typedef std::vector<mst_type::size_type, std::allocator<mst_type::size_type> > index_set_type;
typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_type;

typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;
typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
typedef sentence_type phrase_type;

path_type input_file = "-";
path_type output_file = "-";

std::string goal = "[s]";
std::string non_terminal = "[x]";
bool head_mode = false;
bool pos_mode = false;
bool relation_mode = false;
bool leaf_mode = false;
bool dependency_mode = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (pos_mode && relation_mode)
      throw std::runtime_error("either pos or relation or none");

    if (leaf_mode && dependency_mode)
      throw std::runtime_error("either leaf or dependency");
    
    typedef boost::spirit::istream_iterator iter_type;
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file);
    is.unsetf(std::ios::skipws);
    
    mst_parser<iter_type> parser;
    
    mst_type        mst;
    hypergraph_type hypergraph;
    dependency_type dependency;
    node_map_type   node_map;

    symbol_type __goal(symbol_type(goal).non_terminal());
    symbol_type __non_terminal(symbol_type(non_terminal).non_terminal());

    if (! __goal.is_non_terminal())
      throw std::runtime_error(goal + " is not a non-terminal");
    
    if (! __non_terminal.is_non_terminal())
      throw std::runtime_error(non_terminal + " is not a non-terminal");
    
    tail_set_type   tails;
    phrase_type     phrase;
    phrase_type     non_terminals;
    
    std::string line;
    iter_type iter(is);
    iter_type iter_end;
    
    int num = 0;
    for (;;) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      if (debug)
	std::cerr << "parsing: " << num << std::endl;
      
      mst.clear();
      
      if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, mst))
	throw std::runtime_error("parsing failed");

      if (! mst.verify())
	throw std::runtime_error("invalid mst format");
      
      if (debug >= 2)
	std::cerr << "size: " << mst.words.size() << std::endl;
      
      ++ num;
      
      if (leaf_mode) {
	if (! mst.words.empty()) {
	  if (pos_mode) {
	    for (size_t i = 0; i != mst.words.size() - 1; ++ i)
	      os << mst.words[i] << "|[" << mst.labels[i] << "] ";
	    os << mst.words.back() << "|[" << mst.labels.back() << ']';
	  } else {
	    std::copy(mst.words.begin(), mst.words.end() - 1, std::ostream_iterator<mst_type::label_type>(os, " "));
	    os << mst.words.back();
	  }
	}
	os << '\n';
	
	continue;
      }
      
      if (dependency_mode) {
	if (! mst.words.empty()) {
	  std::copy(mst.positions.begin(), mst.positions.end() - 1, std::ostream_iterator<mst_type::size_type>(os, " "));
	  os << mst.positions.back();
	}
	os << '\n';
	
	continue;
      }

      
      hypergraph.clear();
      
      if (mst.words.empty()) {
	os << hypergraph << '\n';
	continue;
      }
      
      dependency.clear();
      dependency.resize(mst.words.size() + 1);

      node_map.clear();
      node_map.resize(mst.words.size() + 1, hypergraph_type::invalid);
      
      non_terminals.clear();
      non_terminals.push_back(__goal);
      
      for (size_t i = 0; i != mst.words.size(); ++ i) {
	dependency[mst.positions[i]].push_back(i + 1);
	
	if (pos_mode)
	  non_terminals.push_back('[' + mst.poss[i] + ']');
	else if (relation_mode)
	  non_terminals.push_back('[' + mst.labels[i] + ']');
	else
	  non_terminals.push_back(__non_terminal);
      }
      
      if (! dependency[0].empty()) {
	tails.clear();
	phrase.clear();
	
	for (size_t i = 0; i != dependency[0].size(); ++ i) {
	  const size_t antecedent = dependency[0][i];
	  
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
      
      for (size_t id = 1; id != dependency.size(); ++ id) {
	tails.clear();
	phrase.clear();
	
	index_set_type::const_iterator iiter_begin = dependency[id].begin();
	index_set_type::const_iterator iiter_end   = dependency[id].end();
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
	  edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(lhs, hypergraph_type::rule_type::symbol_set_type(1, mst.words[id - 1])));
	  
	  hypergraph.connect_edge(edge.id, tails.back());
	} else
	  phrase.push_back(mst.words[id - 1]);
	
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
    ("relation",   po::bool_switch(&relation_mode),   "use relation as non-terminal")
    ("leaf",       po::bool_switch(&leaf_mode),       "collect leaf nodes only")
    ("dependency", po::bool_switch(&dependency_mode), "collect dependency only")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
