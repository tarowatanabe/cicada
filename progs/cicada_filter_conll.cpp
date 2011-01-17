//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// from conll dependency to hypergraph conversion...
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

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>

#include <cicada/hypergraph.hpp>
#include <cicada/sentence.hpp>
#include <cicada/vocab.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

struct conll_type
{
  typedef size_t size_type;
  typedef boost::variant<size_type, std::string> phead_type;

  struct visitor_phead : public boost::static_visitor<size_type>
  {
    size_type operator()(const size_type& x) const { return x; }
    size_type operator()(const std::string& x) const { return size_type(-1); }
  };

  size_type   id;
  std::string form;
  std::string lemma;
  std::string cpostag;
  std::string postag;
  std::string feats;
  size_type   head;
  std::string deprel;
  phead_type  phead;
  std::string pdeprel;

  conll_type() {}
  conll_type(const size_type&   __id,
	     const std::string& __form,
	     const std::string& __lemma,
	     const std::string& __cpostag,
	     const std::string& __postag,
	     const std::string& __feats,
	     const size_type&   __head,
	     const std::string& __deprel,
	     const phead_type&  __phead,
	     const std::string& __pdeprel)
    : id(__id),
      form(__form),
      lemma(__lemma),
      cpostag(__cpostag),
      feats(__feats),
      head(__head),
      deprel(__deprel),
      phead(__phead),
      pdeprel(__pdeprel) {}
};

BOOST_FUSION_ADAPT_STRUCT(
			  conll_type,
			  (conll_type::size_type,   id)
			  (std::string, form)
			  (std::string, lemma)
			  (std::string, cpostag)
			  (std::string, postag)
			  (std::string, feats)
			  (conll_type::size_type, head)
			  (std::string, deprel)
			  (conll_type::phead_type, phead)
			  (std::string, pdeprel)
			  )

typedef std::vector<conll_type, std::allocator<conll_type> > conll_set_type;

template <typename Iterator>
struct conll_parser : boost::spirit::qi::grammar<Iterator, conll_set_type(), boost::spirit::standard::blank_type>
{
  conll_parser() : conll_parser::base_type(conlls)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    token %= qi::lexeme[+(standard::char_ - standard::space)];
    
    conll  %= size >> token >> token >> token >> token >> token >> size >> token >> (size | token) >> token >> qi::eol;
    conlls %= *conll >> qi::eol;
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<conll_type::size_type, 10, 1, -1>            size;
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>    token;
  boost::spirit::qi::rule<Iterator, conll_type(), blank_type>     conll;
  boost::spirit::qi::rule<Iterator, conll_set_type(), blank_type> conlls;
  
};

typedef boost::filesystem::path path_type;
typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Vocab      vocab_type;
typedef cicada::Symbol     word_type;
typedef cicada::Symbol     symbol_type;

typedef std::vector<conll_type::size_type, std::allocator<conll_type::size_type> > index_set_type;
typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_type;

typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;
typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
typedef sentence_type phrase_type;

path_type input_file = "-";
path_type output_file = "-";

std::string goal = "[s]";
std::string non_terminal = "[x]";
bool pos_mode = false;
bool relation_mode = false;
bool leaf_mode = false;
bool projective_mode = false;
bool split_mode = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (pos_mode && relation_mode)
      throw std::runtime_error("either pos or relation or none");
    
    typedef boost::spirit::istream_iterator iter_type;
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file);
    is.unsetf(std::ios::skipws);
    
    conll_parser<iter_type> parser;
    
    conll_set_type conlls;
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
    while (iter != iter_end) {
      conlls.clear();
      
      if (debug)
	std::cerr << "parsing: " << num << std::endl;
      
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, conlls))
	throw std::runtime_error("parsing failed");
      
      if (debug >= 2)
	std::cerr << "size: " << conlls.size() << std::endl;

      ++ num;
      
      if (leaf_mode) {
	if (! conlls.empty()) {
	  conll_set_type::const_iterator citer_end = conlls.end();
	  for (conll_set_type::const_iterator citer = conlls.begin(); citer != citer_end - 1; ++ citer)
	    os << citer->form << ' ';
	  os << conlls.back().form;
	}
	os << '\n';
	
	continue;
      }

      
      hypergraph.clear();
      
      if (conlls.empty()) {
	os << hypergraph << '\n';
	continue;
      }
      
      dependency.clear();
      dependency.resize(conlls.size() + 1);

      node_map.clear();
      node_map.resize(conlls.size() + 1, hypergraph_type::invalid);

      non_terminals.clear();
      non_terminals.push_back(__goal);
      
      conll_set_type::const_iterator citer_end = conlls.end();
      for (conll_set_type::const_iterator citer = conlls.begin(); citer != citer_end; ++ citer) {
	if (projective_mode) {
	  const conll_type::size_type head = boost::apply_visitor(conll_type::visitor_phead(), citer->phead);
	  if (head == conll_type::size_type(-1))
	    throw std::runtime_error("invalid projective head");
	  
	  dependency[head].push_back(citer->id);
	  
	  if (pos_mode)
	    non_terminals.push_back('[' + citer->cpostag + ']');
	  else if (relation_mode)
	    non_terminals.push_back('[' + citer->pdeprel + ']');
	  else
	    non_terminals.push_back(__non_terminal);
	  
	} else {
	  dependency[citer->head].push_back(citer->id);
	  
	  if (pos_mode)
	    non_terminals.push_back('[' + citer->cpostag + ']');
	  else if (relation_mode)
	    non_terminals.push_back('[' + citer->deprel + ']');
	  else
	    non_terminals.push_back(__non_terminal);
	}
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
	
	if (split_mode && conlls[id - 1].form.size() > 1) {
	  // split multi word expression...!
	  
	  phrase.push_back(conlls[id - 1].form);
	} else
	  phrase.push_back(conlls[id - 1].form);
	
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

    ("pos",        po::bool_switch(&pos_mode),        "use pos as non-terminal")
    ("relation",   po::bool_switch(&relation_mode),   "use relation as non-terminal")
    ("leaf",       po::bool_switch(&leaf_mode),       "collect leaf nodes only")
    ("projective", po::bool_switch(&projective_mode), "use projective filed for dependency")
    ("split",      po::bool_switch(&split_mode),      "split multi word expression")
    
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
