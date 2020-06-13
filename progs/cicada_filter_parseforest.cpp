//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// filter for forest-output from egret and/or parseIt
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>

#include "cicada/hypergraph.hpp"
#include "cicada/vocab.hpp"
#include "cicada/sort.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/chart.hpp"
#include "utils/space_separator.hpp"
#include "utils/unordered_map.hpp"
#include "utils/getline.hpp"

typedef boost::filesystem::path path_type;

struct category_type
{
  std::string cat;
  int first;
  int last;
  
  category_type()
    : cat(), first(-1), last(-1) {}
  category_type(const std::string& __cat)
    : cat(__cat), first(-1), last(-1) {}
  category_type(const std::string& __cat, const int& __first, const int& __last)
    : cat(__cat), first(__first), last(__last) {}
  
  bool is_terminal() const { return first < 0 || last < 0; }

  std::string strip() const
  {
    std::string ::size_type pos = cat.find('_');
    if (pos != std::string::npos)
      return cat.substr(0, pos);
    else {
      std::string ::size_type pos = cat.find('^');
      if (pos != std::string::npos)
	return cat.substr(0, pos);
      else
	return cat;
    }
  }
};

inline
std::string normalize_cat(const std::string& cat)
{
  if (cat.size() == 1) {
    switch (cat[0]) {
    case '.' : return "PERIOD";
    case ',' : return "COMMA";
    case ':' : return "COLON";
    case ';' : return "SEMICOLON";
    default: return cat;
    }
  } else
    return cat;
}

template <typename Iterator>
struct terminal_parser : boost::spirit::qi::grammar<Iterator, std::string()>
{
  terminal_parser() : terminal_parser::base_type(terminal)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    escape_char.add
      ("-LRB-", '(')
      ("-RRB-", ')')
      ("-LSB-", '[')
      ("-RSB-", ']')
      ("-LCB-", '{')
      ("-RCB-", '}')
      ("-PLUS-", '+') // added for ATB
      ("\\/", '/')
      ("\\*", '*');
    
    terminal %= +(escape_char | standard::char_);
  }
  
  boost::spirit::qi::symbols<char, char> escape_char;
  boost::spirit::qi::rule<Iterator, std::string()> terminal;
};

typedef std::vector<category_type, std::allocator<category_type> > category_set_type;

typedef boost::fusion::tuple<category_type, category_set_type, double> item_type;

typedef std::vector<item_type, std::allocator<item_type> > item_set_type;
typedef std::vector<std::string, std::allocator<std::string> > sentence_type;

struct forest_type
{
  int id;
  sentence_type sentence;
  item_set_type items;

  forest_type()
    : id(0), sentence(), items() {}
  forest_type(const sentence_type& __sentence, const item_set_type& __items)
    : id(0), sentence(__sentence), items(__items) {}

  forest_type(const int& __id, const sentence_type& __sentence, const item_set_type& __items)
    : id(__id), sentence(__sentence), items(__items) {}

  void clear()
  {
    id = 0;
    sentence.clear();
    items.clear();
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  category_type,
			  (std::string, cat)
			  (int, first)
			  (int, last)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  forest_type,
			  (int, id)
			  (sentence_type, sentence)
			  (item_set_type, items)
			  )

template <typename Iterator>
struct forest_parser : boost::spirit::qi::grammar<Iterator, forest_type(), boost::spirit::standard::blank_type>
{
  forest_parser() : forest_parser::base_type(forest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    cat %= qi::lexeme[+(standard::char_ - standard::space - '[') - ("|||" >> (standard::space | qi::eoi))];
    
    category %= qi::hold[cat >> '[' >> qi::int_ >> ',' >> qi::int_ >> ']'] | cat;
    
    item %= category >> "=>" >> (+category) >> "|||" >> qi::double_ >> -qi::lit("EXTRAVAL") >> qi::eol;
    
    sentence %= +cat >> qi::eol;
    
    forest %= ("sentence" >> qi::int_ >> ':' >> qi::eol) >> (-sentence) >> (+item | qi::eol) >> qi::eol;
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>   cat;
  boost::spirit::qi::rule<Iterator, category_type(), blank_type> category;
  
  boost::spirit::qi::rule<Iterator, item_type(), blank_type> item;
  
  boost::spirit::qi::rule<Iterator, sentence_type(), blank_type>   sentence;

  boost::spirit::qi::rule<Iterator, forest_type(), blank_type> forest;
};

typedef cicada::HyperGraph hypergraph_type;


typedef utils::unordered_map<std::string, hypergraph_type::id_type, boost::hash<utils::piece>, std::equal_to<std::string>,
			     std::allocator<std::pair<const std::string, hypergraph_type::id_type> > >::type node_map_type;

typedef utils::chart<node_map_type, std::allocator<node_map_type> > node_chart_type;

typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;

typedef cicada::Vocab  vocab_type;
typedef cicada::Symbol word_type;
typedef std::vector<word_type, std::allocator<word_type> > phrase_type;

struct filter_edge
{
  typedef std::vector<bool, std::allocator<bool> > removed_type;
  
  filter_edge(const removed_type& __removed) : removed(__removed) {}
  
  bool operator()(const hypergraph_type::edge_type& edge) const
  {
    return removed[edge.id];
  }

  const removed_type& removed;
};

path_type input_file = "-";
path_type output_file = "-";
path_type map_file;

std::string root;
bool unescape_terminal = false;
bool normalize = false;
bool collapse = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    typedef boost::spirit::istream_iterator iter_type;
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file);
    is.unsetf(std::ios::skipws);

    boost::shared_ptr<utils::compress_istream> ms;
    
    if (! map_file.empty()) {
      if (! boost::filesystem::exists(map_file))
	throw std::runtime_error("no map file: " + map_file.string());
      
      ms.reset(new utils::compress_istream(map_file, 1024 * 1024));
    }

    const bool mapping = ms.get();
    
    forest_parser<iter_type> parser;

    forest_type     forest;
    hypergraph_type hypergraph;
    node_chart_type chart;
    tail_set_type   tails;
    phrase_type     phrase;
    sentence_type   sentence;

    hypergraph_type::feature_set_type::feature_type feature("parse-cost");
    
    std::string line;
    iter_type iter(is);
    iter_type iter_end;
    
    hypergraph_type::symbol_type root_label(root);
    if (! root_label.empty())
      if (! root_label.is_non_terminal())
	throw std::runtime_error("invalid root label: " + root);

    int num = 0;
    while (iter != iter_end) {
      forest.clear();

      if (debug)
	std::cerr << "parsing: " << num << std::endl;

      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, forest))
	throw std::runtime_error("parsing failed");

      if (debug >= 2) {
	std::cerr << "sentence: ";
	std::copy(forest.sentence.begin(), forest.sentence.end(), std::ostream_iterator<std::string>(std::cerr, " "));
	std::cerr << std::endl;
	std::cerr << "forest size: " << forest.items.size() << std::endl;
      }

      hypergraph.clear();

      if (mapping) {	
	if (! utils::getline(*ms, line))
	  throw std::runtime_error("# of lines do not match with map-file");
	
	utils::piece line_piece(line);
	boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer(line_piece);
	sentence.clear();
	sentence.insert(sentence.end(), tokenizer.begin(), tokenizer.end());

	if (sentence.size() != forest.sentence.size())
	  throw std::runtime_error("# of words in map file and parse output differ");
	
	forest.sentence.swap(sentence);
      }
      
      ++ num;
      
      if (forest.sentence.empty() || forest.items.empty()) {
	os << hypergraph << '\n';
	continue;
      }
      
      chart.clear();
      chart.reserve(forest.sentence.size() + 1);
      chart.resize(forest.sentence.size() + 1);
      
      hypergraph_type::id_type node_last = hypergraph_type::invalid;

      // transform into hypergraph...
      item_set_type::const_iterator iiter_end = forest.items.end();
      for (item_set_type::const_iterator iiter = forest.items.begin(); iiter != iiter_end; ++ iiter) {
	const item_type& item = *iiter;
	
	const category_type&     lhs   = boost::fusion::get<0>(item);
	const category_set_type& rhs   = boost::fusion::get<1>(item);
	const double&            score = boost::fusion::get<2>(item);
		
	tails.clear();
	phrase.clear();
	
	category_set_type::const_iterator riter_end = rhs.end();
	for (category_set_type::const_iterator riter = rhs.begin(); riter != riter_end; ++ riter) {
	  const category_type& cat = *riter;
	  
	  if (cat.is_terminal()) {
	    if (mapping)
	      phrase.push_back(forest.sentence[lhs.first]);
	    else if (unescape_terminal) {
	      namespace qi = boost::spirit::qi;
	      
	      static terminal_parser<std::string::const_iterator> parser;
	      
	      std::string::const_iterator iter = cat.cat.begin();
	      std::string::const_iterator iter_end = cat.cat.end();
	      
	      std::string terminal;
	      
	      if (! qi::parse(iter, iter_end, parser, terminal) || iter != iter_end)
		throw std::runtime_error("terminal parsing failed?");
	      
	      phrase.push_back(terminal);
	    } else
	      phrase.push_back(cat.cat);
	  } else {
	    // perform normalization if specified...!
	    if (normalize)
	      phrase.push_back('[' + normalize_cat(cat.strip()) + ']');
	    else
	      phrase.push_back('[' + cat.strip() + ']');

	    node_map_type& node_map = chart(cat.first, cat.last + 1);
	    
	    std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(cat.cat, 0));
	    if (result.second)
	      result.first->second = hypergraph.add_node().id;
	    
	    tails.push_back(result.first->second);
	  }
	}
	
	node_map_type& node_map = chart(lhs.first, lhs.last + 1);
	
	const word_type cat = (normalize
			       ? '[' + normalize_cat(lhs.strip()) + ']'
			       : '[' + lhs.strip() + ']');
	
	std::pair<node_map_type::iterator, bool> result = node_map.insert(std::make_pair(lhs.cat, 0));
	if (result.second)
	  result.first->second = hypergraph.add_node().id;
	
	const hypergraph_type::id_type& node_id = result.first->second;
	
	hypergraph_type::edge_type& edge = hypergraph.add_edge(tails.begin(), tails.end());
	edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(cat, phrase.begin(), phrase.end()));
	
	if (score != 0.0)
	  edge.features[feature] = score;
	
	hypergraph.connect_edge(edge.id, node_id);
	
	node_last = node_id;
      }
      
      // we assume that the last item is always the last rule leading to goal...
      if (node_last != hypergraph_type::invalid) {
	hypergraph.goal = node_last;
	
	if (collapse) {
	  // collapse root label...
	  
	  bool found = false;
	  filter_edge::removed_type removed;
	  
	  do {
	    found = false;
	    removed.clear();
	    removed.resize(hypergraph.edges.size(), false);
	    
	    hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = hypergraph.nodes[hypergraph.goal].edges.end();
	    for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = hypergraph.nodes[hypergraph.goal].edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = hypergraph.edges[*eiter];
	      
	      if (edge.tails.size() != 1 || edge.rule->rhs.size() != 1 || edge.rule->rhs.front() != edge.rule->lhs) continue;
	      
	      // we will collapse this edge...
	      
	      found = true;
	      removed[*eiter] = true;
	      
	      hypergraph_type::node_type::edge_set_type::const_iterator aiter_end = hypergraph.nodes[edge.tails.front()].edges.end();
	      for (hypergraph_type::node_type::edge_set_type::const_iterator aiter = hypergraph.nodes[edge.tails.front()].edges.begin(); aiter != aiter_end; ++ aiter) {
		const hypergraph_type::edge_type& edge_antecedent = hypergraph.edges[*aiter];
		
		removed[*aiter] = true;
		
		hypergraph_type::edge_type& edge_new = hypergraph.add_edge(edge_antecedent);
		edge_new.features += edge.features;
		
		hypergraph.connect_edge(edge_new.id, hypergraph.goal);
	      }
	    }
	    
	    // resize again...
	    removed.resize(hypergraph.edges.size(), false);
	    
	    // topologically sort...
	    if (found) {
	      hypergraph_type sorted;
	      topologically_sort(hypergraph, sorted, filter_edge(removed), true);
	      hypergraph.swap(sorted);
	    }
	    
	  } while (found);
	}
	
	if (! root_label.empty()) {
	  hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = hypergraph.nodes[hypergraph.goal].edges.end();
	  for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = hypergraph.nodes[hypergraph.goal].edges.begin(); eiter != eiter_end; ++ eiter) {
	    hypergraph_type::edge_type& edge = hypergraph.edges[*eiter];
	    
	    edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(root_label, edge.rule->rhs));
	  }
	}
	
	hypergraph.topologically_sort();
      }
      
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
    ("map",       po::value<path_type>(&map_file)->default_value(map_file), "map terminal symbols")
    ("root",      po::value<std::string>(&root), "root label")
    
    ("unescape",  po::bool_switch(&unescape_terminal),  "unescape terminal symbols, such as -LRB-, \\* etc.")
    ("normalize", po::bool_switch(&normalize), "normalize category, such as [,] [.] etc.")
    ("collapse",  po::bool_switch(&collapse),  "collapse root labels")
    
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
