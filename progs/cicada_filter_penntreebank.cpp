//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/xpressive/xpressive.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/shared_ptr.hpp>

#include "cicada/hypergraph.hpp"
#include "cicada/stemmer.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/space_separator.hpp"
#include "utils/chart.hpp"
#include "utils/piece.hpp"
#include "utils/getline.hpp"

typedef boost::filesystem::path path_type;

// tree-bank parser...

struct treebank_type
{
  typedef std::vector<treebank_type> antecedents_type;

  std::string cat;
  antecedents_type antecedents;
  bool removed;
  
  treebank_type() : removed(false) {}
  treebank_type(const std::string& __cat) : cat(__cat), removed(false) {}

  void clear()
  {
    cat.clear();
    antecedents.clear();
    removed = false;
  }
};


BOOST_FUSION_ADAPT_STRUCT(
			  treebank_type,
			  (std::string, cat)
			  (std::vector<treebank_type>, antecedents)
			  )



template <typename Iterator>
struct penntreebank_grammar : boost::spirit::qi::grammar<Iterator, treebank_type(), boost::spirit::standard::space_type>
{
  penntreebank_grammar() : penntreebank_grammar::base_type(root)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    comment_last %= qi::no_skip[!(*standard::blank >> qi::lit('(')) >> *(standard::char_ - qi::eol) >> qi::eoi];
    comment      %= qi::no_skip[*(!(*standard::blank >> qi::lit('(')) >> *(standard::char_ - qi::eol) >> qi::eol) >> -comment_last];
    
    cat %= qi::lexeme[+(standard::char_ - standard::space - '(' - ')')];
    treebank %= qi::hold['(' >> cat >> +treebank >> ')'] | cat;
    root %= (qi::omit[comment]
	     >> (qi::hold['(' >> cat >> +treebank >> ')']
		 | qi::hold['(' >> qi::attr("ROOT") >> +treebank >> ')']
		 | qi::hold[qi::lit('(') >> qi::attr("") >> qi::lit('(') >> qi::lit(')') >> qi::lit(')')])
	     >> qi::omit[comment]);
  }
  
  boost::spirit::qi::rule<Iterator, std::string()> comment;
  boost::spirit::qi::rule<Iterator, std::string()> comment_last;
  
  boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type>   cat;
  boost::spirit::qi::rule<Iterator, treebank_type(), boost::spirit::standard::space_type> treebank;
  boost::spirit::qi::rule<Iterator, treebank_type(), boost::spirit::standard::space_type> root;
};

typedef cicada::HyperGraph hypergraph_type;
typedef std::vector<std::string, std::allocator<std::string> > sentence_type;

typedef cicada::Stemmer stemmer_type;

typedef std::pair<int, int> span_type;
typedef std::pair<span_type, std::string> span_cat_type;
typedef std::vector<span_cat_type, std::allocator<span_cat_type> > span_set_type;

typedef std::set<span_cat_type, std::less<span_cat_type>, std::allocator<span_cat_type> > span_sorted_type;


void transform(const hypergraph_type::id_type node_id,
	       const treebank_type& treebank,
	       hypergraph_type& graph)
{
  typedef hypergraph_type::id_type   id_type;
  typedef hypergraph_type::rule_type rule_type;
  typedef std::vector<id_type, std::allocator<id_type> > node_set_type;
  
  if (treebank.antecedents.empty()) return;
  
  std::string rule = "[" + treebank.cat + "] |||";
  
  node_set_type nodes;
  for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter) {
    if (aiter->antecedents.empty())
      rule += " " + aiter->cat;
    else {
      rule += " [" + aiter->cat + "]";
      nodes.push_back(graph.add_node().id);
    }
  }
  
  hypergraph_type::edge_type& edge = graph.add_edge(nodes.begin(), nodes.end());
  edge.rule = rule_type::create(rule_type(rule));
  graph.connect_edge(edge.id, node_id);
  
  node_set_type::const_iterator niter = nodes.begin();
  for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter) {
    if (! aiter->antecedents.empty()) {
      transform(*niter, *aiter, graph);
      ++ niter;
    }
  }
}

void transform(const treebank_type& treebank, hypergraph_type& graph)
{
  graph.goal = graph.add_node().id;
  
  transform(graph.goal, treebank, graph);
}

void transform(const utils::piece pos, const treebank_type& treebank, sentence_type& sent) 
{
  if (treebank.antecedents.empty())
    sent.push_back(treebank.cat + '_' + pos);
  else
    for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform(treebank.cat, *aiter, sent);
}

void transform(const treebank_type& treebank, sentence_type& sent) 
{
  if (treebank.antecedents.empty())
    sent.push_back(treebank.cat);
  else
    for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform(*aiter, sent);
}

void transform_normalize(treebank_type& treebank)
{
  // no terminal...
  if (treebank.antecedents.empty()) return;
  
  // normalize treebank-category...
  if (treebank.cat.size() == 1) {
    switch (treebank.cat[0]) {
    case '.' : treebank.cat = "PERIOD"; break;
    case ',' : treebank.cat = "COMMA"; break;
    case ':' : treebank.cat = "COLON"; break;
    case ';' : treebank.cat = "SEMICOLON"; break;
    }
  } else {
    namespace xpressive = boost::xpressive;
    
    typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
    typedef xpressive::match_results<utils::piece::const_iterator> pmatch;

    static pregex re = (xpressive::s1= -+(~xpressive::_s)) >> (xpressive::as_xpr('-') | xpressive::as_xpr('=')) >> +(~xpressive::_s);
    
    pmatch what;
    if (xpressive::regex_match(utils::piece(treebank.cat), what, re))
      treebank.cat = what[1];
  }
  
  for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
    transform_normalize(*aiter);
}

void transform_cycle(treebank_type& treebank)
{
  // no terminal...
  if (treebank.antecedents.empty()) return;
  
  for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
    transform_cycle(*aiter);

  // unary rule + the same category...
  if (treebank.antecedents.size() == 1 && treebank.antecedents.front().antecedents.size() == 1 && treebank.cat == treebank.antecedents.front().cat) {
    treebank_type::antecedents_type antecedents;
    
    antecedents.swap(treebank.antecedents.front().antecedents);
    
    treebank.antecedents.swap(antecedents);
  }
}

void transform_collapse(treebank_type& treebank)
{
  // no terminal...
  if (treebank.antecedents.empty()) return;
  
  for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
    transform_collapse(*aiter);
  
  // unary rule + the same category...
  if (treebank.antecedents.size() == 1 && treebank.antecedents.front().antecedents.size() == 1) {
    treebank_type::antecedents_type antecedents;

    antecedents.swap(treebank.antecedents.front().antecedents);
    
    treebank.antecedents.swap(antecedents);
  }
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

void transform_unescape(treebank_type& treebank)
{
  if (treebank.antecedents.empty()) {
    // terminal...
    
    namespace qi = boost::spirit::qi;
    
    static terminal_parser<std::string::const_iterator> parser;

    std::string::const_iterator iter = treebank.cat.begin();
    std::string::const_iterator iter_end = treebank.cat.end();

    std::string terminal;
    
    if (! qi::parse(iter, iter_end, parser, terminal) || iter != iter_end)
      throw std::runtime_error("terminal parsing failed?");
    
    treebank.cat.swap(terminal);
    
  } else
    for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform_unescape(*aiter);
}

void transform_bracket(treebank_type& treebank)
{
  if (treebank.antecedents.empty()) {
    const size_t size = treebank.cat.size();
    
    if (treebank.cat[0] == '[' && treebank.cat[size - 1] == ']')
      treebank.cat = treebank.cat.substr(1, size - 2);
    
  } else
    for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform_bracket(*aiter);
}

void transform_stemmer(treebank_type& treebank, const stemmer_type& stemmer)
{
  if (treebank.antecedents.empty())
    treebank.cat = stemmer(treebank.cat);
  else
    for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform_stemmer(*aiter, stemmer);
}

void transform_terminals(treebank_type& treebank)
{
  if (treebank.antecedents.empty()) return;
  
  treebank_type::antecedents_type antecedents;
  
  std::string cat;
  for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter) {
    if (aiter->antecedents.empty())
      cat += aiter->cat;
    else {
      if (! cat.empty()) {
	antecedents.push_back(treebank_type(cat));
	cat.clear();
      }
      antecedents.push_back(*aiter);
    }
  }
  
  if (! cat.empty()) {
    antecedents.push_back(treebank_type(cat));
    cat.clear();
  }
  
  treebank.antecedents.swap(antecedents);
  for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
    transform_terminals(*aiter);
}

void transform_remove_none(treebank_type& treebank)
{
  if (treebank.cat == "-NONE-") {
    treebank.removed = true;
    return;
  }
  
  if (treebank.antecedents.empty()) return;
  
  treebank_type::antecedents_type antecedents;
  
  for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
    if (aiter->cat != "-NONE-") {
      transform_remove_none(*aiter);
      if (! aiter->removed)
	antecedents.push_back(*aiter);
    }
  
  treebank.removed = antecedents.empty();
  treebank.antecedents.swap(antecedents);
}

void transform_map(treebank_type& treebank, sentence_type& sent)
{
  if (treebank.antecedents.empty()) {
    if (sent.empty())
      throw std::runtime_error("no words for mapping?");
    
    treebank.cat = sent.back();
    sent.pop_back();
  } else
    for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform_map(*aiter, sent);
}

bool is_non_terminal(const std::string& cat)
{
  const size_t size = cat.size();
  return size >= 3 && cat[0] == '[' && cat[size - 1] == ']';
}

struct config_span_type
{
  bool exclude_terminal;
  bool binarize;
  bool unary_top;
  bool unary_bottom;
  bool unary_root;
};

// top down traversal of this treebank structure to collect spans(w/ category)
void transform_span(const treebank_type& treebank, span_set_type& spans, int& terminal, const config_span_type& config, const int level)
{
  if (treebank.antecedents.empty()) {
    if (config.exclude_terminal)
      spans.push_back(std::make_pair(std::make_pair(terminal, terminal + 1), std::string()));
    else
      spans.push_back(std::make_pair(std::make_pair(terminal, terminal + 1), treebank.cat));
    ++ terminal;
  } else {
    span_set_type spans_rule;
    
    for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter) {
      transform_span(*aiter, spans, terminal, config, level + 1);
      spans_rule.push_back(spans.back());
    }
    
    if (spans_rule.size() >= 3 && config.binarize)  {
      typedef utils::chart<std::string, std::allocator<std::string> > chart_type;
      
      // its a chart parsing!
      chart_type chart(spans_rule.size() + 1);
      for (size_t pos = 0; pos < spans_rule.size(); ++ pos) 
	chart(pos, pos + 1) = (is_non_terminal(spans_rule[pos].second)
			       ? spans_rule[pos].second.substr(1, spans_rule[pos].second.size() - 2)
			       : spans_rule[pos].second);
      
      for (size_t length = 2; length != spans_rule.size(); ++ length)
	for (size_t first = 0; first + length <= spans_rule.size(); ++ first) {
	  const size_t last = first + length;
	  
	  chart(first, last) = chart(first, last - 1) + '+' + chart(last - 1, last);
	  
	  spans.push_back(std::make_pair(span_type(spans_rule[first].first.first, spans_rule[last - 1].first.second), '[' + chart(first, last) + ']'));
	}
    }
    
    // when unary rule, and not immediate parent of terminal
    if (spans_rule.size() == 1 && ! treebank.antecedents.front().antecedents.empty()) {
      if (config.unary_top || (config.unary_root && level == 0))
	spans.back().second = '[' + treebank.cat + ']';
      else if (config.unary_bottom)
	;
      else 
	spans.back().second = spans.back().second.substr(0, spans.back().second.size() - 1) + ':' + treebank.cat + ']';

    } else
      spans.push_back(std::make_pair(span_type(spans_rule.front().first.first, spans_rule.back().first.second), '[' + treebank.cat + ']'));
  }
}


void transform_span(const treebank_type& treebank, span_set_type& spans, const config_span_type& config)
{
  int terminal = 0;
  
  transform_span(treebank, spans, terminal, config, 0);
}

std::ostream& treebank_output(const treebank_type& treebank, std::ostream& os)
{
  if (treebank.antecedents.empty())
    os << treebank.cat;
  else {
    os << '(';
    os << treebank.cat;
    os << ' ';
    
    for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      treebank_output(*aiter, os);
    os << ')';
  }
  
  return os;
}


path_type input_file = "-";
path_type output_file = "-";
path_type map_file;

std::string root_symbol;
bool normalize = false;
bool remove_none = false;
bool unescape_terminal = false;
bool remove_cycle = false;
bool collapse = false;
bool remove_bracket = false;
bool add_bos_eos = false;
bool fix_terminal = false;
std::string stemmer;

bool leaf = false;
bool leaf_pos = false;
bool rule = false;
bool treebank_mode = false;

bool span = false;
bool category = false;
bool binarize = false;
bool unary_top = false;
bool unary_bottom = false;
bool unary_root = false;
bool exclude_terminal = false;

bool skip_invalid = false;
bool validate = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (unary_top && unary_bottom)
      throw std::runtime_error("which strategy? unary-[top|bottom]");

    if (int(leaf) + leaf_pos + rule + span + treebank_mode > 1)
      throw std::runtime_error("multiple output options specified: leaf/leaf-pos/rule/span/treebank/(default for cicada-forest)");
      

    typedef boost::spirit::istream_iterator iter_type;

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    boost::shared_ptr<utils::compress_istream> ms;

    if (! map_file.empty()) {
      if (! boost::filesystem::exists(map_file))
	throw std::runtime_error("no map file: " + map_file.string());
      
      ms.reset(new utils::compress_istream(map_file, 1024 * 1024));
    }

    const stemmer_type* __stemmer = (! stemmer.empty() ? &stemmer_type::create(stemmer) : 0);
    
    penntreebank_grammar<iter_type>         grammar;

    treebank_type   parsed;
    hypergraph_type graph;
    sentence_type   sent;
    
    span_set_type    spans;
    span_sorted_type spans_sorted;
    config_span_type config_span;
    config_span.exclude_terminal = exclude_terminal;
    config_span.binarize         = binarize;
    config_span.unary_top        = unary_top;
    config_span.unary_bottom     = unary_bottom;
    config_span.unary_root       = unary_root;

    std::string line;
    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      parsed.clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, grammar, boost::spirit::standard::space, parsed)) {
	std::string buffer;
	for (int i = 0; i != 64 && iter != iter_end; ++ i, ++iter)
	  buffer += *iter;
	
	throw std::runtime_error("parsing failed: " + buffer);
      }

      
      if (ms) {
	if (! utils::getline(*ms, line))
	  throw std::runtime_error("# of lines do not match with map-file");
	
	utils::piece line_piece(line);
	boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer(line_piece);
	
	sent.assign(tokenizer.begin(), tokenizer.end());
	
	std::reverse(sent.begin(), sent.end());
	
	if (! sent.empty()) {
	  transform_map(parsed, sent);
	  if (! sent.empty())
	    throw std::runtime_error("# of words do not match?");
	}
      }

      if (! root_symbol.empty())
	parsed.cat = root_symbol;
      else if (parsed.cat.empty())
	parsed.cat = "ROOT";

      if (fix_terminal)
	transform_terminals(parsed);
      
      if (remove_none)
	transform_remove_none(parsed);
      
      if (normalize)
	transform_normalize(parsed);

      if (remove_cycle)
	transform_cycle(parsed);

      if (collapse)
	transform_collapse(parsed);
      
      if (unescape_terminal)
	transform_unescape(parsed);

      if (__stemmer)
	transform_stemmer(parsed, *__stemmer);

      if (remove_bracket)
	transform_bracket(parsed);

      if (add_bos_eos && ! parsed.antecedents.empty()) {
	parsed.antecedents.insert(parsed.antecedents.begin(), treebank_type("BOS"));
	parsed.antecedents.insert(parsed.antecedents.end(), treebank_type("EOS"));
	parsed.antecedents.front().antecedents.push_back(treebank_type("<s>"));
	parsed.antecedents.back().antecedents.push_back(treebank_type("</s>"));
      }
      
      if (leaf) {
	sent.clear();
	
	transform(parsed, sent);
	
	if (! sent.empty()) {
	  std::copy(sent.begin(), sent.end() - 1, std::ostream_iterator<std::string>(os, " "));
	  os << sent.back();
	  os << '\n';
	} else
	  os << '\n';
	
	if (flush_output)
	  os << std::flush;
      } else if (leaf_pos) {
	sent.clear();
	
	transform("UNK", parsed, sent);
	
	if (! sent.empty()) {
	  std::copy(sent.begin(), sent.end() - 1, std::ostream_iterator<std::string>(os, " "));
	  os << sent.back();
	  os << '\n';
	} else
	  os << '\n';
	
	if (flush_output)
	  os << std::flush;
      } else if (span) {
	spans.clear();
	
	transform_span(parsed, spans, config_span);
	
	if (spans.empty())
	  os << '\n';
	else {
	  
	  spans_sorted.clear();
	  spans_sorted.insert(spans.begin(), spans.end());
	  
	  if (category) {
	    bool initial = true;
	    span_sorted_type::const_iterator siter_end = spans_sorted.end();
	    for (span_sorted_type::const_iterator siter = spans_sorted.begin(); siter != siter_end; ++ siter) {
	      if (! initial)
		os << ' ';
	      
	      if (! siter->second.empty())
		os << siter->first.first << '-' << siter->first.second << ':' << siter->second;
	      
	      initial = false;
	    }
	    
	  } else {
	    span_sorted_type::const_iterator siter_end = spans_sorted.end();
	    span_sorted_type::const_iterator siter_prev = spans_sorted.end();
	    for (span_sorted_type::const_iterator siter = spans_sorted.begin(); siter != siter_end; ++ siter) {
	      
	      if (siter_prev == siter_end)
		os << siter->first.first << '-' << siter->first.second;
	      else if (siter_prev->first != siter->first)
		os << ' ' << siter->first.first << '-' << siter->first.second;
	      
	      siter_prev = siter;
	    }
	  }
	  
	  os << '\n';
	}
	
	if (flush_output)
	  os << std::flush;

      } else if (treebank_mode) {
	if (parsed.antecedents.empty())
	  os << "(())";
	else
	  treebank_output(parsed, os);
	os << '\n';
	if (flush_output)
	  os << std::flush;
      } else {
	graph.clear();
	
	transform(parsed, graph);

	if (debug)
	  std::cerr << "transformed into hypergraph" << std::endl;

	if (! graph.edges.empty()) {
	  graph.topologically_sort();
	  
	  if (validate) {
	    hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
	    for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	      const hypergraph_type::edge_type& edge = *eiter;
	      
	      if (edge.tails.empty()) {
		if (edge.rule->rhs.size() != 1)
		  throw std::runtime_error("terminal rule has more than one terminal?");
		if (edge.rule->rhs.front().is_non_terminal())
		  throw std::runtime_error("terminal rule has non-terminal in rhs?");
	      } else {
		hypergraph_type::rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
		for (hypergraph_type::rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
		  if (riter->is_terminal())
		    throw std::runtime_error("non-terminal rule has terminal!");
	      }
	    }
	  }
	  
	} else
	  graph.clear();
	
	if (rule) {
	  hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
	  for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter)
	    if (eiter->rule)
	      os << *(eiter->rule) << '\n';
	} else {
	  if (! skip_invalid || graph.is_valid())
	    os << graph << '\n';
	}
	
	if (flush_output)
	  os << std::flush;
      }
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

    ("fix-terminal",   po::bool_switch(&fix_terminal),       "fix fragmented terminals")
    ("replace-root",   po::value<std::string>(&root_symbol), "replace root symbol")
    ("unescape",       po::bool_switch(&unescape_terminal),  "unescape terminal symbols, such as -LRB-, \\* etc.")
    ("normalize",      po::bool_switch(&normalize),          "normalize category, such as [,] [.] etc.")
    ("remove-none",    po::bool_switch(&remove_none),        "remove -NONE-")
    ("remove-cycle",   po::bool_switch(&remove_cycle),       "remove cycle unary rules")
    ("remove-bracket", po::bool_switch(&remove_bracket),     "remove bracket terminal (which may be confued with non-terminal!)")
    ("collapse",       po::bool_switch(&collapse),           "collapse unary rules")
    ("add-bos-eos",    po::bool_switch(&add_bos_eos),        "add [BOS]/[EOS] and <s>/</s>")
    ("stemmer",        po::value<std::string>(&stemmer),     "stemming for terminals")
    
    ("leaf",      po::bool_switch(&leaf),          "output leaf nodes")
    ("leaf-pos",  po::bool_switch(&leaf_pos),      "output leaf/pos nodes")
    ("rule",      po::bool_switch(&rule),          "output rules")
    ("treebank",  po::bool_switch(&treebank_mode), "output treebank")
    
    ("span",      po::bool_switch(&span),            "output spans")
    ("binarize",  po::bool_switch(&binarize),        "span: perform binarization")
    ("category",  po::bool_switch(&category),        "span: added category to span")
    
    ("unary-top",    po::bool_switch(&unary_top),    "span: use top-most category for unary rules")
    ("unary-bottom", po::bool_switch(&unary_bottom), "span: use bottom-most category for unary rules")
    ("unary-root",   po::bool_switch(&unary_root),   "span: use single category for root")
    
    ("exclude-terminal", po::bool_switch(&exclude_terminal), "span: no terminal in span")
    
    ("skip",     po::bool_switch(&skip_invalid), "skip invalid penntree")
    ("validate", po::bool_switch(&validate),     "validate penntreebank")
    
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
