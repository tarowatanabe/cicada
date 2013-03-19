//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// dependency filter to project source dependency into target dependency using the word alignment
//
// or, transform various dependency format into cicada one-line format:
//
// sentence ||| POS ||| dependency
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <memory>

// this is important!
#define FUSION_MAX_VECTOR_SIZE 15

#include <boost/variant.hpp>
#include <boost/tokenizer.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <cicada/alignment.hpp>
#include <cicada/dependency.hpp>
#include <cicada/sentence.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/sort_topologically.hpp>
#include <cicada/debinarize.hpp>
#include <cicada/remove_non_terminal.hpp>

#include "utils/bithack.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/space_separator.hpp"
#include "utils/chart.hpp"
#include "utils/piece.hpp"

typedef cicada::Alignment  alignment_type;
typedef cicada::Dependency dependency_type;
typedef cicada::HyperGraph hypergraph_type;

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef std::vector<std::string, std::allocator<std::string> > sentence_type;

template <typename Iterator>
struct sentence_parser : boost::spirit::qi::grammar<Iterator, sentence_type(), boost::spirit::standard::blank_type>
{
  sentence_parser() : sentence_parser::base_type(tokens)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    word    %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];
    tokens  %= *word >> (qi::eoi | qi::eol);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type> word;
  boost::spirit::qi::rule<Iterator, sentence_type(), blank_type> tokens;
};

path_type input_file = "-";
path_type output_file = "-";
path_type map_file;

std::string goal = "[s]";
std::string non_terminal = "[x]";

bool mst_mode = false;
bool conll_mode = false;
bool malt_mode = false;
bool cabocha_mode = false;
bool khayashi_mode = false;
bool khayashi_forest_mode = false;
bool cicada_mode = false;

bool projective_mode = false;
bool relation_mode = false;
bool unescape_mode = false;
bool normalize_mode = false;
bool func_mode = false;
bool forest_mode = false;
bool head_mode = false;
bool span_mode = false;
bool binarize_mode = false;
bool category_mode = false;
bool leaf_mode = false;
bool leaf_pos_mode = false;

void options(int argc, char** argv);

template <typename Dep>
void apply(const path_type& file, const path_type& map, const path_type& output)
{
  Dep dep;
  dep(file, map, output);
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

std::string unescape(const std::string& word)
{
  namespace qi = boost::spirit::qi;
  
  static terminal_parser<std::string::const_iterator> parser;
  
  std::string::const_iterator iter = word.begin();
  std::string::const_iterator iter_end = word.end();
  
  std::string terminal;
  
  if (! qi::parse(iter, iter_end, parser, terminal) || iter != iter_end)
    throw std::runtime_error("terminal parsing failed?");
  
  return terminal;
}

template <typename Iterator>
void unescape(Iterator first, Iterator last)
{
  for (/**/; first != last; ++ first)
    *first = unescape(*first);
}

std::string normalize(const std::string& pos)
{
  if (pos.size() == 1) {
    switch (pos[0]) {
    case '.' : return "PERIOD";
    case ',' : return "COMMA";
    case ':' : return "COLON";
    case ';' : return "SEMICOLON";
    }
    return pos;
  } else
    return pos;
}

template <typename Iterator>
void normalize(Iterator first, Iterator last)
{
  for (/**/; first != last; ++ first)
    *first = normalize(*first);
}

struct MST;
struct CoNLL;
struct Malt;
struct Cabocha;
struct KHayashi;
struct KHayashiForest;
struct Cicada;

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    //default to cicada-native-mode
    if (int(mst_mode) + conll_mode + malt_mode + cabocha_mode + khayashi_mode + khayashi_forest_mode + cicada_mode == 0)
      cicada_mode = true;

    if (int(mst_mode) + conll_mode + malt_mode + cabocha_mode + khayashi_mode + khayashi_forest_mode + cicada_mode > 1)
      throw std::runtime_error("one of mst/conll/malt/khayashi/khayashi-forest/cicada mode");

    if (int(forest_mode) + span_mode + leaf_mode + leaf_pos_mode> 1)
      throw std::runtime_error("one of forest/span/leaf/leaf-pos");
    
    if (mst_mode)
      apply<MST>(input_file, map_file, output_file);
    else if (conll_mode)
      apply<CoNLL>(input_file, map_file, output_file);
    else if (malt_mode)
      apply<Malt>(input_file, map_file, output_file);    
    else if (cabocha_mode)
      apply<Cabocha>(input_file, map_file, output_file);
    else if (khayashi_mode)
      apply<KHayashi>(input_file, map_file, output_file);
    else if (khayashi_forest_mode)
      apply<KHayashiForest>(input_file, map_file, output_file);
    else if (cicada_mode)
      apply<Cicada>(input_file, map_file, output_file);
    else
      throw std::runtime_error("one of mst|conll|cabocha|khayashi|khayashi-forest|cicada?");
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}

struct Dependency
{
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Dependency dependency_type;
  typedef cicada::SpanVector span_set_type;
  typedef cicada::HyperGraph hypergraph_type;

  Dependency() {}

  void clear()
  {
    sentence.clear();
    postag.clear();
    dependency.clear();

    hypergraph.clear();
    spans.clear();
  }

  bool verify() const
  {
    return sentence.size() == postag.size() && sentence.size() == dependency.size();
  }
  
  sentence_type   sentence;
  sentence_type   postag;
  dependency_type dependency;

  hypergraph_type hypergraph;
  span_set_type   spans;
};

struct MapFile
{
  typedef std::vector<std::string, std::allocator<std::string> > sentence_type;
  typedef boost::spirit::istream_iterator iterator_type;
  typedef sentence_parser<iterator_type> parser_type;
  
  MapFile(const path_type& path) : is(0), sentence(), iter(), iter_end()
  {
    if (! path.empty()) {
      if (path != "-" && ! boost::filesystem::exists(path))
	throw std::runtime_error("invalid map file? " + path.string());
      
      is = new utils::compress_istream(path, 1024 * 1024);
      is->unsetf(std::ios::skipws);
      
      iter = iterator_type(*is);
    }
  }
  ~MapFile() { if (is) { delete is; } }

  operator bool() const { return is != 0; }
  
  const sentence_type& operator()()
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    sentence.clear();
    if (iter != iter_end)
      if (! qi::phrase_parse(iter, iter_end, parser, standard::blank, sentence))
	throw std::runtime_error("parsing failed");
    
    return sentence;
  }
  
  std::istream* is;
  sentence_type sentence;
  
  parser_type   parser;
  iterator_type iter;
  iterator_type iter_end;
};

struct TransformSpan
{
  typedef size_t size_type;
  
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Dependency dependency_type;
  typedef cicada::SpanVector span_set_type;
  typedef cicada::Symbol     symbol_type;
  
  typedef span_set_type::span_type span_type;

  typedef std::vector<size_type, std::allocator<size_type> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_map_type;
  typedef span_set_type span_map_type;

  typedef utils::chart<std::string, std::allocator<std::string> > chart_type;
  
  TransformSpan(const symbol_type& __goal, const bool __binarize, const bool __category)
    : goal(__goal), binarize(__binarize), category(__category) {} 
  
  dependency_map_type dependency_map;
  span_map_type       span_map;
  chart_type          chart;

  const symbol_type goal;
  const bool binarize;
  const bool category;
  
  void operator()(const sentence_type& sentence,
		  const sentence_type& postag,
		  const dependency_type& dependency,
		  span_set_type& spans)
  {
    dependency_map.clear();
    dependency_map.resize(dependency.size() + 1);
    
    span_map.clear();
    if (category)
      span_map.resize(dependency.size() + 1, span_type(0, 0, goal));
    else
      span_map.resize(dependency.size() + 1, span_type(0, 0));
    
    for (size_type i = 0; i != dependency.size(); ++ i)
      dependency_map[dependency[i]].push_back(i + 1);
    
    // post order traversal of dependency-map to compute spans...
    spans.clear();
    traverse(0, postag, spans);
  }
  
  void traverse(const size_type node, const sentence_type& postag, span_set_type& spans)
  {
    span_set_type antecedents;
    
    if (category)
      span_map[node].label = (node == 0 ? goal : postag[node - 1]);
    
    int span_pos =  span_map[node].first;
    
    index_set_type::const_iterator iiter_begin = dependency_map[node].begin();
    index_set_type::const_iterator iiter_end   = dependency_map[node].end();
    index_set_type::const_iterator iiter_lex   = std::lower_bound(iiter_begin, iiter_end, node);
    
    for (index_set_type::const_iterator iiter = iiter_begin; iiter != iiter_lex; ++ iiter) {
      span_map[*iiter].first = span_pos;
      
      traverse(*iiter, postag, spans);
      
      span_pos = span_map[*iiter].last;

      antecedents.push_back(span_map[*iiter]);
    }
    
    // increment when non-root!
    // TODO: we need to handle span for head, like head-mode in transform-forest
    
    const int head_first = span_pos;
    
    span_pos += (node != 0);
    
    const int head_last = span_pos;
    
    if (node != 0 && ! dependency_map[node].empty()) {
      if (category)
	antecedents.push_back(span_type(head_first, head_last, postag[node - 1]));
      else
	antecedents.push_back(span_type(head_first, head_last));
      spans.push_back(antecedents.back());
    }
    
    for (index_set_type::const_iterator iiter = iiter_lex; iiter != iiter_end; ++ iiter) {
      span_map[*iiter].first = span_pos;
      
      traverse(*iiter, postag, spans);
      
      span_pos = span_map[*iiter].last;
      
      antecedents.push_back(span_map[*iiter]);
    }
    
    span_map[node].last = span_pos;
    
    // parents...!
    spans.push_back(span_map[node]);
    
    if (! binarize || antecedents.size() <= 2) return;
    
    chart.clear();
    chart.resize(antecedents.size() + 1);
    
    for (size_type pos = 0; pos != antecedents.size(); ++ pos) 
      chart(pos, pos + 1) = antecedents[pos].label.non_terminal_strip();
    
    for (size_type length = 2; length != antecedents.size(); ++ length)
      for (size_type first = 0; first + length <= antecedents.size(); ++ first) {
	const size_type last = first + length;
	
	chart(first, last) = chart(first, last - 1) + '+' + chart(last - 1, last);
	
	if (category)
	  spans.push_back(span_type(antecedents[first].first, antecedents[last - 1].last, '[' + chart(first, last) +']'));
	else
	  spans.push_back(span_type(antecedents[first].first, antecedents[last - 1].last));
      }
  }
};


struct TransformForest
{
  typedef size_t size_type;

  typedef cicada::Sentence   sentence_type;
  typedef cicada::Dependency dependency_type;
  typedef cicada::Rule       rule_type;
  typedef cicada::Symbol     symbol_type;

  typedef std::vector<size_type, std::allocator<size_type> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > dependency_map_type;
  
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > node_map_type;
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > tail_set_type;
  typedef sentence_type rhs_type;
  
  dependency_map_type dependency_map;
  node_map_type       node_map;
  tail_set_type       tails;
  rhs_type            rhs;
  
  const symbol_type goal;
  const bool head_mode;
  
  TransformForest(const symbol_type& __goal,
		  const bool __head_mode)
    : goal(__goal),
      head_mode(__head_mode) {} 

  void operator()(const sentence_type& sentence,
		  const sentence_type& postag,
		  const dependency_type& dependency,
		  hypergraph_type& hypergraph)
  {
    if (sentence.size() != postag.size() || sentence.size() != dependency.size())
      throw std::runtime_error("invalid transformaiton");
    
    hypergraph.clear();
    
    if (sentence.empty()) return;
    
    dependency_map.clear();
    dependency_map.resize(dependency.size() + 1);
    
    node_map.clear();
    node_map.resize(dependency.size() + 1, hypergraph_type::invalid); 
    
    for (size_type i = 0; i != dependency.size(); ++ i)
      dependency_map[dependency[i]].push_back(i + 1);
    
    if (dependency_map.front().empty())
      throw std::runtime_error("invalid dependency structure without root");
    
    tails.clear();
    rhs.clear();
    
    index_set_type::const_iterator iiter_end = dependency_map.front().end();
    for (index_set_type::const_iterator iiter = dependency_map.front().begin(); iiter != iiter_end; ++ iiter) {
      const size_type antecedent = *iiter;
      
      if (node_map[antecedent] == hypergraph_type::invalid)
	node_map[antecedent] = hypergraph.add_node().id;
      
      tails.push_back(node_map[antecedent]);
      rhs.push_back(postag[antecedent - 1]);
    }
    
    if (node_map[0] == hypergraph_type::invalid)
      node_map[0] = hypergraph.add_node().id;

    hypergraph_type::edge_type& edge = hypergraph.add_edge(tails.begin(), tails.end());
    edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(goal, rhs.begin(), rhs.end()));
    
    hypergraph.connect_edge(edge.id, node_map[0]);
    hypergraph.goal = node_map[0];
    
    for (size_t id = 1; id != dependency_map.size(); ++ id) {
      tails.clear();
      rhs.clear();
      
      index_set_type::const_iterator iiter_begin = dependency_map[id].begin();
      index_set_type::const_iterator iiter_end   = dependency_map[id].end();
      index_set_type::const_iterator iiter_lex   = std::lower_bound(iiter_begin, iiter_end, id);
      
      for (index_set_type::const_iterator iiter = iiter_begin; iiter != iiter_lex; ++ iiter) {
	const size_type antecedent = *iiter;
	
	if (node_map[antecedent] == hypergraph_type::invalid)
	  node_map[antecedent] = hypergraph.add_node().id;
	
	tails.push_back(node_map[antecedent]);
	rhs.push_back(postag[antecedent - 1]);
      }
      
      if (head_mode) {
	const symbol_type lhs = '[' + postag[id - 1].non_terminal_strip() + "*]";
	tails.push_back(hypergraph.add_node().id);
	rhs.push_back(lhs);
	
	hypergraph_type::edge_type& edge = hypergraph.add_edge();
	edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(lhs, hypergraph_type::rule_type::symbol_set_type(1, sentence[id - 1])));
	
	hypergraph.connect_edge(edge.id, tails.back());
      } else
	rhs.push_back(sentence[id - 1]);
      
      for (index_set_type::const_iterator iiter = iiter_lex; iiter != iiter_end; ++ iiter) {
	const size_type antecedent = *iiter;
	
	if (node_map[antecedent] == hypergraph_type::invalid)
	  node_map[antecedent] = hypergraph.add_node().id;
	
	tails.push_back(node_map[antecedent]);
	rhs.push_back(postag[antecedent - 1]);
      }
      
      if (node_map[id] == hypergraph_type::invalid)
	node_map[id] = hypergraph.add_node().id;
      
      const symbol_type& lhs = postag[id - 1];
		
      hypergraph_type::edge_type& edge = hypergraph.add_edge(tails.begin(), tails.end());
      edge.rule = hypergraph_type::rule_type::create(hypergraph_type::rule_type(lhs, rhs.begin(), rhs.end()));
      
      hypergraph.connect_edge(edge.id, node_map[id]);
    }

    if (! hypergraph.nodes.empty() && hypergraph.is_valid())
      hypergraph.topologically_sort();
  }
};

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


struct MST
{
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

  
  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_istream is(file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iiter_type iter(is);
    iiter_type iter_end;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    MapFile mapper(map);
    
    mst_parser<iiter_type> parser;
    mst_type mst;
    
    Dependency      deps;
    TransformForest transform_forest(goal, head_mode);
    TransformSpan   transform_span(goal, binarize_mode, category_mode);
    
    while (iter != iter_end) {
      mst.clear();
      
      if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, mst))
	throw std::runtime_error("parsing failed");

      if (mapper) {
	const sentence_type& mapped = mapper();
	
	if (mapped.size() != mst.words.size())
	  throw std::runtime_error("mst size and mapped size differ");
	
	mst.words.assign(mapped.begin(), mapped.end());
      }
      
      if (! mst.verify())
	throw std::runtime_error("invalid mst format");
      
      if (unescape_mode)
	unescape(mst.words.begin(), mst.words.end());
      
      if (normalize_mode) {
	if (relation_mode)
	  normalize(mst.labels.begin(), mst.labels.end());
	else
	  normalize(mst.poss.begin(), mst.poss.end());
      }

      if (forest_mode || span_mode) {
	deps.clear();
	
	deps.sentence.assign(mst.words.begin(), mst.words.end());
	
	if (relation_mode) {
	  mst_type::label_set_type::const_iterator liter_end = mst.labels.end();
	  for (mst_type::label_set_type::const_iterator liter = mst.labels.begin(); liter != liter_end; ++ liter)
	    deps.postag.push_back('[' + *liter + ']');
	} else {
	  mst_type::label_set_type::const_iterator liter_end = mst.poss.end();
	  for (mst_type::label_set_type::const_iterator liter = mst.poss.begin(); liter != liter_end; ++ liter)
	    deps.postag.push_back('[' + *liter + ']');	    
	}
	  
	deps.dependency.assign(mst.positions.begin(), mst.positions.end());
	
	if (forest_mode) {
	  transform_forest(deps.sentence, deps.postag, deps.dependency, deps.hypergraph);
	  
	  os << deps.hypergraph << '\n';
	} else if (span_mode) {
	  transform_span(deps.sentence, deps.postag, deps.dependency, deps.spans);
	  
	  os << deps.spans << '\n';
	}
	
	if (flush_output)
	  os << std::flush;
      } else {
	if (leaf_mode) {
	  if (! karma::generate(oiter, (-(standard::string % ' ') << '\n'), mst.words))
	    throw std::runtime_error("generation failed");
	} else if (leaf_pos_mode) {
	  if (relation_mode) {
	    for (size_t i = 0; i != mst.words.size(); ++ i) {
	      if (i)
		os << ' ';
	      os << mst.words[i] << '_' << mst.labels[i];
	    }
	  } else {
	    for (size_t i = 0; i != mst.words.size(); ++ i) {
	      if (i)
		os << ' ';
	      os << mst.words[i] << '_' << mst.poss[i];
	    }
	  }
	  os << '\n';
	} else if (relation_mode) {
	  if (! karma::generate(oiter, (-(standard::string % ' ')
					<< " ||| " << -(standard::string % ' ')
					<< " ||| " << -(karma::int_ % ' ')
					<< '\n'),
				mst.words, mst.labels, mst.positions))
	    throw std::runtime_error("generation failed");
	} else {
	  if (! karma::generate(oiter, (-(standard::string % ' ')
					<< " ||| " << -(standard::string % ' ')
					<< " ||| " << -(karma::int_ % ' ')
					<< '\n'),
				mst.words, mst.poss, mst.positions))
	    throw std::runtime_error("generation failed");
	}
	
	if (flush_output)
	  os << std::flush;
      }
    }
  }
};

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


struct CoNLL
{
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
    
    boost::spirit::qi::uint_parser<conll_type::size_type, 10, 1, -1> size;
    boost::spirit::qi::rule<Iterator, std::string(), blank_type>    token;
    boost::spirit::qi::rule<Iterator, conll_type(), blank_type>     conll;
    boost::spirit::qi::rule<Iterator, conll_set_type(), blank_type> conlls;
  };

  
  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_istream is(file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iiter_type iter(is);
    iiter_type iter_end;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    MapFile mapper(map);

    conll_parser<iiter_type> parser;
    conll_set_type conll;
    
    Dependency      deps;
    TransformForest transform_forest(goal, head_mode);
    TransformSpan   transform_span(goal, binarize_mode, category_mode);

    while (iter != iter_end) {
      conll.clear();
	
      if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, conll))
	throw std::runtime_error("parsing failed");
      
      if (mapper) {
	const sentence_type& mapped = mapper();

	if (mapped.size() != conll.size())
	  throw std::runtime_error("conll size and mapped size differ");
	
	for (size_t i = 0; i != mapped.size(); ++ i)
	  conll[i].form = mapped[i];
      }


      if (normalize_mode || unescape_mode) {
	conll_set_type::iterator citer_end = conll.end();
	for (conll_set_type::iterator citer = conll.begin(); citer != citer_end; ++ citer) {

	  if (unescape_mode)
	    citer->form = unescape(citer->form);

	  if (normalize_mode) {
	    if (relation_mode)
	      citer->deprel = normalize(citer->deprel);
	    else
	      citer->cpostag = normalize(citer->cpostag);
	  }
	}
      }

      if (forest_mode || span_mode) {
	deps.clear();
	  
	conll_set_type::const_iterator citer_end = conll.end();
	for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer) {
	  deps.sentence.push_back(citer->form);
	    
	  if (relation_mode)
	    deps.postag.push_back('[' + citer->deprel + ']');
	  else
	    deps.postag.push_back('[' + citer->cpostag + ']');
	    
	  if (projective_mode) {
	    const conll_type::size_type head = boost::apply_visitor(conll_type::visitor_phead(), citer->phead);
	    if (head == conll_type::size_type(-1))
	      throw std::runtime_error("invalid projective head");
	      
	    deps.dependency.push_back(head);
	  } else
	    deps.dependency.push_back(citer->head);
	}
	
	if (forest_mode) {
	  transform_forest(deps.sentence, deps.postag, deps.dependency, deps.hypergraph);
	  
	  os << deps.hypergraph << '\n';
	} else if (span_mode) {
	  transform_span(deps.sentence, deps.postag, deps.dependency, deps.spans);
	  
	  os << deps.spans << '\n';
	}
	
	if (flush_output)
	  os << std::flush;
      } else {
	if (leaf_mode) {
	  if (! conll.empty()) {
	    conll_set_type::const_iterator citer_end = conll.end() - 1;
	    for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer)
	      os << citer->form << ' ';
	    os << conll.back().form;
	  }
	  os << '\n';
	} else if (leaf_pos_mode) {
	  
	  conll_set_type::const_iterator citer_begin = conll.begin();
	  conll_set_type::const_iterator citer_end   = conll.end();

	  if (relation_mode) {
	    for (conll_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) {
	      if (citer != citer_begin)
		os << ' ';
	      os << citer->form << '_' << citer->deprel;
	    }
	  } else {
	    for (conll_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) {
	      if (citer != citer_begin)
		os << ' ';
	      os << citer->form << '_' << citer->cpostag;
	    }
	  }
	  os << '\n';
	  
	} else {
	  conll_set_type::const_iterator citer_end = conll.end();
	  for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer)
	    os << citer->form << ' ';
	  os << "||| ";
	  
	  if (relation_mode) {
	    conll_set_type::const_iterator citer_end = conll.end();
	    for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer)
	      os << citer->deprel << ' ';
	  } else {
	    conll_set_type::const_iterator citer_end = conll.end();
	    for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer)
	      os << citer->cpostag << ' ';
	  }
	  os << "|||";
	  
	  if (projective_mode) {
	    conll_set_type::const_iterator citer_end = conll.end();
	    for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer) {
	      const conll_type::size_type head = boost::apply_visitor(conll_type::visitor_phead(), citer->phead);
	      if (head == conll_type::size_type(-1))
		throw std::runtime_error("invalid projective head");
	      
	      os << ' ' << head;
	    }
	  } else {
	    conll_set_type::const_iterator citer_end = conll.end();
	    for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer)
	      os << ' ' << citer->head;
	  }
	  os << '\n';
	}
	if (flush_output)
	  os << std::flush;
      }
    }
  }
};

struct malt_type
{
  typedef int index_type;
  
  index_type   id;
  std::string  word;
  std::string  tag;
  index_type   dep;

  malt_type() {}
  malt_type(const index_type& __id,
	    const std::string& __word,
	    const std::string& __tag,
	    const index_type& __dep)
    : id(__id), word(__word), tag(__tag), dep(__dep) {}
};

BOOST_FUSION_ADAPT_STRUCT(
			  malt_type,
			  (malt_type::index_type, id)
			  (std::string, word)
			  (std::string, tag)
			  (malt_type::index_type, dep)
			  )
struct Malt
{
  typedef std::vector<malt_type, std::allocator<malt_type> > malt_set_type;

  template <typename Iterator>
  struct malt_parser : boost::spirit::qi::grammar<Iterator, malt_set_type()>
  {
    malt_parser() : malt_parser::base_type(malts)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      token %= qi::lexeme[+(standard::char_ - standard::space)];
      
      malt  %= (qi::omit[*standard::blank] >> qi::int_
		>> qi::omit[+standard::blank] >> token
		>> qi::omit[+standard::blank] >> token
		>> qi::omit[+standard::blank] >> qi::int_
		>> qi::omit[*standard::blank] >> qi::eol);
      malts %= *malt >> qi::omit[*standard::blank] >> qi::eol;
    }
    
    typedef boost::spirit::standard::blank_type blank_type;
    
    boost::spirit::qi::rule<Iterator, std::string()>      token;
    boost::spirit::qi::rule<Iterator, malt_type()>     malt;
    boost::spirit::qi::rule<Iterator, malt_set_type()> malts;
  };

  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iiter_type iter(is);
    iiter_type iter_end;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    MapFile mapper(map);

    malt_parser<iiter_type> parser;
    malt_set_type malt;
    
    Dependency      deps;
    TransformForest transform_forest(goal, head_mode);
    TransformSpan   transform_span(goal, binarize_mode, category_mode);
    
    size_t sent_no = 0;
    size_t line_no = 0;
    while (iter != iter_end) {
      malt.clear();
	
      if (! qi::parse(iter, iter_end, parser, malt)) {
	std::string str;
	for (size_t i = 0; i != malt.size(); ++ i)
	  str += ' ' + malt[i].word;
	throw std::runtime_error("parsing failed at # sent: " + utils::lexical_cast<std::string>(sent_no)
				 + " # line: " + utils::lexical_cast<std::string>(line_no) + str);
      }
      
      ++ sent_no;
      line_no += malt.size() + 1;
      
      if (mapper) {
	const sentence_type& mapped = mapper();

	if (mapped.size() != malt.size())
	  throw std::runtime_error("malt size and mapped size differ");
	
	for (size_t i = 0; i != mapped.size(); ++ i)
	  malt[i].word = mapped[i];
      }

      if (normalize_mode || unescape_mode) {
	malt_set_type::iterator citer_end = malt.end();
	for (malt_set_type::iterator citer = malt.begin(); citer != citer_end; ++ citer) {

	  if (unescape_mode)
	    citer->word = unescape(citer->word);

	  if (normalize_mode)
	    citer->tag = normalize(citer->tag);
	}
      }
      
      if (forest_mode || span_mode) {
	deps.clear();
	
	malt_set_type::const_iterator citer_end = malt.end();
	for (malt_set_type::const_iterator citer = malt.begin(); citer != citer_end; ++ citer) {
	  deps.sentence.push_back(citer->word);
	  
	  deps.postag.push_back('[' + citer->tag + ']');
	  
	  // +1....
	  deps.dependency.push_back(citer->dep + 1);
	}
	
	if (forest_mode) {
	  transform_forest(deps.sentence, deps.postag, deps.dependency, deps.hypergraph);
	  
	  os << deps.hypergraph << '\n';
	} else if (span_mode) {
	  transform_span(deps.sentence, deps.postag, deps.dependency, deps.spans);
	  
	  os << deps.spans << '\n';
	}
	
	if (flush_output)
	  os << std::flush;
      } else {
	if (leaf_mode) {
	  if (! malt.empty()) {
	    malt_set_type::const_iterator citer_end = malt.end() - 1;
	    for (malt_set_type::const_iterator citer = malt.begin(); citer != citer_end; ++ citer)
	      os << citer->word << ' ';
	    os << malt.back().word;
	  }
	  os << '\n';
	  
	} else if (leaf_pos_mode) {
	  malt_set_type::const_iterator citer_begin = malt.begin();
	  malt_set_type::const_iterator citer_end   = malt.end();
	  for (malt_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) {
	    if (citer != citer_begin)
	      os << ' ';
	    os << citer->word << '_' << citer->tag;
	  }
	  os << '\n';
	} else {
	  malt_set_type::const_iterator citer_end = malt.end();
	  for (malt_set_type::const_iterator citer = malt.begin(); citer != citer_end; ++ citer)
	    os << citer->word << ' ';
	  os << "||| ";
	  
	  for (malt_set_type::const_iterator citer = malt.begin(); citer != citer_end; ++ citer)
	    os << citer->tag << ' ';
	  os << "|||";
	  
	  for (malt_set_type::const_iterator citer = malt.begin(); citer != citer_end; ++ citer)
	    os << ' ' << (citer->dep + 1);
	  os << '\n';
	}
	if (flush_output)
	  os << std::flush;
      }
    }
  }
};

struct Cabocha
{
  typedef std::pair<std::string, std::string> terminal_type;
  typedef std::vector<terminal_type, std::allocator<terminal_type> > terminal_set_type;
  
  struct node_type
  {
    node_type(int __pos, int __head) : pos(__pos), head(__head) {}
    
    int pos;
    int head;
    
    terminal_set_type terminals;
  };
  typedef std::vector<node_type, std::allocator<node_type> > node_set_type;
  
  typedef std::vector<utils::piece, std::allocator<utils::piece> > tokens_type;
  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
  
  typedef std::vector<int, std::allocator<int> > offset_set_type;

  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);

    MapFile mapper(map);

    node_set_type nodes;
    dependency_type   dependency;
    offset_set_type   offsets;
    terminal_set_type terminals;
    
    std::string line;
    tokens_type tokens;

    Dependency      deps;
    TransformForest transform_forest(goal, head_mode);
    TransformSpan   transform_span(goal, binarize_mode, category_mode);
    
    while (std::getline(is, line)) {
      utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
	
      tokens.clear();
      tokens.insert(tokens.end(), tokenizer.begin(), tokenizer.end());
	
      if (tokens.empty()) continue;
	
      if (tokens.size() == 1) {
	if (tokens.front() != "EOS")
	  throw std::runtime_error("invalid cabocha F1 format: no EOS");
	  
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
	
	if (mapper) {
	  const sentence_type& mapped = mapper();
	  
	  if (terminals.size() != mapped.size())
	    throw std::runtime_error("cabocha size and mapped size differ");
	  
	  for (size_t i = 0; i != mapped.size(); ++ i)
	    terminals[i].first = mapped[i];
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
	  
	if (forest_mode || span_mode) {
	  deps.clear();
	    
	  terminal_set_type::const_iterator titer_end = terminals.end();
	  for (terminal_set_type::const_iterator titer = terminals.begin(); titer != titer_end; ++ titer) {
	    deps.sentence.push_back(titer->first);
	    deps.postag.push_back('[' + titer->second + ']');
	  }
	  
	  deps.dependency.assign(dependency.begin(), dependency.end());
	  
	  if (forest_mode) {
	    transform_forest(deps.sentence, deps.postag, deps.dependency, deps.hypergraph);
	    
	    os << deps.hypergraph << '\n';
	  } else if (span_mode) {
	    transform_span(deps.sentence, deps.postag, deps.dependency, deps.spans);
	    
	    os << deps.spans << '\n';
	  }
	  
	  if (flush_output)
	    os << std::flush;
	} else {
	  if (leaf_mode) {
	    if (! terminals.empty()) {
	      terminal_set_type::const_iterator titer_end = terminals.end() - 1;
	      for (terminal_set_type::const_iterator titer = terminals.begin(); titer != titer_end; ++ titer)
		os << titer->first << ' ';
	      os << terminals.back().first;
	    }
	    os << '\n';
	  } else if (leaf_pos_mode) {
	    terminal_set_type::const_iterator titer_begin = terminals.begin();
	    terminal_set_type::const_iterator titer_end = terminals.end();
	    for (terminal_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	      if (titer != titer_begin)
		os << ' ';
	      os << titer->first << '_' << titer->second;
	    }
	    os << '\n';
	  } else {
	    terminal_set_type::const_iterator titer_end = terminals.end();
	    for (terminal_set_type::const_iterator titer = terminals.begin(); titer != titer_end; ++ titer)
	      os << titer->first << ' ';
	    os << "||| ";
	    for (terminal_set_type::const_iterator titer = terminals.begin(); titer != titer_end; ++ titer)
	      os << titer->second << ' ';
	    os << "||| ";
	    
	    if (! karma::generate(std::ostream_iterator<char>(os), -(karma::int_ % ' ') << '\n', dependency))
	      throw std::runtime_error("generation failed");
	  }
	  
	  if (flush_output)
	    os << std::flush;
	}

	nodes.clear();
	  
      } else if (tokens.size() == 5) {
	if (tokens.front() != "*")
	  throw std::runtime_error("invalid cabocha F1 format: no star");
	  
	const int index = utils::lexical_cast<int>(tokens[1]);
	if (index != static_cast<int>(nodes.size()))
	  throw std::runtime_error("invalid cabocha F1 format: node size do not match: "
				   + tokens[1]
				   + " " + utils::lexical_cast<std::string>(nodes.size()));
	
	int dep;
	int head;
	int func;
	
	{
	  utils::piece::const_iterator iter = tokens[2].begin();
	  utils::piece::const_iterator iter_end = tokens[2].end();
	  
	  if (! qi::parse(iter, iter_end, qi::int_ >> qi::lit('D'), dep) || iter != iter_end)
	    throw std::runtime_error("dep parsing failed?");
	}
	
	{
	  utils::piece::const_iterator iter = tokens[3].begin();
	  utils::piece::const_iterator iter_end = tokens[3].end();
	  
	  if (! qi::parse(iter, iter_end, qi::int_ >> qi::lit('/') >> qi::int_, head, func) || iter != iter_end)
	    throw std::runtime_error("head/func parsing failed?");
	}
	
	nodes.push_back(node_type(dep, utils::bithack::branch(func_mode, func, head)));
	
      } else if (tokens.size() == 3 || tokens.size() == 2) {
	boost::tokenizer<boost::char_separator<char>, utils::piece::const_iterator, utils::piece> tokenizer(tokens[1], boost::char_separator<char>(","));
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
	throw std::runtime_error("invalid cabocha F1 format: # of columns do not match");
    }
  }
};

struct khayashi_type
{
  typedef size_t size_type;
  typedef std::string label_type;
  
  typedef std::vector<label_type, std::allocator<label_type> > label_set_type;
  typedef std::vector<size_type, std::allocator<size_type> >   position_set_type;
  
  
  label_set_type words;
  label_set_type poss;
  position_set_type positions;
  
  
  bool verify() const
  {
    return words.size() == poss.size() && words.size() == positions.size();
  }
  
  void clear()
  {
    words.clear();
    poss.clear();
    positions.clear();
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  khayashi_type,
			  (khayashi_type::label_set_type,    words)
			  (khayashi_type::label_set_type,    poss)
			  (khayashi_type::position_set_type, positions)
			  )


struct KHayashi
{

  template <typename Iterator>
  struct khayashi_parser : boost::spirit::qi::grammar<Iterator, khayashi_type(), boost::spirit::standard::blank_type>
  {
    khayashi_parser() : khayashi_parser::base_type(khayashi)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      word %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];

      words %= +word;
      poss  %= +word;
      positions %= +qi::int_;
      
      khayashi %= (words >> "|||" >> poss >> qi::eol
		   >> positions >> qi::eol
		   >> qi::omit[*(positions >> qi::eol)]
		   >> (qi::eol || qi::eoi));
    }
    
    typedef boost::spirit::standard::blank_type blank_type;
    
    boost::spirit::qi::rule<Iterator, std::string(), blank_type> word;
    boost::spirit::qi::rule<Iterator, khayashi_type::label_set_type(), blank_type> words;
    boost::spirit::qi::rule<Iterator, khayashi_type::label_set_type(), blank_type> poss;
    boost::spirit::qi::rule<Iterator, khayashi_type::position_set_type(), blank_type> positions;
    
    boost::spirit::qi::rule<Iterator, khayashi_type(), blank_type> khayashi;
  };
  
  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_istream is(file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iiter_type iter(is);
    iiter_type iter_end;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    MapFile mapper(map);

    khayashi_parser<iiter_type> parser;
    khayashi_type khayashi;

    Dependency      deps;
    TransformForest transform_forest(goal, head_mode);
    TransformSpan   transform_span(goal, binarize_mode, category_mode);
    
    size_t num = 0;
    while (iter != iter_end) {
      khayashi.clear();
	
      if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, khayashi))
	throw std::runtime_error("parsing failed");
      
      if (mapper) {
	const sentence_type& mapped = mapper();
	
	if (mapped.size() != khayashi.words.size())
	  throw std::runtime_error("khayashi size and mapped size differ: " + utils::lexical_cast<std::string>(num));
	
	khayashi.words.assign(mapped.begin(), mapped.end());
      }
	
      if (! khayashi.verify())
	throw std::runtime_error("invalid khayashi format");
      
      if (unescape_mode)
	unescape(khayashi.words.begin(), khayashi.words.end());
	
      if (normalize_mode)
	normalize(khayashi.poss.begin(), khayashi.poss.end());

      if (forest_mode || span_mode) {
	deps.clear();
	
	deps.sentence.assign(khayashi.words.begin(), khayashi.words.end());
	  
	khayashi_type::label_set_type::const_iterator liter_end = khayashi.poss.end();
	for (khayashi_type::label_set_type::const_iterator liter = khayashi.poss.begin(); liter != liter_end; ++ liter)
	  deps.postag.push_back('[' + *liter + ']');	    
	  
	deps.dependency.assign(khayashi.positions.begin(), khayashi.positions.end());
	
	if (forest_mode) {
	  transform_forest(deps.sentence, deps.postag, deps.dependency, deps.hypergraph);
	  
	  os << deps.hypergraph << '\n';
	} else if (span_mode) {
	  transform_span(deps.sentence, deps.postag, deps.dependency, deps.spans);
	  
	  os << deps.spans << '\n';
	}
	
	if (flush_output)
	  os << std::flush;
      } else {
	
	if (leaf_mode) {
	  if (! karma::generate(oiter, (-(standard::string % ' ') << '\n'), khayashi.words))
	    throw std::runtime_error("generation failed");
	} else if (leaf_pos_mode) {
	  for (size_t i = 0; i != khayashi.words.size(); ++ i) {
	    if (i)
	      os << ' ';
	    os << khayashi.words[i] << '_' << khayashi.poss[i];
	  }
	  os << '\n';
	} else {
	  if (! karma::generate(oiter, (-(standard::string % ' ')
					<< " ||| " << -(standard::string % ' ')
					<< " ||| " << -(karma::int_ % ' ')
					<< '\n'),
				khayashi.words, khayashi.poss, khayashi.positions))
	    throw std::runtime_error("generation failed");
	}
	  
	if (flush_output)
	  os << std::flush;
      }

      ++ num;
    }
  }
};

struct khayashi_forest_type
{
  typedef size_t  size_type;
  typedef int32_t id_type;
  typedef int32_t head_type;
  typedef double  score_type;
  
  typedef std::string label_type;
  typedef std::vector<label_type, std::allocator<label_type> > label_set_type;
  
  
  struct edge_type
  {
    id_type    node1;
    id_type    node2;
    head_type  head;
    score_type score;
  };

  typedef std::vector<edge_type, std::allocator<edge_type> > edge_set_type;
  
  struct node_type
  {
    id_type id;
    id_type first;
    id_type last;
    id_type parent;
    id_type child_left;
    id_type child_right;
    id_type child_left2;
    id_type child_right2;
    
    edge_set_type edges;

    bool is_terminal() const { return edges.empty(); }
  };

  typedef std::vector<node_type, std::allocator<node_type> > node_set_type;

  bool verify() const
  {
    return words.size() == poss.size();
  }

  void clear()
  {
    words.clear();
    poss.clear();
    nodes.clear();
  }
  
  label_set_type words;
  label_set_type poss;
  node_set_type  nodes;
};

BOOST_FUSION_ADAPT_STRUCT(
			  khayashi_forest_type::edge_type,
			  (khayashi_forest_type::id_type,    node1)
			  (khayashi_forest_type::id_type,    node2)
			  (khayashi_forest_type::head_type,  head)
			  (khayashi_forest_type::score_type, score)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  khayashi_forest_type::node_type,
			  (khayashi_forest_type::id_type,     id)
			  (khayashi_forest_type::id_type,     first)
			  (khayashi_forest_type::id_type,     last)
			  (khayashi_forest_type::id_type,     parent)
			  (khayashi_forest_type::id_type,     child_left)
			  (khayashi_forest_type::id_type,     child_right)
			  (khayashi_forest_type::id_type,     child_left2)
			  (khayashi_forest_type::id_type,     child_right2)
			  (khayashi_forest_type::edge_set_type, edges)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  khayashi_forest_type,
			  (khayashi_forest_type::label_set_type, words)
			  (khayashi_forest_type::label_set_type, poss)
			  (khayashi_forest_type::node_set_type, nodes)
			  )
struct KHayashiForest
{
  template <typename Iterator>
  struct khayashi_parser : boost::spirit::qi::grammar<Iterator, khayashi_forest_type()>
  {
    khayashi_parser() : khayashi_parser::base_type(khayashi)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      token %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];
      
      words %= token % (+standard::blank);
      poss  %= token % (+standard::blank);
      
      edge %= (qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::lit("|||")
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::double_
	       >> qi::omit[*standard::blank >> qi::eol]);
      
      node %= (qi::int_
	       >> qi::omit[+standard::blank] >> qi::int_ >> qi::lit('-') >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::lit("|||")
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[+standard::blank] >> qi::int_
	       >> qi::omit[*qi::blank >> qi::eol]
	       >> *edge);
      
      khayashi %= (qi::omit[*standard::blank] >> words >> qi::omit[+standard::blank]
		   >> qi::lit("|||")
		   >>  qi::omit[+standard::blank] >> poss >> qi::omit[*standard::blank >> qi::eol]
		   >> *node
		   >> qi::omit[qi::eol]);
    }
    
    boost::spirit::qi::rule<Iterator, std::string()> token;
    
    boost::spirit::qi::rule<Iterator, khayashi_forest_type::label_set_type()> words;
    boost::spirit::qi::rule<Iterator, khayashi_forest_type::label_set_type()> poss;
    
    boost::spirit::qi::rule<Iterator, khayashi_forest_type::edge_type()> edge;
    boost::spirit::qi::rule<Iterator, khayashi_forest_type::node_type()> node;
    
    boost::spirit::qi::rule<Iterator, khayashi_forest_type()> khayashi;
  };

  struct remove_nodes
  {
    typedef cicada::Symbol     symbol_type;

    bool operator()(const symbol_type& x) const
    {
      if (! x.is_non_terminal()) return false;
      
      const symbol_type non_terminal = x.non_terminal();
      
      return (non_terminal[non_terminal.size() - 2] == '*') || x.binarized();
    }
  };

  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    typedef size_t size_type;
    
    typedef cicada::Sentence   sentence_type;
    typedef cicada::Rule       rule_type;
    typedef cicada::Symbol     symbol_type;
    typedef cicada::Feature    feature_type;
    typedef cicada::Attribute  attribute_type;
    
    typedef std::pair<hypergraph_type::id_type, hypergraph_type::id_type> node_id_type;
    typedef std::vector<node_id_type, std::allocator<node_id_type> > node_id_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > label_set_type;
    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > goal_id_set_type;

    const feature_type feature("depdency-parse");

    //const attribute_type attr_dependency_head("dependency-head");
    //const attribute_type attr_dependency_dependent("dependency-dependent");
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_istream is(file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iiter_type iter(is);
    iiter_type iter_end;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);

    MapFile mapper(map);
    
    khayashi_parser<iiter_type> parser;
    khayashi_forest_type khayashi;

    hypergraph_type hypergraph;
    node_id_set_type nodes;
    node_id_set_type terminals;
    label_set_type   labels;
    goal_id_set_type goals;

    while (iter != iter_end) {
      khayashi.clear();
      hypergraph.clear();
	
      if (! qi::parse(iter, iter_end, parser, khayashi))
	throw std::runtime_error("parsing failed");
      
      if (mapper) {
	const MapFile::sentence_type& mapped = mapper();
	
	if (mapped.size() != khayashi.words.size())
	  throw std::runtime_error("khayashi size and mapped size differ");
	
	khayashi.words.assign(mapped.begin(), mapped.end());
      }

      if (! khayashi.verify())
	throw std::runtime_error("invalid khayashi forest format");
	
      if (unescape_mode)
	unescape(khayashi.words.begin(), khayashi.words.end());
	
      if (normalize_mode)
	normalize(khayashi.poss.begin(), khayashi.poss.end());
	
      if (khayashi.nodes.empty()) {
	os << hypergraph << '\n';
	if (flush_output)
	  os << std::flush;
      } else {
	nodes.clear();
	nodes.resize(khayashi.nodes.size() + 1, std::make_pair(hypergraph_type::invalid, hypergraph_type::invalid));
	  
	terminals.clear();
	terminals.resize(khayashi.words.size(), std::make_pair(hypergraph_type::invalid, hypergraph_type::invalid));
	  
	labels.clear();
	goals.clear();
	  
	rule_type::symbol_set_type rhs(2);
	hypergraph_type::edge_type::node_set_type tails(2);
	  
	khayashi_forest_type::node_set_type::const_iterator niter_end = khayashi.nodes.end();
	for (khayashi_forest_type::node_set_type::const_iterator niter = khayashi.nodes.begin(); niter != niter_end; ++ niter) {
	  const khayashi_forest_type::node_type& node = *niter;
	    
	  if (node.is_terminal()) {
	    // we will remove terminal sharing code, since we need to differentiate head and dependent...??
	      
	    if (terminals[node.parent].first == hypergraph_type::invalid) {
	      const hypergraph_type::id_type head_id = hypergraph.add_node().id;
	      const hypergraph_type::id_type node_id = hypergraph.add_node().id;
		
	      terminals[node.parent].first  = head_id;
	      terminals[node.parent].second = node_id;
		
	      if (head_id >= labels.size())
		labels.resize(head_id + 1);		
	      if (node_id >= labels.size())
		labels.resize(node_id + 1);
		
	      labels[head_id] = '[' + khayashi.poss[node.parent] + "*]";
	      labels[node_id] = '[' + khayashi.poss[node.parent] + ']';
		
	      hypergraph_type::edge_type& edge = hypergraph.add_edge();
		
	      edge.rule = rule_type::create(rule_type(labels[head_id], rule_type::symbol_set_type(1, khayashi.words[node.parent])));
		
	      hypergraph.connect_edge(edge.id, head_id);
		
	      hypergraph_type::edge_type& edge2 = hypergraph.add_edge(&head_id, (&head_id) + 1);
		
	      edge2.rule = rule_type::create(rule_type(labels[node_id], rule_type::symbol_set_type(1, labels[head_id])));
		
	      hypergraph.connect_edge(edge2.id, node_id);
	    }
	      
	    nodes[node.id] = terminals[node.parent];
	  } else {
	    const hypergraph_type::id_type node_id   = hypergraph.add_node().id;
	    const hypergraph_type::id_type binary_id = hypergraph.add_node().id;
	      
	    nodes[node.id].first = node_id;
	    nodes[node.id].second = binary_id;
	      
	    if (node_id >= labels.size())
	      labels.resize(node_id + 1);
	    if (binary_id >= labels.size())
	      labels.resize(binary_id + 1);
	      
	    labels[node_id]   = '[' + khayashi.poss[node.parent] + ']';
	    labels[binary_id] = '[' + khayashi.poss[node.parent] + "^]";

	    if (node.first == 0 && node.last == static_cast<int>(khayashi.words.size()))
	      goals.push_back(node_id);
	      
	    khayashi_forest_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	    for (khayashi_forest_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	      const khayashi_forest_type::node_type& node1 = khayashi.nodes[eiter->node1 - 1];
	      const khayashi_forest_type::node_type& node2 = khayashi.nodes[eiter->node2 - 1];

	      if (node1.is_terminal())
		tails[0] = utils::bithack::branch(node.parent == node1.parent, nodes[node1.id].first, nodes[node1.id].second);
	      else
		tails[0] = utils::bithack::branch(node.parent != node1.parent, nodes[node1.id].first, nodes[node1.id].second);
		
	      if (node2.is_terminal())
		tails[1] = utils::bithack::branch(node.parent == node2.parent, nodes[node2.id].first, nodes[node2.id].second);
	      else
		tails[1] = utils::bithack::branch(node.parent != node2.parent, nodes[node2.id].first, nodes[node2.id].second);
		
	      rhs[0] = labels[tails[0]];
	      rhs[1] = labels[tails[1]];
		
	      hypergraph_type::edge_type& edge1 = hypergraph.add_edge(tails.begin(), tails.end());
	      hypergraph_type::edge_type& edge2 = hypergraph.add_edge(tails.begin(), tails.end());
		
	      edge1.rule = rule_type::create(rule_type(labels[node_id], rhs));
	      edge2.rule = rule_type::create(rule_type(labels[binary_id], rhs));
		
	      edge1.features[feature] = eiter->score;
	      edge2.features[feature] = eiter->score;

	      hypergraph.connect_edge(edge1.id, node_id);
	      hypergraph.connect_edge(edge2.id, binary_id);
	    }
	  }
	}
	  
	if (! goals.empty()) {
	  hypergraph.goal = hypergraph.add_node().id;

	  goal_id_set_type::const_iterator giter_end = goals.end();
	  for (goal_id_set_type::const_iterator giter = goals.begin(); giter != giter_end; ++ giter) {
	    const hypergraph_type::id_type tail = *giter;
	      
	    hypergraph_type::edge_type& edge = hypergraph.add_edge(&tail, (&tail) + 1);
	      
	    edge.rule = rule_type::create(rule_type(goal, rule_type::symbol_set_type(1, labels[tail])));
	      
	    hypergraph.connect_edge(edge.id, hypergraph.goal);
	  }
	    
	  cicada::topologically_sort(hypergraph);
	    
	  if (head_mode)
	    cicada::debinarize(hypergraph);
	  else
	    cicada::remove_non_terminal(hypergraph, remove_nodes());
	}
	  
	os << hypergraph << '\n';
	if (flush_output)
	  os << std::flush;
      }
    }
  }
};

struct cicada_type
{
  typedef int index_type;
  typedef std::string label_type;
  
  typedef std::vector<label_type, std::allocator<label_type> > label_set_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;
  
  label_set_type tok;
  label_set_type pos;
  index_set_type dep;
  
  cicada_type() : tok(), pos(), dep() {}
  
  void clear()
  {
    tok.clear();
    pos.clear();
    dep.clear();
  }

  bool verify() const
  {
    return tok.size() == pos.size() && tok.size() == dep.size();
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  cicada_type,
			  (cicada_type::label_set_type, tok)
			  (cicada_type::label_set_type, pos)
			  (cicada_type::index_set_type, dep)
			  )

struct Cicada
{
  
  template <typename Iterator>
  struct cicada_parser : boost::spirit::qi::grammar<Iterator, cicada_type(), boost::spirit::standard::blank_type>
  {
    cicada_parser() : cicada_parser::base_type(cicada)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      word %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];
      
      tok %= +word;
      pos %= +word;
      dep %= +qi::int_;
      
      cicada %= tok >> "|||" >> pos >> "|||" >> dep >> (qi::eol || qi::eoi);
    }
    
    typedef boost::spirit::standard::blank_type blank_type;
    
    boost::spirit::qi::rule<Iterator, std::string(), blank_type> word;
    boost::spirit::qi::rule<Iterator, cicada_type::label_set_type(), blank_type> tok;
    boost::spirit::qi::rule<Iterator, cicada_type::label_set_type(), blank_type> pos;
    boost::spirit::qi::rule<Iterator, cicada_type::index_set_type(), blank_type> dep;
    
    boost::spirit::qi::rule<Iterator, cicada_type(), blank_type> cicada;
  };

  
  void operator()(const path_type& file, const path_type& map, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_istream is(file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iiter_type iter(is);
    iiter_type iter_end;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    MapFile mapper(map);

    cicada_parser<iiter_type> parser;
    cicada_type cicada;

    Dependency      deps;
    TransformForest transform_forest(goal, head_mode);
    TransformSpan   transform_span(goal, binarize_mode, category_mode);
    
    size_t num = 0;
    while (iter != iter_end) {
      cicada.clear();
      
      if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, cicada))
	throw std::runtime_error("parsing failed");
      
      if (mapper) {
	const sentence_type& mapped = mapper();
	
	if (mapped.size() != cicada.tok.size())
	  throw std::runtime_error("cicada size and mapped size differ: " + utils::lexical_cast<std::string>(num));
	
	cicada.tok.assign(mapped.begin(), mapped.end());
      }
	
      if (! cicada.verify())
	throw std::runtime_error("invalid cicada format");
      
      if (unescape_mode)
	unescape(cicada.tok.begin(), cicada.tok.end());

      if (normalize_mode)
	normalize(cicada.pos.begin(), cicada.pos.end());

      if (forest_mode || span_mode) {
	deps.clear();
	deps.sentence.assign(cicada.tok.begin(), cicada.tok.end());
	  
	cicada_type::label_set_type::const_iterator liter_end = cicada.pos.end();
	for (cicada_type::label_set_type::const_iterator liter = cicada.pos.begin(); liter != liter_end; ++ liter)
	  deps.postag.push_back('[' + *liter + ']');	    
	
	deps.dependency.assign(cicada.dep.begin(), cicada.dep.end());
	
	if (forest_mode) {
	  transform_forest(deps.sentence, deps.postag, deps.dependency, deps.hypergraph);
	  
	  os << deps.hypergraph << '\n';
	} else if (span_mode) {
	  transform_span(deps.sentence, deps.postag, deps.dependency, deps.spans);
	  
	  os << deps.spans << '\n';
	}
	
	if (flush_output)
	  os << std::flush;
      } else {
	
	if (leaf_mode) {
	  if (! karma::generate(oiter, (-(standard::string % ' ') << '\n'), cicada.tok))
	    throw std::runtime_error("generation failed");
	} else if (leaf_pos_mode) {
	  for (size_t i = 0; i != cicada.tok.size(); ++ i) {
	    if (i)
	      os << ' ';
	    os << cicada.tok[i] << '_' << cicada.pos[i];
	  }
	  os << '\n';
	} else {
	  if (! karma::generate(oiter, (-(standard::string % ' ')
					<< " ||| " << -(standard::string % ' ')
					<< " ||| " << -(karma::int_ % ' ')
					<< '\n'),
				cicada.tok, cicada.pos, cicada.dep))
	    throw std::runtime_error("generation failed");
	}
	  
	if (flush_output)
	  os << std::flush;
      }

      ++ num;
    }    
  }
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",      po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",     po::value<path_type>(&output_file)->default_value(output_file), "output file")
    ("map",        po::value<path_type>(&map_file),                                "map file")
    
    ("goal",         po::value<std::string>(&goal)->default_value(goal),                 "goal symbol")
    ("non-terminal", po::value<std::string>(&non_terminal)->default_value(non_terminal), "non-terminal symbol")
    
    ("mst",             po::bool_switch(&mst_mode),             "tranform MST dependency")
    ("conll",           po::bool_switch(&conll_mode),           "tranform CoNLL dependency")
    ("malt",            po::bool_switch(&malt_mode),            "tranform Malt dependency")
    ("cabocha",         po::bool_switch(&cabocha_mode),         "tranform Cabocha dependency")
    ("khayashi",        po::bool_switch(&khayashi_mode),        "tranform KHayashi dependency")
    ("khayashi-forest", po::bool_switch(&khayashi_forest_mode), "tranform KHayashi Forest dependency")
    ("cicada",          po::bool_switch(&cicada_mode),          "tranform CICADA dependency")
    
    ("projective", po::bool_switch(&projective_mode), "project into projective dependency")
    ("relation",   po::bool_switch(&relation_mode),   "assing relation to POS")
    ("unescape",   po::bool_switch(&unescape_mode),   "unescape terminal symbols, such as -LRB-, \\* etc.")
    ("normalize",  po::bool_switch(&normalize_mode),  "normalize category, such as [,] [.] etc.")
    ("func",       po::bool_switch(&func_mode),       "use function word as a head in bunsetsu")
    ("forest",     po::bool_switch(&forest_mode),     "output as a forest")
    ("head",       po::bool_switch(&head_mode),       "output hypergraph with explicit head")
    ("span",       po::bool_switch(&span_mode),       "output as a set of spans")
    ("binarize",   po::bool_switch(&binarize_mode),   "output binarized spans")
    ("category",   po::bool_switch(&category_mode),   "output spans with category labels")
    ("leaf",       po::bool_switch(&leaf_mode),       "output leafs (temrinals)")
    ("leaf-pos",   po::bool_switch(&leaf_pos_mode),   "output leaf-pos (temrinal+POS)")
            
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
