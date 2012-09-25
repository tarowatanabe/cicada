//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include <cicada/hypergraph.hpp>

#include "utils/bithack.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/space_separator.hpp"

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
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"] >> (qi::eoi | qi::eol);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, sentence_type(), blank_type> tokens;
};

path_set_type input_files;
path_set_type source_files;
path_set_type target_files;
path_set_type alignment_files;
path_set_type dependency_files;

path_type list_file;
path_type list_source_file;
path_type list_target_file;
path_type list_alignment_file;
path_type list_dependency_file;

std::string goal = "[s]";
std::string non_terminal = "[x]";

bool bilingual_mode = false;
bool mst_mode = false;
bool conll_mode = false;
bool cabocha_mode = false;
bool forest_mode = false;

bool projective_mode = false;
bool relation_mode = false;
bool hypergraph_mode = false;
bool head_mode = false;

// dependency output
path_type output_file = "-";

void read_list(const path_type& path, path_set_type& files);
void options(int argc, char** argv);


template <typename Dep>
void apply(const path_set_type& files, const path_type& output)
{
  Dep dep;
  dep(files, output);
}

struct MST;
struct CoNLL;
struct Cabocha;
struct Forest;

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (int(bilingual_mode) + mst_mode + conll_mode + cabocha_mode + forest_mode == 0)
      throw std::runtime_error("one of bilingual/mst/conll mode");

    if (int(bilingual_mode) + mst_mode + conll_mode + cabocha_mode + forest_mode > 1)
      throw std::runtime_error("one of bilingual/mst/conll mode");
    
    if (mst_mode) {
      read_list(list_file, input_files);
      if (input_files.empty())
	input_files.push_back("-");
      
      apply<MST>(input_files, output_file);
    } else if (conll_mode) {
      read_list(list_file, input_files);
      if (input_files.empty())
	input_files.push_back("-");
      
      apply<CoNLL>(input_files, output_file);
    } else if (cabocha_mode) {
      read_list(list_file, input_files);
      if (input_files.empty())
	input_files.push_back("-");
      
      apply<Cabocha>(input_files, output_file);
    } else if (forest_mode) {
      read_list(list_file, input_files);
      if (input_files.empty())
	input_files.push_back("-");
      
      apply<Forest>(input_files, output_file);
    } else if (bilingual_mode) {
      read_list(list_source_file, source_files);
      read_list(list_target_file, target_files);
      read_list(list_alignment_file, alignment_files);
      read_list(list_dependency_file, dependency_files);
      
      if (source_files.empty())
	source_files.push_back("-");
      if (target_files.empty())
	target_files.push_back("-");
      if (alignment_files.empty())
	alignment_files.push_back("-");
      if (dependency_files.empty())
	dependency_files.push_back("-");
      
      if (source_files.size() != target_files.size())
	throw std::runtime_error("# of files do not match");
      
      if (source_files.size() != alignment_files.size())
	throw std::runtime_error("# of alignment files do not match");
      if (target_files.size() != alignment_files.size())
	throw std::runtime_error("# of alignment files do not match");
      if (source_files.size() != dependency_files.size())
	throw std::runtime_error("# of dependency files do not match");
      if (target_files.size() != dependency_files.size())
	throw std::runtime_error("# of dependency files do not match");
      
      utils::compress_ostream os(output_file, 1024 * 1024);
      
      typedef boost::spirit::istream_iterator iiter_type;
      
      sentence_parser<iiter_type>    parser;
      
      for (size_t i = 0; i != source_files.size(); ++ i) {
	utils::compress_istream is_src(source_files[i], 1024 * 1024);
	utils::compress_istream is_trg(target_files[i], 1024 * 1024);
	utils::compress_istream is_align(alignment_files[i], 1024 * 1024);
	utils::compress_istream is_dep(dependency_files[i], 1024 * 1024);
	
	is_src.unsetf(std::ios::skipws);
	is_trg.unsetf(std::ios::skipws);
	
	iiter_type siter(is_src);
	iiter_type titer(is_trg);
	iiter_type siter_end;
	iiter_type titer_end;
	
	sentence_type   source;
	sentence_type   target;
	alignment_type  alignment;
	dependency_type dependency;
	dependency_type projected;
	
	for (size_t line_no = 0; siter != siter_end && titer != titer_end; ++ line_no) {
	  source.clear();
	  target.clear();
	  alignment.clear();
	  dependency.clear();
	  
	  if (! boost::spirit::qi::phrase_parse(siter, siter_end, parser, boost::spirit::standard::blank, source))
	    throw std::runtime_error("source sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	  if (! boost::spirit::qi::phrase_parse(titer, titer_end, parser, boost::spirit::standard::blank, target))
	    throw std::runtime_error("target sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	  
	  if (! (is_align >> alignment))
	    throw std::runtime_error("no alignment?");
	  if (! (is_dep >> dependency))
	    throw std::runtime_error("no dependency?");
	  
	  const int source_size = source.size();
	  const int target_size = target.size();
	  
	  if (source_size != static_cast<int>(dependency.size()))
	    throw std::runtime_error("source size do not match with dependency size at # " + utils::lexical_cast<std::string>(line_no));
	  
	  if (source_size == 0 || target_size == 0) {
	    os << '\n';
	    continue;
	  };
	  
	  projected.clear();
	  projected.resize(target_size, -1);
	  
	  // perform projection....
	  
	  
	  os << projected << '\n';
	}
      }
    } else
      throw std::runtime_error("one of bilingual|mst|conll|cabocha|forest?");
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}


void read_list(const path_type& path, path_set_type& files)
{
  if (path.empty()) return;
  if (path != "-" && ! boost::filesystem::exists(path))
    throw std::runtime_error("no file? " + path.string());
  
  utils::compress_istream is(path);
  std::string file;
  while (std::getline(is, file)) {
    if (file.empty()) continue;
    if (! boost::filesystem::exists(file))
      throw std::runtime_error("no file? " + file);
    files.push_back(file);
  }
}

struct Transform
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

  // assign here for faster memory access!
  sentence_type   sentence;
  sentence_type   pos;
  dependency_type dependency;
  hypergraph_type hypergraph;

  Transform(const symbol_type& __goal,
	    const bool __head_mode)
    : goal(__goal),
      head_mode(__head_mode) {} 

  void clear()
  {
    sentence.clear();
    pos.clear();
    dependency.clear();
    hypergraph.clear();
  }
  
  void operator()()
  {
    if (sentence.size() != pos.size() || sentence.size() != dependency.size())
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
      rhs.push_back(pos[antecedent - 1]);
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
	rhs.push_back(pos[antecedent - 1]);
      }
      
      if (head_mode) {
	const symbol_type lhs = '[' + pos[id - 1].non_terminal_strip() + "*]";
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
	rhs.push_back(pos[antecedent - 1]);
      }
      
      if (node_map[id] == hypergraph_type::invalid)
	node_map[id] = hypergraph.add_node().id;
      
      const symbol_type& lhs = pos[id - 1];
		
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

  
  void operator()(const path_set_type& files, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    mst_parser<iiter_type> parser;
    mst_type mst;

    Transform transform(goal, head_mode);
    
    path_set_type::const_iterator fiter_end = files.end();
    for (path_set_type::const_iterator fiter = files.begin(); fiter != fiter_end; ++ fiter) {
      utils::compress_istream is(*fiter, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iiter_type iter(is);
      iiter_type iter_end;

      while (iter != iter_end) {
	mst.clear();
	
	if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, mst))
	  throw std::runtime_error("parsing failed");
	
	if (! mst.verify())
	  throw std::runtime_error("invalid mst format");

	if (hypergraph_mode) {
	  transform.clear();
	  transform.sentence.assign(mst.words.begin(), mst.words.end());
	  
	  if (relation_mode) {
	    mst_type::label_set_type::const_iterator liter_end = mst.labels.end();
	    for (mst_type::label_set_type::const_iterator liter = mst.labels.begin(); liter != liter_end; ++ liter)
	      transform.pos.push_back('[' + *liter + ']');
	  } else {
	    mst_type::label_set_type::const_iterator liter_end = mst.poss.end();
	    for (mst_type::label_set_type::const_iterator liter = mst.poss.begin(); liter != liter_end; ++ liter)
	      transform.pos.push_back('[' + *liter + ']');	    
	  }
	  
	  transform.dependency.assign(mst.positions.begin(), mst.positions.end());
	  
	  transform();
	  
	  os << transform.hypergraph << '\n';
	} else {
	  if (relation_mode) {
	    if (! karma::generate(oiter, (-(standard::string % ' ')
					  << " ||| " << -(standard::string % ' ')
					  << " ||| " << -(karma::int_ % ' ')
					  << '\n'),
				  mst.words, mst.labels, mst.positions))
	      throw std::runtime_error("generation failed");
	    else
	      if (! karma::generate(oiter, (-(standard::string % ' ')
					    << " ||| " << -(standard::string % ' ')
					    << " ||| " << -(karma::int_ % ' ')
					    << '\n'),
				    mst.words, mst.poss, mst.positions))
		throw std::runtime_error("generation failed");
	  }
	}
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

  
  void operator()(const path_set_type& files, const path_type& output)
  {
    typedef boost::spirit::istream_iterator iiter_type;
    
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    conll_parser<iiter_type> parser;
    conll_set_type conll;
    
    Transform transform(goal, head_mode);

    path_set_type::const_iterator fiter_end = files.end();
    for (path_set_type::const_iterator fiter = files.begin(); fiter != fiter_end; ++ fiter) {
      utils::compress_istream is(*fiter, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iiter_type iter(is);
      iiter_type iter_end;
      
      while (iter != iter_end) {
	conll.clear();
	
	if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, conll))
	  throw std::runtime_error("parsing failed");

	if (hypergraph_mode) {
	  transform.clear();
	  
	  conll_set_type::const_iterator citer_end = conll.end();
	  for (conll_set_type::const_iterator citer = conll.begin(); citer != citer_end; ++ citer) {
	    transform.sentence.push_back(citer->form);
	    
	    if (relation_mode)
	      transform.pos.push_back('[' + citer->deprel + ']');
	    else
	      transform.pos.push_back('[' + citer->cpostag + ']');
	    
	    if (projective_mode) {
	      const conll_type::size_type head = boost::apply_visitor(conll_type::visitor_phead(), citer->phead);
	      if (head == conll_type::size_type(-1))
		throw std::runtime_error("invalid projective head");
	      
	      transform.dependency.push_back(head);
	    } else
	      transform.dependency.push_back(citer->head);
	  }
	  
	  transform();
	  
	  os << transform.hypergraph << '\n';
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
  
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;
  typedef std::vector<std::string, std::allocator<std::string> > tokens_type;
  
  typedef std::vector<int, std::allocator<int> > offset_set_type;

  void operator()(const path_set_type& files, const path_type& output)
  {
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    utils::compress_ostream os(output_file, 1024 * 1024);

    node_set_type nodes;
    dependency_type   dependency;
    offset_set_type   offsets;
    terminal_set_type terminals;
    
    std::string line;
    tokens_type tokens;

    Transform transform(goal, head_mode);
    
    path_set_type::const_iterator fiter_end = files.end();
    for (path_set_type::const_iterator fiter = files.begin(); fiter != fiter_end; ++ fiter) {
      utils::compress_istream is(*fiter, 1024 * 1024);
      
      while (std::getline(is, line)) {
	tokenizer_type tokenizer(line);
	
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
	  
	  if (hypergraph_mode) {
	    transform.clear();
	    
	    terminal_set_type::const_iterator titer_end = terminals.end();
	    for (terminal_set_type::const_iterator titer = terminals.begin(); titer != titer_end; ++ titer) {
	      transform.sentence.push_back(titer->first);
	      transform.pos.push_back('[' + titer->second + ']');
	    }
	    
	    transform.dependency.assign(dependency.begin(), dependency.end());
	    
	    transform();
	    
	    os << transform.hypergraph << '\n';
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
	  
	} else if (tokens.size() == 5) {
	  if (tokens.front() != "*")
	    throw std::runtime_error("invalid cabocha F1 format: no star");
	  
	  const int index = utils::lexical_cast<int>(tokens[1]);
	  if (index != static_cast<int>(nodes.size()))
	    throw std::runtime_error("invalid cabocha F1 format: node size do not match");
	  
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
	  throw std::runtime_error("invalid cabocha F1 format: # of columns do not match");
      }
    }
  }
};

struct Forest
{
  
  void operator()(const path_set_type& files, const path_type& output)
  {
    
    
  }
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",      po::value<path_set_type>(&input_files)->multitoken(),           "input file(s)")
    ("output",     po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("source",     po::value<path_set_type>(&source_files)->multitoken(),     "source file(s)")
    ("target",     po::value<path_set_type>(&target_files)->multitoken(),     "target file(s)")
    ("alignment",  po::value<path_set_type>(&alignment_files)->multitoken(),  "alignment file(s)")
    ("dependency", po::value<path_set_type>(&dependency_files)->multitoken(), "dependency file(s)")
    
    ("list",            po::value<path_type>(&list_file),            "list file")
    ("list-source",     po::value<path_type>(&list_source_file),     "source list file")
    ("list-target",     po::value<path_type>(&list_target_file),     "target list file")
    ("list-alignment",  po::value<path_type>(&list_alignment_file),  "alignment list file")
    ("list-dependency", po::value<path_type>(&list_dependency_file), "dependency list file")

    ("goal",         po::value<std::string>(&goal)->default_value(goal),                 "goal symbol")
    ("non-terminal", po::value<std::string>(&non_terminal)->default_value(non_terminal), "non-terminal symbol")
    
    ("bilingual",  po::bool_switch(&bilingual_mode),  "project source dependency into target dependency")
    ("mst",        po::bool_switch(&mst_mode),        "tranform MST dependency")
    ("conll",      po::bool_switch(&conll_mode),      "tranform CoNLL dependency")
    ("cabocha",    po::bool_switch(&cabocha_mode),    "tranform Cabocha dependency")
    ("forest",     po::bool_switch(&forest_mode),     "tranform FOREST dependency")
    ("projective", po::bool_switch(&projective_mode), "project into projective dependency")
    ("relation",   po::bool_switch(&relation_mode),   "assing relation to POS")
    ("hypergraph", po::bool_switch(&hypergraph_mode), "output as a hypergraph")
    ("head",       po::bool_switch(&head_mode),       "output hypergraph with explicit head")
            
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
