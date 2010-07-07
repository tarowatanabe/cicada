#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

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

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>

#include "cicada/hypergraph.hpp"

#include "utils/program_options.hpp"
#include "utils/icu_filter.hpp"
#include "utils/compress_stream.hpp"
#include "utils/space_separator.hpp"


#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>
#include <unicode/translit.h>
#include <unicode/regex.h>

struct TransLit
{
  TransLit(const std::string& name, 
	   const std::string& pattern) { initialize(name, pattern); }
  
  TransLit(const std::string& name) { initialize(name); }
  
  void operator()(UnicodeString& data) { trans->transliterate(data); }
  
  void initialize(const std::string& name)
  {
    UErrorCode status = U_ZERO_ERROR;
    trans.reset(Transliterator::createInstance(UnicodeString::fromUTF8(name),
					       UTRANS_FORWARD, status));
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
  }
  
  void initialize(const std::string& name, const std::string& pattern)
  {
    UErrorCode status = U_ZERO_ERROR;
    UParseError status_parse;
    trans.reset(Transliterator::createFromRules(UnicodeString::fromUTF8(name), UnicodeString::fromUTF8(pattern),
						UTRANS_FORWARD, status_parse, status));
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
  }
  
  boost::shared_ptr<Transliterator> trans;
};

class Replace
{
  
public:
  Replace(const char* pattern, const char* subst) : matcher() { initialize(pattern, subst); }
  Replace(const UnicodeString& pattern, const UnicodeString& subst) : matcher() { initialize(pattern, subst); }
  ~Replace() { clear(); }


  const UnicodeString& operator()(UnicodeString& uline)
  {
    UErrorCode status = U_ZERO_ERROR;
    matcher->reset(uline);
    uline = matcher->replaceAll(substitute, status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
    return uline;
  }

  void clear()
  {
    if (matcher) delete matcher;
    matcher = 0;
  }
  
  void initialize(const char* pattern, const char* subst)
  {
    initialize(UnicodeString(pattern, "utf-8"), UnicodeString(subst, "utf-8"));
  }
  
  
  void initialize(const UnicodeString& pattern, const UnicodeString& subst)
  {
    clear();
    
    UErrorCode status = U_ZERO_ERROR;
    matcher = new RegexMatcher(pattern, 0, status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
    substitute = subst;
  }
  

private:
  RegexMatcher* matcher;
  UnicodeString substitute;
};

class ReplaceAll
{
  
public:
  ReplaceAll(const char* pattern, const char* subst) : matcher() { initialize(pattern, subst); }
  ReplaceAll(const UnicodeString& pattern, const UnicodeString& subst) : matcher() { initialize(pattern, subst); }
  ~ReplaceAll() { clear(); }


  const UnicodeString& operator()(UnicodeString& uline)
  {
    
    while (1) {
      matcher->reset(uline);
      if (! matcher->find()) break;
      
      UErrorCode status = U_ZERO_ERROR;
      uline = matcher->replaceAll(substitute, status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
    }
    return uline;
  }

  void clear()
  {
    if (matcher) delete matcher;
    matcher = 0;
  }
  
  void initialize(const char* pattern, const char* subst)
  {
    initialize(UnicodeString(pattern, "utf-8"), UnicodeString(subst, "utf-8"));
  }
  
  
  void initialize(const UnicodeString& pattern, const UnicodeString& subst)
  {
    clear();
    
    UErrorCode status = U_ZERO_ERROR;
    matcher = new RegexMatcher(pattern, 0, status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
    substitute = subst;
  }
  

private:
  RegexMatcher* matcher;
  UnicodeString substitute;
};

typedef boost::filesystem::path path_type;

// tree-bank parser...

struct treebank_type
{
  typedef std::vector<treebank_type> antecedents_type;

  std::string cat;
  antecedents_type antecedents;
  
  treebank_type() {}
  treebank_type(const std::string& __cat) : cat(__cat) {}

  void clear()
  {
    cat.clear();
    antecedents.clear();
  }
};


BOOST_FUSION_ADAPT_STRUCT(
			  treebank_type,
			  (std::string, cat)
			  (std::vector<treebank_type>, antecedents)
			  )


template <typename Iterator>
struct penntreebank_escaped_grammar : boost::spirit::qi::grammar<Iterator, treebank_type(), boost::spirit::standard::space_type>
{
  penntreebank_escaped_grammar() : penntreebank_escaped_grammar::base_type(treebank)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    using qi::lit;
    using qi::lexeme;
    using qi::hold;
    using qi::repeat;
    using qi::attr;
    using qi::on_error;
    using qi::int_;
    using qi::double_;
    using standard::char_;
    using standard::space;

    escaped_char.add
      ("\\/",   '/')
      ("\\*",   '*');
    
    escaped_word.add
      ("-LRB-", "(")
      ("-RRB-", ")")
      ("-LSB-", "[")
      ("-RSB-", "]")
      ("-LCB-", "{")
      ("-RCB-", "}");
    
    cat %= lexeme[escaped_word | +(escaped_char | (char_ - space - '(' - ')'))];
    treebank %= hold['(' >> cat >> +treebank >> ')'] | cat;
  }
  
  boost::spirit::qi::symbols<char, char>        escaped_char;
  boost::spirit::qi::symbols<char, const char*> escaped_word;
  
  boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type>   cat;
  boost::spirit::qi::rule<Iterator, treebank_type(), boost::spirit::standard::space_type> treebank;
};

template <typename Iterator>
struct penntreebank_grammar : boost::spirit::qi::grammar<Iterator, treebank_type(), boost::spirit::standard::space_type>
{
  penntreebank_grammar() : penntreebank_grammar::base_type(treebank)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    using qi::lit;
    using qi::lexeme;
    using qi::hold;
    using qi::repeat;
    using qi::attr;
    using qi::on_error;
    using qi::int_;
    using qi::double_;
    using standard::char_;
    using standard::space;
    
    cat %= lexeme[+(char_ - space - '(' - ')')];
    treebank %= hold['(' >> cat >> +treebank >> ')'] | cat;
  }
  
  boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type>   cat;
  boost::spirit::qi::rule<Iterator, treebank_type(), boost::spirit::standard::space_type> treebank;
};

typedef cicada::HyperGraph hypergraph_type;
typedef std::vector<std::string, std::allocator<std::string> > sentence_type;


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
  rule += " |||";
  
  hypergraph_type::edge_type& edge = graph.add_edge(nodes.begin(), nodes.end());
  edge.rule.reset(new rule_type(rule));
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

void transform(const treebank_type& treebank, sentence_type& sent) 
{
  if (treebank.antecedents.empty())
    sent.push_back(treebank.cat);
  else
    for (treebank_type::antecedents_type::const_iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter)
      transform(*aiter, sent);
}

void transform(treebank_type& treebank)
{
  if (treebank.antecedents.empty()) {
    
    
    
  } else {
    treebank_type::antecedents_type antecedents;
  
    for (treebank_type::antecedents_type::iterator aiter = treebank.antecedents.begin(); aiter != treebank.antecedents.end(); ++ aiter) {
      if (aiter->antecedents.empty()) {
	UnicodeString ucat = UnicodeString::fromUTF8(aiter->cat);
    
	ucat = ' ' + ucat + ' ';
	
	static ReplaceAll replace("(?<=[[:White_Space:]])([[:^Numeric_Type=None:][\\u4E07][\\u842C][\\u4EBF][\\u5104][\\u5146]])(?=[[:^Numeric_Type=None:][\\u4E07][\\u842C][\\u4EBF][\\u5104][\\u5146]]+[[:White_Space:]])", "$1 ");
	
	replace(ucat);
	ucat.trim();
	
	std::string categories;
	StringByteSink<std::string> sink(&categories);
	ucat.toUTF8(sink);
	
	typedef boost::tokenizer<utils::space_separator> tokenizer_type;
	typedef std::vector<std::string> tokens_type;
	
	tokenizer_type tokenizer(categories);
	tokens_type tokens(tokenizer.begin(), tokenizer.end());
	
	if (tokens.size() == 1) {
	  antecedents.push_back(*aiter);
	  antecedents.back().cat = tokens.front();
	} else {
	  for (tokens_type::iterator titer = tokens.begin(); titer != tokens.end(); ++ titer) {
	    antecedents.push_back(*aiter);
	    antecedents.back().cat = treebank.cat + "-SPLIT";
	    antecedents.back().antecedents.push_back(*titer);
	  }
	}
      } else {
	antecedents.push_back(*aiter);
	transform(antecedents.back());
      }
    }
    
    treebank.antecedents = antecedents;
  }
}

path_type input_file = "-";
path_type output_file = "-";

std::string codepage = "utf-8";

bool escaped = false;
bool leaf = false;

bool split_digits = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);


    typedef boost::spirit::istream_iterator iter_type;

    boost::iostreams::filtering_istream is;
    is.push(utils::icu_filter(codepage, "utf-8", utils::icu_filter::stop));
    utils::push_compress_istream(is, input_file, 1024 * 1024);
    
    is.unsetf(std::ios::skipws);

    boost::iostreams::filtering_ostream os;
    os.push(utils::icu_filter("utf-8", codepage, utils::icu_filter::stop));
    utils::push_compress_ostream(os, output_file, 1024 * 1024);
    
    penntreebank_grammar<iter_type>         grammar;
    penntreebank_escaped_grammar<iter_type> grammar_escaped;

    treebank_type   parsed;
    hypergraph_type graph;
    sentence_type   sent;

    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      parsed.clear();
      
      if (escaped) {
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, grammar_escaped, boost::spirit::standard::space, parsed))
	  throw std::runtime_error("parsing failed");
      } else {
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, grammar, boost::spirit::standard::space, parsed))
	  throw std::runtime_error("parsing failed");
      }

      if (split_digits)
	transform(parsed);
      
      if (leaf) {
	sent.clear();
	
	transform(parsed, sent);
	
	if (! sent.empty()) {
	  std::copy(sent.begin(), sent.end() - 1, std::ostream_iterator<std::string>(os, " "));
	  os << sent.back();
	  os << '\n';
	} else
	  os << '\n';
	
      } else {
	graph.clear();
	
	transform(parsed, graph);

	graph.topologically_sort();
	
	os << graph << '\n';
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
    
    ("codepage",  po::value<std::string>(&codepage)->default_value(codepage),     "codepage")
    
    ("escape",    po::bool_switch(&escaped), "escape English penntreebank")
    ("leaf",      po::bool_switch(&leaf),    "collect leaf nodes only")
    
    ("split-digits", po::bool_switch(&split_digits), "split digits")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
