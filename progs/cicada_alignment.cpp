// cicada alignment tool:
//
// we can produce intersection/union/grow-{,diag}-{final,final-and}/source/target/itg/max-match from GIZA++ alingment
//        invert alignment
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/spirit/include/support_istream_iterator.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>

#include <cicada/alignment.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;

typedef cicada::Alignment alignment_type;

struct BitextGiza
{
  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::pair<std::string, point_set_type> word_align_type;
  typedef std::vector<word_align_type, std::allocator<word_align_type> > word_align_set_type;
  typedef std::vector<std::string, std::allocator<std::string> > sent_type;
  
  std::string         comment;
  sent_type           target;
  word_align_set_type source;
  
  BitextGiza() {}

  void clear()
  {
    comment.clear();
    target.clear();
    source.clear();
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  BitextGiza,
			  (std::string,                     comment)
			  (BitextGiza::sent_type,           target)
			  (BitextGiza::word_align_set_type, source)
			  )

typedef BitextGiza bitext_giza_type;


path_type source_target_file;
path_type target_source_file;
path_type span_source_file;
path_type span_target_file;
path_type input_file;
path_type output_file = "-";

bool source_target_mode = false;
bool target_source_mode = false;

bool itg_mode = false;
bool max_match_mode = false;

bool intersection_mode = false;
bool union_mode = false;
bool grow_mode = false;
bool final_mode = false;
bool diag_mode = false;
bool final_and_mode = false;
bool invert_mode = false;

int threads = 1;
int debug = 0;

void process(std::istream& is, std::ostream& os);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    process(std::cin, std::cout);
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Iterator>
struct bitext_giza_parser : boost::spirit::qi::grammar<Iterator, bitext_giza_type(), boost::spirit::standard::blank_type>
{
  bitext_giza_parser() : bitext_giza_parser::base_type(bitext)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    comment %= qi::lexeme[*(standard::char_ - qi::eol)] >> qi::eol;
    target  %= *qi::lexeme[+(standard::char_ - standard::space - qi::eol)] >> qi::eol;
    
    points     %= "({" >> *qi::int_ >> "})";
    word_align %= qi::lexeme[+(standard::char_ - standard::space - qi::eol)] >> points;
    source     %= *word_align >> (qi::eol | qi::eoi);
    
    bitext %= comment >> target >> source;
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>                           comment;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::sent_type(), blank_type>           target;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::point_set_type(), blank_type>      points;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::word_align_type(), blank_type>     word_align;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::word_align_set_type(), blank_type> source;
  
  boost::spirit::qi::rule<Iterator, bitext_giza_type(), blank_type> bitext;
};

void process(std::istream& is, 
	     std::ostream& os)
{
  typedef boost::spirit::istream_iterator iter_type;
  
  bitext_giza_parser<iter_type> parser;
  
  bitext_giza_type bitext_giza;
  
  is.unsetf(std::ios::skipws);
  
  iter_type iter(is);
  iter_type iter_end;
  
  while (iter != iter_end) {
    bitext_giza.clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, bitext_giza))
      if (iter != iter_end)
	throw std::runtime_error("parsing failed");
    
    os << "comment: " << bitext_giza.comment << '\n';
    os << "target: ";
    std::copy(bitext_giza.target.begin(), bitext_giza.target.end(), std::ostream_iterator<std::string>(os, " "));
    os << '\n';
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source-target", po::value<path_type>(&source_target_file), "P(target | source) viterbi output")
    ("target-source", po::value<path_type>(&target_source_file), "P(source | target) viterbi output")
    ("span-source",   po::value<path_type>(&span_source_file),   "source span data")
    ("span-target",   po::value<path_type>(&span_target_file),   "target span data")
    ("input",         po::value<path_type>(&input_file),                      "input alignement")
    ("output",        po::value<path_type>(&output_file)->default_value("-"), "output alignment")
    
    ("f2e", po::bool_switch(&source_target_mode), "source target")
    ("e2f", po::bool_switch(&target_source_mode), "target source")
    
    ("itg",          po::bool_switch(&itg_mode),          "itg")
    ("max-match",    po::bool_switch(&max_match_mode),    "max-match")
    ("intersection", po::bool_switch(&intersection_mode), "intersection")
    ("union",        po::bool_switch(&union_mode),        "union")
    ("grow",         po::bool_switch(&grow_mode),         "grow")
    ("diag",         po::bool_switch(&diag_mode),         "diag")
    ("final",        po::bool_switch(&final_mode),        "final")
    ("final-and",    po::bool_switch(&final_and_mode),    "final-and")
    ("invert",       po::bool_switch(&invert_mode),       "invert alignment")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
