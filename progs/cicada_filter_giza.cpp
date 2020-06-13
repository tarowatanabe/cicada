//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// giza alignment conversion

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

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

#include <utils/compress_stream.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;

struct BitextGiza
{
  typedef std::set<int, std::less<int>, std::allocator<int> > point_set_type;
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

  void swap(BitextGiza& x)
  {
    comment.swap(x.comment);
    target.swap(x.target);
    source.swap(x.source);
  }
};

namespace std
{
  inline
  void swap(BitextGiza& x, BitextGiza& y)
  {
    x.swap(y);
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  BitextGiza,
			  (std::string,                     comment)
			  (BitextGiza::sent_type,           target)
			  (BitextGiza::word_align_set_type, source)
			  )

typedef BitextGiza bitext_giza_type;

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

path_type input_file = "-";
path_type output_file = "-";
path_type output_source_file;
path_type output_target_file;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    typedef boost::spirit::istream_iterator iter_type;
    
    bitext_giza_parser<iter_type> parser_giza;
    
    bitext_giza_type bitext;
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    std::unique_ptr<std::ostream> os_source(! output_source_file.empty() ?
					  new utils::compress_ostream(output_source_file, 1024 * 1024) : 0);
    std::unique_ptr<std::ostream> os_target(! output_target_file.empty() ?
					  new utils::compress_ostream(output_target_file, 1024 * 1024) : 0);
    
    
    utils::compress_istream is(input_file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      bitext.clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser_giza, boost::spirit::standard::blank, bitext))
	if (iter != iter_end)
	  throw std::runtime_error("parsing failed");
      
      bool initial = true;
      for (size_t src = 1; src < bitext.source.size(); ++ src) {
	const bitext_giza_type::point_set_type& aligns = bitext.source[src].second;
	
	bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
	for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) {
	  if (! initial)
	    os << ' ';
	  initial = false;

	  os << (*titer - 1) << '-' << (src - 1);
	}
      }
      os << '\n';
      
      if (os_source.get()) {
	bool initial = true;
	for (size_t src = 1; src < bitext.source.size(); ++ src) {
	  if (! initial)
	    *os_source << ' ';
	  initial = false;
	  
	  *os_source << bitext.source[src].first;
	}
	*os_source << '\n';
      }
      
      if (os_target.get()) {
	if (! bitext.target.empty()) {
	  std::copy(bitext.target.begin(), bitext.target.end() - 1, std::ostream_iterator<std::string>(*os_target, " "));
	  *os_target << bitext.target.back();
	}
	*os_target << '\n';
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
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",         po::value<path_type>(&input_file)->default_value(input_file),   "input giza alignment")
    ("output",        po::value<path_type>(&output_file)->default_value(output_file), "output alignment")
    ("output-source", po::value<path_type>(&output_source_file),                      "output source")
    ("output-target", po::value<path_type>(&output_target_file),                      "output target")
    
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
