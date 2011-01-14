//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// filter for kbest-output of parseIt
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

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/space_separator.hpp"
#include "utils/chart.hpp"

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
};

typedef std::vector<category_type, std::allocator<category_type> > category_set_type;

typedef boost::fusion::tuple<category_type, category_set_type, double> item_type;

typedef std::vector<item_type, std::allocator<item_type> > item_set_type;
typedef std::vector<std::string, std::allocator<std::string> > sentence_type;

struct forest_type
{
  sentence_type sentence;
  item_set_type items;

  forest_type()
    : sentence(), items() {}
  forest_type(const sentence_type& __sentence, const item_set_type& __items)
    : sentence(__sentence), items(__items) {}

  void clear()
  {
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
    
    cat %= qi::lexeme[+(standard::char_ - standard::space - '[') - "|||"];
    
    category %= qi::hold[cat >> '[' >> qi::int_ >> ',' >> qi::int_ >> ']'] | cat;
    
    item %= category >> "=>" >> (+category) >> "|||" >> qi::double_ >> -qi::lit("EXTRAVAL") >> qi::eol;
    sentence %= *cat >> qi::eol;
    
    forest %= (-sentence) >> (*item) >> qi::eol;
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>   cat;
  boost::spirit::qi::rule<Iterator, category_type(), blank_type> category;
  
  boost::spirit::qi::rule<Iterator, item_type(), blank_type> item;
  boost::spirit::qi::rule<Iterator, sentence_type(), blank_type>   sentence;

  boost::spirit::qi::rule<Iterator, forest_type(), blank_type> forest;
};


path_type input_file = "-";
path_type output_file = "-";
path_type map_file;

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
    
    forest_parser<iter_type> parser;

    forest_type forest;

    iter_type iter(is);
    iter_type iter_end;

    int num = 0;
    while (iter != iter_end) {
      forest.clear();

      if (debug)
	std::cerr << "parsing: " << num << std::endl;

      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, forest))
	throw std::runtime_error("parsing failed");
      
      ++ num;
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
