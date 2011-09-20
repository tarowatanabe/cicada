//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// kbest filter which transforms moses style features into cicada-style features.
//
// id ||| sentence ||| features etc.
//
// TODO: add subprocess filter so that we can uncover the kbest-format again...
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <deque>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"

typedef boost::filesystem::path path_type;

typedef size_t size_type;
typedef std::vector<std::string, std::allocator<std::string> > tokens_type;

typedef boost::fusion::tuple<size_type, tokens_type, tokens_type, tokens_type> kbest_type;

template <typename Iterator>
struct kbest_parser : boost::spirit::qi::grammar<Iterator, kbest_type(), boost::spirit::standard::blank_type>
{
  kbest_parser() : kbest_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens   %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    features %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    remains  %= *qi::lexeme[+(standard::char_ - standard::space)];
    
    kbest %= size >> "|||" >> tokens >> -("|||" >> features) >> -("|||" >> remains) >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> features;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  
  boost::spirit::qi::rule<Iterator, kbest_type(), blank_type>  kbest;
};

path_type input_file = "-";
path_type output_file = "-";

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    typedef boost::spirit::istream_iterator iter_type;
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    kbest_parser<iter_type> parser;
    iter_type iter(is);
    iter_type iter_end;
    
    kbest_type kbest;
    tokens_type features_new;
    std::string feature_name;
    
    while (iter != iter_end) {
      boost::fusion::get<1>(kbest).clear();
      boost::fusion::get<2>(kbest).clear();
      boost::fusion::get<3>(kbest).clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	if (iter != iter_end)
	  throw std::runtime_error("kbest parsing failed");
      
      os << boost::fusion::get<0>(kbest) << " ||| ";
      
      const tokens_type& tokens = boost::fusion::get<1>(kbest);
      
      if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -(boost::spirit::standard::string % ' '), tokens))
	throw std::runtime_error("tokens generation failed...?");
      
      const tokens_type& features = boost::fusion::get<2>(kbest);
      
      
      if (! features.empty()) {
	features_new.clear();
	feature_name = "";
	int id = 0;

	tokens_type::const_iterator fiter_end = features.end();
	for (tokens_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter) {
	  const size_t feature_size = fiter->size();
	  
	  if (fiter->operator[](feature_size - 1) == ':') {
	    feature_name = *fiter;
	    id = 0;
	  } else {
	    features_new.push_back(feature_name + utils::lexical_cast<std::string>(id) + "=" + *fiter);
	    ++ id;
	  }
	}
	
	if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), " ||| " << -(boost::spirit::standard::string % ' '), features_new))
	  throw std::runtime_error("tokens generation failed...?");
      }
      
      const tokens_type& remains = boost::fusion::get<3>(kbest);
      
      if (! remains.empty())
	if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), " ||| " << -(boost::spirit::standard::string % ' '), remains))
	  throw std::runtime_error("tokens generation failed...?");
      
      os << '\n';
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
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output")
    
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
