//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// transform moses, cdec style reference translations into cicada indexed format

#include <iostream>
#include <vector>
#include <string>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>

#include <boost/fusion/adapted.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef std::string sentence_type;
typedef std::vector<sentence_type, std::allocator<sentence_type> >         sentence_set_type;
typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_map_type;

path_set_type input_files;
path_type     output_file = "-";
std::string   prefix_file;

bool cicada_mode = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (input_files.empty())
      input_files.push_back("-");

    if (cicada_mode) {
      if (prefix_file.empty())
	throw std::runtime_error("prefix must be non-empty");
      
      sentence_map_type refsets;
      
      for (path_set_type::const_iterator iter = input_files.begin(); iter != input_files.end(); ++ iter) {
	typedef boost::spirit::istream_iterator iter_type;
	  
	utils::compress_istream is(*iter, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	int seg;
	std::string sentence;
	
	while (iter != iter_end) {
	  namespace qi = boost::spirit::qi;
	  namespace standard = boost::spirit::standard;
	  
	  sentence.clear();
	  if (! qi::phrase_parse(iter, iter_end, qi::int_ >> "|||" >> qi::lexeme[*(standard::char_ - qi::eol - qi::eoi)] >> (qi::eol | qi::eoi), standard::blank, seg, sentence))
	    if (iter != iter_end)
	      throw std::runtime_error("refset parsing failed");
	  
	  if (seg >= refsets.size())
	    refsets.resize(seg + 1);
	  refsets[seg].push_back(sentence);
	}
      }
      
      size_t num_ref = 0;
      for (size_t seg = 0; seg != refsets.size(); ++ seg)
	num_ref = utils::bithack::max(num_ref, refsets[seg].size());
      
      for (size_t i = 0; i != num_ref; ++ i) {
	utils::compress_ostream os(prefix_file + "." + utils::lexical_cast<std::string>(i));
	
	for (size_t seg = 0; seg != refsets.size(); ++ seg) {
	  if (i < refsets[seg].size())
	    os << refsets[seg][i] << '\n';
	  else
	    os << '\n';
	}
	  
      }
      
    } else {
      sentence_map_type refsets;
      
      for (path_set_type::const_iterator iter = input_files.begin(); iter != input_files.end(); ++ iter) {
	utils::compress_istream is(*iter, 1024 * 1024);
	std::string line;
	
	for (size_t seg = 0; std::getline(is, line); ++ seg) {
	  if (seg >= refsets.size())
	    refsets.resize(seg + 1);
	  
	  refsets[seg].push_back(line);
	}
      }
      
      utils::compress_ostream os(output_file, 1024 * 1024);
      
      for (size_t seg = 0; seg != refsets.size(); ++ seg) {
	sentence_set_type::const_iterator iter_end = refsets[seg].end();
	for (sentence_set_type::const_iterator iter = refsets[seg].begin(); iter != iter_end; ++ iter)
	  if (! iter->empty())
	    os << seg << " ||| " << *iter << '\n';
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
    ("input",  po::value<path_set_type>(&input_files)->multitoken(),   "input file(s)")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output")
    ("prefix", po::value<std::string>(&prefix_file), "prefix")
    ("cicada", po::bool_switch(&cicada_mode),      "refset in cicada foramt and dump in moses/cdec format")
    
    ("help", "help message");
  
  po::positional_options_description pos;
  pos.add("input", -1); // all the files

  po::command_line_parser parser(argc, argv);
  parser.style(po::command_line_style::unix_style & (~po::command_line_style::allow_guessing));
  parser.options(desc);
  parser.positional(pos);
  
  po::variables_map vm;
  po::store(parser.run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options] reference translations" << '\n' << desc << '\n';
    exit(0);
  }
}
