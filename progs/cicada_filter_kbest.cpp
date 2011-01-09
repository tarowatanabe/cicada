//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// kbest filter:
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

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/subprocess.hpp"
#include "utils/lockfree_list_queue.hpp"

typedef boost::filesystem::path path_type;

typedef size_t size_type;
typedef std::vector<std::string, std::allocator<std::string> > tokens_type;
typedef boost::fusion::tuple<size_type, tokens_type, tokens_type> kbest_type;

template <typename Iterator>
struct kbest_parser : boost::spirit::qi::grammar<Iterator, kbest_type(), boost::spirit::standard::blank_type>
{
  kbest_parser() : kbest_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    remains %= *qi::lexeme[+(standard::char_ - standard::space)];
    
    kbest %= size >> "|||" >> tokens >> -("|||" >> remains) >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  boost::spirit::qi::rule<Iterator, kbest_type(), blank_type>  kbest;
};

struct Task
{
  typedef utils::lockfree_list_queue<kbest_type, std::allocator<kbest_type> > queue_type;
  typedef utils::subprocess subprocess_type;

  Task(queue_type&      __queue,
       subprocess_type& __subprocess,
       std::ostream&    __os,
       const bool __id_mode,
       const bool __features_mode)
    : queue(__queue),
      subprocess(__subprocess),
      os(__os),
      id_mode(__id_mode),
      features_mode(__features_mode) {}
  
  queue_type& queue;
  subprocess_type& subprocess;
  std::ostream& os;
  
  bool id_mode;
  bool features_mode;

  void operator()()
  {
    boost::iostreams::filtering_istream is;
#if BOOST_VERSION >= 104400
    is.push(boost::iostreams::file_descriptor_source(subprocess.desc_read(), boost::iostreams::close_handle));
#else
    is.push(boost::iostreams::file_descriptor_source(subprocess.desc_read(), true));
#endif
    
    kbest_type kbest;
    std::string line;
    while (queue.pop(kbest)) {
      if (boost::fusion::get<0>(kbest) == size_type(-1)) break;
      
      if (id_mode)
	os << boost::fusion::get<0>(kbest) << " ||| ";
      
      // read from is...
      std::getline(is, line);
      os << line;
      
      if (features_mode) {
	const tokens_type& features = boost::fusion::get<2>(kbest);
	
	os << " ||| ";
	if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -((+boost::spirit::standard::char_) % ' '), features))
	  throw std::runtime_error("tokens generation failed...?");
      }
      
      os << '\n';
    }
  }
};

path_type input_file = "-";
path_type output_file = "-";

std::string filter;

bool id_mode = false;
bool features_mode = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    typedef boost::spirit::istream_iterator iter_type;
    
    kbest_parser<iter_type> parser;
    iter_type iter(is);
    iter_type iter_end;
    
    kbest_type kbest;
    
    if (! filter.empty()) {
      typedef Task task_type;
      
      task_type::subprocess_type subprocess(filter);
      task_type::queue_type      queue;
      
      boost::thread thread(task_type(queue, subprocess, os, id_mode, features_mode));

      {
	boost::iostreams::filtering_ostream os_filter;
#if BOOST_VERSION >= 104400
	os_filter.push(boost::iostreams::file_descriptor_sink(subprocess.desc_write(), boost::iostreams::close_handle));
#else
	os_filter.push(boost::iostreams::file_descriptor_sink(subprocess.desc_write(), true));
#endif
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  const tokens_type& tokens = boost::fusion::get<1>(kbest);
	  
	  if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os_filter), -((+boost::spirit::standard::char_) % ' '), tokens))
	    throw std::runtime_error("tokens generation failed...?");
	  os_filter << '\n';
	  
	  queue.push(kbest);
	}
      }
      
      boost::fusion::get<0>(kbest) = size_type(-1);
      queue.push(kbest);
      
      thread.join();
    } else {
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
      
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
      
	if (id_mode)
	  os << boost::fusion::get<0>(kbest) << " ||| ";
      
	const tokens_type& tokens = boost::fusion::get<1>(kbest);
      
	if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -((+boost::spirit::standard::char_) % ' '), tokens))
	  throw std::runtime_error("tokens generation failed...?");
      
	if (features_mode) {
	  const tokens_type& features = boost::fusion::get<2>(kbest);
	
	  os << " ||| ";
	  if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -((+boost::spirit::standard::char_) % ' '), features))
	    throw std::runtime_error("tokens generation failed...?");
	}

	os << '\n';
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
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output")
    
    ("filter", po::value<std::string>(&filter), "filter for sentences")

    ("id",       po::bool_switch(&id_mode),       "output id")
    ("features", po::bool_switch(&features_mode), "output features")
    
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
