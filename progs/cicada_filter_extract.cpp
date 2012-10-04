//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada_filter_extract_impl.hpp"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>
#include <cfloat>
#include <cmath>
#include <queue>

#include <boost/program_options.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>

typedef boost::filesystem::path path_type;

typedef PhrasePair phrase_pair_type;

struct score_phrase_pair_type
{
  double score;
  std::string line;

  score_phrase_pair_type(const double& __score, const std::string& __line)
    : score(__score), line(__line) {}

  friend
  bool operator<(const score_phrase_pair_type& x, const score_phrase_pair_type& y)
  {
    return x.score < y.score;
  }
  
  friend
  bool operator>(const score_phrase_pair_type& x, const score_phrase_pair_type& y)
  {
    return x.score > y.score;
  }
  
  void swap(score_phrase_pair_type& x)
  {
    std::swap(score, x.score);
    line.swap(x.line);
  }
};

namespace std
{
  inline
  void swap(score_phrase_pair_type& x, score_phrase_pair_type& y)
  {
    x.swap(y);
  }
};

typedef std::vector<score_phrase_pair_type, std::allocator<score_phrase_pair_type> > heap_type;

path_type input_file = "-";
path_type output_file = "-";

int buffer_size = 1024 * 1024;
int nbest = 100;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (nbest <= 0)
      throw std::runtime_error("nbest must be positive...");
    
    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, buffer_size);
    
    std::string          source_prev;
    phrase_pair_type     phrase_pair;
    
    heap_type heap;
    
    PhrasePairParser    parser;
    std::string line;
    
    while (std::getline(is, line)) {
      if (! parser(line, phrase_pair)) continue;
      if (phrase_pair.counts.empty()) continue;
      
      if (phrase_pair.source != source_prev) {
	if (! heap.empty()) {
	  heap_type::iterator iter_begin = heap.begin();
	  heap_type::iterator iter       = heap.end();
	  
	  for (int k = 0; k != nbest && iter_begin != iter; ++ k, -- iter) {
	    os << iter_begin->line << '\n';
	    std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
	  }
	  
	  if (iter != iter_begin && iter != heap.end()) {
	    const double threshold = iter->score;
	    
	    for (/**/; iter_begin != iter && iter_begin->score == threshold; -- iter) {
	      os << iter_begin->line << '\n';
	      std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
	    }
	  }
	}
	
	heap.clear();
	source_prev = phrase_pair.source;
      }
      
      heap.push_back(score_phrase_pair_type(phrase_pair.counts.front(), line));
      std::push_heap(heap.begin(), heap.end(), std::less<score_phrase_pair_type>());
    }
    
    if (! heap.empty()) {
      heap_type::iterator iter_begin = heap.begin();
      heap_type::iterator iter       = heap.end();
      
      for (int k = 0; k != nbest && iter_begin != iter; ++ k, -- iter) {
	os << iter_begin->line << '\n';
	std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
      }
      
      if (iter != iter_begin && iter != heap.end()) {
	const double threshold = iter->score;
	
	for (/**/; iter_begin != iter && iter_begin->score == threshold; -- iter) {
	  os << iter_begin->line << '\n';
	  std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
	}
      }
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("nbest", po::value<int>(&nbest)->default_value(nbest), "nbest of pairs (wrt to joint-count)")
    
    ("buffer", po::value<int>(&buffer_size)->default_value(buffer_size), "buffer size")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
