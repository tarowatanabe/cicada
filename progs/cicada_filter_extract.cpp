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
#include <map>

#include <boost/program_options.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/sgi_hash_set.hpp>

typedef boost::filesystem::path path_type;

typedef PhrasePair phrase_pair_type;
typedef std::multimap<double, std::string, std::greater<double>,
		      std::allocator<std::pair<const double, std::string> > > phrase_pair_set_type;

path_type input_file = "-";
path_type output_file = "-";

int nbest = 100;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (nbest <= 0)
      throw std::runtime_error("nbest must be positive...");
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    std::string          source_prev;
    phrase_pair_type     phrase_pair;
    phrase_pair_set_type phrase_pairs;
    
    PhrasePairParser    parser;
    std::string line;
    
    while (std::getline(is, line)) {
      if (! parser(line, phrase_pair)) continue;
      if (phrase_pair.counts.empty()) continue;
      
      if (phrase_pair.source != source_prev) {
	if (! phrase_pairs.empty()) {
	  phrase_pair_set_type::const_iterator iter = phrase_pairs.begin();
	  phrase_pair_set_type::const_iterator iter_end = phrase_pairs.end();
	  phrase_pair_set_type::const_iterator iter_prev = iter_end;
	  for (int k = 0; k < nbest && iter != iter_end; ++ k, ++ iter) {
	    os << iter->second << '\n';
	    iter_prev = iter;
	  }
	  
	  if (iter != iter_end && iter_prev != iter_end) {
	    iter_end = phrase_pairs.upper_bound(iter_prev->first);
	    for (/**/; iter != iter_end; ++ iter)
	      os << iter->second << '\n';
	  }
	}
	phrase_pairs.clear();
	source_prev = phrase_pair.source;
      }
      
      const double& count = phrase_pair.counts.front();
      if (static_cast<int>(phrase_pairs.size()) < nbest || count >= (-- phrase_pairs.end())->first)
	phrase_pairs.insert(std::make_pair(count, line));
    }
    
    if (! phrase_pairs.empty()) {
      phrase_pair_set_type::const_iterator iter = phrase_pairs.begin();
      phrase_pair_set_type::const_iterator iter_end = phrase_pairs.end();
      phrase_pair_set_type::const_iterator iter_prev = iter_end;
      for (int k = 0; k < nbest && iter != iter_end; ++ k, ++ iter) {
	os << iter->second << '\n';
	iter_prev = iter;
      }
      
      if (iter != iter_end && iter_prev != iter_end) {
	iter_end = phrase_pairs.upper_bound(iter_prev->first);
	for (/**/; iter != iter_end; ++ iter)
	  os << iter->second << '\n';
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
