//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada_filter_extract_impl.hpp"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>
#include <cfloat>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>

typedef boost::filesystem::path path_type;

typedef RootCount  root_count_type;
typedef PhrasePair phrase_pair_type;

path_type input_file = "-";
path_type output_file = "-";
path_type root_source_file;
path_type root_target_file;

double dirichlet_prior = 0.1;

bool mode_cicada = false;
bool mode_moses = false;

bool mode_reordering = false;
bool mode_monotonicity = false;
bool mode_bidirectional = false;
bool mode_source_only = false;
bool mode_target_only = false;

int debug = 0;

template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_type& root_source,
	     const root_count_type& root_target);

struct ScorerCICADA;
struct ScorerCICADAReordering;
struct ScorerMOSES;
struct ScorerMOSESReordering;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (! boost::filesystem::exists(root_source_file))
      throw std::runtime_error("no root count file for source side");
    if (! boost::filesystem::exists(root_target_file))
      throw std::runtime_error("no root count file for target side");

    if (mode_cicada && mode_moses)
      throw std::runtime_error("specify either --{cicada,moses} (default to cicada)");
    if (int(mode_cicada) + mode_moses == 0)
      mode_cicada = true;

    if (mode_reordering)
      if (mode_source_only && mode_target_only)
	throw std::runtime_error("you cannot specify both of source-only and target-only");
    
    dirichlet_prior = std::max(dirichlet_prior, 0.0);
    
    root_count_type root_source;
    root_count_type root_target;
    
    {
      RootCountParser parser;
      std::string line;
      
      utils::compress_istream is_source(root_source_file);
      while (std::getline(is_source, line))
	if (parser(line, root_source)) break;
      
      utils::compress_istream is_target(root_target_file);
      while (std::getline(is_target, line))
	if (parser(line, root_target)) break;
    }

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    if (mode_cicada) {
      if (mode_reordering)
	process<ScorerCICADA>(is, os, root_source, root_target);
      else
	process<ScorerCICADAReordering>(is, os, root_source, root_target);
    } else if (mode_moses) {
      if (mode_reordering)
	process<ScorerMOSESReordering>(is, os, root_source, root_target);
      else
	process<ScorerMOSES>(is, os, root_source, root_target);
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_type& root_source,
	     const root_count_type& root_target)
{
  Scorer scorer;

  phrase_pair_type phrase_pair;
  PhrasePairParser parser;
  std::string line;
  
  while (std::getline(is, line)) {
    phrase_pair.clear();
    
    if (! parser(line, phrase_pair)) continue;
    
    scorer(phrase_pair, root_source, root_target, os);
  }
}

struct ScorerCICADA
{
  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  std::ostream& os)
  {
    if (phrase_pair.counts.size() != 5)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 5)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 5)
      throw std::runtime_error("target counts size do not match");
    
    const double& count = phrase_pair.counts.front();
    const double& count_source = phrase_pair.counts_source.front();
    const double& count_target = phrase_pair.counts_target.front();
    
    const double prob_source_target = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_source + count_source);
    const double prob_target_source = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_target + count_target);

    os << phrase_pair.source
       << " ||| " << phrase_pair.target
       << " |||"
       << ' ' << std::log(prob_source_target) << ' ' << std::log(phrase_pair.lexicon_source_target)
       << ' ' << std::log(prob_target_source) << ' ' << std::log(phrase_pair.lexicon_target_source)
       << '\n';
  }
};

struct ScorerCICADAReordering
{
  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  std::ostream& os)
  {
    if (phrase_pair.counts.size() != 5)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 5)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 5)
      throw std::runtime_error("target counts size do not match");
    
    if (mode_source_only)
      os << phrase_pair.source << " |||";
    else if (mode_target_only)
      os << phrase_pair.target << " |||";
    else
      os << phrase_pair.source << " ||| " << phrase_pair.target << " |||";

    const phrase_pair_type::counts_type& counts = (mode_source_only 
						   ? phrase_pair.counts_source
						   : (mode_target_only
						      ? phrase_pair.counts_target
						      : phrase_pair.counts));
    const double count = counts[0];
    const double count_prev_mono   = counts[1];
    const double count_prev_swap   = counts[2];
    const double count_prev_others = count - (count_prev_mono + count_prev_swap);
    const double count_next_mono   = counts[3];
    const double count_next_swap   = counts[4];
    const double count_next_others = count - (count_next_mono + count_next_swap);
    
    if (mode_monotonicity) {
      if (mode_bidirectional) {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	const double prob_next_mono   = (dirichlet_prior + count_next_mono) / (dirichlet_prior * 2 + count);
	const double prob_next_others = (dirichlet_prior + count_next_swap + count_next_others) / (dirichlet_prior * 2 + count);
	
	os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_others)
	   << ' ' << std::log(prob_next_mono) << ' ' << std::log(prob_next_others)
	   << '\n';
      } else {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	
	os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_others) << '\n';
      }
    } else {
      if (mode_bidirectional) {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	const double prob_next_mono   = (dirichlet_prior + count_next_mono)   / (dirichlet_prior * 3 + count);
	const double prob_next_swap   = (dirichlet_prior + count_next_swap)   / (dirichlet_prior * 3 + count);
	const double prob_next_others = (dirichlet_prior + count_next_others) / (dirichlet_prior * 3 + count);
	
	os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_swap) << ' ' << std::log(prob_prev_others)
	   << ' ' << std::log(prob_next_mono) << ' ' << std::log(prob_next_swap) << ' ' << std::log(prob_next_others)
	   << '\n';
      } else {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_swap) << ' ' << std::log(prob_prev_others) << '\n';
      }
    }
  }
};

struct ScorerMOSES
{
  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  std::ostream& os)
  {
    if (phrase_pair.counts.size() != 5)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 5)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 5)
      throw std::runtime_error("target counts size do not match");
    
    const double& count = phrase_pair.counts.front();
    const double& count_source = phrase_pair.counts_source.front();
    const double& count_target = phrase_pair.counts_target.front();
    
    const double prob_source_target = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_source + count_source);
    const double prob_target_source = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_target + count_target);

    os << phrase_pair.source
       << " ||| " << phrase_pair.target
       << " |||"
       << ' ' << prob_source_target << ' ' << phrase_pair.lexicon_source_target
       << ' ' << prob_target_source << ' ' << phrase_pair.lexicon_target_source
       << ' ' << boost::math::constants::e<double>()
       << '\n';
  }
};

struct ScorerMOSESReordering
{
  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  std::ostream& os)
  {
    if (phrase_pair.counts.size() != 5)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 5)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 5)
      throw std::runtime_error("target counts size do not match");
    
    if (mode_source_only)
      os << phrase_pair.source << " |||";
    else if (mode_target_only)
      os << phrase_pair.target << " |||";
    else
      os << phrase_pair.source << " ||| " << phrase_pair.target << " |||";

    const phrase_pair_type::counts_type& counts = (mode_source_only 
						   ? phrase_pair.counts_source
						   : (mode_target_only
						      ? phrase_pair.counts_target
						      : phrase_pair.counts));
    const double count = counts[0];
    const double count_prev_mono   = counts[1];
    const double count_prev_swap   = counts[2];
    const double count_prev_others = count - (count_prev_mono + count_prev_swap);
    const double count_next_mono   = counts[3];
    const double count_next_swap   = counts[4];
    const double count_next_others = count - (count_next_mono + count_next_swap);
    
    if (mode_monotonicity) {
      if (mode_bidirectional) {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	const double prob_next_mono   = (dirichlet_prior + count_next_mono) / (dirichlet_prior * 2 + count);
	const double prob_next_others = (dirichlet_prior + count_next_swap + count_next_others) / (dirichlet_prior * 2 + count);
	
	os << ' ' << prob_prev_mono << ' ' << prob_prev_others
	   << ' ' << prob_next_mono << ' ' << prob_next_others
	   << '\n';
      } else {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	
	os << ' ' << prob_prev_mono << ' ' << prob_prev_others << '\n';
      }
    } else {
      if (mode_bidirectional) {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	const double prob_next_mono   = (dirichlet_prior + count_next_mono)   / (dirichlet_prior * 3 + count);
	const double prob_next_swap   = (dirichlet_prior + count_next_swap)   / (dirichlet_prior * 3 + count);
	const double prob_next_others = (dirichlet_prior + count_next_others) / (dirichlet_prior * 3 + count);
	
	os << ' ' << prob_prev_mono << ' ' << prob_prev_swap << ' ' << prob_prev_others
	   << ' ' << prob_next_mono << ' ' << prob_next_swap << ' ' << prob_next_others
	   << '\n';
      } else {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	os << ' ' << prob_prev_mono << ' ' << prob_prev_swap << ' ' << prob_prev_others << '\n';
      }
    }
  }  
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("root-source", po::value<path_type>(&root_source_file), "root source file")
    ("root-target", po::value<path_type>(&root_target_file), "root target file")
    
    ("dirichlet-prior", po::value<double>(&dirichlet_prior)->default_value(dirichlet_prior), "dirichlet prior weight")
    
    ("cicada",   po::bool_switch(&mode_cicada),   "output in cicada format")
    ("moses",    po::bool_switch(&mode_moses),    "output in moses format")
    
    ("reordering",    po::bool_switch(&mode_reordering),    "reordering table")
    ("monotonicity",  po::bool_switch(&mode_monotonicity),  "monotonicity")
    ("bidirectional", po::bool_switch(&mode_bidirectional), "bidirectional")
    ("source-only",   po::bool_switch(&mode_source_only),   "source only")
    ("target-only",   po::bool_switch(&mode_target_only),   "target only")
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
