
#include "cicada_filter_extract_impl.hpp"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>
#include <cfloat>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>
#include <boost/math/constants/constants.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>

#include <google/dense_hash_set>

typedef boost::filesystem::path path_type;

typedef RootCount  root_count_type;
typedef PhrasePair phrase_pair_type;

typedef google::dense_hash_set<root_count_type, boost::hash<root_count_type>, std::equal_to<root_count_type> > root_count_set_type;

path_type input_file = "-";
path_type output_file = "-";
path_type root_source_file;
path_type root_target_file;

double dirichlet_prior = 0.1;

bool mode_cicada = false;
bool feature_root = false;

int debug = 0;

template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_set_type& root_source,
	     const root_count_set_type& root_target);

struct ScorerCICADA;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (! boost::filesystem::exists(root_source_file))
      throw std::runtime_error("no root count file for source side");
    if (! boost::filesystem::exists(root_target_file))
      throw std::runtime_error("no root count file for target side");

    mode_cicada = true;
    dirichlet_prior = std::max(dirichlet_prior, 0.0);

    root_count_set_type root_source;
    root_count_set_type root_target;
    root_source.set_empty_key(root_count_type());
    root_target.set_empty_key(root_count_type());
    
    {
      root_count_type root_count;
      RootCountParser parser;
      std::string line;
      
      utils::compress_istream is_source(root_source_file);
      while (std::getline(is_source, line))
	if (parser(line, root_count))
	  root_source.insert(root_count);
      
      utils::compress_istream is_target(root_target_file);
      while (std::getline(is_target, line))
	if (parser(line, root_count))
	  root_target.insert(root_count);
    }

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    if (mode_cicada)
      process<ScorerCICADA>(is, os, root_source, root_target);
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
	     const root_count_set_type& root_source,
	     const root_count_set_type& root_target)
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
  ExtractPhraseSCFG phrase_extractor;

  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_set_type& root_count_source,
		  const root_count_set_type& root_count_target,
		  std::ostream& os)
  {
    if (phrase_pair.counts.size() != 1)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 1)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 1)
      throw std::runtime_error("target counts size do not match");
    
    const double& count = phrase_pair.counts.front();
    const double& count_source = phrase_pair.counts_source.front();
    const double& count_target = phrase_pair.counts_target.front();
    
    const double prob_source_target = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_source + count_source);
    const double prob_target_source = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_target + count_target);

    const std::pair<std::string, std::string> phrase_source = phrase_extractor(phrase_pair.source);
    const std::pair<std::string, std::string> phrase_target = phrase_extractor(phrase_pair.target);

    if (phrase_source.first != phrase_target.first)
      throw std::runtime_error("synchronous-CFG, but different lh?");

    if (feature_root) {
      root_count_set_type::const_iterator siter = root_count_source.find(phrase_source.first);
      root_count_set_type::const_iterator titer = root_count_target.find(phrase_target.first);
      
      if (siter == root_count_source.end())
	throw std::runtime_error("no root count for " + phrase_source.first);
      if (titer == root_count_target.end())
	throw std::runtime_error("no root count for " + phrase_target.first);
      
      if (siter->counts.size() != 1)
	throw std::runtime_error("invalid root count for source: " + phrase_source.first);
      if (titer->counts.size() != 1)
	throw std::runtime_error("invalid root count for target: " + phrase_target.first);
      
      const double prob_root_source = (dirichlet_prior + count_source) / (dirichlet_prior * siter->observed + siter->counts.front());
      const double prob_root_target = (dirichlet_prior + count_target) / (dirichlet_prior * titer->observed + titer->counts.front());
      
      os << phrase_source.first
	 << " ||| " << phrase_source.second
	 << " ||| " << phrase_target.second
	 << " |||"
	 << ' ' << std::log(prob_source_target) << ' ' << std::log(phrase_pair.lexicon_source_target)
	 << ' ' << std::log(prob_target_source) << ' ' << std::log(phrase_pair.lexicon_target_source)
	 << ' ' << std::log(prob_root_source)
	 << ' ' << std::log(prob_root_target)
	 << '\n';
    } else 
      os << phrase_source.first
	 << " ||| " << phrase_source.second
	 << " ||| " << phrase_target.second
	 << " |||"
	 << ' ' << std::log(prob_source_target) << ' ' << std::log(phrase_pair.lexicon_source_target)
	 << ' ' << std::log(prob_target_source) << ' ' << std::log(phrase_pair.lexicon_target_source)
	 << '\n';
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
    
    ("cicada",       po::bool_switch(&mode_cicada),   "output in cicada format")
    ("feature-root", po::bool_switch(&feature_root),  "output feature of p(lhs | root(lhs)) and p(rhs | root(rhs))")
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
