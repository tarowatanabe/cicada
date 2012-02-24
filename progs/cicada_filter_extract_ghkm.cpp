//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

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
#include <boost/spirit/include/karma.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/dense_hash_set.hpp>

typedef boost::filesystem::path path_type;

typedef RootCount  root_count_type;
typedef PhrasePair phrase_pair_type;

typedef utils::dense_hash_set<root_count_type, boost::hash<root_count_type>, std::equal_to<root_count_type> >::type root_count_set_type;

typedef LexiconModel lexicon_model_type;

path_type input_file = "-";
path_type output_file = "-";
path_type root_source_file;
path_type root_target_file;

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

bool model1_mode = false;
bool noisy_or_mode = false;
bool insertion_deletion_mode = false;

double threshold_insertion = 0.01;
double threshold_deletion = 0.01;

double dirichlet_prior = 0.1;

int buffer_size = 1024 * 1024;
int debug = 0;

template <typename Scorer, typename Lexicon>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_set_type& root_source,
	     const root_count_set_type& root_target,
	     const Lexicon& lexicon);

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

    bool read_lexicon = false;
    
    if (model1_mode || noisy_or_mode || insertion_deletion_mode) {
      if (lexicon_source_target_file.empty() || ! boost::filesystem::exists(lexicon_source_target_file))
	throw std::runtime_error("no lexicon model for lex(target | source): " + lexicon_source_target_file.string());
      if (lexicon_target_source_file.empty() || ! boost::filesystem::exists(lexicon_target_source_file))
	throw std::runtime_error("no lexicon model for lex(source | target): " + lexicon_target_source_file.string());

      read_lexicon = true;
    }
    
    dirichlet_prior = std::max(dirichlet_prior, 0.0);

    root_count_set_type root_source;
    root_count_set_type root_target;
    root_source.set_empty_key(root_count_type());
    root_target.set_empty_key(root_count_type());
    
    lexicon_model_type lexicon_source_target;
    lexicon_model_type lexicon_target_source;

    if (read_lexicon) {
      lexicon_source_target.open(lexicon_source_target_file);
      lexicon_target_source.open(lexicon_target_source_file);
    }

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

    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, buffer_size);
    
    process<ScorerCICADA, LexiconGHKM>(is, os, root_source, root_target, LexiconGHKM(lexicon_source_target, lexicon_target_source));
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Scorer, typename Lexicon>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_set_type& root_source,
	     const root_count_set_type& root_target,
	     const Lexicon& lexicon)
{
  Scorer scorer;

  phrase_pair_type phrase_pair;
  PhrasePairParser parser;
  std::string line;
  
  while (std::getline(is, line)) {
    phrase_pair.clear();
    
    if (! parser(line, phrase_pair)) continue;
    
    scorer(phrase_pair, root_source, root_target, lexicon, os);
  }
}

struct ScorerCICADA
{
  ExtractRootGHKM root_extractor;

  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 10;
    }
  };
  
  boost::spirit::karma::real_generator<double, real_precision> double10;
  
  template <typename Lexicon>
  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_set_type& root_count_source,
		  const root_count_set_type& root_count_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

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
    
    const std::string root_source = root_extractor(phrase_pair.source);
    const std::string root_target = root_extractor(phrase_pair.target);
    
    root_count_set_type::const_iterator siter = root_count_source.find(root_source);
    root_count_set_type::const_iterator titer = root_count_target.find(root_target);
    
    if (siter == root_count_source.end())
      throw std::runtime_error("no root count for source: " + root_source);
    if (titer == root_count_target.end())
      throw std::runtime_error("no root count for target: " + root_target);
    
    if (siter->counts.size() != 1)
      throw std::runtime_error("invalid root count for source: " + root_source);
    if (titer->counts.size() != 1)
      throw std::runtime_error("invalid root count for target: " + root_target);
    
    const double prob_root_source = (dirichlet_prior + count_source) / (dirichlet_prior * siter->observed + siter->counts.front());
    const double prob_root_target = (dirichlet_prior + count_target) / (dirichlet_prior * titer->observed + titer->counts.front());

    std::ostream_iterator<char> iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string << " |||"
			  << ' ' << double10 << ' ' << double10
			  << ' ' << double10 << ' ' << double10
			  << ' ' << double10
			  << ' ' << double10,
			  phrase_pair.source, phrase_pair.target,
			  std::log(prob_source_target), std::log(phrase_pair.lexicon_source_target),
			  std::log(prob_target_source), std::log(phrase_pair.lexicon_target_source),
			  std::log(prob_root_source),
			  std::log(prob_root_target)))
      throw std::runtime_error("failed generation");
    
    if (model1_mode || noisy_or_mode || insertion_deletion_mode) {
      const_cast<Lexicon&>(lexicon).assign_source(phrase_pair.source);
      const_cast<Lexicon&>(lexicon).assign_target(phrase_pair.target);
      
      if (model1_mode) {
	const std::pair<double, double> scores = lexicon.model1();
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }
      
      if (noisy_or_mode) {
	const std::pair<double, double> scores = lexicon.noisy_or();
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }
      
      if (insertion_deletion_mode) {
	const std::pair<double, double> scores = lexicon.insertion_deletion(threshold_insertion, threshold_deletion);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }
    }
    
    os << '\n';
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

    ("lexicon-source-target",  po::value<path_type>(&lexicon_source_target_file),     "lexicon model for lex(target | source)")
    ("lexicon-target-source",  po::value<path_type>(&lexicon_target_source_file),     "lexicon model for lex(source | target)")
    
    ("dirichlet-prior", po::value<double>(&dirichlet_prior)->default_value(dirichlet_prior), "dirichlet prior weight")

    ("model1",             po::bool_switch(&model1_mode),             "Model1 feature (requires lexicon models)")
    ("noisy-or",           po::bool_switch(&noisy_or_mode),           "noisy-or feature (requires lexicon models)")
    ("insertion-deletion", po::bool_switch(&insertion_deletion_mode), "insertion/deletion feature (requires lexicon models)")
    
    ("threshold-insertion", po::value<double>(&threshold_insertion)->default_value(threshold_insertion), "threshold for insertion")
    ("threshold-deletion",  po::value<double>(&threshold_deletion)->default_value(threshold_deletion),   "threshold for deletion")
    
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
