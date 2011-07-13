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

#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/spirit/include/karma.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>

typedef boost::filesystem::path path_type;

typedef RootCount  root_count_type;
typedef PhrasePair phrase_pair_type;

typedef LexiconModel lexicon_model_type;

path_type input_file = "-";
path_type output_file = "-";
path_type root_source_file;
path_type root_target_file;

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

double dirichlet_prior = 0.1;

bool mode_cicada = false;
bool mode_moses = false;

bool mode_reordering = false;
bool mode_monotonicity = false;
bool mode_bidirectional = false;
bool mode_source_only = false;
bool mode_target_only = false;

bool model1_mode = false;
bool noisy_or_mode = false;
bool insertion_deletion_mode = false;

double threshold_insertion = 0.01;
double threshold_deletion = 0.01;

int buffer_size = 1024 * 1024;
int debug = 0;

template <typename Scorer, typename Lexicon>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_type& root_source,
	     const root_count_type& root_target,
	     const Lexicon& lexicon);

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

    bool read_lexicon = false;
    
    if (model1_mode || noisy_or_mode || insertion_deletion_mode) {
      if (lexicon_source_target_file.empty() || ! boost::filesystem::exists(lexicon_source_target_file))
	throw std::runtime_error("no lexicon model for lex(target | source): " + lexicon_source_target_file.string());
      if (lexicon_target_source_file.empty() || ! boost::filesystem::exists(lexicon_target_source_file))
	throw std::runtime_error("no lexicon model for lex(source | target): " + lexicon_target_source_file.string());

      read_lexicon = true;
    }
    
    dirichlet_prior = std::max(dirichlet_prior, 0.0);
    
    root_count_type root_source;
    root_count_type root_target;
    
    lexicon_model_type lexicon_source_target;
    lexicon_model_type lexicon_target_source;

    if (read_lexicon) {
      lexicon_source_target.open(lexicon_source_target_file);
      lexicon_target_source.open(lexicon_target_source_file);
    }

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

    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, buffer_size);
    
    if (mode_cicada) {
      if (mode_reordering)
	process<ScorerCICADAReordering, LexiconPhrase>(is, os, root_source, root_target, LexiconPhrase(lexicon_source_target, lexicon_target_source));
      else
	process<ScorerCICADA, LexiconPhrase>(is, os, root_source, root_target, LexiconPhrase(lexicon_source_target, lexicon_target_source));
    } else if (mode_moses) {
      if (mode_reordering)
	process<ScorerMOSESReordering, LexiconPhrase>(is, os, root_source, root_target, LexiconPhrase(lexicon_source_target, lexicon_target_source));
      else
	process<ScorerMOSES, LexiconPhrase>(is, os, root_source, root_target, LexiconPhrase(lexicon_source_target, lexicon_target_source));
    }
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
	     const root_count_type& root_source,
	     const root_count_type& root_target,
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
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
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
    
    std::ostream_iterator<char> iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string << " |||"
			  << ' ' << double10 << ' ' << double10
			  << ' ' << double10 << ' ' << double10,
			  phrase_pair.source, phrase_pair.target,
			  std::log(prob_source_target), std::log(phrase_pair.lexicon_source_target),
			  std::log(prob_target_source), std::log(phrase_pair.lexicon_target_source)))
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

struct ScorerCICADAReordering
{
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
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    if (phrase_pair.counts.size() != 5)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 5)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 5)
      throw std::runtime_error("target counts size do not match");

    std::ostream_iterator<char> iter(os);

    if (! karma::generate(iter, standard::string << " ||| " << standard::string << " |||", phrase_pair.source, phrase_pair.target))
      throw std::runtime_error("failed generation");
    
    //os << phrase_pair.source << " ||| " << phrase_pair.target << " |||";
    
    const double& count = phrase_pair.counts.front();
    const double& count_source = phrase_pair.counts_source.front();
    const double& count_target = phrase_pair.counts_target.front();
    
    const double prob_source_target = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_source + count_source);
    const double prob_target_source = (dirichlet_prior + count) / (dirichlet_prior * phrase_pair.observed_target + count_target);

    if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10,
			  std::log(prob_source_target), std::log(phrase_pair.lexicon_source_target),
			  std::log(prob_target_source), std::log(phrase_pair.lexicon_target_source)))
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
    
    os << " |||";

    // we will dump reordering table as "attributes"
    {

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
	  
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_others),
				std::log(prob_next_mono), std::log(prob_next_others)))
	    throw std::runtime_error("failed generation");
#if 0
	  os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_others)
	     << ' ' << std::log(prob_next_mono) << ' ' << std::log(prob_next_others)
	     << '\n';
#endif
	} else {
	  const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	  const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_others)))
	    throw std::runtime_error("failed generation");
	  
	  //os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_others) << '\n';
	}
      } else {
	if (mode_bidirectional) {
	  const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	  const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	  const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	  const double prob_next_mono   = (dirichlet_prior + count_next_mono)   / (dirichlet_prior * 3 + count);
	  const double prob_next_swap   = (dirichlet_prior + count_next_swap)   / (dirichlet_prior * 3 + count);
	  const double prob_next_others = (dirichlet_prior + count_next_others) / (dirichlet_prior * 3 + count);
	
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_swap), std::log(prob_prev_others),
				std::log(prob_next_mono), std::log(prob_next_swap), std::log(prob_next_others)))
	    throw std::runtime_error("failed generation");
#if 0
	  os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_swap) << ' ' << std::log(prob_prev_others)
	     << ' ' << std::log(prob_next_mono) << ' ' << std::log(prob_next_swap) << ' ' << std::log(prob_next_others)
	     << '\n';
#endif
	} else {
	  const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	  const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	  const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_swap), std::log(prob_prev_others)))
	    throw std::runtime_error("failed generation");
	  
	  //os << ' ' << std::log(prob_prev_mono) << ' ' << std::log(prob_prev_swap) << ' ' << std::log(prob_prev_others) << '\n';
	}
      }
    }
  }
};

struct ScorerMOSES
{
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
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
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
    
    std::ostream_iterator<char> iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string << " |||"
			  << ' ' << double10 << ' ' << double10
			  << ' ' << double10 << ' ' << double10
			  << ' ' << double10
			  << '\n',
			  phrase_pair.source, phrase_pair.target,
			  prob_source_target, phrase_pair.lexicon_source_target,
			  prob_target_source, phrase_pair.lexicon_target_source,
			  boost::math::constants::e<double>()))
      throw std::runtime_error("failed generation");
  }
};

struct ScorerMOSESReordering
{
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
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    if (phrase_pair.counts.size() != 5)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() != 5)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() != 5)
      throw std::runtime_error("target counts size do not match");
    
    std::ostream_iterator<char> iter(os);

    if (mode_source_only) {
      if (! karma::generate(iter, standard::string << " |||", phrase_pair.source))
	throw std::runtime_error("failed generation");
    } else if (mode_target_only) {
      if (! karma::generate(iter, standard::string << " |||", phrase_pair.target))
	throw std::runtime_error("failed generation");
    } else {
      if (! karma::generate(iter, standard::string << " ||| " << standard::string << " |||", phrase_pair.source, phrase_pair.target))
	throw std::runtime_error("failed generation");
    }

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

	if (! karma::generate(iter,
			      ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
			      prob_prev_mono, prob_prev_others,
			      prob_next_mono, prob_next_others))
	  throw std::runtime_error("failed generation");
      } else {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	
	if (! karma::generate(iter,
			      ' ' << double10 << ' ' << double10 << '\n',
			      prob_prev_mono, prob_prev_others))
	  throw std::runtime_error("failed generation");
      }
    } else {
      if (mode_bidirectional) {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	const double prob_next_mono   = (dirichlet_prior + count_next_mono)   / (dirichlet_prior * 3 + count);
	const double prob_next_swap   = (dirichlet_prior + count_next_swap)   / (dirichlet_prior * 3 + count);
	const double prob_next_others = (dirichlet_prior + count_next_others) / (dirichlet_prior * 3 + count);
	
	if (! karma::generate(iter,
			      ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
			      prob_prev_mono, prob_prev_swap, prob_prev_others,
			      prob_next_mono, prob_next_swap, prob_next_others))
	  throw std::runtime_error("failed generation");
      } else {
	const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	
	if (! karma::generate(iter,
			      ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
			      prob_prev_mono, prob_prev_swap, prob_prev_others))
	  throw std::runtime_error("failed generation");
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

    ("lexicon-source-target",  po::value<path_type>(&lexicon_source_target_file),     "lexicon model for lex(target | source)")
    ("lexicon-target-source",  po::value<path_type>(&lexicon_target_source_file),     "lexicon model for lex(source | target)")
    
    ("dirichlet-prior", po::value<double>(&dirichlet_prior)->default_value(dirichlet_prior), "dirichlet prior weight")
    
    ("cicada",   po::bool_switch(&mode_cicada),   "output in cicada format")
    ("moses",    po::bool_switch(&mode_moses),    "output in moses format")
    
    ("reordering",    po::bool_switch(&mode_reordering),    "reordering table")
    ("monotonicity",  po::bool_switch(&mode_monotonicity),  "monotonicity")
    ("bidirectional", po::bool_switch(&mode_bidirectional), "bidirectional")
    ("source-only",   po::bool_switch(&mode_source_only),   "source only")
    ("target-only",   po::bool_switch(&mode_target_only),   "target only")
    
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
