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
path_type root_joint_file;
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

bool feature_root_mode = false;
bool feature_type_mode = false;
bool feature_singleton_mode = false;
bool feature_cross_mode = false;
bool feature_unaligned_mode = false;

bool feature_lexicon_mode = false;
bool feature_model1_mode = false;
bool feature_noisy_or_mode = false;
bool feature_insertion_deletion_mode = false;

double threshold_insertion = 0.01;
double threshold_deletion = 0.01;

int buffer_size = 1024 * 1024;
int debug = 0;

template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const root_count_type& root_joint,
	     const root_count_type& root_source,
	     const root_count_type& root_target,
	     const Lexicon& lexicon);

struct ScorerCICADA;
struct ScorerMOSES;
struct ScorerMOSESReordering;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (! boost::filesystem::exists(root_joint_file))
      throw std::runtime_error("no root count file");
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
    
    if (feature_lexicon_mode || feature_model1_mode || feature_noisy_or_mode || feature_insertion_deletion_mode) {
      if (lexicon_source_target_file.empty() || ! boost::filesystem::exists(lexicon_source_target_file))
	throw std::runtime_error("no lexicon model for lex(target | source): " + lexicon_source_target_file.string());
      if (lexicon_target_source_file.empty() || ! boost::filesystem::exists(lexicon_target_source_file))
	throw std::runtime_error("no lexicon model for lex(source | target): " + lexicon_target_source_file.string());
      
      read_lexicon = true;
    }
    
    dirichlet_prior = std::max(dirichlet_prior, 0.0);
    
    root_count_type root_joint;
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

      utils::compress_istream is_joint(root_joint_file);
      while (std::getline(is_joint, line))
        if (parser(line, root_joint)) break;
      
      utils::compress_istream is_source(root_source_file);
      while (std::getline(is_source, line))
        if (parser(line, root_source)) break;
      
      utils::compress_istream is_target(root_target_file);
      while (std::getline(is_target, line))
        if (parser(line, root_target)) break;

      if (root_joint.counts.empty())
	throw std::runtime_error("invalid counts");
      if (root_source.counts.empty())
	throw std::runtime_error("invalid source counts");
      if (root_target.counts.empty())
	throw std::runtime_error("invalid target counts");
    }

    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, buffer_size);
    
    if (mode_cicada)
      process<ScorerCICADA>(is, os, root_joint, root_source, root_target, Lexicon(lexicon_source_target, lexicon_target_source));
    else if (mode_moses) {
      if (mode_reordering)
	process<ScorerMOSESReordering>(is, os, root_joint, root_source, root_target, Lexicon(lexicon_source_target, lexicon_target_source));
      else
	process<ScorerMOSES>(is, os, root_joint, root_source, root_target, Lexicon(lexicon_source_target, lexicon_target_source));
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
	     const root_count_type& root_joint,
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
    
    scorer(phrase_pair, root_joint, root_source, root_target, lexicon, os);
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
  
  ExtractPhrase    extract_phrase;
  ExtractAlignment extract_alignment;
  Unaligned        unaligned;
  Cross            cross;
  
  typedef ExtractPhrase::sentence_type         sentence_type;
  typedef ExtractAlignment::alignment_type     alignment_type;
  typedef ExtractAlignment::alignment_set_type alignment_set_type;

  sentence_type      source;
  sentence_type      target;
  alignment_set_type alignments;
  
  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_joint,
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
			  << ' ' << double10 << ' ' << double10,
			  phrase_pair.source, phrase_pair.target,
			  std::log(prob_source_target), 
			  std::log(prob_target_source)))
      throw std::runtime_error("failed generation");
    
    if (feature_root_mode) {
      const double logprob_root = (std::log(dirichlet_prior + count)
				   - std::log(dirichlet_prior * root_joint.observed + root_joint.counts.front()));
      const double logprob_root_source = (std::log(dirichlet_prior + count_source)
					  - std::log(dirichlet_prior * root_source.observed + root_source.counts.front()));
      const double logprob_root_target = (std::log(dirichlet_prior + count_target)
					  - std::log(dirichlet_prior * root_target.observed + root_target.counts.front()));
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10,
			    logprob_root, logprob_root_source, logprob_root_target))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_cross_mode || feature_unaligned_mode || feature_lexicon_mode || feature_model1_mode || feature_noisy_or_mode || feature_insertion_deletion_mode) {
      extract_phrase(phrase_pair.source, source);
      extract_phrase(phrase_pair.target, target);
    }
    
    if (feature_cross_mode || feature_lexicon_mode || feature_unaligned_mode)
      extract_alignment(phrase_pair.alignments, alignments);
    
    if (feature_lexicon_mode || feature_model1_mode || feature_noisy_or_mode || feature_insertion_deletion_mode) {
      if (feature_lexicon_mode) {
	const std::pair<double, double> scores = lexicon.lexicon(source, target, alignments);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }

      if (feature_model1_mode) {
	const std::pair<double, double> scores = lexicon.model1(source, target);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }
      
      if (feature_noisy_or_mode) {
	const std::pair<double, double> scores = lexicon.noisy_or(source, target);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }
      
      if (feature_insertion_deletion_mode) {
	const std::pair<double, double> scores = lexicon.insertion_deletion(source, target, threshold_insertion, threshold_deletion);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, scores.first, scores.second))
	  throw std::runtime_error("failed generation");
      }
    }

    if (feature_unaligned_mode) {
      const std::pair<size_t, size_t> scores = unaligned(source, target, alignments);
      
      if (! karma::generate(iter, ' ' << karma::uint_ << ' ' << karma::uint_, scores.first, scores.second))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_type_mode) {
      const double prob_type_source_target = 1.0 / phrase_pair.observed_source;
      const double prob_type_target_source = 1.0 / phrase_pair.observed_target;
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10, std::log(prob_type_source_target), std::log(prob_type_target_source)))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_singleton_mode) {
      const int singleton_source = phrase_pair.observed_source == 1;
      const int singleton_target = phrase_pair.observed_target == 1;
      const int singleton        = singleton_source && singleton_target;
      
      if (! karma::generate(iter, ' ' << karma::int_ << ' ' << karma::int_ << ' ' << karma::int_,
			    singleton, singleton_source, singleton_target))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_cross_mode)
      if (! karma::generate(iter, ' ' << karma::int_, cross(source, target, alignments)))
	throw std::runtime_error("failed generation");
    
    if (! mode_reordering) 
      os << '\n';
    else {
      // we will dump reordering table as "attributes"
      
      os << " |||";
      
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
	} else {
	  const double prob_prev_mono   = (dirichlet_prior + count_prev_mono) / (dirichlet_prior * 2 + count);
	  const double prob_prev_others = (dirichlet_prior + count_prev_swap + count_prev_others) / (dirichlet_prior * 2 + count);
	
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_others)))
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
	
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_swap), std::log(prob_prev_others),
				std::log(prob_next_mono), std::log(prob_next_swap), std::log(prob_next_others)))
	    throw std::runtime_error("failed generation");
	} else {
	  const double prob_prev_mono   = (dirichlet_prior + count_prev_mono)   / (dirichlet_prior * 3 + count);
	  const double prob_prev_swap   = (dirichlet_prior + count_prev_swap)   / (dirichlet_prior * 3 + count);
	  const double prob_prev_others = (dirichlet_prior + count_prev_others) / (dirichlet_prior * 3 + count);
	  
	  if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10 << '\n',
				std::log(prob_prev_mono), std::log(prob_prev_swap), std::log(prob_prev_others)))
	    throw std::runtime_error("failed generation");
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
  
  ExtractPhrase    extract_phrase;
  ExtractAlignment extract_alignment;
  Unaligned        unaligned;
  Cross            cross;
  
  typedef ExtractPhrase::sentence_type         sentence_type;
  typedef ExtractAlignment::alignment_type     alignment_type;
  typedef ExtractAlignment::alignment_set_type alignment_set_type;
  
  sentence_type      source;
  sentence_type      target;
  alignment_set_type alignments;

  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_joint,
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
			  << ' ' << double10,
			  phrase_pair.source, phrase_pair.target,
			  prob_source_target, prob_target_source, 
			  boost::math::constants::e<double>()))
      throw std::runtime_error("failed generation");
    
    
    if (feature_root_mode) {
      const double prob_root = ((dirichlet_prior + count)
				/ (dirichlet_prior * root_joint.observed + root_joint.counts.front()));
      const double prob_root_source = ((dirichlet_prior + count_source)
				       / (dirichlet_prior * root_source.observed + root_source.counts.front()));
      const double prob_root_target = ((dirichlet_prior + count_target)
				       / (dirichlet_prior * root_target.observed + root_target.counts.front()));
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10,
			    prob_root, prob_root_source, prob_root_target))
	throw std::runtime_error("failed generation");
    }
        
    if (feature_cross_mode || feature_unaligned_mode || feature_lexicon_mode || feature_model1_mode || feature_noisy_or_mode || feature_insertion_deletion_mode) {
      extract_phrase(phrase_pair.source, source);
      extract_phrase(phrase_pair.target, target);
    }
    
    if (feature_cross_mode || feature_lexicon_mode || feature_unaligned_mode)
      extract_alignment(phrase_pair.alignments, alignments);
    
    if (feature_lexicon_mode || feature_model1_mode || feature_noisy_or_mode || feature_insertion_deletion_mode) {
      if (feature_lexicon_mode) {
	const std::pair<double, double> scores = lexicon.lexicon(source, target, alignments);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, utils::mathop::exp(scores.first), utils::mathop::exp(scores.second)))
	  throw std::runtime_error("failed generation");
      }

      if (feature_model1_mode) {
	const std::pair<double, double> scores = lexicon.model1(source, target);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, utils::mathop::exp(scores.first), utils::mathop::exp(scores.second)))
	  throw std::runtime_error("failed generation");
      }
      
      if (feature_noisy_or_mode) {
	const std::pair<double, double> scores = lexicon.noisy_or(source, target);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, utils::mathop::exp(scores.first), utils::mathop::exp(scores.second)))
	  throw std::runtime_error("failed generation");
      }
      
      if (feature_insertion_deletion_mode) {
	const std::pair<double, double> scores = lexicon.insertion_deletion(source, target, threshold_insertion, threshold_deletion);
	
	if (! karma::generate(iter, ' ' << double10 << ' ' << double10, utils::mathop::exp(scores.first), utils::mathop::exp(scores.second)))
	  throw std::runtime_error("failed generation");
      }
    }

    if (feature_unaligned_mode) {
      const std::pair<size_t, size_t> scores = unaligned(source, target, alignments);
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10, utils::mathop::exp(scores.first), utils::mathop::exp(scores.second)))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_type_mode) {
      const double prob_type_source_target = 1.0 / phrase_pair.observed_source;
      const double prob_type_target_source = 1.0 / phrase_pair.observed_target;
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10, prob_type_source_target, prob_type_target_source))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_singleton_mode) {
      const int singleton_source = phrase_pair.observed_source == 1;
      const int singleton_target = phrase_pair.observed_target == 1;
      const int singleton        = singleton_source && singleton_target;
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10,
			    utils::mathop::exp(singleton), utils::mathop::exp(singleton_source), utils::mathop::exp(singleton_target)))
	throw std::runtime_error("failed generation");
    }
    
    if (feature_cross_mode)
      if (! karma::generate(iter, ' ' << double10, utils::mathop::exp(cross(source, target, alignments))))
	throw std::runtime_error("failed generation");
    
    os << '\n';
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

  void operator()(const phrase_pair_type& phrase_pair,
		  const root_count_type& root_joint,
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
    
    ("root-joint",  po::value<path_type>(&root_joint_file),  "root count file")
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

    ("feature-root",       po::bool_switch(&feature_root_mode),       "feature by generative probability")
    ("feature-type",       po::bool_switch(&feature_type_mode),       "feature by obesrved types")
    ("feature-singleton",  po::bool_switch(&feature_singleton_mode),  "singleton features")
    ("feature-cross",      po::bool_switch(&feature_cross_mode),      "crossing features")
    ("feature-unaligned",  po::bool_switch(&feature_unaligned_mode),  "unaligned features")
    
    ("feature-lexicon",            po::bool_switch(&feature_lexicon_mode),            "lexical weight feature (requires lexicon models)")
    ("feature-model1",             po::bool_switch(&feature_model1_mode),             "Model1 feature (requires lexicon models)")
    ("feature-noisy-or",           po::bool_switch(&feature_noisy_or_mode),           "noisy-or feature (requires lexicon models)")
    ("feature-insertion-deletion", po::bool_switch(&feature_insertion_deletion_mode), "insertion/deletion feature (requires lexicon models)")
        
    
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
