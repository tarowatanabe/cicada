//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada_extract_impl.hpp"
#include "cicada_filter_extract_impl.hpp"

#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>
#include <cfloat>
#include <cmath>

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/spirit/include/karma.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/unordered_set.hpp>
#include "utils/getline.hpp"

typedef boost::filesystem::path path_type;

typedef RootCount  root_count_type;
typedef PhrasePair phrase_pair_type;

typedef utils::unordered_set<root_count_type, boost::hash<root_count_type>, std::equal_to<root_count_type>,
			     std::allocator<root_count_type> >::type root_count_set_type;

typedef LexiconModel lexicon_model_type;

typedef Statistic statistic_type;

path_type input_file = "-";
path_type output_file = "-";
path_type root_joint_file;
path_type root_source_file;
path_type root_target_file;

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

path_type statistic_file;

double dirichlet_prior = 0.1;

bool mode_cicada = false;
bool mode_moses = false;

bool feature_root_mode = false;
bool feature_fisher_mode = false;
bool feature_type_mode = false;
bool feature_singleton_mode = false;
bool feature_cross_mode = false;
bool feature_unaligned_mode = false;

bool feature_lexicon_mode = false;
bool feature_model1_mode = false;
bool feature_noisy_or_mode = false;
bool feature_insertion_deletion_mode = false;

double threshold_insertion = 0.5;
double threshold_deletion = 0.5;

int buffer_size = 1024 * 1024;
int debug = 0;

template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const statistic_type& statistic,
	     const root_count_set_type& root_joint,
	     const root_count_set_type& root_source,
	     const root_count_set_type& root_target,
	     const Lexicon& lexicon);
template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const statistic_type& statistic,
	     const root_count_type& root_joint,
	     const root_count_type& root_source,
	     const root_count_type& root_target,
	     const Lexicon& lexicon);

struct ScorerCICADA;

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

    bool read_lexicon = false;
    
    if (feature_lexicon_mode || feature_model1_mode || feature_noisy_or_mode || feature_insertion_deletion_mode) {
      if (lexicon_source_target_file.empty() || ! boost::filesystem::exists(lexicon_source_target_file))
	throw std::runtime_error("no lexicon model for lex(target | source): " + lexicon_source_target_file.string());
      if (lexicon_target_source_file.empty() || ! boost::filesystem::exists(lexicon_target_source_file))
	throw std::runtime_error("no lexicon model for lex(source | target): " + lexicon_target_source_file.string());

      read_lexicon = true;
    }
    
    if (feature_fisher_mode) {
      if (statistic_file.empty() || ! boost::filesystem::exists(statistic_file))
	throw std::runtime_error("no statistic for sigtest?");
    }
    
    statistic_type statistic;
    if (! statistic_file.empty()) {
      utils::compress_istream is(statistic_file);
      is >> statistic;
    }
    
    dirichlet_prior = std::max(dirichlet_prior, 0.0);

    root_count_set_type root_joint;
    root_count_set_type root_source;
    root_count_set_type root_target;

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

      utils::compress_istream is_joint(root_joint_file);
      while (utils::getline(is_joint, line))
	if (parser(line, root_count))
	  root_joint.insert(root_count);
      
      utils::compress_istream is_source(root_source_file);
      while (utils::getline(is_source, line))
	if (parser(line, root_count))
	  root_source.insert(root_count);
      
      utils::compress_istream is_target(root_target_file);
      while (utils::getline(is_target, line))
	if (parser(line, root_count))
	  root_target.insert(root_count);
    }

    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, buffer_size);

    if (root_joint.size() == 1 && root_source.size() == 1 && root_target.size() == 1)
      process<ScorerCICADA>(is, os, statistic, *root_joint.begin(), *root_source.begin(), *root_target.begin(), Lexicon(lexicon_source_target, lexicon_target_source));
    else
      process<ScorerCICADA>(is, os, statistic, root_joint, root_source, root_target, Lexicon(lexicon_source_target, lexicon_target_source));
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
	     const statistic_type& statistic,
	     const root_count_set_type& root_joint,
	     const root_count_set_type& root_source,
	     const root_count_set_type& root_target,
	     const Lexicon& lexicon)
{
  Scorer scorer;
  
  if (debug)
    std::cerr << scorer;
  
  phrase_pair_type phrase_pair;
  PhrasePairParser parser;
  std::string line;
  
  while (utils::getline(is, line)) {
    phrase_pair.clear();
    
    if (! parser(line, phrase_pair)) continue;
    
    scorer(phrase_pair, statistic, root_joint, root_source, root_target, lexicon, os);
  }
}

template <typename Scorer>
void process(std::istream& is,
	     std::ostream& os,
	     const statistic_type& statistic,
	     const root_count_type& root_joint,
	     const root_count_type& root_source,
	     const root_count_type& root_target,
	     const Lexicon& lexicon)
{
  Scorer scorer;

  if (debug)
    std::cerr << scorer;

  phrase_pair_type phrase_pair;
  PhrasePairParser parser;
  std::string line;
  
  while (utils::getline(is, line)) {
    phrase_pair.clear();
    
    if (! parser(line, phrase_pair)) continue;
    
    scorer(phrase_pair, statistic, root_joint, root_source, root_target, lexicon, os);
  }
}

struct ScorerCICADA
{
  ExtractPhraseSCFG phrase_extractor;
  
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 10;
    }
  };
  
  boost::spirit::karma::real_generator<double, real_precision> double10;

  ExtractSCFG      extract_phrase;
  ExtractAlignment extract_alignment;
  Unaligned        unaligned;
  Cross            cross;
  Fisher           fisher;
  
  typedef ExtractSCFG::sentence_type           sentence_type;
  typedef ExtractAlignment::alignment_type     alignment_type;
  typedef ExtractAlignment::alignment_set_type alignment_set_type;

  sentence_type      source;
  sentence_type      target;
  alignment_set_type alignments;
  
  friend
  std::ostream& operator<<(std::ostream& os, const ScorerCICADA& x) 
  {
    os << "log p(target | source)" << '\n';
    os << "log p(source | target)" << '\n';

    if (feature_root_mode) {
      os << "log p(source, target | root(source), root(target))" << '\n';
      os << "log p(source | root(source))" << '\n';
      os << "log p(target | root(target))" << '\n';
    }
    
    if (feature_fisher_mode)
      os << "fisher's exact test score" << '\n';
    
    if (feature_lexicon_mode) {
      os << "log lex(target | source)" << '\n';
      os << "log lex(source | target)" << '\n';
    }

    if (feature_model1_mode) {
      os << "log model1(target | source)" << '\n';
      os << "log model1(source | target)" << '\n';
    }

    if (feature_noisy_or_mode) {
      os << "log noisyor(target | source)" << '\n';
      os << "log noisyor(source | target)" << '\n';
    }

    if (feature_insertion_deletion_mode) {
      os << "# of inserted words" << '\n';
      os << "# of deleted words" << '\n';
    }
    
    if (feature_unaligned_mode) {
      os << "# of aligned source words" << '\n';
      os << "# of aligned target words" << '\n';
      os << "# of unaligned source words" << '\n';
      os << "# of unaligned target words" << '\n';
    }

    if (feature_type_mode) {
      os << "log 1 / type(source)" << '\n';
      os << "log 1 / type(target)" << '\n';
    }
    
    if (feature_singleton_mode) {
      os << "singleton for source" << '\n';
      os << "singleton for target" << '\n';
      os << "singleton for source/target" << '\n';
    }

    if (feature_cross_mode) {
      os << "# of crossed non-terminal alignment" << '\n';
      os << "# of crossed word alignment" << '\n';
    }
    
    return os;
  }

  void operator()(const phrase_pair_type& phrase_pair,
		  const statistic_type& statistic,
		  const root_count_set_type& root_count_joint,
		  const root_count_set_type& root_count_source,
		  const root_count_set_type& root_count_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    if (phrase_pair.counts.size() < 1)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() < 1)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() < 1)
      throw std::runtime_error("target counts size do not match");
    
    const double& count = phrase_pair.counts.front();
    const double& count_source = phrase_pair.counts_source.front();
    const double& count_target = phrase_pair.counts_target.front();
    
    const double logprob_source_target = (std::log(dirichlet_prior + count)
					  - std::log(dirichlet_prior * phrase_pair.observed_source + count_source));
    const double logprob_target_source = (std::log(dirichlet_prior + count)
					  - std::log(dirichlet_prior * phrase_pair.observed_target + count_target));

    const std::pair<std::string, std::string> phrase_source = phrase_extractor(phrase_pair.source);
    const std::pair<std::string, std::string> phrase_target = phrase_extractor(phrase_pair.target);

    if (phrase_source.first != phrase_target.first)
      throw std::runtime_error("synchronous-CFG, but different lhs?" + phrase_source.first + ' ' + phrase_target.first
			       + " original: " + phrase_pair.source + ' ' + phrase_pair.target);
    
    std::ostream_iterator<char> iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string << " ||| " << standard::string << " |||"
			  << ' ' << double10 << ' ' << double10,
			  phrase_source.first, phrase_source.second, phrase_target.second,
			  logprob_source_target,
			  logprob_target_source))
      throw std::runtime_error("failed generation");
    
    if (feature_root_mode) {
      root_count_set_type::const_iterator jiter = root_count_joint.find(phrase_source.first + phrase_target.first);
      root_count_set_type::const_iterator siter = root_count_source.find(phrase_source.first);
      root_count_set_type::const_iterator titer = root_count_target.find(phrase_target.first);
      
      if (jiter == root_count_joint.end())
	throw std::runtime_error("no root count: " + phrase_source.first + phrase_target.first);
      if (siter == root_count_source.end())
	throw std::runtime_error("no root count for source: " + phrase_source.first);
      if (titer == root_count_target.end())
	throw std::runtime_error("no root count for target: " + phrase_target.first);
      
      if (jiter->counts.size() < 1)
	throw std::runtime_error("invalid root count: " + phrase_source.first + phrase_target.first);
      if (siter->counts.size() < 1)
	throw std::runtime_error("invalid root count for source: " + phrase_source.first);
      if (titer->counts.size() < 1)
	throw std::runtime_error("invalid root count for target: " + phrase_target.first);
      
      const double logprob_root = (std::log(dirichlet_prior + count)
				   - std::log(dirichlet_prior * jiter->observed + jiter->counts.front()));
      const double logprob_root_source = (std::log(dirichlet_prior + count_source)
					  - std::log(dirichlet_prior * siter->observed + siter->counts.front()));
      const double logprob_root_target = (std::log(dirichlet_prior + count_target)
					  - std::log(dirichlet_prior * titer->observed + titer->counts.front()));
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10,
			    logprob_root,
			    logprob_root_source,
			    logprob_root_target))
	throw std::runtime_error("failed generation");
    }
    
    operator()(phrase_pair, statistic, lexicon, os);

    os << '\n';
  }

  void operator()(const phrase_pair_type& phrase_pair,
		  const statistic_type& statistic,
		  const root_count_type& root_joint,
		  const root_count_type& root_source,
		  const root_count_type& root_target,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    if (phrase_pair.counts.size() < 1)
      throw std::runtime_error("counts size do not match");
    if (phrase_pair.counts_source.size() < 1)
      throw std::runtime_error("source counts size do not match");
    if (phrase_pair.counts_target.size() < 1)
      throw std::runtime_error("target counts size do not match");
    
    const double& count = phrase_pair.counts.front();
    const double& count_source = phrase_pair.counts_source.front();
    const double& count_target = phrase_pair.counts_target.front();
    
    const double logprob_source_target = (std::log(dirichlet_prior + count)
					  - std::log(dirichlet_prior * phrase_pair.observed_source + count_source));
    const double logprob_target_source = (std::log(dirichlet_prior + count)
					  - std::log(dirichlet_prior * phrase_pair.observed_target + count_target));

    const std::pair<std::string, std::string> phrase_source = phrase_extractor(phrase_pair.source);
    const std::pair<std::string, std::string> phrase_target = phrase_extractor(phrase_pair.target);

    if (phrase_source.first != phrase_target.first)
      throw std::runtime_error("synchronous-CFG, but different lhs?" + phrase_source.first + ' ' + phrase_target.first
			       + " original: " + phrase_pair.source + ' ' + phrase_pair.target);
    
    std::ostream_iterator<char> iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string << " ||| " << standard::string << " |||"
			  << ' ' << double10 << ' ' << double10,
			  phrase_source.first, phrase_source.second, phrase_target.second,
			  logprob_source_target,
			  logprob_target_source))
      throw std::runtime_error("failed generation");
    
    if (feature_root_mode) {
      if (root_source.label != phrase_source.first)
	throw std::runtime_error("invalid root count for source: " + phrase_source.first);
      if (root_target.label != phrase_target.first)
	throw std::runtime_error("invalid root count for target: " + phrase_target.first);

      if (root_joint.counts.size() < 1)
	throw std::runtime_error("invalid root count: " + phrase_source.first + phrase_target.first);
      if (root_source.counts.size() < 1)
	throw std::runtime_error("invalid root count for source: " + phrase_source.first);
      if (root_target.counts.size() < 1)
	throw std::runtime_error("invalid root count for target: " + phrase_target.first);
      
      const double logprob_root = (std::log(dirichlet_prior + count)
				   - std::log(dirichlet_prior * root_joint.observed + root_joint.counts.front()));
      const double logprob_root_source = (std::log(dirichlet_prior + count_source)
					  - std::log(dirichlet_prior * root_source.observed + root_source.counts.front()));
      const double logprob_root_target = (std::log(dirichlet_prior + count_target)
					  - std::log(dirichlet_prior * root_target.observed + root_target.counts.front()));
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10 << ' ' << double10,
			    logprob_root,
			    logprob_root_source,
			    logprob_root_target))
	throw std::runtime_error("failed generation");
    }
        
    operator()(phrase_pair, statistic, lexicon, os);
    
    os << '\n';
  }

  void operator()(const phrase_pair_type& phrase_pair,
		  const statistic_type& statistic,
		  const Lexicon& lexicon,
		  std::ostream& os)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    std::ostream_iterator<char> iter(os);

    if (feature_fisher_mode) {
      if (phrase_pair.counts.size() < 1 + 3)
	throw std::runtime_error("invalid counts for Fisher's exact test");
      if (phrase_pair.counts_source.size() < 1 + 3)
	throw std::runtime_error("invalid source counts for Fisher's exact test");
      if (phrase_pair.counts_target.size() < 1 + 3)
	throw std::runtime_error("invalid target counts for Fisher's exact test");

      const Fisher::count_type cfe = phrase_pair.counts[phrase_pair.counts.size() - 3];
      const Fisher::count_type cf  = phrase_pair.counts_source[phrase_pair.counts_source.size() - 2];
      const Fisher::count_type ce  = phrase_pair.counts_target[phrase_pair.counts_target.size() - 1];
      
      const double score = fisher(cfe, cf, ce, statistic.bitext);
      
      if (! karma::generate(iter, ' ' << double10, score))
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
      const Unaligned::scores_type scores = unaligned(source, target, alignments);
      
      if (! karma::generate(iter, ' ' << karma::uint_ << ' ' << karma::uint_ << ' ' << karma::uint_ << ' ' << karma::uint_,
			    scores.aligned_source,
			    scores.aligned_target,
			    scores.unaligned_source,
			    scores.unaligned_target))
	throw std::runtime_error("failed generation");
    }

    if (feature_type_mode) {
      const double logprob_type_source_target = - std::log(phrase_pair.observed_source);
      const double logprob_type_target_source = - std::log(phrase_pair.observed_target);
      
      if (! karma::generate(iter, ' ' << double10 << ' ' << double10, logprob_type_source_target, logprob_type_target_source))
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
      if (! karma::generate(iter, ' ' << karma::uint_ << ' ' << karma::uint_, cross(source, target), cross(source, target, alignments)))
	throw std::runtime_error("failed generation");
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
    
    ("statistic", po::value<path_type>(&statistic_file),                              "significant test statistic")
    
    ("dirichlet-prior", po::value<double>(&dirichlet_prior)->default_value(dirichlet_prior), "dirichlet prior weight")

    ("cicada",   po::bool_switch(&mode_cicada),   "output in cicada format")
    ("moses",    po::bool_switch(&mode_moses),    "output in moses format")
        
    ("feature-root",       po::bool_switch(&feature_root_mode),       "feature by generative probability")
    ("feature-fisher",     po::bool_switch(&feature_fisher_mode),     "Fisher's exact test")
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
