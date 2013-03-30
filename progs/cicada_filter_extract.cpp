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
#include <boost/math/distributions/hypergeometric.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/unordered_set.hpp>

typedef boost::filesystem::path path_type;

typedef RootCount  root_count_type;
typedef PhrasePair phrase_pair_type;

typedef utils::unordered_set<root_count_type, boost::hash<root_count_type>, std::equal_to<root_count_type>,
			     std::allocator<root_count_type> >::type root_count_set_type;

struct score_phrase_pair_type
{
  double score;
  std::string line;

  score_phrase_pair_type(const double& __score, const std::string& __line)
    : score(__score), line(__line) {}

  score_phrase_pair_type(const double& __score)
    : score(__score), line() {}

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
size_t nbest = 100;
double cutoff = 0.0;
int    types = 0;
double sigtest = 0.0;

bool sigtest_phrase = false;
bool sigtest_scfg   = false;
bool sigtest_ghkm   = false;

path_type root_joint_file;
path_type root_source_file;
path_type root_target_file;

int debug = 0;

template <typename Filter>
void process(const Filter& filter,
	     std::istream& is,
	     std::ostream& os);

struct FilterNone
{
  bool operator()(const phrase_pair_type& phrase_pair) const
  {
    return false;
  }
};

struct FilterCutoff
{
  FilterCutoff(const double& __cutoff, const int __types) : cutoff(__cutoff), types(__types) {}
  
  bool operator()(const phrase_pair_type& phrase_pair) const
  {
    return (phrase_pair.counts.front() < cutoff
	    && (types <= 0
		|| phrase_pair.observed_source > types
		|| phrase_pair.observed_target > types));
  }
  
  const double cutoff;
  const int types;
};

template <typename Extractor>
struct FilterSigtest
{
  FilterSigtest(const root_count_set_type& __root_joint,
		const root_count_set_type& __root_source,
		const root_count_set_type& __root_target,
		const double& __cutoff,
		const double& __sigtest)
    : root_joint(__root_joint),
      root_source(__root_source),
      root_target(__root_target),
      cutoff(__cutoff),
      sigtest(__sigtest) {}

  Extractor extractor;

  bool operator()(const phrase_pair_type& phrase_pair) const
  {
    if (phrase_pair.counts.front() < cutoff) return true;
    
    const std::string source = extractor(phrase_pair.source);
    const std::string target = extractor(phrase_pair.target);
    
    root_count_set_type::const_iterator jiter = root_joint.find(source + target);
    root_count_set_type::const_iterator siter = root_source.find(source);
    root_count_set_type::const_iterator titer = root_target.find(target);
    
    if (jiter == root_joint.end())
      throw std::runtime_error("no root count: " + source + target);
    if (siter == root_source.end())
      throw std::runtime_error("no root count for source: " + source);
    if (titer == root_target.end())
      throw std::runtime_error("no root count for target: " + target);

    if (jiter->counts.size() != 1)
      throw std::runtime_error("invalid root count: " + source + target);
    if (siter->counts.size() != 1)
      throw std::runtime_error("invalid root count for source: " + source);
    if (titer->counts.size() != 1)
      throw std::runtime_error("invalid root count for target: " + target);
    
    const unsigned int n = phrase_pair.observed_source;
    const unsigned int r = phrase_pair.observed_target;
    const unsigned int N = jiter->counts.front();
    
    const double density = boost::math::pdf(boost::math::hypergeometric(r, n, N), 1);
    
    if (debug >= 2)
      std::cerr << "density: " << density
		<< " ||| " << phrase_pair.source << " ||| " << phrase_pair.target
		<< " ||| " << phrase_pair.observed_source << ' ' << phrase_pair.observed_target
		<< std::endl;
	
    
    return density < sigtest;
  }
  
  const root_count_set_type& root_joint;
  const root_count_set_type& root_source;
  const root_count_set_type& root_target;

  const double cutoff;
  const double sigtest;
};

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (nbest <= 0)
      throw std::runtime_error("nbest must be positive...");

    if (sigtest > 0.0) {
      if (int(sigtest_phrase) + sigtest_scfg + sigtest_ghkm != 1)
	throw std::runtime_error("specify either one of --sigtest-phrase|scfg|ghkm");
      
      if (! boost::filesystem::exists(root_joint_file))
	throw std::runtime_error("no root count file");
      if (! boost::filesystem::exists(root_source_file))
	throw std::runtime_error("no root count file for source side");
      if (! boost::filesystem::exists(root_target_file))
	throw std::runtime_error("no root count file for target side");
    }
    
    root_count_set_type root_joint;
    root_count_set_type root_source;
    root_count_set_type root_target;

    {
      root_count_type root_count;
      RootCountParser parser;
      std::string line;
      
      utils::compress_istream is_joint(root_joint_file);
      while (std::getline(is_joint, line))
	if (parser(line, root_count))
	  root_joint.insert(root_count);
      
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

    if (sigtest > 0.0) {
      if (sigtest_ghkm)
	process(FilterSigtest<ExtractRootGHKM>(root_joint, root_source, root_target, cutoff, sigtest), is, os);
      else if (sigtest_scfg)
	process(FilterSigtest<ExtractRootSCFG>(root_joint, root_source, root_target, cutoff, sigtest), is, os);
      else
	process(FilterSigtest<ExtractRootPhrase>(root_joint, root_source, root_target, cutoff, sigtest), is, os);
    } else if (cutoff > 0.0)
      process(FilterCutoff(cutoff, types), is, os);
    else
      process(FilterNone(), is, os);
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Filter>
void process(const Filter& filter,
	     std::istream& is,
	     std::ostream& os)
{
  std::string          source_prev;
  phrase_pair_type     phrase_pair;
  
  heap_type heap;
  double score_min = std::numeric_limits<double>::infinity();
  
  PhrasePairParser  parser;
  std::string line;
  
  while (std::getline(is, line)) {
    if (! parser(line, phrase_pair)) continue;
    if (phrase_pair.counts.empty()) continue;
    if (filter(phrase_pair)) continue;
    
    if (phrase_pair.source != source_prev) {
      if (! heap.empty()) {
	if (heap.size() <= nbest) {
	  heap_type::iterator iter_end = heap.end();
	  for (heap_type::iterator iter = heap.begin(); iter != iter_end; ++ iter)
	    os << iter->line << '\n';
	} else {
	  heap_type::iterator iter_begin = heap.begin();
	  heap_type::iterator iter_kbest = heap.end() - nbest;
	  heap_type::iterator iter       = heap.end();
	    
	  for (/**/; iter_kbest != iter; -- iter) {
	    os << iter_begin->line << '\n';
	    std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
	  }
	    
	  const double threshold = iter->score;
	    
	  for (/**/; iter_begin != iter && iter_begin->score == threshold; -- iter) {
	    os << iter_begin->line << '\n';
	    std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
	  }
	}
      }
	
      heap.clear();
      score_min = std::numeric_limits<double>::infinity();
      source_prev = phrase_pair.source;
    }
      
    // push-heap and swap line
    // we memorize the temporary score_min and perform pruning
    const double& score = phrase_pair.counts.front();
      
    if (heap.size() >= nbest) {
      if (score < score_min)
	continue;
    } else
      score_min = std::min(score_min, score);
      
    heap.push_back(score_phrase_pair_type(score));
    heap.back().line.swap(line);
    std::push_heap(heap.begin(), heap.end(), std::less<score_phrase_pair_type>());
  }
    
  if (! heap.empty()) {
    if (heap.size() <= nbest) {
      heap_type::iterator iter_end = heap.end();
      for (heap_type::iterator iter = heap.begin(); iter != iter_end; ++ iter)
	os << iter->line << '\n';
    } else {
      heap_type::iterator iter_begin = heap.begin();
      heap_type::iterator iter_kbest = heap.end() - nbest;
      heap_type::iterator iter       = heap.end();
	    
      for (/**/; iter_kbest != iter; -- iter) {
	os << iter_begin->line << '\n';
	std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
      }
	
      const double threshold = iter->score;
	
      for (/**/; iter_begin != iter && iter_begin->score == threshold; -- iter) {
	os << iter_begin->line << '\n';
	std::pop_heap(iter_begin, iter, std::less<score_phrase_pair_type>());
      }
    }
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("nbest",   po::value<size_t>(&nbest)->default_value(nbest),     "nbest of pairs (wrt to joint-count)")
    ("cutoff",  po::value<double>(&cutoff)->default_value(cutoff),   "cutoff count")
    ("types",   po::value<int>(&types)->default_value(types),        "cutoff variation")
    ("sigtest", po::value<double>(&sigtest)->default_value(sigtest), "significant test threshold")
    
    ("sigtest-phrase", po::bool_switch(&sigtest_phrase), "significant test for phrase")
    ("sigtest-scfg",   po::bool_switch(&sigtest_scfg),   "significant test for synchronous-CFG")
    ("sigtest-ghkm",   po::bool_switch(&sigtest_ghkm),   "significant test for ghkm")
    
    ("root-joint",  po::value<path_type>(&root_joint_file),  "root count file")
    ("root-source", po::value<path_type>(&root_source_file), "root source file")
    ("root-target", po::value<path_type>(&root_target_file), "root target file")
    
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
