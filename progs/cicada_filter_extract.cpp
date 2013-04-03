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
#include <queue>

#include <boost/program_options.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/unordered_set.hpp>
#include <utils/mathop.hpp>

typedef boost::filesystem::path path_type;

typedef Statistic statistic_type;

typedef PhrasePair phrase_pair_type;

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
size_t nbest = 0;
double cutoff = 0.0;
double threshold = 0.0;
int    types = 0;
double sigtest = 0.0;
path_type statistic_file;

int debug = 0;

template <typename Filter>
void process(const Filter& filter,
	     std::istream& is,
	     std::ostream& os);

template <typename Filter>
void process_kbest(const Filter& filter,
		   std::istream& is,
		   std::ostream& os);

struct FilterNone
{
  bool operator()(const phrase_pair_type& phrase_pair) const
  {
    return false;
  }
};

struct FilterThreshold
{
  FilterThreshold(const double& __cutoff) : cutoff(__cutoff) {}
  
  bool operator()(const phrase_pair_type& phrase_pair) const
  {
    return ((phrase_pair.counts.front() / phrase_pair.counts_source.front()) < cutoff
	    && (phrase_pair.counts.front() / phrase_pair.counts_target.front()) < cutoff);
  }
  
  const double cutoff;
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

struct FilterSigtest
{
  typedef int64_t count_type;
  
  FilterSigtest(const statistic_type& __statistic,
		const double& alpha)
    : statistic(__statistic),
      threshold()
  {
    threshold = - logfisher(1, 1, 1) + alpha;
  }
  
  bool operator()(const phrase_pair_type& phrase_pair) const
  {
    if (phrase_pair.counts.size() <= 3)
      throw std::runtime_error("invalid counts");
    if (phrase_pair.counts_source.size() <= 3)
      throw std::runtime_error("invalid source counts");
    if (phrase_pair.counts_target.size() <= 3)
      throw std::runtime_error("invalid target counts");

    const count_type cfe = phrase_pair.counts[phrase_pair.counts.size() - 3];
    const count_type cf  = phrase_pair.counts_source[phrase_pair.counts_source.size() - 2];
    const count_type ce  = phrase_pair.counts_target[phrase_pair.counts_target.size() - 1];
    
    const double score = - logfisher(cfe, cf, ce);

    if (debug >= 2)
      std::cerr << phrase_pair.source << " ||| " << phrase_pair.target
		<< " ||| "
		<< (score < threshold ? "true" : "false")
		<< ' ' << cfe << ' ' << cf << ' ' << ce << ' ' << score << ' ' << threshold
		<< std::endl;
    
    return score < threshold;
  }

  double logfisher(const count_type cfe, const count_type cf, const count_type ce) const
  {
    count_type a = cfe;
    count_type b = cf - cfe;
    count_type c = ce - cfe;
    count_type d = statistic.bitext - ce - cf + cfe;
    const count_type n = statistic.bitext;
    
    double log_p = (utils::mathop::lgamma<double>(1+a+c)
		    + utils::mathop::lgamma<double>(1+b+d)
		    + utils::mathop::lgamma<double>(1+a+b)
		    + utils::mathop::lgamma<double>(1+c+d)
		    - utils::mathop::lgamma<double>(1+n)
		    - utils::mathop::lgamma<double>(1+a)
		    - utils::mathop::lgamma<double>(1+b)
		    - utils::mathop::lgamma<double>(1+c)
		    - utils::mathop::lgamma<double>(1+d));
    
    double log_total_p = boost::numeric::bounds<double>::lowest();
    
    const count_type tc = utils::bithack::min(b, c);
    for (count_type i = 0; i <= tc; ++ i) {
      log_total_p = utils::mathop::logsum(log_total_p, log_p);
      
      log_p += (utils::mathop::log<double>(b)
		+ utils::mathop::log<double>(c)
		- utils::mathop::log<double>(a + 1)
		- utils::mathop::log<double>(d + 1));
      
      ++ a;
      -- b;
      -- c;
      ++ d;
    }
    
    return log_total_p;
  }
  
  const statistic_type statistic;
  double threshold;
};

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (sigtest != 0.0) {
      if (statistic_file.empty() || ! boost::filesystem::exists(statistic_file))
	throw std::runtime_error("no statistic for sigtest?");
    }
    
    statistic_type statistic;
    if (! statistic_file.empty()) {
      utils::compress_istream is(statistic_file);
      is >> statistic;
    }
    
    utils::compress_istream is(input_file,  1024 * 1024);
    utils::compress_ostream os(output_file, buffer_size);

    if (nbest > 0) {
      if (sigtest != 0.0)
	process_kbest(FilterSigtest(statistic, sigtest), is, os);
      else if (threshold > 0.0)
	process_kbest(FilterThreshold(threshold), is, os);
      else if (cutoff > 0.0)
	process_kbest(FilterCutoff(cutoff, types), is, os);
      else
	process(FilterNone(), is, os);
    } else {
      if (sigtest != 0.0)
	process(FilterSigtest(statistic, sigtest), is, os);
      else if (threshold > 0.0)
	process(FilterThreshold(threshold), is, os);
      else if (cutoff > 0.0)
	process(FilterCutoff(cutoff, types), is, os);
      else
	process(FilterNone(), is, os);
    }
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
  phrase_pair_type     phrase_pair;
  
  PhrasePairParser  parser;
  std::string line;
  
  while (std::getline(is, line)) {
    if (! parser(line, phrase_pair)) continue;
    if (phrase_pair.counts.empty()) continue;
    if (filter(phrase_pair)) continue;
    
    os << line << '\n';
  }
}

template <typename Filter>
void process_kbest(const Filter& filter,
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
    
    ("nbest",     po::value<size_t>(&nbest)->default_value(nbest),         "nbest of pairs (wrt to joint-count)")
    ("cutoff",    po::value<double>(&cutoff)->default_value(cutoff),       "cutoff count")
    ("threshold", po::value<double>(&threshold)->default_value(threshold), "probability threshold")
    ("types",     po::value<int>(&types)->default_value(types),            "cutoff variation")
    ("sigtest",   po::value<double>(&sigtest)->default_value(sigtest),     "significant test threshold")
    
    ("statistic", po::value<path_type>(&statistic_file),                   "significant test statistic")
    
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
