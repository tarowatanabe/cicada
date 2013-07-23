//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <map>
#include <cmath>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/compact_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/mathop.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef cicada::Symbol   word_type;
typedef cicada::Sentence sentence_type;
typedef boost::filesystem::path path_type;

typedef uint64_t count_type;
typedef std::vector<count_type, std::allocator<count_type> > count_unigram_type;

struct count_map_type
{
  typedef utils::compact_map<word_type, count_type,
			     utils::unassigned<word_type>, utils::unassigned<word_type>,
			     boost::hash<word_type>, std::equal_to<word_type>,
			     std::allocator<std::pair<const word_type, count_type> > > counts_type;
  
  typedef counts_type::value_type      value_type;
  typedef counts_type::size_type       size_type;
  typedef counts_type::difference_type difference_type;
      
  typedef counts_type::mapped_type     mapped_type;
  typedef counts_type::key_type        key_type;

  typedef counts_type::const_iterator const_iterator;
  typedef counts_type::iterator       iterator;

  typedef counts_type::const_reference const_reference;
  typedef counts_type::reference       reference;
  
  count_map_type() {  }

  inline const_iterator begin() const { return counts.begin(); }
  inline       iterator begin()       { return counts.begin(); }
  inline const_iterator end() const { return counts.end(); }
  inline       iterator end()       { return counts.end(); }

  mapped_type& operator[](const key_type& key) { return counts[key]; }
      
  size_type size() const { return counts.size(); }
  bool empty() const { return counts.empty(); }

  void swap(count_map_type& x) { counts.swap(x.counts); }
  void clear() { counts.clear(); }
  
  counts_type counts;
};

typedef utils::alloc_vector<count_map_type, std::allocator<count_map_type> > count_dict_type;

typedef utils::unordered_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_set_type;


double fisher(const count_type cfe, const count_type cf, const count_type ce, const count_type n)
{
  if (cfe > cf || cfe > ce || cfe <= 0 || cf <= 0 || ce <= 0 || n <= 0 || n < cfe || n < cf || n < ce)
    throw std::runtime_error(std::string("invalid count:")
			     + ' ' + boost::lexical_cast<std::string>(cfe)
			     + ' ' + boost::lexical_cast<std::string>(cf)
			     + ' ' + boost::lexical_cast<std::string>(ce)
			     + ' ' + boost::lexical_cast<std::string>(n));

  count_type a = cfe;
  count_type b = cf - cfe;
  count_type c = ce - cfe;
  count_type d = n - ce - cf + cfe;
    
  double log_p = (utils::mathop::lgamma<double>(1+a+c)
		  + utils::mathop::lgamma<double>(1+a+b)
		  - utils::mathop::lgamma<double>(1+a)
		  - utils::mathop::lgamma<double>(1+b)
		  - utils::mathop::lgamma<double>(1+c)
		  + utils::mathop::lgamma<double>(1+b+d)
		  - utils::mathop::lgamma<double>(1+n)
		  + utils::mathop::lgamma<double>(1+c+d)
		  - utils::mathop::lgamma<double>(1+d));

  if (! std::isfinite(log_p))
    return std::numeric_limits<double>::infinity();
    
  double log_total_p = log_p;
    
  const count_type total_count = utils::bithack::min(b, c);
  for (count_type i = 0; i < total_count; ++ i, ++ a, -- b, -- c, ++ d) {
    log_p += std::log(b) + std::log(c) - std::log(a + 1) - std::log(d + 1);
      
    if (! std::isfinite(log_p))
      return - log_total_p;
      
    log_total_p = utils::mathop::logsum(log_total_p, log_p);
  }
    
  return - log_total_p; 
}

path_type source_file = "-";
path_type target_file = "-";
path_type output_file = "-";
bool inverse = false;
bool normalize = false;
bool cutoff = false;

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    count_unigram_type sources;
    count_unigram_type targets;
    count_dict_type    dict;

    utils::compress_istream is_src(! inverse ? source_file : target_file, 1024 * 1024);
    utils::compress_istream is_trg(! inverse ? target_file : source_file, 1024 * 1024);

    sentence_type source;
    sentence_type target;
    
    word_set_type words_source;
    word_set_type words_target;

    count_type n = 0;
    
    while (1) {
      is_src >> source;
      is_trg >> target;

      if (! is_src || ! is_trg) break;

      if (source.empty() || target.empty()) continue;

      ++ n;
      
      words_source.clear();
      words_source.insert(source.begin(), source.end());

      words_target.clear();
      words_target.insert(target.begin(), target.end());
      
      word_set_type::const_iterator siter_begin = words_source.begin();
      word_set_type::const_iterator siter_end   = words_source.end();
      word_set_type::const_iterator titer_begin = words_target.begin();
      word_set_type::const_iterator titer_end   = words_target.end();

      for (word_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	if (titer->id() >= targets.size())
	  targets.resize(titer->id() + 1, 0);
	++ targets[titer->id()];
      }
      
      for (word_set_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	if (siter->id() >= sources.size())
	  sources.resize(siter->id() + 1, 0);
	++ sources[siter->id()];
	
	for (word_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	  ++ dict[siter->id()][*titer];
      }
    }
    if (is_src || is_trg)
      throw std::runtime_error("# of lines do not match");
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    os.precision(10);
    
    typedef std::multimap<double, word_type, std::greater<double>, std::allocator<std::pair<const double, word_type> > > sorted_type;
    
    sorted_type sorted;

    const double threshold = fisher(1, 1, 1, n) - 0.001;

    if (normalize) {
      // dump..
      count_dict_type::const_iterator siter_begin = dict.begin();
      count_dict_type::const_iterator siter_end   = dict.end();
      for (count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	if (*siter) {
	  const word_type source(word_type::id_type(siter - siter_begin));
	  const count_map_type& target = *(*siter);

	  sorted.clear();

	  double logsum = boost::numeric::bounds<double>::lowest();
	
	  count_map_type::const_iterator titer_end = target.end();
	  for (count_map_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) {
	    const double score = fisher(titer->second, sources[source.id()], targets[titer->first.id()], n);

	    if (cutoff && score < threshold) continue;
	  
	    sorted.insert(std::make_pair(score, titer->first));
	    
	    logsum = utils::mathop::logsum(logsum, score);
	  }
	
	  sorted_type::const_iterator iter_end = sorted.end();
	  for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	    os << iter->second << ' ' << source << ' ' << utils::mathop::exp(iter->first - logsum) << '\n';
	}
    } else {
      // dump..
      count_dict_type::const_iterator siter_begin = dict.begin();
      count_dict_type::const_iterator siter_end   = dict.end();
      for (count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	if (*siter) {
	  const word_type source(word_type::id_type(siter - siter_begin));
	  const count_map_type& target = *(*siter);

	  sorted.clear();
	
	  count_map_type::const_iterator titer_end = target.end();
	  for (count_map_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) {
	    const double score = fisher(titer->second, sources[source.id()], targets[titer->first.id()], n);

	    if (cutoff && score < threshold) continue;
	  
	    sorted.insert(std::make_pair(score, titer->first));
	  }
	
	  sorted_type::const_iterator iter_end = sorted.end();
	  for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	    os << iter->second << ' ' << source << ' ' << iter->first << '\n';
	}
    }    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    ("output", po::value<path_type>(&output_file)->default_value("-"), "dice: for target source")
    ("inverse",   po::bool_switch(&inverse),   "inverse source/target")
    ("normalize", po::bool_switch(&normalize), "normalize as probability")
    ("cutoff",    po::bool_switch(&cutoff),    "cutoff by 1-1-1-N")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
