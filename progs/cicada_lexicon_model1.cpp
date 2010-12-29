//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/vector2.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <google/dense_hash_map>


typedef cicada::Symbol   word_type;
typedef cicada::Sentence sentence_type;
typedef cicada::Vocab    vocab_type;
typedef boost::filesystem::path path_type;

typedef double count_type;
typedef double prob_type;

struct ttable_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  struct count_map_type
  {
    typedef google::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type> > counts_type;

    typedef counts_type::value_type      value_type;
    typedef counts_type::size_type       size_type;
    typedef counts_type::difference_type difference_type;
      
    typedef counts_type::mapped_type     mapped_type;
    typedef counts_type::key_type        key_type;
  
    typedef counts_type::const_iterator const_iterator;
    typedef counts_type::iterator       iterator;
  
    typedef counts_type::const_reference const_reference;
    typedef counts_type::reference       reference;
  
    count_map_type() { counts.set_empty_key(word_type()); }

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
  
  ttable_type(const double __smooth=1e-20) : smooth(__smooth) {}
  
  count_map_type& operator[](const word_type& word)
  {
    return ttable[word.id()];
  }
  
  double operator()(const word_type& source, const word_type& target) const
  {
    if (! ttable.exists(source.id())) return smooth;
    
    const count_map_type& counts = ttable[source.id()];
    count_map_type::const_iterator citer = counts.find(target);
    
    return (citer == counts.end() ? smooth : citer->second);
  }
  
  void clear() { ttable.clear(); }
  void swap(ttable_type& x)
  {
    ttable.swap(x.ttable);
    std::swap(smooth, x.smooth);
  }

  size_type size() const { return ttable.size(); }
  bool empty() const { return ttable.empty(); }

  ttable_type& operator+=(const ttable_type& x)
  {
    for (size_type i = 0; i != x.ttable.size(); ++ i) 
      if (x.ttable.exists(i)) {
	count_map_type& counts = ttable[i];
	
	count_map_type::const_iterator titer_end = x.ttable[i].end();
	for (count_map_type::const_iterator titer = x.ttable[i].begin(); titer != titer_end; ++ titer)
	  counts[titer->first] += titer->second;
      }
    
    return *this;
  }
  
  count_dict_type ttable;
  double smooth;
}

path_type source_file = "-";
path_type target_file = "-";
path_type output_source_target_file = "-";
path_type output_target_source_file = "-";

bool individual_mode = false;
bool symmetric_mode = false;
bool posterior_mode = false;

// parameter...
double p0    = 1e-4;
double prior = 1e-4;
double smooth = 1e-20;

int threads = 2;

int debug = 0;

struct LearnIndividual;
struct LearnSymmetric;
struct LearnPosterior;

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(individual_mode) + symmetric_mode + posterior_mode > 1)
      throw std::runtime_error("specify either individual|symmetric|posterior");
    
    if (int(individual_mode) + symmetric_mode + posterior_mode == 0)
      individual_mode = true;
    
    threads = std::max(threads, 1);
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct LearnBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  LearnBase(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source,
	    const double smooth=1e-20)
    : ttable_source_target(__ttable_target_source),
      ttable_target_source(__ttable_source_target),
      objective_source_target(0),
      objective_target_source(0)
  {}
  
  const ttable_type& ttable_source_target;
  const ttable_type& ttable_target_source;
  ttable_type counts_source_target;
  ttable_type counts_target_source;
  
  double objective_source_target;
  double objective_target_source;
};

struct LearnIndividual : public LearnBase
{
  typedef std::vector<double, std::allocator<double> > prob_set_type;

  void operator()(const sentence_type& source, const sentence_type& target)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    probs.resize(source_size + 1);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter = probs.begin();
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      piter = probs.begin();
      counts_source_target[vocab_type::NONE][target[trg]] += (*piter) * factor;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter)
	counts_source_target[source[src]][target[trg]] += (*piter) * factor;
    }
    
    // inverse direction...
    probs.resize(target_size + 1);
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter = probs.begin();
      *piter = ttable_target_source(vocab_type::NONE, source[src]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      piter = probs.begin();
      counts_target_source[vocab_type::NONE][source[src]] += (*piter) * factor;
      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter)
	counts_target_source[target[trg]][source[src]] += (*piter) * factor;
    }
    
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
  }

  prob_set_type probs;
};

struct LearnSymmetric : public LearnBase
{
  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!
    
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1);
    posterior_target_source.resize(source_size + 1, target_size + 1);
    
    prob_source_target.reserve(source_size + 1);
    prob_target_source.reserve(target_size + 1);
    
    prob_source_target.resize(source_size + 1);
    prob_target_source.resize(target_size + 1);

    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_source_target.begin();
      prob_set_type::iterator piter_end = prob_source_target.end();
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin();
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_target_source.begin();
      prob_set_type::iterator piter_end = prob_target_source.end();
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin();
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
    }
    
    // accumulate!
    for (size_type src = 0; src != source_size + 1; ++ src)
      for (size_type trg = (src == 0); trg != target_size + 1; ++ trg) {
	double count = ((trg == 0 ? 1.0 : posterior_source_target(trg, src))
			* (src == 0 ? 1.0 : posterior_target_source(src, trg)));
	
	if (src != 0 && trg != 0)
	  count = utils::mathop::sqrt(count);
	
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  counts_source_target[source_word][target_word] += count;
	
	if (src != 0)
	  counts_target_source[target_word][source_word] += count;
      }
    
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
  }

  prob_set_type      prob_source_target;
  prob_set_type      prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;
};

struct LearnPosterior : public LearnBase
{
  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!
    
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1);
    posterior_target_source.resize(source_size + 1, target_size + 1);
    
    prob_source_target.reserve(source_size + 1);
    prob_target_source.reserve(target_size + 1);
    
    prob_source_target.resize(source_size + 1);
    prob_target_source.resize(target_size + 1);

    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_source_target.begin();
      prob_set_type::iterator piter_end = prob_source_target.end();
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin();
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_target_source.begin();
      prob_set_type::iterator piter_end = prob_target_source.end();
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin();
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
    }
    
    // accumulate!
    for (size_type src = 0; src != source_size + 1; ++ src)
      for (size_type trg = (src == 0); trg != target_size + 1; ++ trg) {
	double count = ((trg == 0 ? 1.0 : posterior_source_target(trg, src))
			* (src == 0 ? 1.0 : posterior_target_source(src, trg)));
	
	if (src != 0 && trg != 0)
	  count = utils::mathop::sqrt(count);
	
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  counts_source_target[source_word][target_word] += count;
	
	if (src != 0)
	  counts_target_source[target_word][source_word] += count;
      }
    
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
  }

  prob_set_type      prob_source_target;
  prob_set_type      prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    
    ("output_source_targert", po::value<path_type>(&output_source_target_file), "output for P(target | source)")
    ("output_targert_source", po::value<path_type>(&output_target_source_file), "output for P(source | target)")
    
    ("individual", po::bool_switch(&individual_mode), "individual model1 training")
    ("symmetric",  po::bool_switch(&symmetric_mode),  "symmetric model1 training")
    ("posterior",  po::bool_switch(&posterior_mode),  "posterior constrained model1 training")
    
    ("p0",     po::value<double>(&p0),     "parameter for NULL alignment")
    ("prior",  po::value<double>(&prior),  "prior for Dirichlet smoothing")
    ("smooth", po::value<double>(&smooth), "smoothing parameter for cutoff probability")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
