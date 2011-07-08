//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// learn lexical weight table from alignment data...
//
//
// the code is based on cicada-lexicon-model1, but without count accumulation....
//


#include <cicada/sentence.hpp>
#include <cicada/alignment.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include <google/dense_hash_map>

typedef cicada::Symbol    word_type;
typedef cicada::Sentence  sentence_type;
typedef cicada::Alignment alignment_type;
typedef cicada::Vocab     vocab_type;
typedef boost::filesystem::path path_type;

typedef double count_type;

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

    inline const_iterator find(const key_type& x) const { return counts.find(x); }
    inline       iterator find(const key_type& x)       { return counts.find(x); }
    
    mapped_type& operator[](const key_type& key) { return counts[key]; }
    
    size_type size() const { return counts.size(); }
    bool empty() const { return counts.empty(); }

    void swap(count_map_type& x) { counts.swap(x.counts); }
    void clear() { counts.clear(); }

    count_map_type& operator+=(const count_map_type& x)
    {
      const_iterator citer_end = x.counts.end();
      for (const_iterator citer = x.counts.begin(); citer != citer_end; ++ citer)
	counts[citer->first] += citer->second;
      return *this;
    }
    
    counts_type counts;
  };
  
  typedef utils::alloc_vector<count_map_type, std::allocator<count_map_type> > count_dict_type;
  
  ttable_type() {}
  
  count_map_type& operator[](const word_type& word)
  {
    return ttable[word.id()];
  }

  const count_map_type& operator[](const word_type& word) const
  {
    return ttable[word.id()];
  }  
  
  void clear() { ttable.clear(); }
  void swap(ttable_type& x)
  {
    ttable.swap(x.ttable);
  }
  
  size_type size() const { return ttable.size(); }
  bool empty() const { return ttable.empty(); }
  bool exists(size_type pos) const { return ttable.exists(pos); }
  bool exists(const word_type& word) const { return ttable.exists(word.id()); }
  
  void resize(size_type __size) { ttable.resize(__size); }
  
  
  ttable_type& operator+=(const ttable_type& x)
  {
    for (size_type i = 0; i != x.ttable.size(); ++ i) 
      if (x.ttable.exists(i))
	ttable[i] += x.ttable[i];
    
    return *this;
  }
  
  count_dict_type ttable;
};


path_type source_file = "-";
path_type target_file = "-";
path_type alignment_file = "-";
path_type output_source_target_file = "-";
path_type output_target_source_file = "-";

bool variational_bayes_mode = false;
double prior = 0.1;
bool inverse_mode = false;
bool logprob_mode = false;

int threads = 2;

int debug = 0;

struct Maximize;
struct MaximizeBayes;


void dump(const path_type& path, const ttable_type& lexicon);

template <typename Maximizer>
void learn(ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    ttable_type ttable_source_target;
    ttable_type ttable_target_source;

    if (variational_bayes_mode)
      learn<MaximizeBayes>(ttable_source_target, ttable_target_source);
    else
      learn<Maximize>(ttable_source_target, ttable_target_source);
    
    
    // final dumping...
    boost::thread_group workers_dump;

    if (! output_source_target_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_source_target_file),
							    boost::cref(ttable_source_target))));
    
    if (! output_target_source_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_target_source_file),
							    boost::cref(ttable_target_source))));
    
    workers_dump.join_all();
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct greater_second
{
  template <typename Tp>
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

void dump(const path_type& path, const ttable_type& lexicon)
{
  typedef ttable_type::count_map_type::value_type value_type;
  typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;

  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(10);
  
  sorted_type sorted;
  
  ttable_type::count_dict_type::const_iterator siter_begin = lexicon.ttable.begin();
  ttable_type::count_dict_type::const_iterator siter_end   = lexicon.ttable.end();
  for (ttable_type::count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) 
    if (*siter) {
      const word_type source(word_type::id_type(siter - siter_begin));
      const ttable_type::count_map_type& dict = *(*siter);
      
      if (dict.empty()) continue;
      
      sorted.clear();
      sorted.reserve(dict.size());
      
      ttable_type::count_map_type::const_iterator titer_end = dict.end();
      for (ttable_type::count_map_type::const_iterator titer = dict.begin(); titer != titer_end; ++ titer)
	if (titer->second >= 0.0)
	  sorted.push_back(&(*titer));
      
      std::sort(sorted.begin(), sorted.end(), greater_second());
      
      if (logprob_mode) {
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  os << (*iter)->first << ' ' << source << ' '  << std::log((*iter)->second) << '\n';
      } else {
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  os << (*iter)->first << ' ' << source << ' '  << (*iter)->second << '\n';
      }
    }
}


template <typename LearnerSet, typename Maximizer>
struct TaskMaximize : public Maximizer
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  TaskMaximize(const LearnerSet& __learners,
	       const int __id,
	       ttable_type& __ttable_source_target,
	       ttable_type& __ttable_target_source)
    : learners(__learners),
      id(__id),
      ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source) {}
  
  void operator()()
  {
    for (word_type::id_type source_id = id; source_id < ttable_source_target.size(); source_id += learners.size()) {
      for (size_t i = 0; i != learners.size(); ++ i)
	if (learners[i].counts_source_target.exists(source_id))
	  ttable_source_target[source_id] += learners[i].counts_source_target[source_id];
      
      if (ttable_source_target.exists(source_id))
	Maximizer::operator()(ttable_source_target[source_id]);
    }
    
    for (word_type::id_type target_id = id; target_id < ttable_target_source.size(); target_id += learners.size()) {
      for (size_t i = 0; i != learners.size(); ++ i)
	if (learners[i].counts_target_source.exists(target_id))
	  ttable_target_source[target_id] += learners[i].counts_target_source[target_id];
      
      if (ttable_target_source.exists(target_id))
	Maximizer::operator()(ttable_target_source[target_id]);
    }
  }
  
  const LearnerSet& learners;
  const int id;

  ttable_type& ttable_source_target;
  ttable_type& ttable_target_source;
};

struct TaskLearn
{
  struct bitext_type
  {
    bitext_type() : source(), target(), alignment() {}
    bitext_type(const sentence_type& __source, const sentence_type& __target, const alignment_type& __alignment)
      : source(__source), target(__target), alignment(__alignment) {}
    
    sentence_type  source;
    sentence_type  target;
    alignment_type alignment;
  };
  
  typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
  typedef utils::lockfree_list_queue<bitext_set_type, std::allocator<bitext_set_type> > queue_type;
  
  TaskLearn(queue_type& __queue) : queue(__queue) {}
  
  void operator()()
  {
    typedef std::vector<int, std::allocator<int> > aligned_type;
    
    bitext_set_type bitexts;
    aligned_type aligned_source;
    aligned_type aligned_target;

    alignment_type align;
    
    for (;;) {
      bitexts.clear();
      queue.pop_swap(bitexts);
      if (bitexts.empty()) break;
      
      bitext_set_type::iterator biter_end = bitexts.end();
      for (bitext_set_type::iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {
	
	if (inverse_mode) {
	  align.clear();
	  alignment_type::const_iterator aiter_end = biter->alignment.end();
	  for (alignment_type::const_iterator aiter = biter->alignment.begin(); aiter != aiter_end; ++ aiter)
	    align.push_back(std::make_pair(aiter->target, aiter->source));
	  biter->alignment.swap(align);
	}

	std::sort(biter->alignment.begin(), biter->alignment.end());
	
	aligned_source.clear();
	aligned_target.clear();
	aligned_source.resize(biter->source.size(), 0);
	aligned_target.resize(biter->target.size(), 0);
	
	alignment_type::const_iterator aiter_prev = biter->alignment.end();
	alignment_type::const_iterator aiter_end = biter->alignment.end();
	for (alignment_type::const_iterator aiter = biter->alignment.begin(); aiter != aiter_end; ++ aiter) {
	  if (aiter_prev == aiter_end || *aiter != *aiter_prev) {
	    if (aiter->source >= static_cast<int>(biter->source.size()))
	      throw std::runtime_error("invalid word alignment");
	    if (aiter->target >= static_cast<int>(biter->target.size()))
	      throw std::runtime_error("invalid word alignment");
	    
	    const word_type& source_word = biter->source[aiter->source];
	    const word_type& target_word = biter->target[aiter->target];
	    
	    ++ aligned_source[aiter->source];
	    ++ aligned_target[aiter->target];
	    
	    ++ counts_source_target[source_word][target_word];
	    ++ counts_target_source[target_word][source_word];
	  }
	  
	  aiter_prev = aiter;
	}
	
	for (size_t trg = 0; trg != biter->target.size(); ++ trg)
	  if (! aligned_target[trg]) {
	    const word_type& target_word = biter->target[trg];
	    
	    ++ counts_source_target[vocab_type::EPSILON][target_word];
	    ++ counts_target_source[target_word][vocab_type::EPSILON];
	  }
	
	for (size_t src = 0; src != biter->source.size(); ++ src)
	  if (! aligned_source[src]) {
	    const word_type& source_word = biter->source[src];
	    
	    ++ counts_source_target[source_word][vocab_type::EPSILON];
	    ++ counts_target_source[vocab_type::EPSILON][source_word];
	  }
      }
    }
  }
  
  queue_type& queue;
  ttable_type counts_source_target;
  ttable_type counts_target_source;
};

template <typename Maximizer>
void learn(ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source)
{
  typedef TaskLearn learner_type;
  
  typedef typename learner_type::bitext_type     bitext_type;
  typedef typename learner_type::bitext_set_type bitext_set_type;
  typedef typename learner_type::queue_type      queue_type;
  
  typedef std::vector<learner_type, std::allocator<learner_type> > learner_set_type;
  
  typedef TaskMaximize<learner_set_type, Maximizer> maximizer_type;
  
  queue_type       queue;
  learner_set_type learners(threads, learner_type(queue));
  
  utils::resource accumulate_start;
  
  boost::thread_group workers_learn;
  for (size_t i = 0; i != learners.size(); ++ i)
    workers_learn.add_thread(new boost::thread(boost::ref(learners[i])));
  
  utils::compress_istream is_src(source_file, 1024 * 1024);
  utils::compress_istream is_trg(target_file, 1024 * 1024);
  utils::compress_istream is_align(alignment_file, 1024 * 1024);
  
  bitext_type     bitext;
  bitext_set_type bitexts;
  
  size_t num_bitext = 0;
  
  for (;;) {
    is_src   >> bitext.source;
    is_trg   >> bitext.target;
    is_align >> bitext.alignment;
    
    if (! is_src || ! is_trg || ! is_align) break;
    
    if (bitext.source.empty() || bitext.target.empty()) continue;
    
    bitexts.push_back(bitext);

    ++ num_bitext;
    if (debug) {
      if (num_bitext % 10000 == 0)
	std::cerr << '.';
      if (num_bitext % 1000000 == 0)
	std::cerr << '\n';
    }
      
    if (bitexts.size() == 64) {
      queue.push_swap(bitexts);
      bitexts.clear();
    }
  }
  
  if (! bitexts.empty())
    queue.push_swap(bitexts);
  
  if (debug && num_bitext >= 10000)
    std::cerr << std::endl;
  if (debug)
    std::cerr << "# of bitexts: " << num_bitext << std::endl;
  
  for (size_t i = 0; i != learners.size(); ++ i) {
    bitexts.clear();
    queue.push_swap(bitexts);
  }
  
  workers_learn.join_all();
  
  // merge and normalize...! 
  ttable_source_target.resize(word_type::allocated());
  ttable_target_source.resize(word_type::allocated());
  
  boost::thread_group workers_maximize;
  for (size_t i = 0; i != learners.size(); ++ i)
    workers_maximize.add_thread(new boost::thread(maximizer_type(learners, i,
								 ttable_source_target, ttable_target_source)));
  workers_maximize.join_all();
  
  utils::resource accumulate_end;
  
  if (debug)
    std::cerr << "cpu time:  " << accumulate_end.cpu_time() - accumulate_start.cpu_time() << std::endl
	      << "user time: " << accumulate_end.user_time() - accumulate_start.user_time() << std::endl;
}

struct Maximize
{
  template <typename Counts>
  void operator()(Counts& counts)
  {
    double sum = 0.0;
    typename Counts::iterator iter_end = counts.end();
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior;
    
    const double factor = 1.0 / sum;
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      iter->second = (iter->second + prior) * factor;
  }
};

struct MaximizeBayes
{
  template <typename Counts>
  void operator()(Counts& counts)
  {
    double sum = 0.0;
    typename Counts::iterator iter_end = counts.end();
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior;
    
    const double sum_digamma = utils::mathop::digamma(sum);
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      iter->second = utils::mathop::exp(utils::mathop::digamma(iter->second + prior) - sum_digamma);
  }
};


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source",    po::value<path_type>(&source_file), "source file")
    ("target",    po::value<path_type>(&target_file), "target file")
    ("alignment", po::value<path_type>(&alignment_file), "alignment file")
    
    ("output-source-target", po::value<path_type>(&output_source_target_file), "output for P(target | source)")
    ("output-target-source", po::value<path_type>(&output_target_source_file), "output for P(source | target)")
    
    ("variational-bayes", po::bool_switch(&variational_bayes_mode), "variational Bayes estimates")
    
    ("prior",  po::value<double>(&prior)->default_value(prior),   "Dirichlet prior for variational Bayes")
    
    ("inverse",   po::bool_switch(&inverse_mode),                          "inverse alignment")
    ("logprob",   po::bool_switch(&logprob_mode),                          "dump in log-domain")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
