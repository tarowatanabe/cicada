//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/vector2.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

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
  
  ttable_type(const double __smooth=1e-20) : smooth(__smooth) {}
  
  count_map_type& operator[](const word_type& word)
  {
    return ttable[word.id()];
  }

  const count_map_type& operator[](const word_type& word) const
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
  bool exists(size_type pos) const { return ttable.exists(pos); }
  
  void resize(size_type __size) { ttable.resize(__size); }

  void initialize()
  {
    for (size_type i = 0; i != ttable.size(); ++ i)
      if (ttable.exists(i))
	ttable[i].clear();
  }

  ttable_type& operator+=(const ttable_type& x)
  {
    for (size_type i = 0; i != x.ttable.size(); ++ i) 
      if (x.ttable.exists(i))
	ttable[i] += x.ttable[i];
    
    return *this;
  }
  
  count_dict_type ttable;
  double smooth;
};

struct aligned_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  struct aligned_map_type
  {
    typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > map_type;

    typedef map_type::value_type      value_type;
    typedef map_type::size_type       size_type;
    typedef map_type::difference_type difference_type;
    
    typedef map_type::const_iterator const_iterator;
    typedef map_type::iterator       iterator;
  
    typedef map_type::const_reference const_reference;
    typedef map_type::reference       reference;

    aligned_map_type() { aligned.set_empty_key(word_type()); }
    
    inline const_iterator begin() const { return aligned.begin(); }
    inline       iterator begin()       { return aligned.begin(); }
    inline const_iterator end() const { return aligned.end(); }
    inline       iterator end()       { return aligned.end(); }
    
    inline const_iterator find(const value_type& x) const { return aligned.find(x); }
    inline       iterator find(const value_type& x)       { return aligned.find(x); }
    
    std::pair<iterator, bool> insert(const value_type& x) { return aligned.insert(x); }
    
    size_type size() const { return aligned.size(); }
    bool empty() const { return aligned.empty(); }

    void clear() { aligned.clear(); }
    
    aligned_map_type& operator+=(const aligned_map_type& x)
    {
      const_iterator citer_end = x.aligned.end();
      for (const_iterator citer = x.aligned.begin(); citer != citer_end; ++ citer)
	aligned.insert(*citer);
      
      return *this;
    }
    
    map_type aligned;
  };
  
  typedef utils::alloc_vector<aligned_map_type, std::allocator<aligned_map_type> > aligned_set_type;
  
  aligned_map_type& operator[](const word_type& word)
  {
    return aligned[word.id()];
  }
  
  const aligned_map_type& operator[](const word_type& word) const
  {
    return aligned[word.id()];
  }
  
  size_type size() const { return aligned.size(); }
  bool empty() const { return aligned.empty(); }
  bool exists(size_type pos) const { return aligned.exists(pos); }
  bool exists(const word_type& word) const { return aligned.exists(word.id()); }
  
  void resize(size_type __size) { aligned.resize(__size); }
  void clear() { aligned.clear(); }

  void initialize()
  {
    for (size_type i = 0; i != aligned.size(); ++ i)
      if (aligned.exists(i))
	aligned[i].clear();
  }
  
  aligned_set_type aligned;
};

path_type source_file = "-";
path_type target_file = "-";
path_type output_source_target_file = "-";
path_type output_target_source_file = "-";

int iteration = 5;

bool symmetric_mode = false;
bool posterior_mode = false;
bool variational_bayes_mode = false;

// parameter...
double p0    = 1e-4;
double prior = 0.1;
double smooth = 1e-7;

double threshold = 0.0;
bool   logprob_mode = false;

int threads = 2;

int debug = 0;

struct LearnIndividual;
struct LearnIndividualPosterior;
struct LearnSymmetric;
struct LearnSymmetricPosterior;

struct Maximize;
struct MaximizeBayes;

void dump(const path_type& path, const ttable_type& lexicon, const aligned_type& aligned);
template <typename Learner, typename Maximizer>
void learn(ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source,
	   aligned_type& aligned_source_target,
	   aligned_type& aligned_target_source);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    ttable_type ttable_source_target(smooth);
    ttable_type ttable_target_source(smooth);
    
    aligned_type aligned_source_target;
    aligned_type aligned_target_source;

    if (iteration > 0) {
      if (variational_bayes_mode) {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnSymmetricPosterior, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnSymmetric, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnIndividualPosterior, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnIndividual, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	}
	
      } else {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnSymmetricPosterior, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnSymmetric, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnIndividualPosterior, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnIndividual, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	}
      }
    }
      
    // final dumping...
    boost::thread_group workers_dump;

    if (! output_source_target_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_source_target_file),
							    boost::cref(ttable_source_target),
							    boost::cref(aligned_source_target))));
    
    if (! output_target_source_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_target_source_file),
							    boost::cref(ttable_target_source),
							    boost::cref(aligned_target_source))));
    
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

void dump(const path_type& path, const ttable_type& lexicon, const aligned_type& aligned)
{
  typedef ttable_type::count_map_type::value_type value_type;
  typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;

  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(10);

  const aligned_type::aligned_map_type __empty;
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
      
      if (threshold > 0.0) {
	const double prob_max       = sorted.front()->second;
	const double prob_threshold = prob_max * threshold;
	
	const aligned_type::aligned_map_type& viterbi = (aligned.exists(source) ? aligned[source] : __empty);
	
	// TODO: extra checking to keep Viterbi alignemnt in the final output!
	
	if (logprob_mode) {
	  sorted_type::const_iterator iter_end = sorted.end();
	  for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	    if ((*iter)->second >= prob_threshold || viterbi.find((*iter)->first) != viterbi.end())
	    os << (*iter)->first << ' ' << source << ' '  << std::log((*iter)->second) << '\n';
	} else {
	  sorted_type::const_iterator iter_end = sorted.end();
	  for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	    if ((*iter)->second >= prob_threshold || viterbi.find((*iter)->first) != viterbi.end())
	      os << (*iter)->first << ' ' << source << ' '  << (*iter)->second << '\n';
	}
      } else {
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
}

template <typename LearnerSet, typename Maximizer>
struct TaskMaximize : public Maximizer
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  TaskMaximize(const LearnerSet& __learners,
	       const int __id,
	       ttable_type& __ttable_source_target,
	       ttable_type& __ttable_target_source,
	       aligned_type& __aligned_source_target,
	       aligned_type& __aligned_target_source)
    : learners(__learners),
      id(__id),
      ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      aligned_source_target(__aligned_source_target),
      aligned_target_source(__aligned_target_source) {}
  
  void operator()()
  {
    for (word_type::id_type source_id = id; source_id < ttable_source_target.size(); source_id += learners.size()) {
      for (size_t i = 0; i != learners.size(); ++ i) {
	if (learners[i].counts_source_target.exists(source_id))
	  ttable_source_target[source_id] += learners[i].counts_source_target[source_id];
	
	if (learners[i].aligned_source_target.exists(source_id))
	  aligned_source_target[source_id] += learners[i].aligned_source_target[source_id];
      }
      
      if (ttable_source_target.exists(source_id))
	Maximizer::operator()(ttable_source_target[source_id]);
    }
    
    for (word_type::id_type target_id = id; target_id < ttable_target_source.size(); target_id += learners.size()) {
      for (size_t i = 0; i != learners.size(); ++ i) {
	if (learners[i].counts_target_source.exists(target_id))
	  ttable_target_source[target_id] += learners[i].counts_target_source[target_id];
	
	if (learners[i].aligned_target_source.exists(target_id))
	  aligned_target_source[target_id] += learners[i].aligned_target_source[target_id];
      }
      
      if (ttable_target_source.exists(target_id))
	Maximizer::operator()(ttable_target_source[target_id]);
    }
  }
  
  const LearnerSet& learners;
  const int id;

  ttable_type& ttable_source_target;
  ttable_type& ttable_target_source;
  
  aligned_type& aligned_source_target;
  aligned_type& aligned_target_source;
};

template <typename Learner>
struct TaskLearn : public Learner
{
  struct bitext_type
  {
    bitext_type() : source(), target() {}
    bitext_type(const sentence_type& __source, const sentence_type& __target)
      : source(__source), target(__target) {}
    
    sentence_type source;
    sentence_type target;
  };
  
  typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
  typedef utils::lockfree_list_queue<bitext_set_type, std::allocator<bitext_set_type> > queue_type;
  
  TaskLearn(queue_type& __queue,
	    const ttable_type& ttable_source_target,
	    const ttable_type& ttable_target_source)
    : Learner(ttable_source_target, ttable_target_source), queue(__queue) {}
  
  void operator()()
  {
    Learner::initialize();

    bitext_set_type bitexts;
    
    for (;;) {
      bitexts.clear();
      queue.pop_swap(bitexts);
      if (bitexts.empty()) break;
      
      typename bitext_set_type::const_iterator biter_end = bitexts.end();
      for (typename bitext_set_type::const_iterator biter = bitexts.begin(); biter != biter_end; ++ biter)
	Learner::operator()(biter->source, biter->target);
    }
  }
  
  queue_type& queue;
};

template <typename Learner, typename Maximizer>
void learn(ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source,
	   aligned_type& aligned_source_target,
	   aligned_type& aligned_target_source)
{
  typedef TaskLearn<Learner> learner_type;
  
  typedef typename learner_type::bitext_type     bitext_type;
  typedef typename learner_type::bitext_set_type bitext_set_type;
  typedef typename learner_type::queue_type      queue_type;
  
  typedef std::vector<learner_type, std::allocator<learner_type> > learner_set_type;
  
  typedef TaskMaximize<learner_set_type, Maximizer> maximizer_type;

  queue_type       queue;
  learner_set_type learners(threads, learner_type(queue, ttable_source_target, ttable_target_source));
  
  for (int iter = 0; iter < iteration; ++ iter) {
    if (debug)
      std::cerr << "iteration: " << iter << std::endl;

    utils::resource accumulate_start;
    
    boost::thread_group workers_learn;
    for (size_t i = 0; i != learners.size(); ++ i)
      workers_learn.add_thread(new boost::thread(boost::ref(learners[i])));
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    
    bitext_type     bitext;
    bitext_set_type bitexts;
    
    size_t num_bitext = 0;
    
    for (;;) {
      is_src >> bitext.source;
      is_trg >> bitext.target;
      
      if (! is_src || ! is_trg) break;      
      
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
    aligned_source_target.resize(word_type::allocated());
    aligned_target_source.resize(word_type::allocated());
    
    ttable_source_target.initialize();
    ttable_target_source.initialize();
    aligned_source_target.initialize();
    aligned_target_source.initialize();

    double objective_source_target = 0;
    double objective_target_source = 0;
    
    boost::thread_group workers_maximize;
    for (size_t i = 0; i != learners.size(); ++ i) {
      objective_source_target += learners[i].objective_source_target;
      objective_target_source += learners[i].objective_target_source;
      
      workers_maximize.add_thread(new boost::thread(maximizer_type(learners, i,
								   ttable_source_target, ttable_target_source,
								   aligned_source_target, aligned_target_source)));
    }
    
    if (debug)
      std::cerr << "perplexity for P(target | source): " << objective_source_target << '\n'
		<< "perplexity for P(source | target): " << objective_target_source << '\n';
    
    workers_maximize.join_all();
    
    utils::resource accumulate_end;
    
    if (debug)
      std::cerr << "cpu time:  " << accumulate_end.cpu_time() - accumulate_start.cpu_time() << std::endl
		<< "user time: " << accumulate_end.user_time() - accumulate_start.user_time() << std::endl;
    
  }
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

struct LearnBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  LearnBase(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      objective_source_target(0),
      objective_target_source(0)
  {}

  void initialize()
  {
    counts_source_target.initialize();
    counts_target_source.initialize();
    
    aligned_source_target.initialize();
    aligned_target_source.initialize();
    
    objective_source_target = 0.0;
    objective_target_source = 0.0;
  }
  
  const ttable_type& ttable_source_target;
  const ttable_type& ttable_target_source;
  ttable_type counts_source_target;
  ttable_type counts_target_source;
  
  aligned_type aligned_source_target;
  aligned_type aligned_target_source;
  
  double objective_source_target;
  double objective_target_source;
};

struct LearnIndividual : public LearnBase
{
  LearnIndividual(const ttable_type& __ttable_source_target,
		  const ttable_type& __ttable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source) {}
  
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const ttable_type& ttable,
	     ttable_type& counts,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum = 0.0;
    
    probs.resize(source_size + 1);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      
      prob_set_type::iterator piter = probs.begin();
      *piter = ttable(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;
      
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      piter = probs.begin();
      counts[vocab_type::NONE][target[trg]] += (*piter) * factor;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter)
	counts[source[src]][target[trg]] += (*piter) * factor;
      
      aligned[word_max].insert(target[trg]);
    }
    
    objective += logsum / target_size;
  }
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    learn(source, target, ttable_source_target, counts_source_target, aligned_source_target, objective_source_target);
    learn(target, source, ttable_target_source, counts_target_source, aligned_target_source, objective_target_source);
  }

  prob_set_type probs;
};

struct LearnIndividualPosterior : public LearnBase
{
  LearnIndividualPosterior(const ttable_type& __ttable_source_target,
			   const ttable_type& __ttable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source) {}

  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  typedef std::vector<prob_type, std::allocator<prob_type> > prob_set_type;
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const ttable_type& ttable,
	     ttable_type& counts,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum = 0.0;

    posterior.reserve(target_size + 1, source_size + 1);
    posterior.resize(target_size + 1, source_size + 1, 0.0);
    
    probs.reserve(target_size + 1, source_size + 1);
    probs.resize(target_size + 1, source_size + 1, 0.0);
    
    phi.clear();
    phi.resize(source_size + 1, 0.0);
    
    exp_phi.clear();
    exp_phi.resize(source_size + 1, 1.0);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      posterior_set_type::iterator piter     = probs.begin(trg + 1);
      posterior_set_type::iterator piter_end = probs.end(trg + 1);
      *piter = ttable(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;
      
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      piter = probs.begin(trg + 1);
      posterior_set_type::iterator siter = posterior.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
      
      aligned[word_max].insert(target[trg]);
    }
    
    objective += logsum / target_size;
    
    for (int iter = 0; iter < 5; ++ iter) {
      // update phi.. but ignore NULL...
      
      bool updated = false;
      for (size_type src = 1; src <= source_size; ++ src) {
	double sum = 0.0;
	for (size_type trg = 1; trg <= target_size; ++ trg)
	  sum += posterior(trg, src);
	
	phi[src] += 1.0 - sum;
	if (phi[src] > 0.0)
	  phi[src] = 0.0;
	
	updated |= (phi[src] != 0.0);
	exp_phi[src] = utils::mathop::exp(phi[src]);
      }
      
      if (! updated) break;
      
      for (size_type trg = 1; trg <= target_size; ++ trg) {
	double sum = 0.0;
	for (size_type src = 0; src <= source_size; ++ src)
	  sum += probs(trg, src) * exp_phi[src];
	
	const double factor = 1.0 / sum;
	for (size_type src = 0; src <= source_size; ++ src)
	  posterior(trg, src) = probs(trg, src) * factor * exp_phi[src];
      }
    }
    
    // update...
    for (size_type trg = 1; trg <= target_size; ++ trg)
      for (size_type src = 0; src <= source_size; ++ src)
	counts[src == 0 ? vocab_type::NONE : source[src - 1]][target[trg - 1]] += posterior(trg, src);
  }

  void operator()(const sentence_type& source, const sentence_type& target)
  {
    learn(source, target, ttable_source_target, counts_source_target, aligned_source_target, objective_source_target);
    learn(target, source, ttable_target_source, counts_target_source, aligned_target_source, objective_target_source);
  }

  posterior_set_type posterior;
  posterior_set_type probs;
  
  prob_set_type      phi;
  prob_set_type      exp_phi;
};

struct LearnSymmetric : public LearnBase
{
  LearnSymmetric(const ttable_type& __ttable_source_target,
		 const ttable_type& __ttable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source) {}
      
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
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin();
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
      
      aligned_source_target[word_max].insert(target[trg]);
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_target_source.begin();
      prob_set_type::iterator piter_end = prob_target_source.end();
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;

	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = target[trg];
	}
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin();
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
      
      aligned_target_source[word_max].insert(source[src]);
    }
    
    // accumulate!
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
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

struct LearnSymmetricPosterior : public LearnBase
{
  LearnSymmetricPosterior(const ttable_type& __ttable_source_target,
			  const ttable_type& __ttable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source) {}

  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  
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

    prob_source_target.reserve(target_size + 1, source_size + 1);
    prob_target_source.reserve(source_size + 1, target_size + 1);
    
    prob_source_target.resize(target_size + 1, source_size + 1);
    prob_target_source.resize(source_size + 1, target_size + 1);
    
    phi.clear();
    exp_phi.clear();
    
    phi.resize(target_size + 1, source_size + 1, 0.0);
    exp_phi.resize(target_size + 1, source_size + 1, 1.0);
    
    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      posterior_set_type::iterator piter     = prob_source_target.begin(trg + 1);
      posterior_set_type::iterator piter_end = prob_source_target.end(trg + 1);
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin(trg + 1);
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;

      aligned_source_target[word_max].insert(target[trg]);
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      posterior_set_type::iterator piter     = prob_target_source.begin(src + 1);
      posterior_set_type::iterator piter_end = prob_target_source.end(src + 1);
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = target[trg];
	}
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin(src + 1);
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
      
      aligned_target_source[word_max].insert(source[src]);
    }
    
    // perplexity..
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
    
    // now we will adjust posterior..
    
    for (int iter = 0; iter != 5; ++ iter) {
      
      bool updated = false;
      
      // update phi... we do not consider null alignment!
      for (size_type src = 1; src <= source_size; ++ src)
	for (size_type trg = 1; trg <= target_size; ++ trg) {
	  const double epsi = posterior_source_target(trg, src) - posterior_target_source(src, trg);
	  const double update = - epsi;
	  
	  phi(trg, src) += update;
	  
	  updated |= (phi(trg, src) != 0.0);
	  exp_phi(trg, src) = utils::mathop::exp(phi(trg, src));
	}
      
      if (! updated) break;
      
      // recompute...
      for (size_type trg = 1; trg <= target_size; ++ trg) {
	double prob_sum = 0.0;
	for (size_type src = 0; src <= source_size; ++ src)
	  prob_sum += prob_source_target(trg, src) * exp_phi(trg, src);
	
	const double factor = 1.0 / prob_sum;
	for (size_type src = 0; src <= source_size; ++ src)
	  posterior_source_target(trg, src) = prob_source_target(trg, src) * factor *  exp_phi(trg, src);
      }
      
      for (size_type src = 1; src <= source_size; ++ src) {
	double prob_sum = 0.0;
	for (size_type trg = 0; trg <= target_size; ++ trg)
	  prob_sum += prob_target_source(src, trg) / exp_phi(trg, src);
	
	const double factor = 1.0 / prob_sum;
	for (size_type trg = 0; trg <= target_size; ++ trg)
	  posterior_target_source(src, trg) = prob_target_source(src, trg) * factor / exp_phi(trg, src);
      }
    }
    
    // since we have already adjusted posterior, we simply accumulate individual counts...
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  counts_source_target[source_word][target_word] += posterior_source_target(trg, src);
	
	if (src != 0)
	  counts_target_source[target_word][source_word] += posterior_target_source(src, trg);
      }
  }

  posterior_set_type prob_source_target;
  posterior_set_type prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;
  
  posterior_set_type phi;
  posterior_set_type exp_phi;
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    
    ("output-source-target", po::value<path_type>(&output_source_target_file), "output for P(target | source)")
    ("output-target-source", po::value<path_type>(&output_target_source_file), "output for P(source | target)")

    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max iteration")
    
    ("symmetric",  po::bool_switch(&symmetric_mode),  "symmetric model1 training")
    ("posterior",  po::bool_switch(&posterior_mode),  "posterior constrained model1 training")
    ("variational-bayes", po::bool_switch(&variational_bayes_mode), "variational Bayes estimates")
    
    ("p0",     po::value<double>(&p0)->default_value(p0),         "parameter for NULL alignment")
    ("prior",  po::value<double>(&prior)->default_value(prior),   "Dirichlet prior for variational Bayes")
    ("smooth", po::value<double>(&smooth)->default_value(smooth), "smoothing parameter for uniform distribution")

    ("threshold", po::value<double>(&threshold)->default_value(threshold), "dump with beam-threshold (<= 0.0 implies no beam)")
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
