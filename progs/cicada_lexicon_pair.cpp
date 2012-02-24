//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// learn lexical pair table from corpus
//
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
#include "utils/dense_hash_set.hpp"
#include "utils/dense_hash_map.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

typedef cicada::Symbol    word_type;
typedef cicada::Sentence  sentence_type;
typedef cicada::Alignment alignment_type;
typedef cicada::Vocab     vocab_type;
typedef boost::filesystem::path path_type;

typedef uint64_t count_type;

struct ttable_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  struct count_map_type
  {
    typedef utils::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type> >::type counts_type;

    typedef counts_type::value_type      value_type;
    typedef counts_type::size_type       size_type;
    typedef counts_type::difference_type difference_type;
      
    typedef counts_type::key_type        key_type;
    typedef counts_type::mapped_type     mapped_type;
  
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


struct ttable_pair_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef std::pair<word_type, word_type> word_pair_type;
  
  struct count_map_type
  {
    typedef utils::dense_hash_map<word_pair_type, count_type, utils::hashmurmur<size_t>, std::equal_to<word_pair_type> >::type counts_type;

    typedef counts_type::value_type      value_type;
    typedef counts_type::size_type       size_type;
    typedef counts_type::difference_type difference_type;
      
    typedef counts_type::key_type        key_type;
    typedef counts_type::mapped_type     mapped_type;
  
    typedef counts_type::const_iterator const_iterator;
    typedef counts_type::iterator       iterator;
  
    typedef counts_type::const_reference const_reference;
    typedef counts_type::reference       reference;
  
    count_map_type() { counts.set_empty_key(word_pair_type()); }

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
  
  ttable_pair_type() {}
  
  count_map_type& operator[](const word_type& word)
  {
    return ttable[word.id()];
  }

  const count_map_type& operator[](const word_type& word) const
  {
    return ttable[word.id()];
  }  
  
  void clear() { ttable.clear(); }
  void swap(ttable_pair_type& x)
  {
    ttable.swap(x.ttable);
  }
  
  size_type size() const { return ttable.size(); }
  bool empty() const { return ttable.empty(); }
  bool exists(size_type pos) const { return ttable.exists(pos); }
  bool exists(const word_type& word) const { return ttable.exists(word.id()); }
  
  void resize(size_type __size) { ttable.resize(__size); }
  
  
  ttable_pair_type& operator+=(const ttable_pair_type& x)
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
path_type output_file;
path_type output_prefix_file;
path_type output_suffix_file;

int kbest = 0;

bool inverse_mode = false;

int threads = 2;

int debug = 0;

void dump_pair(const path_type& path, const ttable_pair_type& lexicon);
void dump(const path_type& path, const ttable_type& lexicon);

void learn(ttable_type& ttable, ttable_pair_type& ttable_prefix, ttable_pair_type& ttable_suffix);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    ttable_type      ttable;
    ttable_pair_type ttable_prefix;
    ttable_pair_type ttable_suffix;
      
    learn(ttable, ttable_prefix, ttable_suffix);
      
    // final dumping...
    boost::thread_group workers_dump;

    if (! output_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_file),
							    boost::cref(ttable))));
       
    if (! output_prefix_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump_pair,
							    boost::cref(output_prefix_file),
							    boost::cref(ttable_prefix))));
      
    if (! output_suffix_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump_pair,
							    boost::cref(output_suffix_file),
							    boost::cref(ttable_suffix))));
      
    workers_dump.join_all();
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Tp>
struct compare_pvalue
{
  bool operator()(const Tp& x, const Tp& y) const
  {
    return x.first->second > y.first->second;
  }
};

void dump(const path_type& path, const ttable_type& lexicon)
{
  // we do sorting...
  typedef std::pair<const ttable_type::count_map_type::value_type*, word_type> pvalue_type;
  typedef std::vector<pvalue_type, std::allocator<pvalue_type> > sorted_type;
  
  sorted_type sorted;
  
  ttable_type::count_dict_type::const_iterator siter_begin = lexicon.ttable.begin();
  ttable_type::count_dict_type::const_iterator siter_end   = lexicon.ttable.end();
  for (ttable_type::count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) 
    if (*siter) {
      const word_type source(word_type::id_type(siter - siter_begin));
      const ttable_type::count_map_type& dict = *(*siter);
      
      if (dict.empty()) continue;
      
      ttable_type::count_map_type::const_iterator titer_end = dict.end();
      for (ttable_type::count_map_type::const_iterator titer = dict.begin(); titer != titer_end; ++ titer)
	sorted.push_back(std::make_pair(&(*titer), source));
    }
  
  std::sort(sorted.begin(), sorted.end(), compare_pvalue<pvalue_type>());
  
  utils::compress_ostream os(path, 1024 * 1024);
  
  if (kbest <= 0 || sorted.size() <= static_cast<size_t>(kbest)) {
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << siter->first->first << ' ' << siter->second << ' ' << siter->first->second << '\n';
  } else {
    sorted_type::const_iterator siter_end  = sorted.end();
    sorted_type::const_iterator siter_last = sorted.begin() + kbest;
    
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_last; ++ siter)
      os << siter->first->first << ' ' << siter->second << ' ' << siter->first->second << '\n';
    for (sorted_type::const_iterator siter = siter_last; siter != siter_end && siter->first->second == siter_last->first->second; ++ siter)
      os << siter->first->first << ' ' << siter->second << ' ' << siter->first->second << '\n';
    
  }
}

void dump_pair(const path_type& path, const ttable_pair_type& lexicon)
{
  // we do sorting...
  typedef std::pair<const ttable_pair_type::count_map_type::value_type*, word_type> pvalue_type;
  typedef std::vector<pvalue_type, std::allocator<pvalue_type> > sorted_type;

  sorted_type sorted;
  
  ttable_pair_type::count_dict_type::const_iterator siter_begin = lexicon.ttable.begin();
  ttable_pair_type::count_dict_type::const_iterator siter_end   = lexicon.ttable.end();
  for (ttable_pair_type::count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) 
    if (*siter) {
      const word_type source(word_type::id_type(siter - siter_begin));
      const ttable_pair_type::count_map_type& dict = *(*siter);
      
      if (dict.empty()) continue;
      
      ttable_pair_type::count_map_type::const_iterator titer_end = dict.end();
      for (ttable_pair_type::count_map_type::const_iterator titer = dict.begin(); titer != titer_end; ++ titer)
	sorted.push_back(std::make_pair(&(*titer), source));
    }
  
  std::sort(sorted.begin(), sorted.end(), compare_pvalue<pvalue_type>());
    
  utils::compress_ostream os(path, 1024 * 1024);

  if (kbest <= 0 || sorted.size() <= static_cast<size_t>(kbest)) {
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << siter->first->first.first << ' ' << siter->first->first.second << ' ' << siter->second << ' ' << siter->first->second << '\n';
  } else {
    sorted_type::const_iterator siter_end = sorted.end();
    sorted_type::const_iterator siter_last = sorted.begin() + kbest;
    
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_last; ++ siter)
      os << siter->first->first.first << ' ' << siter->first->first.second << ' ' << siter->second << ' ' << siter->first->second << '\n';
    
    for (sorted_type::const_iterator siter = siter_last; siter != siter_end && siter->first->second == siter_last->first->second; ++ siter)
      os << siter->first->first.first << ' ' << siter->first->first.second << ' ' << siter->second << ' ' << siter->first->second << '\n';
  }
}

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
    bitext_set_type bitexts;

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
		
	alignment_type::const_iterator aiter_prefix = biter->alignment.end();
	alignment_type::const_iterator aiter_end = biter->alignment.end();
	for (alignment_type::const_iterator aiter = biter->alignment.begin(); aiter != aiter_end; ++ aiter) {
	  if (aiter_prefix == aiter_end || *aiter != *aiter_prefix) {
	    if (aiter->source >= static_cast<int>(biter->source.size()))
	      throw std::runtime_error("invalid word alignment");
	    if (aiter->target >= static_cast<int>(biter->target.size()))
	      throw std::runtime_error("invalid word alignment");
	    
	    const word_type& source_word = biter->source[aiter->source];
	    const word_type& source_prefix = (aiter->source == 0 ? vocab_type::BOS : biter->source[aiter->source - 1]);
	    const word_type& source_suffix = (aiter->source + 1 == biter->source.size() ? vocab_type::EOS : biter->source[aiter->source + 1]);
	    const word_type& target_word = biter->target[aiter->target];
	    
	    ++ ttable[target_word][source_word];
	    ++ ttable_prefix[target_word][std::make_pair(source_prefix, source_word)];
	    ++ ttable_suffix[target_word][std::make_pair(source_word, source_suffix)];
	  }
	  
	  aiter_prefix = aiter;
	}
      }
    }
  }
  
  queue_type& queue;
  ttable_type      ttable;
  ttable_pair_type ttable_prefix;
  ttable_pair_type ttable_suffix;
};

void learn(ttable_type& ttable, ttable_pair_type& ttable_prefix, ttable_pair_type& ttable_suffix)
{
  typedef TaskLearn learner_type;
  
  typedef learner_type::bitext_type     bitext_type;
  typedef learner_type::bitext_set_type bitext_set_type;
  typedef learner_type::queue_type      queue_type;
  
  typedef std::vector<learner_type, std::allocator<learner_type> > learner_set_type;
  
  ttable.clear();
  ttable_prefix.clear();
  ttable_suffix.clear();
  
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
  
  // merge!  
  for (size_t i = 0; i != learners.size(); ++ i) {
    if (ttable.empty())
      learners[i].ttable.swap(ttable);
    else
      ttable += learners[i].ttable;
    learners[i].ttable.clear();

    if (ttable_prefix.empty())
      learners[i].ttable_prefix.swap(ttable_prefix);
    else
      ttable_prefix += learners[i].ttable_prefix;
    learners[i].ttable_prefix.clear();

    if (ttable_suffix.empty())
      learners[i].ttable_suffix.swap(ttable_suffix);
    else
      ttable_suffix += learners[i].ttable_suffix;
    learners[i].ttable_suffix.clear();    
  }
  utils::resource accumulate_end;
  
  if (debug)
    std::cerr << "cpu time:  " << accumulate_end.cpu_time() - accumulate_start.cpu_time() << std::endl
	      << "user time: " << accumulate_end.user_time() - accumulate_start.user_time() << std::endl;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source",    po::value<path_type>(&source_file), "source file")
    ("target",    po::value<path_type>(&target_file), "target file")
    ("alignment", po::value<path_type>(&alignment_file), "alignment file")
    
    ("output",      po::value<path_type>(&output_file),          "output for source, target")
    ("output-prefix", po::value<path_type>(&output_prefix_file), "output for source-1, source, target")
    ("output-suffix", po::value<path_type>(&output_suffix_file), "output for source, source+1, target")
    
    ("kbest", po::value<int>(&kbest), "kbest output")
    
    ("inverse", po::bool_switch(&inverse_mode), "inverse alignment")

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
