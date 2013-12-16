//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_ALIGNMENT_IMPL__HPP__
#define __CICADA_ALIGNMENT_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <numeric>
#include <limits>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/filesystem.hpp>
#include <boost/range.hpp>

#include <cicada/sentence.hpp>
#include <cicada/alignment.hpp>
#include <cicada/dependency.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/bithack.hpp>
#include <utils/compact_map.hpp>
#include <utils/compact_set.hpp>
#include <utils/alloc_vector.hpp>
#include <utils/bithack.hpp>
#include <utils/unordered_map.hpp>
#include <utils/hashmurmur3.hpp>
#include <utils/mathop.hpp>
#include <utils/simple_vector.hpp>
#include <utils/spinlock.hpp>
#include <utils/vector2.hpp>
#include <utils/chart.hpp>

typedef cicada::Symbol     word_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Alignment  alignment_type;
typedef cicada::Dependency dependency_type;
typedef cicada::SpanVector span_set_type;
typedef cicada::Vocab      vocab_type;
typedef boost::filesystem::path path_type;

typedef double count_type;
typedef double prob_type;
typedef double logprob_type;

struct log_likelihood_type
{
  typedef double   logprob_type;
  typedef uint64_t count_type;
  
  log_likelihood_type() : average_(0), count_(0) {}
  log_likelihood_type(const logprob_type& x) : average_(x), count_(1) {}
  
  log_likelihood_type& operator+=(const logprob_type& x)
  {
    average_ += (x - average_) / (++ count_);
    return *this;
  }
  
  log_likelihood_type& operator+=(const log_likelihood_type& x)
  {
    const count_type total = count_ + x.count_;
    
    average_ = average_ * (logprob_type(count_) / total) + x.average_ * (logprob_type(x.count_) / total);
    count_ = total;
    
    return *this;
  }
  
  operator const logprob_type&() const { return average_; }
  
  logprob_type average_;
  count_type   count_;
};

struct classes_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef std::vector<word_type, std::allocator<word_type> > map_type;
  
  classes_type(const bool __unknown=false, const bool __surface=false) : classes(), unknown(__unknown), surface(__surface) {}
  
  void clear() { classes.clear(); }
  void shrink() { map_type(classes).swap(classes); }
  
  word_type operator[](const word_type& word) const
  {
    if (unknown)
      return (word == vocab_type::BOS || word == vocab_type::EOS || word == vocab_type::EPSILON ? word : vocab_type::UNK);
    else if (surface)
      return word;
    else
      return (word.id() >= classes.size() ? vocab_type::UNK : classes[word.id()]);
  }
  
  word_type& operator[](const word_type& word)
  {
    if (word.id() >= classes.size())
      classes.resize(word.id() + 1, vocab_type::UNK);
    
    return classes[word.id()];
  }
  
  map_type classes;
  
  bool unknown;
  bool surface;
};

struct atable_counts_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef int       index_type;

  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };


  class difference_map_type
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef int       index_type;
    
    typedef utils::simple_vector<count_type, std::allocator<count_type> > difference_set_type;

    difference_map_type() : counts(), offset(0) {}
    
    difference_map_type& operator+=(const difference_map_type& x)
    {
      reserve(x.min(), x.max());
      
      for (index_type i = x.min(); i <= x.max(); ++ i)
	operator[](i) += x[i];
      return *this;
    }
    
    count_type& operator[](const index_type& diff)
    {
      const index_type pos = diff + offset;
      
      if (pos < 0) {
	difference_set_type counts_new(counts.size() - pos, 0.0);
	std::copy(counts.begin(), counts.end(), counts_new.begin() - pos);
	
	counts.swap(counts_new);
	offset -= pos;
	
	return counts[0];
      } else {
	if (pos >= static_cast<index_type>(counts.size()))
	  counts.resize(pos + 1, 0.0);
	
	return counts[pos];
      }
    }
    
    count_type operator[](const index_type& diff) const
    {
      const index_type pos = diff + offset;
      
      return (pos < 0 || pos >= static_cast<int>(counts.size()) ? 0.0 : counts[pos]);
    }
    
    void reserve(const index_type& min, const index_type& max)
    {
      operator[](min);
      operator[](max);
    }
    
    void initialize()
    {
      std::fill(counts.begin(), counts.end(), 0.0);
    }
    
    bool empty() const { return counts.empty(); }
    
    index_type min() const { return - offset; }
    index_type max() const { return static_cast<index_type>(counts.size()) - 1 - offset; }
    
    void shrink()
    {
      
    }
    
    void swap(difference_map_type& x)
    {
      counts.swap(x.counts);
      std::swap(offset, x.offset);
    }
    
    difference_set_type counts;
    index_type offset;
  };
  
  typedef std::pair<word_type, word_type> class_pair_type;
  
  typedef utils::unordered_map<class_pair_type, difference_map_type, utils::hashmurmur3<size_t>, std::equal_to<class_pair_type>,
			       std::allocator<std::pair<const class_pair_type, difference_map_type> > >::type counts_type;

  typedef difference_map_type mapped_type;
  
  typedef counts_type::const_iterator const_iterator;
  typedef counts_type::iterator       iterator;

  inline const_iterator begin() const { return counts.begin(); }
  inline       iterator begin()       { return counts.begin(); }

  inline const_iterator end() const { return counts.end(); }
  inline       iterator end()       { return counts.end(); }

  inline const_iterator find(const class_pair_type& x) const { return counts.find(x); }
  inline       iterator find(const class_pair_type& x)       { return counts.find(x); }
  
  difference_map_type& operator[](const class_pair_type& x)
  {
    return counts[x];
  }
  
  difference_map_type& operator()(const word_type& source, const word_type& target)
  {
    return counts[class_pair_type(source, target)];
  }
  
  count_type& operator()(const word_type& source, const word_type& target, const index_type& diff)
  {
    return counts[class_pair_type(source, target)][diff];
  }

  bool empty() const { return counts.empty(); }
  
  void clear() { counts.clear(); }
  
  void swap(atable_counts_type& x) { counts.swap(x.counts); }
  
  void estimate_unk()
  {
    if (counts.empty()) return;

    difference_map_type counts_source_target;
    counts_type counts_source;
    counts_type counts_target;
    
    counts_type::const_iterator aiter_end = counts.end();
    for (counts_type::const_iterator aiter = counts.begin(); aiter != aiter_end; ++ aiter) {
      const class_pair_type& pair = aiter->first;
      
      if (pair.first != vocab_type::BOS
	  && pair.first != vocab_type::EOS
	  && pair.first != vocab_type::EPSILON
	  && pair.first != vocab_type::NONE)
	counts_source[class_pair_type(vocab_type::UNK, pair.second)] += aiter->second;
      
      if (pair.second != vocab_type::BOS
	  && pair.second != vocab_type::EOS
	  && pair.second != vocab_type::EPSILON
	  && pair.second != vocab_type::NONE)
	counts_target[class_pair_type(pair.first, vocab_type::UNK)] += aiter->second;
      
      if (pair.first != vocab_type::BOS
	  && pair.first != vocab_type::EOS
	  && pair.first != vocab_type::EPSILON
	  && pair.first != vocab_type::NONE
	  && pair.second != vocab_type::BOS
	  && pair.second != vocab_type::EOS
	  && pair.second != vocab_type::EPSILON
	  && pair.second != vocab_type::NONE)
	counts_source_target += aiter->second;
    }
    
    counts_type::const_iterator siter_end = counts_source.end();
    for (counts_type::const_iterator siter = counts_source.begin(); siter != siter_end; ++ siter)
      counts[siter->first] = siter->second;
    
    counts_type::const_iterator titer_end = counts_target.end();
    for (counts_type::const_iterator titer = counts_target.begin(); titer != titer_end; ++ titer)
      counts[titer->first] = titer->second;
    
    counts[class_pair_type(vocab_type::UNK, vocab_type::UNK)] = counts_source_target;
  }

  void initialize()
  {
    counts_type::iterator aiter_end = counts.end();
    for (counts_type::iterator aiter = counts.begin(); aiter != aiter_end; ++ aiter)
      aiter->second.initialize();
  }
  
  atable_counts_type& operator+=(const atable_counts_type& x)
  {
    counts_type::const_iterator aiter_end = x.counts.end();
    for (counts_type::const_iterator aiter = x.counts.begin(); aiter != aiter_end; ++ aiter)
      counts[aiter->first] += aiter->second;
    
    return *this;
  }
  
private:
  counts_type counts;
};

struct atable_type
{
  typedef atable_counts_type::size_type       size_type;
  typedef atable_counts_type::difference_type difference_type;
  typedef atable_counts_type::index_type      index_type;
  
  typedef atable_counts_type::difference_map_type difference_map_type;
  
  typedef atable_counts_type::class_pair_type class_pair_type;
  
  typedef std::pair<index_type, index_type> range_type;
  
  struct cache_type
  {
    typedef utils::spinlock            spinlock_type;
    typedef spinlock_type::scoped_lock lock_type;
    
    struct counts_type
    {
      typedef utils::spinlock            spinlock_type;
      typedef spinlock_type::scoped_lock lock_type;
      typedef utils::simple_vector<difference_map_type, std::allocator<difference_map_type> > count_set_type;
      
      counts_type() : counts() {}
      counts_type(const counts_type& x) : counts(x.counts) {}
      counts_type& operator=(const counts_type& x)
      {
	counts = x.counts;
	return *this;
      }
      
      void clear() { counts.clear(); }
      void resize(size_type x) { counts.resize(x); }
      bool empty() const { return counts.empty(); }
      
      count_set_type counts;
      spinlock_type  mutex;
    };

    typedef counts_type value_type;
    
    typedef utils::array_power2<counts_type, 64, std::allocator<counts_type> > counts_static_type;
    typedef std::deque<counts_type, std::allocator<counts_type> > counts_mutable_type;
    
    cache_type() { clear(); }
    cache_type(const cache_type& x): counts_static(x.counts_static), counts_mutable(x.counts_mutable) {}
    cache_type& operator=(const cache_type& x)
    {
      counts_static  = x.counts_static;
      counts_mutable = x.counts_mutable;
      return *this;
    }
    
    void clear()
    {
      counts_static.clear();
      counts_mutable.clear();
      
      for (size_type i = 0; i != counts_static.size(); ++ i)
	counts_static[i].resize(i + 2);
    }
    
    counts_type& operator[](const range_type& range) const
    {
      const size_type length = range.second - range.first;
      
      if (length >= counts_static.size()) {
	const size_type pos = length - counts_static.size();
	
	lock_type lock(const_cast<spinlock_type&>(mutex));
	
	if (pos >= counts_mutable.size())
	  const_cast<counts_mutable_type&>(counts_mutable).resize(pos + 1);
	
	if (counts_mutable[pos].empty())
	  const_cast<counts_type&>(counts_mutable[pos]).resize(length + 2);
	
	return const_cast<counts_type&>(counts_mutable[pos]);
      } else
	return const_cast<counts_type&>(counts_static[length]);
    }
    
    counts_static_type  counts_static;
    counts_mutable_type counts_mutable;
    spinlock_type       mutex;
  };
  
  typedef utils::unordered_map<class_pair_type, cache_type, utils::hashmurmur3<size_t>, std::equal_to<class_pair_type>,
			       std::allocator<std::pair<const class_pair_type, cache_type> > >::type cache_set_type;
  
  
  atable_type(const double __prior=0.1, const double __smooth=1e-20)
    : table(), prior(__prior), smooth(__smooth) { initialize_cache(); }
  
  atable_type(const atable_type& x)
    : table(x.table), prior(x.prior), smooth(x.smooth) { initialize_cache(); }
  atable_type& operator=(const atable_type& x)
  {
    clear();
    
    table = x.table;
    prior  = x.prior;
    smooth = x.smooth;

    initialize_cache();
    
    return *this;
  }
  
  prob_type operator()(const word_type& source,
		       const word_type& target,
		       const index_type& source_size,
		       const index_type& target_size,
		       const index_type& i_prev,
		       const index_type& i) const
  {
    if (table.empty()) return 1.0 / source_size;
    
    // i_prev < 0 implies BOS
    // i >= souce_size implies EOS
    //
    // we will cache wrt class_pair_type and diff's range
    //
    
    if (source == vocab_type::BOS) {
      // 0 <= i < source_size
      // which implies: 1 <= diff < source_size + 1
      
      return estimate(class_pair_type(source, target), range_type(1, source_size + 1))[i - i_prev];
    } else if (target == vocab_type::EOS) {
      // which implies: 1 <= diff < source_size - i_prev + 1
      
      return estimate(class_pair_type(source, target), range_type(1, source_size - i_prev + 1))[i - i_prev];
    } else {
      // 0 <= i < source_size
      // which implies: 0 - i_prev <= diff < source_size - i_prev
      
      return estimate(class_pair_type(source, target), range_type(0 - i_prev, source_size - i_prev))[i - i_prev];
    }
  }
  
  const difference_map_type& estimate(const class_pair_type& classes, const range_type& range) const
  {
    //
    // range-second is always positive, the range-second can range from 0 to (range.second - range.first) + 2, including BOS/EOS
    //

    cache_type* cache = &const_cast<cache_type&>(cache_unk);
    cache_set_type::const_iterator citer = caches.find(classes);
    if (citer != caches.end())
      cache = &const_cast<cache_type&>(citer->second);
    
    cache_type::value_type& value = cache->operator[](range);
    
    cache_type::value_type::lock_type lock(value.mutex);
    
    difference_map_type& diffs = value.counts[range.second];
    if (diffs.empty()) {
      diffs.reserve(range.first, range.second - 1);
      
      double sum = 0.0;
      
      atable_counts_type::const_iterator aiter = table.find(classes);
      
      for (index_type i = range.first; i != range.second; ++ i) {
	const double count = (aiter != table.end() ? aiter->second[i] + prior : prior);
	
	diffs[i] = count;
	sum += count;
      }
      
      const double sum_digamma = utils::mathop::digamma(sum);
      for (index_type i = range.first; i != range.second; ++ i)
	diffs[i] = std::max(utils::mathop::exp(utils::mathop::digamma(diffs[i]) - sum_digamma), smooth);
    }
    
    return diffs;
  }
  
  
  difference_map_type& operator[](const class_pair_type& x)
  {
    return table[x];
  }
  
  difference_map_type& operator()(const word_type& source, const word_type& target)
  {
    return table[class_pair_type(source, target)];
  }

  count_type& operator()(const word_type& source, const word_type& target, const index_type& diff)
  {
    return table[class_pair_type(source, target)][diff];
  }
  
  void clear()
  {
    table.clear();

    initialize_cache();
  }
  void swap(atable_type& x)
  {
    table.swap(x.table);
    std::swap(prior,  x.prior);
    std::swap(smooth, x.smooth);
    
    initialize_cache();
    x.initialize_cache();
  }
  
  void estimate_unk()
  {
    table.estimate_unk();
  }
  
  void shrink() {}
  
  void initialize()
  {
    table.initialize();
    
    initialize_cache();
  }

  void initialize_cache()
  {
    caches.clear();
    atable_counts_type::const_iterator aiter_end = table.end();
    for (atable_counts_type::const_iterator aiter = table.begin(); aiter != aiter_end; ++ aiter) 
      caches[aiter->first].clear();
    
    cache_unk.clear();
  };
  
  atable_type& operator+=(const atable_type& x)
  {
    table += x.table;

    return *this;
  }

  atable_type& operator+=(const atable_counts_type& x)
  {
    table += x;
    
    return *this;
  }
  
  bool empty() const { return table.empty(); }
  
  atable_counts_type table;
  double prior;
  double smooth;
  
  // caching....
  cache_set_type caches;
  cache_type     cache_unk;
};

struct ttable_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };

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

    inline const_iterator find(const key_type& x) const { return counts.find(x); }
    inline       iterator find(const key_type& x)       { return counts.find(x); }
    
    mapped_type& operator[](const key_type& key) { return counts[key]; }
    
    size_type size() const { return counts.size(); }
    bool empty() const { return counts.empty(); }

    void swap(count_map_type& x) { counts.swap(x.counts); }
    void clear()
    {
      counts.clear();
    }

    void rehash(size_t hint) { counts.rehash(hint); }

    count_map_type& operator+=(const count_map_type& x)
    {
      counts.rehash(counts.size() + x.counts.size());

      const_iterator citer_end = x.counts.end();
      for (const_iterator citer = x.counts.begin(); citer != citer_end; ++ citer)
	counts[citer->first] += citer->second;
      return *this;
    }
    
    count_map_type& operator|=(const count_map_type& x)
    {
      counts.rehash(counts.size() + x.counts.size());
      
      counts.insert(x.begin(), x.end());
      return *this;
    }
    
    counts_type counts;
  };
  
  typedef utils::alloc_vector<count_map_type, std::allocator<count_map_type> > count_dict_type;
  
  ttable_type(const double __prior=0.1, const double __smooth=1e-20) : ttable(), prior(__prior), smooth(__smooth) {}
  
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
    if (source == vocab_type::BOS || source == vocab_type::EOS || target == vocab_type::BOS || target == vocab_type::EOS)
      return source == target;
    
    if (! ttable.exists(source.id())) return smooth;
    
    const count_map_type& counts = ttable[source.id()];
    count_map_type::const_iterator citer = counts.find(target);
    
    return (citer == counts.end() ? smooth : std::max(citer->second, smooth));
  }
  
  
  void shrink() { ttable.shrink(); }
  void clear() { ttable.clear(); }
  
  void clear(const word_type& word)
  {
    if (! ttable.exists(word.id())) return;
    
    ttable.erase(ttable.begin() + word.id());
  }

  void swap(ttable_type& x)
  {
    ttable.swap(x.ttable);
    std::swap(prior,  x.prior);
    std::swap(smooth, x.smooth);
  }

  size_type size() const { return ttable.size(); }
  bool empty() const { return ttable.empty(); }
  bool exists(size_type pos) const { return ttable.exists(pos); }
  
  void resize(size_type __size) { ttable.resize(__size); }
  void reserve(size_type __size) { ttable.reserve(__size); }

  void initialize()
  {
    for (size_type i = 0; i != ttable.size(); ++ i)
      if (ttable.exists(i))
        ttable[i].clear();

    shrink();
  }

  ttable_type& operator+=(const ttable_type& x)
  {
    for (size_type i = 0; i != x.ttable.size(); ++ i) 
      if (x.ttable.exists(i))
	ttable[i] += x.ttable[i];
    
    return *this;
  }

  ttable_type& operator|=(const ttable_type& x)
  {
    for (size_type i = 0; i != x.ttable.size(); ++ i) 
      if (x.ttable.exists(i))
	ttable[i] |= x.ttable[i];
    
    return *this;
  }
  
  count_dict_type ttable;
  double prior;
  double smooth;
};

struct aligned_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  struct aligned_map_type
  {
    typedef utils::compact_set<word_type,
			       utils::unassigned<word_type>, utils::unassigned<word_type>,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<word_type> > map_type;

    typedef map_type::value_type      value_type;
    typedef map_type::size_type       size_type;
    typedef map_type::difference_type difference_type;
    
    typedef map_type::const_iterator const_iterator;
    typedef map_type::iterator       iterator;
  
    typedef map_type::const_reference const_reference;
    typedef map_type::reference       reference;

    aligned_map_type() { }
    
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
    void swap(aligned_map_type& x)
    {
      aligned.swap(x.aligned);
    }

    aligned_map_type& operator=(const aligned_map_type& x)
    {
      aligned = x.aligned;
      return *this;
    }
    
    aligned_map_type& operator+=(const aligned_map_type& x)
    {
      aligned.rehash(aligned.size() + x.aligned.size());
      
      const_iterator citer_end = x.aligned.end();
      for (const_iterator citer = x.aligned.begin(); citer != citer_end; ++ citer)
	aligned.insert(*citer);
      
      return *this;
    }
    
    map_type aligned;
  };
  
  typedef utils::alloc_vector<aligned_map_type, std::allocator<aligned_map_type> > aligned_set_type;

  aligned_type() : aligned() {}
  
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
  void reserve(size_type __size) { aligned.reserve(__size); }

  void swap(aligned_type& x) { aligned.swap(x.aligned); }
  
  void shrink() { aligned.shrink(); }
  void clear() { aligned.clear(); }
  void clear(const word_type& word)
  {
    if (! aligned.exists(word.id())) return;
    
    aligned.erase(aligned.begin() + word.id());
  }
  
  void initialize()
  {
    for (size_type i = 0; i != aligned.size(); ++ i)
      if (aligned.exists(i))
        aligned[i].clear();
    
    shrink();
  }
  
  aligned_set_type aligned;
};

//
// e-table for length modeling
//
struct etable_type
{
  typedef utils::vector2<double, std::allocator<double> > table_type;

  static const int max_length = 128;

  etable_type(const double __lambda=1.0) : lambda(__lambda) { assign(lambda); }
  
  double operator()(const int m, const int l) const
  {
    if (m < max_length && l < max_length)
      return table(m, l);
    else
      return utils::mathop::pow(lambda * l, double(m)) * utils::mathop::exp(- lambda * l) / utils::mathop::factorial<double>(m);
  }

  void assign(const double __lambda)
  {
    lambda = __lambda;
    
    table.clear();
    table.reserve(max_length, max_length);
    table.resize(max_length, max_length, 0.0);
    
    for (int m = 1; m != max_length; ++ m)
      for (int l = 1; l != max_length; ++ l)
	table(m, l) = utils::mathop::pow(lambda * l, double(m)) * utils::mathop::exp(- lambda * l) / utils::mathop::factorial<double>(m);
  }
  
  table_type table;
  double lambda;
};

//
// p-table for NULL alignment
//
struct ptable_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  //typedef utils::chart<double, std::allocator<double> > table_type;
  typedef utils::vector2<double, std::allocator<double> > table_type;

  static const int max_length = 128;
  
  ptable_type(const double __p1=0.01) : p1(__p1) { initialize(); }
  
  double operator()(const int m, const int phi0) const
  {
    if (m < max_length && phi0 < max_length)
      return table(phi0, m);
    else
      return (std::pow(1.0 - p1, utils::bithack::max(m - phi0 * 2, 0))
	      * std::pow(p1, phi0)
	      * binomial(m - phi0, phi0));
  }
  
  void initialize()
  {
    table.clear();
    table.reserve(max_length, max_length);
    table.resize(max_length, max_length, 0.0);
    
    for (int m = 1; m != max_length; ++ m)
      for (int phi0 = 0; phi0 <= m; ++ phi0)
	table(phi0, m) = (std::pow(1.0 - p1, utils::bithack::max(m - phi0 * 2, 0))
			  * std::pow(p1, phi0)
			  * binomial(m - phi0, phi0));
  }

  double binomial(const int n, const int k) const
  {
    double result = 1.0;
    if (k < n)
      for (int i = 1; i <= k; ++ i)
	result *= double(n - i + 1) / double(i);
    return result;
  }
  
  table_type table;
  double p1;
};

// n-model table...
struct ntable_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef uint32_t  index_type;
  
  typedef utils::vector2<count_type, std::allocator<count_type> > table_type;
  typedef std::vector<index_type, std::allocator<index_type> > map_type;
  
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };

  typedef table_type::const_iterator const_iterator;
  typedef table_type::iterator       iterator;
  
  static const size_type fertility_size = 16;
  
  ntable_type(const double __prior=0.1, const double __smooth=1e-20)
    : map(), table(0, fertility_size), prior(__prior), smooth(__smooth) {}
  
  count_type& operator()(const word_type& word, const size_type fertility)
  {
    if (word.id() >= map.size())
      map.resize(word.id() + 1, index_type(-1));
    
    if (map[word.id()] == index_type(-1)) {
      table.resize(table.size1() + 1, fertility_size, count_type(0));
      map[word.id()] = table.size1() - 1;
    }
    
    return table(map[word.id()], utils::bithack::min(fertility, fertility_size - 1));
  }
  
  count_type operator()(word_type word, const size_type length, const size_type fertility) const
  {
    if (word.id() >= map.size() || map[word.id()] == index_type(-1))
      word = vocab_type::UNK;
    
    // uniform... 0 <= fertility <= length
    if (word.id() >= map.size() || map[word.id()] == index_type(-1))
      return 1.0 / (length + 1);
    
    // the last value (15, or fertility_size - 1) is reserved for further smoothing...
    if (fertility >= fertility_size - 1)
      return table(map[word.id()], fertility_size - 1) / (length  + 1 - (fertility_size - 1));
    else
      return table(map[word.id()], fertility);
  }

  void assign(const word_type& word)
  {
    if (word.id() >= map.size())
      map.resize(word.id() + 1, index_type(-1));
    
    if (map[word.id()] == index_type(-1)) {
      table.resize(table.size1() + 1, fertility_size, count_type(0));
      map[word.id()] = table.size1() - 1;
    }
  }

  void estimate()
  {
    // first, estimate for UNK...
    std::vector<count_type, std::allocator<count_type> > counts(fertility_size, 0.0);
    
    for (size_type i = 0; i != table.size1(); ++ i)
      std::transform(table.begin(i), table.end(i), counts.begin(), counts.begin(), std::plus<count_type>());
    
    assign(vocab_type::UNK);
    std::copy(counts.begin(), counts.end(), begin(vocab_type::UNK));
    
    // then, perform estimation...
    for (size_type i = 0; i != table.size1(); ++ i) {
      const double sum = std::accumulate(table.begin(i), table.end(i), prior * fertility_size);
      const double sum_digamma = utils::mathop::digamma(sum);
      
      table_type::iterator iter_end = table.end(i);
      for (table_type::iterator iter = table.begin(i); iter != iter_end; ++ iter)
	*iter = std::max(utils::mathop::exp(utils::mathop::digamma(*iter + prior) - sum_digamma), smooth);
    }
  }

  inline const_iterator begin(const word_type& x) const { return table.begin(map[x.id()]); }
  inline       iterator begin(const word_type& x)       { return table.begin(map[x.id()]); }

  inline const_iterator end(const word_type& x) const { return table.end(map[x.id()]); }
  inline       iterator end(const word_type& x)       { return table.end(map[x.id()]); }
  
  void swap(ntable_type& x)
  {
    map.swap(x.map);
    table.swap(x.table);
    std::swap(prior,  x.prior);
    std::swap(smooth, x.smooth);
  }

  void clear()
  {
    map.clear();
    table.clear();
  }
  
  void initialize()
  {
    map.clear();
    table.clear();
  }
  
  ntable_type& operator+=(const ntable_type& x)
  {
    for (word_type::id_type id = 0; id != x.map.size(); ++ id)
      if (x.map[id] != index_type(-1)) {
	const word_type word(id);
	
	if (word.id() >= map.size())
	  map.resize(word.id() + 1, index_type(-1));
	
	if (map[word.id()] == index_type(-1)) {
	  table.resize(table.size1() + 1, fertility_size, count_type(0));
	  map[word.id()] = table.size1() - 1;
	}
	
	std::transform(x.table.begin(x.map[word.id()]), x.table.end(x.map[word.id()]),
		       table.begin(map[word.id()]),
		       table.begin(map[word.id()]),
		       std::plus<count_type>());
      }
    
    return *this;
  }
  
  map_type   map;
  table_type table;
  
  double prior;
  double smooth;
};


// d-table type...
// it is almost identical to atable_type..
typedef atable_counts_type dtable_counts_type;

struct dtable_type
{
  typedef dtable_counts_type::size_type       size_type;
  typedef dtable_counts_type::difference_type difference_type;
  typedef dtable_counts_type::index_type      index_type;
  
  typedef dtable_counts_type::difference_map_type difference_map_type;
  
  typedef dtable_counts_type::class_pair_type class_pair_type;
  
  typedef std::pair<index_type, index_type> range_type;

  
  struct cache_type
  {
    typedef utils::spinlock            spinlock_type;
    typedef spinlock_type::scoped_lock lock_type;
    
    struct counts_type
    {
      typedef utils::spinlock            spinlock_type;
      typedef spinlock_type::scoped_lock lock_type;
      typedef utils::simple_vector<difference_map_type, std::allocator<difference_map_type> > count_set_type;
      
      counts_type() : counts() {}
      counts_type(const counts_type& x) : counts(x.counts) {}
      counts_type& operator=(const counts_type& x)
      {
	counts = x.counts;
	return *this;
      }
      
      void clear() { counts.clear(); }
      void resize(size_type x) { counts.resize(x); }
      bool empty() const { return counts.empty(); }
      
      count_set_type counts;
      spinlock_type  mutex;
    };

    typedef counts_type value_type;
    
    typedef utils::array_power2<counts_type, 64, std::allocator<counts_type> > counts_static_type;
    typedef std::deque<counts_type, std::allocator<counts_type> > counts_mutable_type;
    
    cache_type() { clear(); }
    cache_type(const cache_type& x): counts_static(x.counts_static), counts_mutable(x.counts_mutable) {}
    cache_type& operator=(const cache_type& x)
    {
      counts_static  = x.counts_static;
      counts_mutable = x.counts_mutable;
      return *this;
    }
    
    void clear()
    {
      counts_static.clear();
      counts_mutable.clear();
      
      for (size_type i = 0; i != counts_static.size(); ++ i)
	counts_static[i].resize(i + 2);
    }
    
    counts_type& operator[](const range_type& range) const
    {
      const size_type length = range.second - range.first;
      
      if (length >= counts_static.size()) {
	const size_type pos = length - counts_static.size();
	
	lock_type lock(const_cast<spinlock_type&>(mutex));
	
	if (pos >= counts_mutable.size())
	  const_cast<counts_mutable_type&>(counts_mutable).resize(pos + 1);
	
	if (counts_mutable[pos].empty())
	  const_cast<counts_type&>(counts_mutable[pos]).resize(length + 2);
	
	return const_cast<counts_type&>(counts_mutable[pos]);
      } else
	return const_cast<counts_type&>(counts_static[length]);
    }
    
    counts_static_type  counts_static;
    counts_mutable_type counts_mutable;
    spinlock_type       mutex;
  };
  
  typedef utils::unordered_map<class_pair_type, cache_type, utils::hashmurmur3<size_t>, std::equal_to<class_pair_type>,
			       std::allocator<std::pair<const class_pair_type, cache_type> > >::type cache_set_type;
  
  
  dtable_type(const double __prior=0.1, const double __smooth=1e-20)
    : table(), prior(__prior), smooth(__smooth) { initialize_cache(); }
  
  dtable_type(const dtable_type& x)
    : table(x.table), prior(x.prior), smooth(x.smooth) { initialize_cache(); }
  dtable_type& operator=(const dtable_type& x)
  {
    clear();
    
    table = x.table;
    prior  = x.prior;
    smooth = x.smooth;

    initialize_cache();
    
    return *this;
  }
  
  // head....
  prob_type operator()(const word_type& source,
		       const word_type& target,
		       const index_type& source_size,
		       const index_type& target_size,
		       const index_type& j_prev,
		       const index_type& j) const
  {
    if (table.empty()) return 1.0 / target_size;
    
    // 1 <= j <= target_size
    
    return estimate(class_pair_type(source, target), range_type(1 - j_prev, target_size - j_prev + 1))[j - j_prev];
  }
  
  
  // non-head...
  prob_type operator()(const word_type& target,
		       const index_type& source_size,
		       const index_type& target_size,
		       const index_type& j_prev,
		       const index_type& j) const
  {
    if (table.empty()) return 1.0 / target_size;
    
    // j_prev < j <= target_size, thus, range starts from 1!
    
    return estimate(class_pair_type(vocab_type::NONE, target), range_type(1, target_size - j_prev + 1))[j - j_prev];
  }
  
  const difference_map_type& estimate(const class_pair_type& classes, const range_type& range) const
  {
    //
    // range-second is always positive, the range-second can range from 0 to (range.second - range.first) + 2, including BOS/EOS
    //
    
    cache_type* cache = &const_cast<cache_type&>(classes.first == vocab_type::NONE
						 ? cache_none_unk
						 : cache_unk);
    cache_set_type::const_iterator citer = caches.find(classes);
    if (citer != caches.end())
      cache = &const_cast<cache_type&>(citer->second);
    
    cache_type::value_type& value = cache->operator[](range);
    
    cache_type::value_type::lock_type lock(value.mutex);
    
    difference_map_type& diffs = value.counts[range.second];
    if (diffs.empty()) {
      diffs.reserve(range.first, range.second - 1);
      
      double sum = 0.0;
      
      dtable_counts_type::const_iterator aiter = table.find(classes);
      
      for (index_type i = range.first; i != range.second; ++ i) {
	const double count = (aiter != table.end() ? aiter->second[i] + prior : prior);
	
	diffs[i] = count;
	sum += count;
      }
      
      const double sum_digamma = utils::mathop::digamma(sum);
      for (index_type i = range.first; i != range.second; ++ i)
	diffs[i] = std::max(utils::mathop::exp(utils::mathop::digamma(diffs[i]) - sum_digamma), smooth);
    }
    
    return diffs;
  }
  
  
  difference_map_type& operator[](const class_pair_type& x)
  {
    return table[x];
  }
  
  difference_map_type& operator()(const word_type& source, const word_type& target)
  {
    return table[class_pair_type(source, target)];
  }

  count_type& operator()(const word_type& source, const word_type& target, const index_type& diff)
  {
    return table[class_pair_type(source, target)][diff];
  }
  
  void clear()
  {
    table.clear();

    initialize_cache();
  }
  void swap(dtable_type& x)
  {
    table.swap(x.table);
    std::swap(prior,  x.prior);
    std::swap(smooth, x.smooth);
    
    initialize_cache();
    x.initialize_cache();
  }
  
  void estimate_unk()
  {
    table.estimate_unk();
  }
  
  void shrink() {}
  
  void initialize()
  {
    table.initialize();
    
    initialize_cache();
  }

  void initialize_cache()
  {
    caches.clear();
    dtable_counts_type::const_iterator aiter_end = table.end();
    for (dtable_counts_type::const_iterator aiter = table.begin(); aiter != aiter_end; ++ aiter) 
      caches[aiter->first].clear();
    
    cache_unk.clear();
    cache_none_unk.clear();
  };
  
  dtable_type& operator+=(const dtable_type& x)
  {
    table += x.table;

    return *this;
  }

  dtable_type& operator+=(const dtable_counts_type& x)
  {
    table += x;
    
    return *this;
  }
  
  bool empty() const { return table.empty(); }
  
  dtable_counts_type table;
  double prior;
  double smooth;
  
  // caching....
  cache_set_type caches;
  cache_type     cache_unk;
  cache_type     cache_none_unk;
};

struct LearnBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef int       index_type;

  static const classes_type& __classes()
  {
    static classes_type __tmp;
    return __tmp;
  }

  static const atable_type& __atable()
  {
    static atable_type __tmp;
    return __tmp;
  }

  static const dtable_type& __dtable()
  {
    static dtable_type __tmp;
    return __tmp;
  }

  static const ntable_type& __ntable()
  {
    static ntable_type __tmp;
    return __tmp;
  }

  static const ptable_type& __ptable()
  {
    static ptable_type __tmp;
    return __tmp;
  }
  
  LearnBase(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      ttable_counts_source_target(__ttable_source_target.prior, __ttable_source_target.smooth),
      ttable_counts_target_source(__ttable_target_source.prior, __ttable_target_source.smooth),
      atable_source_target(__atable()),
      atable_target_source(__atable()),
      dtable_source_target(__dtable()),
      dtable_target_source(__dtable()),
      ntable_source_target(__ntable()),
      ntable_target_source(__ntable()),
      ptable_source_target(__ptable()),
      ptable_target_source(__ptable()),
      classes_source(__classes()),
      classes_target(__classes()),
      objective_source_target(0),
      objective_target_source(0)
  {}
  
  LearnBase(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source,
	    const atable_type& __atable_source_target,
	    const atable_type& __atable_target_source,
	    const classes_type& __classes_source,
	    const classes_type& __classes_target)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      ttable_counts_source_target(__ttable_source_target.prior, __ttable_source_target.smooth),
      ttable_counts_target_source(__ttable_target_source.prior, __ttable_target_source.smooth),
      atable_source_target(__atable_source_target),
      atable_target_source(__atable_target_source),
      dtable_source_target(__dtable()),
      dtable_target_source(__dtable()),
      ntable_source_target(__ntable()),
      ntable_target_source(__ntable()),
      ptable_source_target(__ptable()),
      ptable_target_source(__ptable()),
      classes_source(__classes_source),
      classes_target(__classes_target),
      objective_source_target(0),
      objective_target_source(0)
  {}


  LearnBase(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source,
	    const atable_type& __atable_source_target,
	    const atable_type& __atable_target_source,
	    const dtable_type& __dtable_source_target,
	    const dtable_type& __dtable_target_source,
	    const ntable_type& __ntable_source_target,
	    const ntable_type& __ntable_target_source,
	    const ptable_type& __ptable_source_target,
	    const ptable_type& __ptable_target_source,
	    const classes_type& __classes_source,
	    const classes_type& __classes_target)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      ttable_counts_source_target(__ttable_source_target.prior, __ttable_source_target.smooth),
      ttable_counts_target_source(__ttable_target_source.prior, __ttable_target_source.smooth),
      atable_source_target(__atable_source_target),
      atable_target_source(__atable_target_source),
      dtable_source_target(__dtable_source_target),
      dtable_target_source(__dtable_target_source),
      ntable_source_target(__ntable_source_target),
      ntable_target_source(__ntable_target_source),
      ntable_counts_source_target(__ntable_source_target.prior, __ntable_source_target.smooth),
      ntable_counts_target_source(__ntable_target_source.prior, __ntable_target_source.smooth),
      ptable_source_target(__ptable_source_target),
      ptable_target_source(__ptable_target_source),
      classes_source(__classes_source),
      classes_target(__classes_target),
      objective_source_target(0),
      objective_target_source(0)
  {}

  void initialize()
  {
    ttable_counts_source_target.clear();
    ttable_counts_target_source.clear();
    ttable_counts_source_target.reserve(word_type::allocated());
    ttable_counts_target_source.reserve(word_type::allocated());
    ttable_counts_source_target.resize(word_type::allocated());
    ttable_counts_target_source.resize(word_type::allocated());
    
    aligned_source_target.clear();
    aligned_target_source.clear();
    aligned_source_target.reserve(word_type::allocated());
    aligned_target_source.reserve(word_type::allocated());
    aligned_source_target.resize(word_type::allocated());
    aligned_target_source.resize(word_type::allocated());
    
    atable_counts_source_target.initialize();
    atable_counts_target_source.initialize();
    
    dtable_counts_source_target.initialize();
    dtable_counts_target_source.initialize();
    
    ntable_counts_source_target.initialize();
    ntable_counts_target_source.initialize();
    
    objective_source_target = 0.0;
    objective_target_source = 0.0;
  }

  void shrink()
  {
    //atable_source_target.shrink();
    //atable_target_source.shrink();
  };
  
  const ttable_type& ttable_source_target;
  const ttable_type& ttable_target_source;
  ttable_type ttable_counts_source_target;
  ttable_type ttable_counts_target_source;

  aligned_type aligned_source_target;
  aligned_type aligned_target_source;

  const atable_type& atable_source_target;
  const atable_type& atable_target_source;
  atable_counts_type atable_counts_source_target;
  atable_counts_type atable_counts_target_source;

  const dtable_type& dtable_source_target;
  const dtable_type& dtable_target_source;
  dtable_counts_type dtable_counts_source_target;
  dtable_counts_type dtable_counts_target_source;

  const ntable_type& ntable_source_target;
  const ntable_type& ntable_target_source;
  ntable_type ntable_counts_source_target;
  ntable_type ntable_counts_target_source;

  const ptable_type& ptable_source_target;
  const ptable_type& ptable_target_source;

  const classes_type& classes_source;
  const classes_type& classes_target;
    
  double objective_source_target;
  double objective_target_source;
};

struct ViterbiBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef int       index_type;

  static const classes_type& __classes()
  {
    static classes_type __tmp;
    return __tmp;
  }
  
  static const atable_type& __atable()
  {
    static atable_type __tmp;
    return __tmp;
  }

  static const dtable_type& __dtable()
  {
    static dtable_type __tmp;
    return __tmp;
  }

  static const ntable_type& __ntable()
  {
    static ntable_type __tmp;
    return __tmp;
  }

  static const ptable_type& __ptable()
  {
    static ptable_type __tmp;
    return __tmp;
  }


  ViterbiBase(const ttable_type& __ttable_source_target,
	      const ttable_type& __ttable_target_source)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      atable_source_target(__atable()),
      atable_target_source(__atable()),
      dtable_source_target(__dtable()),
      dtable_target_source(__dtable()),
      ntable_source_target(__ntable()),
      ntable_target_source(__ntable()),
      ptable_source_target(__ptable()),
      ptable_target_source(__ptable()),
      classes_source(__classes()),
      classes_target(__classes())
  {}
  ViterbiBase(const ttable_type& __ttable_source_target,
	      const ttable_type& __ttable_target_source,
	      const atable_type& __atable_source_target,
	      const atable_type& __atable_target_source,
	      const classes_type& __classes_source,
	      const classes_type& __classes_target)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      atable_source_target(__atable_source_target),
      atable_target_source(__atable_target_source),
      dtable_source_target(__dtable()),
      dtable_target_source(__dtable()),
      ntable_source_target(__ntable()),
      ntable_target_source(__ntable()),
      ptable_source_target(__ptable()),
      ptable_target_source(__ptable()),
      classes_source(__classes_source),
      classes_target(__classes_target)
  {}

  ViterbiBase(const ttable_type& __ttable_source_target,
	      const ttable_type& __ttable_target_source,
	      const atable_type& __atable_source_target,
	      const atable_type& __atable_target_source,
	      const dtable_type& __dtable_source_target,
	      const dtable_type& __dtable_target_source,
	      const ntable_type& __ntable_source_target,
	      const ntable_type& __ntable_target_source,
	      const ptable_type& __ptable_source_target,
	      const ptable_type& __ptable_target_source,
	      const classes_type& __classes_source,
	      const classes_type& __classes_target)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      atable_source_target(__atable_source_target),
      atable_target_source(__atable_target_source),
      dtable_source_target(__dtable_source_target),
      dtable_target_source(__dtable_target_source),
      ntable_source_target(__ntable_source_target),
      ntable_target_source(__ntable_target_source),
      ptable_source_target(__ptable_source_target),
      ptable_target_source(__ptable_target_source),
      classes_source(__classes_source),
      classes_target(__classes_target)
  {}
  
  void shrink()
  {
    //atable_source_target.shrink();
    //atable_target_source.shrink();
  };

  const ttable_type& ttable_source_target;
  const ttable_type& ttable_target_source;
  
  const atable_type& atable_source_target;
  const atable_type& atable_target_source;

  const dtable_type& dtable_source_target;
  const dtable_type& dtable_target_source;

  const ntable_type& ntable_source_target;
  const ntable_type& ntable_target_source;

  const ptable_type& ptable_source_target;
  const ptable_type& ptable_target_source;
  
  const classes_type& classes_source;
  const classes_type& classes_target;
};

void read_lexicon(const path_type& path, ttable_type& lexicon)
{
  typedef boost::fusion::tuple<std::string, std::string, double > lexicon_parsed_type;
  typedef boost::spirit::istream_iterator iterator_type;

  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  
  qi::rule<iterator_type, std::string(), standard::blank_type>         word;
  qi::rule<iterator_type, lexicon_parsed_type(), standard::blank_type> parser; 
  
  word   %= qi::lexeme[+(standard::char_ - standard::space)];
  parser %= word >> word >> qi::double_ >> (qi::eol | qi::eoi);
  
  lexicon.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iterator_type iter(is);
  iterator_type iter_end;
  
  lexicon_parsed_type lexicon_parsed;
  
  while (iter != iter_end) {
    boost::fusion::get<0>(lexicon_parsed).clear();
    boost::fusion::get<1>(lexicon_parsed).clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, lexicon_parsed))
      if (iter != iter_end)
	throw std::runtime_error("global lexicon parsing failed");
    
    const word_type target(boost::fusion::get<0>(lexicon_parsed));
    const word_type source(boost::fusion::get<1>(lexicon_parsed));
    const double&   prob(boost::fusion::get<2>(lexicon_parsed));
    
    lexicon[source][target] = prob;
  }

  lexicon.shrink();
}

struct ttable_greater_second
{
  template <typename Tp>
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

void write_lexicon(const path_type& path, const ttable_type& lexicon, const aligned_type& aligned, const double threshold)
{
  typedef ttable_type::count_map_type::value_type value_type;
  typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;
  typedef std::ostream_iterator<char> oiterator_type;
  
  namespace karma = boost::spirit::karma;
  namespace standard = boost::spirit::standard;

  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(20);
  
  oiterator_type oiter(os);
  karma::real_generator<double, ttable_type::real_precision> real;

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
      
      std::sort(sorted.begin(), sorted.end(), ttable_greater_second());
      
      if (threshold > 0.0) {
	const double prob_max       = sorted.front()->second;
	const double prob_threshold = prob_max * threshold;
	
	const aligned_type::aligned_map_type& viterbi = (aligned.exists(source) ? aligned[source] : __empty);
	
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  if ((*iter)->second >= prob_threshold || viterbi.find((*iter)->first) != viterbi.end())
	    if (! karma::generate(oiter, standard::string << ' ' << standard::string << ' ' << real << '\n',
				  (*iter)->first,
				  source,
				  (*iter)->second))
	      throw std::runtime_error("generation failed");
      } else {
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  if (! karma::generate(oiter, standard::string << ' ' << standard::string << ' ' << real << '\n',
				(*iter)->first,
				source,
				(*iter)->second))
	    throw std::runtime_error("generation failed");
      }
    }
}

void read_alignment(const path_type& path, atable_type& align)
{
  typedef boost::fusion::tuple<std::string, std::string, int, double > align_parsed_type;
  typedef boost::spirit::istream_iterator iterator_type;

  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  
  qi::rule<iterator_type, std::string(), standard::blank_type>       word;
  qi::rule<iterator_type, align_parsed_type(), standard::blank_type> parser; 
  
  word   %= qi::lexeme[+(standard::char_ - standard::space)];
  parser %= word >> word >> qi::int_ >> qi::double_ >> (qi::eol | qi::eoi);
  
  align.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iterator_type iter(is);
  iterator_type iter_end;
  
  align_parsed_type align_parsed;
  
  while (iter != iter_end) {
    boost::fusion::get<0>(align_parsed).clear();
    boost::fusion::get<1>(align_parsed).clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, align_parsed))
      if (iter != iter_end)
	throw std::runtime_error("global lexicon parsing failed");
    
    const word_type source(boost::fusion::get<0>(align_parsed));
    const word_type target(boost::fusion::get<1>(align_parsed));
    const int&      index(boost::fusion::get<2>(align_parsed));
    const double&   prob(boost::fusion::get<3>(align_parsed));
    
    align[std::make_pair(source, target)][index] = prob;
  }

  align.initialize_cache();
}


void write_alignment(const path_type& path, const atable_type& align)
{
  namespace karma = boost::spirit::karma;
  namespace standard = boost::spirit::standard;

  typedef std::ostream_iterator<char> iterator_type;

  karma::real_generator<double, atable_counts_type::real_precision> real;

  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(20);

  iterator_type iter(os);
  
  atable_counts_type::const_iterator citer_end = align.table.end();
  for (atable_counts_type::const_iterator citer = align.table.begin(); citer != citer_end; ++ citer)
    for (int diff = citer->second.min(); diff <= citer->second.max(); ++ diff)
      if (citer->second[diff] > 0.0)
	if (! karma::generate(iter, standard::string << ' ' << standard::string << ' ' << karma::int_ << ' '<< real << '\n',
			      citer->first.first,
			      citer->first.second,
			      diff,
			      std::min(std::max(citer->second[diff], std::numeric_limits<double>::min()),
				       std::numeric_limits<double>::max())))
	throw std::runtime_error("generation failed");
}

void read_distortion(const path_type& path, dtable_type& distortion)
{
  typedef boost::fusion::tuple<std::string, std::string, int, double > distortion_parsed_type;
  typedef boost::spirit::istream_iterator iterator_type;

  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  
  qi::rule<iterator_type, std::string(), standard::blank_type>       word;
  qi::rule<iterator_type, distortion_parsed_type(), standard::blank_type> parser; 
  
  word   %= qi::lexeme[+(standard::char_ - standard::space)];
  parser %= word >> word >> qi::int_ >> qi::double_ >> (qi::eol | qi::eoi);
  
  distortion.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iterator_type iter(is);
  iterator_type iter_end;
  
  distortion_parsed_type distortion_parsed;
  
  while (iter != iter_end) {
    boost::fusion::get<0>(distortion_parsed).clear();
    boost::fusion::get<1>(distortion_parsed).clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, distortion_parsed))
      if (iter != iter_end)
	throw std::runtime_error("global lexicon parsing failed");
    
    const word_type source(boost::fusion::get<0>(distortion_parsed));
    const word_type target(boost::fusion::get<1>(distortion_parsed));
    const int&      index(boost::fusion::get<2>(distortion_parsed));
    const double&   prob(boost::fusion::get<3>(distortion_parsed));
    
    distortion[std::make_pair(source, target)][index] = prob;
  }

  distortion.initialize_cache();
}

void write_distortion(const path_type& path, const dtable_type& distortion)
{
  namespace karma = boost::spirit::karma;
  namespace standard = boost::spirit::standard;

  typedef std::ostream_iterator<char> iterator_type;

  karma::real_generator<double, dtable_counts_type::real_precision> real;

  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(20);

  iterator_type iter(os);
  
  dtable_counts_type::const_iterator citer_end = distortion.table.end();
  for (dtable_counts_type::const_iterator citer = distortion.table.begin(); citer != citer_end; ++ citer)
    for (int diff = citer->second.min(); diff <= citer->second.max(); ++ diff)
      if (citer->second[diff] > 0.0)
	if (! karma::generate(iter, standard::string << ' ' << standard::string << ' ' << karma::int_ << ' '<< real << '\n',
			      citer->first.first,
			      citer->first.second,
			      diff,
			      std::min(std::max(citer->second[diff], std::numeric_limits<double>::min()),
				       std::numeric_limits<double>::max())))
	throw std::runtime_error("generation failed");
}

void read_fertility(const path_type& path, ntable_type& fertility)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  typedef boost::fusion::tuple<std::string, prob_set_type> fertility_parsed_type;
  typedef boost::spirit::istream_iterator iterator_type;
    
  qi::rule<iterator_type, std::string(), standard::blank_type>         word;
  qi::rule<iterator_type, fertility_parsed_type(), standard::blank_type> parser; 
  
  // fertility_size of ntable_type!
  word   %= qi::lexeme[+(standard::char_ - standard::space)];
  parser %= word >> qi::repeat(16)[qi::double_] >> (qi::eol | qi::eoi);
  
  fertility.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iterator_type iter(is);
  iterator_type iter_end;

  fertility_parsed_type fertility_parsed;
  
  while (iter != iter_end) {
    boost::fusion::get<0>(fertility_parsed).clear();
    boost::fusion::get<1>(fertility_parsed).clear();
    
    if (! qi::phrase_parse(iter, iter_end, parser, standard::blank, fertility_parsed))
      if (iter != iter_end)
	throw std::runtime_error("fertility parsing failed");

    const word_type word(boost::fusion::get<0>(fertility_parsed));
    
    fertility.assign(word);
    
    std::copy(boost::fusion::get<1>(fertility_parsed).begin(),
	      boost::fusion::get<1>(fertility_parsed).end(),
	      fertility.begin(word));
  }
  
}

void write_fertility(const path_type& path, const ntable_type& fertility)
{
  namespace karma = boost::spirit::karma;
  namespace standard = boost::spirit::standard;

  typedef std::ostream_iterator<char> iterator_type;
  
  karma::real_generator<double, ntable_type::real_precision> real;
  
  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(20);

  iterator_type iter(os);
  
  for (word_type::id_type id = 0; id != fertility.map.size(); ++ id)
    if (fertility.map[id] != ntable_type::index_type(-1))
      karma::generate(iter, standard::string << ' ' << (real % ' ') << '\n',
		      word_type(id),
		      boost::make_iterator_range(fertility.begin(word_type(id)),
						 fertility.end(word_type(id))));
}


void read_classes(const path_type& path, classes_type& classes)
{
  typedef boost::spirit::istream_iterator iterator_type;
  typedef std::pair<std::string, std::string> word_pair_type;
  
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  
  qi::rule<iterator_type, std::string(), standard::blank_type>    word;
  qi::rule<iterator_type, word_pair_type(), standard::blank_type> parser; 
  
  word   %= qi::lexeme[+(standard::char_ - standard::space)];
  parser %= word >> word >> (qi::eol | qi::eoi);

  classes.clear();
  
  utils::compress_istream is(path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iterator_type iter(is);
  iterator_type iter_end;
  
  word_pair_type word_pair;
  
  while (iter != iter_end) {
    word_pair.first.clear();
    word_pair.second.clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, word_pair))
      if (iter != iter_end)
	throw std::runtime_error("cluster parsing failed");
    
    const word_type cluster(word_pair.first);
    const word_type word(word_pair.second);
    
    classes[word] = cluster;
  }
  
  classes[vocab_type::BOS] = vocab_type::BOS;
  classes[vocab_type::EOS] = vocab_type::EOS;
  classes[vocab_type::UNK] = vocab_type::UNK;
  classes[vocab_type::EPSILON] = vocab_type::EPSILON;
  classes[vocab_type::NONE] = vocab_type::NONE;

  classes.shrink();
}

#endif
