//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_IMPL__HPP__
#define __CICADA_LEXICON_IMPL__HPP__ 1

#include <cicada/sentence.hpp>
#include <cicada/alignment.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <boost/filesystem.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

#include <utils/alloc_vector.hpp>
#include <utils/bithack.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>

typedef cicada::Symbol     word_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Alignment  alignment_type;
typedef cicada::SpanVector span_set_type;
typedef cicada::Vocab      vocab_type;
typedef boost::filesystem::path path_type;

typedef double count_type;
typedef double prob_type;

struct classes_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef std::vector<word_type, std::allocator<word_type> > map_type;
  
  classes_type() : classes() {}

  void clear() { classes.clear(); }
  
  word_type operator[](const word_type& word) const
  {
    return (word.id() >= classes.size() ? vocab_type::UNK : classes[word.id()]);
  }
  
  word_type& operator[](const word_type& word)
  {
    if (word.id() >= classes.size())
      classes.resize(word.id() + 1, vocab_type::UNK);
    
    return classes[word.id()];
  }
  
  map_type classes;
};

struct atable_type
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef int       index_type;

  struct difference_map_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef int       index_type;
    
    typedef std::vector<count_type, std::allocator<count_type> > difference_set_type;

    difference_map_type() : positives(), negatives() {}

    count_type& operator[](const index_type& diff)
    {
      const size_type pos = utils::bithack::branch(diff >= 0, diff, - diff - 1);
      difference_set_type& diffs = (diff >= 0 ? positives : negatives);
      
      if (pos >= diffs.size())
	diffs.resize(pos + 1, 0.0);
      return diffs[pos];
    }
    
    count_type operator[](const index_type& diff) const
    {
      const size_type pos = utils::bithack::branch(diff >= 0, diff, - diff - 1);
      const difference_set_type& diffs = (diff >= 0 ? positives : negatives);
      
      return (pos >= diffs.size() ? 0.0 : diffs[pos]);
    }

    void initialize()
    {
      std::fill(positives.begin(), positives.end(), 0.0);
      std::fill(negatives.begin(), negatives.end(), 0.0);
    }

    bool empty() const { return positives.empty() && negatives.empty(); }
    
    difference_set_type positives;
    difference_set_type negatives;
  };

  typedef std::pair<word_type, word_type> class_pair_type;
  typedef std::pair<index_type, index_type> range_type;
  
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<class_pair_type, difference_map_type, utils::hashmurmur<size_t>, std::equal_to<class_pair_type>,
				  std::allocator<std::pair<const class_pair_type, difference_map_type> > > count_dict_type;
#else
  typedef sgi::hash_map<class_pair_type, difference_map_type, utils::hashmurmur<size_t>, std::equal_to<class_pair_type>,
			std::allocator<std::pair<const class_pair_type, difference_map_type> > > count_dict_type;
#endif

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<range_type, difference_map_type, utils::hashmurmur<size_t>, std::equal_to<range_type>,
				  std::allocator<std::pair<const range_type, difference_map_type> > > cache_type;
  typedef std::tr1::unordered_map<class_pair_type, cache_type, utils::hashmurmur<size_t>, std::equal_to<class_pair_type>,
				  std::allocator<std::pair<const class_pair_type, cache_type> > > cache_set_type;
#else
  typedef sgi::hash_map<range_type, difference_map_type, utils::hashmurmur<size_t>, std::equal_to<range_type>,
		      std::allocator<std::pair<const range_type, difference_map_type> > > cache_type;
  typedef sgi::hash_map<class_pair_type, cache_type, utils::hashmurmur<size_t>, std::equal_to<class_pair_type>,
			std::allocator<std::pair<const class_pair_type, cache_type> > > cache_set_type;
#endif

  atable_type(const double __prior=0.1) : prior(__prior) {}
  
  prob_type operator()(const word_type& source,
		       const word_type& target,
		       const index_type& source_size,
		       const index_type& target_size,
		       const index_type& i_prev,
		       const index_type& i) const
  {
    if (atable.empty()) return 1.0 / source_size;
    
    // i_prev < 0 implies BOS
    // i >= souce_size implies EOS
    //
    // we will cache wrt class_pair_type and diff's range
    //
    
    if (i_prev < 0) {
      // 0 <= i < source_size
      // which implies: 1 <= diff < source_size + 1
      //
      
      return estimate(class_pair_type(source, target), range_type(1, source_size + 1))[i + 1];
    } else if (i >= source_size) {
      // which implies: 1 <= diff < source_size - i_prev + 1
      // 
      
      return estimate(class_pair_type(source, target), range_type(1, source_size - i_prev + 1))[source_size - i_prev];
    } else {
      // 0 <= i < source_size
      // which implies: 0 - i_prev <= diff < source_size - i_prev
      //
      
      return estimate(class_pair_type(source, target), range_type(0 - i_prev, source_size - i_prev))[i - i_prev];
    }
  }
  
  const difference_map_type& estimate(const class_pair_type& classes, const range_type& range) const
  {
    difference_map_type& diffs = const_cast<cache_set_type&>(caches)[classes][range];
    if (diffs.empty()) {
      double sum = 0.0;
      for (index_type i = range.first; i != range.second; ++ i) {
	count_dict_type::const_iterator aiter = atable.find(classes);
	const double count = (aiter != atable.end() ? aiter->second[i] + prior : prior);
	
	diffs[i] = count;
	sum += count;
      }
      
      const double sum_digamma = utils::mathop::digamma(sum);
      for (index_type i = range.first; i != range.second; ++ i)
	diffs[i] = utils::mathop::exp(utils::mathop::digamma(diffs[i]) - sum_digamma);
    }
    
    return diffs;
  }

  difference_map_type& operator[](const class_pair_type& x)
  {
    return atable[x];
  }

  difference_map_type& operator()(const word_type& source, const word_type& target)
  {
    return atable[class_pair_type(source, target)];
  }

  count_type& operator()(const word_type& source, const word_type& target, const index_type& diff)
  {
    return atable[class_pair_type(source, target)][diff];
  }
  
  void clear() { atable.clear(); caches.clear(); }
  
  void initialize()
  {
    count_dict_type::iterator aiter_end = atable.end();
    for (count_dict_type::iterator aiter = atable.begin(); aiter != aiter_end; ++ aiter)
      aiter->second.initialize();
    
    cache_set_type::iterator citer_end = caches.end();
    for (cache_set_type::iterator citer = caches.begin(); citer != citer_end; ++ citer)
      citer->second.clear();
  }
  
  count_dict_type atable;
  cache_set_type  caches;
  double prior;
};

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

  LearnBase(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source,
	    const atable_type& __atable_source_target,
	    const atable_type& __atable_target_source)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      atable_source_target(__atable_source_target),
      atable_target_source(__atable_target_source),
      objective_source_target(0),
      objective_target_source(0)
  {}

  void initialize()
  {
    ttable_counts_source_target.initialize();
    ttable_counts_target_source.initialize();

    atable_counts_source_target.initialize();
    atable_counts_target_source.initialize();
    
    aligned_source_target.initialize();
    aligned_target_source.initialize();
    
    objective_source_target = 0.0;
    objective_target_source = 0.0;
  }
  
  const ttable_type& ttable_source_target;
  const ttable_type& ttable_target_source;
  ttable_type ttable_counts_source_target;
  ttable_type ttable_counts_target_source;

  atable_type atable_source_target;
  atable_type atable_target_source;
  atable_type atable_counts_source_target;
  atable_type atable_counts_target_source;
  
  aligned_type aligned_source_target;
  aligned_type aligned_target_source;
  
  double objective_source_target;
  double objective_target_source;
};

struct ViterbiBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  ViterbiBase(const ttable_type& __ttable_source_target,
	      const ttable_type& __ttable_target_source)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source)
  {}
  ViterbiBase(const ttable_type& __ttable_source_target,
	      const ttable_type& __ttable_target_source,
	      const atable_type& __atable_source_target,
	      const atable_type& __atable_target_source)
    : ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      atable_source_target(__atable_source_target),
      atable_target_source(__atable_target_source)
  {}
  
  const ttable_type& ttable_source_target;
  const ttable_type& ttable_target_source;
  
  atable_type atable_source_target;
  atable_type atable_target_source;
};


#endif
