// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__COMPACT_HASHTABLE__HPP__
#define __UTILS__COMPACT_HASHTABLE__HPP__

#include <stdint.h>
#include <climits>
#include <stdexcept>
#include <algorithm>
#include <iterator>

#include <boost/numeric/conversion/bounds.hpp>

#include <utils/bithack.hpp>
#include <utils/memory.hpp>

namespace utils
{
  // bucket impelemntation
  template <typename _Tp, typename _Alloc>
  struct __compact_hashtable_bucket : public _Alloc
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef _Tp    value_type;
    typedef _Alloc allocator_type;
    typedef _Tp*   pointer;
    
    typedef       _Tp* iterator;
    typedef const _Tp* const_iterator;
    typedef       _Tp& reference;
    typedef const _Tp& const_reference;
    
    __compact_hashtable_bucket() : __first(0), __last(0) {}
    __compact_hashtable_bucket(size_type __n, const value_type& x)
      : __first(0), __last(0)
    {
      if (__n) {
	try {
	  __first = allocator().allocate(__n);
	  __last  = __first + __n;
	}
	catch(...) {
	  throw;
	}
	
	std::uninitialized_fill(__first, __last, x);
      }
    }
    
    template <typename Iterator>
    __compact_hashtable_bucket(Iterator first, Iterator last)
      : __first(0), __last(0)
    {
      typedef typename boost::is_integral<Iterator>::type __integral;
	
      __initialize_dispatch(first, last, __integral());
    }
    
    ~__compact_hashtable_bucket() 
    {
      if (! empty()) {
	utils::destroy_range(__first, __last);
	allocator().deallocate(__first, std::distance(__first, __last));
      }
    }
      
    bool empty() const { return __first == __last; }
    size_type size() const { return std::distance(__first, __last); }

    inline const_iterator begin() const { return reinterpret_cast<const_iterator>(__first); }
    inline       iterator begin()       { return reinterpret_cast<iterator>(__first); }
    inline const_iterator end() const { return reinterpret_cast<const_iterator>(__last); }
    inline       iterator end()       { return reinterpret_cast<iterator>(__last); }

    const_reference operator[](size_type x) const { return __first[x]; }
    reference operator[](size_type x) { return __first[x]; }
    
    void swap(__compact_hashtable_bucket& x)
    {
      std::swap(__first, x.__first);
      std::swap(__last,  x.__last);
    }
    
    const allocator_type& allocator() const { return static_cast<allocator_type&>(*this); }
    allocator_type& allocator() { return static_cast<allocator_type&>(*this); }

    template <typename Integer>
    void __initialize_dispatch(Integer __n, Integer x, boost::true_type)
    {
      if (__n) {
	try {
	  __first = allocator().allocate(__n);
	  __last  = __first + __n;
	}
	catch(...) {
	  throw;
	}
	
	std::uninitialized_fill(__first, __last, x);
      }
    }
    
    template <typename _InputIterator>
    void __initialize_dispatch(_InputIterator first, _InputIterator last, boost::false_type)
    {
      size_type __n = std::distance(first, last);

      if (__n) {
	try {
	  __first = allocator().allocate(__n);
	  __last  = __first + __n;
	}
	catch(...) {
	  throw;
	}
	
	std::uninitialized_copy(first, last, __first);
      }
    }
    
  private:
    pointer __first;
    pointer __last;
  };
  
  template <typename Table, typename Iterator, typename Reference, typename Tp>
  struct __compact_hashtable_iterator
  {
    typedef std::forward_iterator_tag iterator_category;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef Tp        value_type;
    typedef Reference reference;
    typedef Tp*       pointer;

    typedef __compact_hashtable_iterator<Table, Iterator, Reference, Tp> iterator;

    __compact_hashtable_iterator() : table(0), pos(), last() {}
    __compact_hashtable_iterator(const __compact_hashtable_iterator& x)
      : table(x.table), pos(x.pos), last(x.last) {}
    
    template <typename T, typename I, typename R, typename V>
    __compact_hashtable_iterator(const __compact_hashtable_iterator<T,I,R,V>& x)
      : table(x.table), pos(const_cast<pointer>(x.pos)), last(const_cast<pointer>(x.last)) {}
    
    __compact_hashtable_iterator(const Table& __table,
				 Iterator __pos,
				 Iterator __last,
				 bool forward)
      : table(&__table), pos(const_cast<pointer>(__pos)), last(const_cast<pointer>(__last))
    {
      if (forward)
	advance();
    }
    
    reference operator*() const { return *pos; }
    Iterator operator->() const { return &(operator*()); }
    
    iterator& operator++()
    {
      ++ pos;
      advance();
      return *this;
    }
    
    iterator operator++(int)
    {
      iterator tmp(*this);
      ++ *this;
      return tmp;
    }

    void advance()
    {
      for (/**/; pos != last && (table->is_deleted(*this) || table->is_empty(*this)); ++ pos);
    }
    
    const Table* table;
    pointer pos;
    pointer last;
  };

  template <typename T1, typename I1, typename R1, typename V1,
	    typename T2, typename I2, typename R2, typename V2>
  inline
  bool operator==(const __compact_hashtable_iterator<T1,I1,R1,V1>& x,
		  const __compact_hashtable_iterator<T2,I2,R2,V2>& y)
  {
    return x.pos == y.pos;
  }

  template <typename T1, typename I1, typename R1, typename V1,
	    typename T2, typename I2, typename R2, typename V2>
  inline
  bool operator!=(const __compact_hashtable_iterator<T1,I1,R1,V1>& x,
		  const __compact_hashtable_iterator<T2,I2,R2,V2>& y)
  {
    return x.pos != y.pos;
  }
  
  template <typename Key, typename Value, typename ExtractKey, typename Hash, typename Pred, typename Alloc>
  class compact_hashtable : public ExtractKey,
                            public Hash,
                            public Pred
  {
  public:
    typedef Key        key_type;
    typedef Value      value_type;
    typedef ExtractKey extract_key_type;
    typedef Hash       hash_type;
    typedef Pred       pred_type;
    typedef Alloc      allocator_type;
    
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    
  private:
    typedef compact_hashtable<Key, Value, ExtractKey, Hash, Pred, Alloc> self_type;
    typedef __compact_hashtable_bucket<Value, Alloc> bucket_type;

  public:
    typedef typename bucket_type::reference       reference;
    typedef typename bucket_type::const_reference const_reference;
    typedef typename bucket_type::pointer         pointer;

  public:
    typedef __compact_hashtable_iterator<self_type, typename bucket_type::iterator, reference, value_type> iterator;
    typedef __compact_hashtable_iterator<self_type, typename bucket_type::const_iterator, const_reference, value_type> const_iterator;
    
  private:
    static const size_type __cache_bytes = 128;
    static const size_type __cache_size = __cache_bytes / sizeof(value_type);
    static const size_type __cache_max = (__cache_size == 0 ? size_type(1) : __cache_size);
    
  public:
    compact_hashtable(size_type hint=0)
      : __bucket(),
	__value_empty(),
	__value_deleted(),
	__size_element(0),
	__size_deleted(0) {  }
    compact_hashtable(const compact_hashtable& x,
		      size_type minimum_size=0)
      : __bucket(),
	__value_empty(x.__value_empty),
	__value_deleted(x.__value_deleted),
	__size_element(0),
	__size_deleted(0) 
    {
      if (minimum_size == 0)
	assign(x);
      else
	initialize(x, minimum_size);
    }
    
  public:
    compact_hashtable& operator=(const compact_hashtable& x)
    {
      assign(x);
      return *this;
    }
    
    void assign(const compact_hashtable& x)
    {
      if (this == &x) return;
      
      // new bucket, then, swap
      if (x.empty()) {
	bucket_type __bucket_new;
	__bucket.swap(__bucket_new);
      } else {
	bucket_type __bucket_new(x.__bucket.begin(), x.__bucket.end());
	__bucket.swap(__bucket_new);
      }
      
      set_value(__value_empty, x.__value_empty);
      set_value(__value_deleted, x.__value_deleted);
      __size_element = x.__size_element;
      __size_deleted = x.__size_deleted;

      extract_key() = x.extract_key();
      hash() = x.hash();
      pred() = x.pred();
    }

    void swap(compact_hashtable& x)
    {
      __bucket.swap(x.__bucket);
      
      // very strange swapping
      {
	const value_type tmptmp(x.__value_empty);
	set_value(x.__value_empty, __value_empty);
	set_value(__value_empty, tmptmp);
      }
      
      {
	const value_type tmptmp(x.__value_deleted);
	set_value(x.__value_deleted, __value_deleted);
	set_value(__value_deleted, tmptmp);
      }
      
      std::swap(__size_element, x.__size_element);
      std::swap(__size_deleted, x.__size_deleted);
      
      std::swap(extract_key(), x.extract_key());
      std::swap(hash(), x.hash());
      std::swap(pred(), x.pred());
    }
    
    void clear()
    {
      typename bucket_type::iterator biter_end = __bucket.end();
      for (typename bucket_type::iterator biter = __bucket.begin(); biter != biter_end; ++ biter)
	set_value(*biter, __value_empty);
      
      __size_element = 0;
      __size_deleted = 0;
    }
    
    bool empty() const { return size() == 0; }
    size_type size() const { return __size_element - __size_deleted; }
    size_type bucket_count() const { return __bucket.size(); }
    
    void resize(size_type __n) { rehash(__n); }
    void rehash(size_type minimum_size)
    {
      if (minimum_size <= __size_element) return;
      
      size_type target = 0;
      
      if (__bucket.size() <= __cache_max) { // we are linear-table mode...
	const size_type size_power2 = capacity_power2(minimum_size);
	
	if (size_power2 <= __cache_max) // still, keep linear-table
	  target = size_power2;
	else // migrate to hash-table mode...
	  target = capacity_power2(minimum_size + (minimum_size >> 2));
      } else // we are hash-table mode...
	target = capacity_power2(minimum_size + (minimum_size >> 2));

      if (target > __bucket.size()) {
	compact_hashtable table_new(*this, minimum_size);
	
	swap(table_new);
      }
    }
    
    const_iterator begin() const { return const_iterator(*this, __bucket.begin(), __bucket.end(), true); }
    iterator begin() { return iterator(*this, __bucket.begin(), __bucket.end(), true); }
      
    const_iterator end() const { return const_iterator(*this, __bucket.end(), __bucket.end(), true); }
    iterator end() { return iterator(*this, __bucket.end(), __bucket.end(), true); }

    bool is_deleted(const_iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(__value_deleted));
    }
    bool is_deleted(iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(__value_deleted));
    }
      
    bool is_empty(const_iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(__value_empty));
    }
    bool is_empty(iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(__value_empty));
    }
    
    const_iterator find(const key_type& key) const
    {
      if (empty()) return end();
      
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_max
						   ? find_linear(key)
						   : find_bucket(key));
      
      if (pos.first == size_type(-1))
	return end();
      else
	return const_iterator(*this, __bucket.begin() + pos.first, __bucket.end(), false);
    }

    iterator find(const key_type& key)
    {
      if (empty()) return end();
      
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_max
						   ? find_linear(key)
						   : find_bucket(key));
      
      if (pos.first == size_type(-1))
	return end();
      else
	return iterator(*this, __bucket.begin() + pos.first, __bucket.end(), false);
    }
    
    size_type erase(const key_type& key)
    {
      iterator iter = find(key);
      
      if (iter != end()) {
	set_value(*iter, __value_deleted);
	++ __size_deleted;
	return 1;
      } else
	return 0;
    }
    
    void erase(iterator iter)
    {
      if (iter == end()) return;
      
      const key_type& key = extract_key()(*iter.pos);
      
      if (! pred()(key, extract_key()(__value_empty)) && ! pred()(key, extract_key()(__value_deleted))) {
	set_value(*iter.pos, __value_deleted);
	++ __size_deleted;
      }
    }
    
    void erase(iterator first, iterator last)
    {
      const key_type& key_empty   = extract_key()(__value_empty);
      const key_type& key_deleted = extract_key()(__value_deleted);
      
      for (/**/; first != last; ++ first) {
	const key_type& key = extract_key()(*first.pos);
	
	if (! pred()(key, key_empty) && ! pred()(key, key_deleted)) {
	  set_value(*first.pos, __value_deleted);
	  ++ __size_deleted;
	}
      }
    }

    void erase(const_iterator iter)
    {
      if (iter == end()) return;
      
      const key_type& key = extract_key()(*iter.pos);
      
      if (! pred()(key, extract_key()(__value_empty)) && ! pred()(key, extract_key()(__value_deleted))) {
	set_value(*iter.pos, __value_deleted);
	++ __size_deleted;
      }
    }
    
    void erase(const_iterator first, const_iterator last)
    {
      const key_type& key_empty   = extract_key()(__value_empty);
      const key_type& key_deleted = extract_key()(__value_deleted);
      
      for (/**/; first != last; ++ first) {
	const key_type& key = extract_key()(*first.pos);
	
	if (! pred()(key, key_empty) && ! pred()(key, key_deleted)) {
	  set_value(*first.pos, __value_deleted);
	  ++ __size_deleted;
	}
      }
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      insert(first, last, typename std::iterator_traits<Iterator>::iterator_category());
    }
    
    std::pair<iterator, bool> insert(const value_type& x)
    {
      rehash(__size_element + 1);
      
      return insert_noresize(x);
    }

    
  private:

    std::pair<size_type, size_type> find_linear(const key_type& key) const
    {
      const key_type& key_empty   = extract_key()(__value_empty);
      const key_type& key_deleted = extract_key()(__value_deleted);

      size_type pos_insert = size_type(-1);
      for (size_type pos_buck = 0; pos_buck != __bucket.size(); ++ pos_buck) {
	const key_type& key_buck = extract_key()(__bucket[pos_buck]);
	
	if (pred()(key_buck, key_empty) || pred()(key_buck, key_deleted))
	  pos_insert = utils::bithack::branch(pos_insert == size_type(-1), pos_buck, pos_insert);
	else if (pred()(key_buck, key))
	  return std::make_pair(pos_buck, size_type(-1));
      }
      
      return std::make_pair(size_type(-1), pos_insert);
    }
    
    std::pair<size_type, size_type> find_bucket(const key_type& key) const
    {
      size_type num_probes = 0;
      
      size_type pos_buck = hash()(key) & (__bucket.size() - 1);
      size_type pos_insert = size_type(-1);
      
      const key_type& key_empty   = extract_key()(__value_empty);
      const key_type& key_deleted = extract_key()(__value_deleted);
      
      for (;;) {
	const key_type& key_buck = extract_key()(__bucket[pos_buck]);
	
	if (pred()(key_buck, key_empty)) { // no searching further
	  if (pos_insert == size_type(-1))
	    return std::make_pair(size_type(-1), pos_buck);
	  else
	    return std::make_pair(size_type(-1), pos_insert);
	} else if (pred()(key_buck, key_deleted)) // searching...
	  pos_insert = utils::bithack::branch(pos_insert == size_type(-1), pos_buck, pos_insert);
	else if (pred()(key_buck, key))
	  return std::make_pair(pos_buck, size_type(-1));
	
	// quadratic probing...
	++ num_probes;
	pos_buck = (pos_buck + num_probes) & (__bucket.size() - 1);
      }
    }

    size_type capacity_power2(size_type n)
    {
      return bithack::branch(bithack::is_power2(n), n, static_cast<size_type>(bithack::next_largest_power2(n)));
    }
    
    void initialize(const compact_hashtable& table, size_type minimum_size = 0)
    {
      minimum_size = utils::bithack::max(table.size(), minimum_size);
      
      if (minimum_size == 0) return;
      
      const size_type size_power2 = capacity_power2(minimum_size);
      
      if (size_power2 <= __cache_max) {
	bucket_type bucket_new(size_power2, __value_empty);

	__bucket.swap(bucket_new);
	
	initialize_linear(table);
      } else {
	// allocate enough space...
	const size_type size_power2 = capacity_power2(minimum_size + (minimum_size >> 2));

	bucket_type bucket_new(size_power2, __value_empty);
	__bucket.swap(bucket_new);
	
	initialize_bucket(table);
      }
    }
    
    void initialize_linear(const compact_hashtable& table)
    {
      // we assume that: we are empty, and the old bucket has no duplicates
      // enough bucket allocated

      if (table.empty()) return;

      const key_type& key_empty   = extract_key()(__value_empty);
      const key_type& key_deleted = extract_key()(__value_deleted);

      typename bucket_type::iterator iter = __bucket.begin();
      
      typename bucket_type::const_iterator biter_end = table.__bucket.end();
      for (typename bucket_type::const_iterator biter = table.__bucket.begin(); biter != biter_end; ++ biter) {
	const key_type& key = extract_key()(*biter);
	
	if (pred()(key, key_empty) || pred()(key, key_deleted)) continue;
	
	set_value(*iter, *biter);
	++ iter;
	++ __size_element;
      }
    }
    
    void initialize_bucket(const compact_hashtable& table)
    {
      // we assume that: we are empty, and the old bucket has no duplicates
      // enough bucket allocated

      if (table.empty()) return;
      
      const key_type& key_empty   = extract_key()(__value_empty);
      const key_type& key_deleted = extract_key()(__value_deleted);
      
      typename bucket_type::const_iterator biter_end = table.__bucket.end();
      for (typename bucket_type::const_iterator biter = table.__bucket.begin(); biter != biter_end; ++ biter) {
	const key_type& key = extract_key()(*biter);
	
	if (pred()(key, key_empty) || pred()(key, key_deleted)) continue;
	
	size_type num_probes = 0;
	size_type pos_buck = hash()(key) & (__bucket.size() - 1);
	
	for (;;) {
	  const key_type& key_buck = extract_key()(__bucket[pos_buck]);
	  
	  if (pred()(key_buck, key_empty)) break;
	  
	  ++ num_probes;
	  pos_buck = (pos_buck + num_probes) & (__bucket.size() - 1);
	}
	
	set_value(__bucket[pos_buck], *biter);
	++ __size_element;
      }
    }
    
    // insert when we already know that the storage is big enough
    std::pair<iterator, bool> insert_noresize(const value_type& x)
    {
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_max
						   ? find_linear(extract_key()(x))
						   : find_bucket(extract_key()(x)));
      if (pos.first != size_type(-1))
	return std::make_pair(iterator(*this, __bucket.begin() + pos.first, __bucket.end(), false), false);
      else {
	if (pred()(extract_key()(__bucket[pos.second]), extract_key()(__value_deleted)))
	  -- __size_deleted;
	else
	  ++ __size_element;
	
	set_value(__bucket[pos.second], x);
	return std::make_pair(iterator(*this, __bucket.begin() + pos.second, __bucket.end(), false), true);
      }
    }

    template <typename ForwardIterator>
    void insert(ForwardIterator first, ForwardIterator last, std::forward_iterator_tag)
    {
      size_type __n = std::distance(first, last);
      
      rehash(__size_element + __n);
      
      for (/**/; __n != 0; -- __n, ++ first)
	insert_noresize(*first);
    }
    
    template <typename InputIterator>
    void insert(InputIterator first, InputIterator last, std::input_iterator_tag)
    {
      for (/**/; first != last; ++ first)
	insert(*first);
    }
    
  public:
    void set_empty_key(const value_type& x)
    {
      set_value(__value_empty, x);
    }
    
    void set_deleted_key(const value_type& x)
    {
      set_value(__value_deleted, x);
    }

    extract_key_type& extract_key() { return static_cast<extract_key_type&>(*this); }
    const extract_key_type& extract_key() const { return static_cast<const extract_key_type&>(*this); }

    hash_type& hash() { return static_cast<hash_type&>(*this); }
    const hash_type& hash() const { return static_cast<const hash_type&>(*this); }

    pred_type& pred() { return static_cast<pred_type&>(*this); }
    const pred_type& pred() const { return static_cast<const pred_type&>(*this); }
    
  private:
    void set_value(value_type& dest, const value_type& x) 
    {
      utils::destroy_object(&dest);
      ::new(&dest) value_type(x);
    }
    
  private:
    bucket_type __bucket;
    value_type  __value_empty;
    value_type  __value_deleted;
    size_type   __size_element;
    size_type   __size_deleted;
  };
  
};

namespace std
{
  template <typename K, typename V, typename E, typename H, typename P, typename A>
  inline
  void swap(utils::compact_hashtable<K,V,E,H,P,A>& x,
	    utils::compact_hashtable<K,V,E,H,P,A>& y)
  {
    x.swap(y);
  }
  
};

#endif