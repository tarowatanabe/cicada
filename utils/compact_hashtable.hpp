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
#include <cstring>

#include <boost/type_traits.hpp>

#include <utils/bithack.hpp>
#include <utils/memory.hpp>
#include <utils/compact_func.hpp>

namespace utils
{
  template <typename Key, typename Value, typename Empty, typename Deleted,
	    typename ExtractKey, typename Hash, typename Pred, typename Alloc>
  class compact_hashtable;
  
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

    typedef __compact_hashtable_bucket<_Tp,_Alloc> self_type;
    
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
    
    const allocator_type& allocator() const { return static_cast<const allocator_type&>(*this); }
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

    void assign(const self_type& x)
    {
      if (x.empty()) {
	self_type __bucket;
	swap(__bucket);
      } else if (size() != x.size()) {
	self_type __bucket(x.begin(), x.end());
	swap(__bucket);
      } else
	__assign_dispatch(x, boost::has_trivial_assign<value_type>());
    }
    
    void __assign_dispatch(const self_type& x, boost::true_type)
    {
      std::memcpy(__first, x.__first, sizeof(value_type) * size());
    }
    
    void __assign_dispatch(const self_type& x, boost::false_type)
    {
      utils::destroy_range(__first, __last);
      std::uninitialized_copy(x.__first, x.__last, __first);
    }
    
  private:
    pointer __first;
    pointer __last;
  };
  
  template <typename Table, typename Iterator, typename Reference, typename Tp>
  struct __compact_hashtable_iterator
  {
    template <typename K, typename V, typename _E, typename _D, typename E,typename H,typename P, typename A>
    friend class compact_hashtable;

    template <typename T, typename I, typename R, typename V>
    friend struct __compact_hashtable_iterator;

    typedef std::forward_iterator_tag iterator_category;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef Tp        value_type;
    typedef Reference reference;
    typedef Tp*       pointer;

    typedef __compact_hashtable_iterator<Table, Iterator, Reference, Tp> iterator;

    __compact_hashtable_iterator() : table(0), pos() {}
    __compact_hashtable_iterator(const __compact_hashtable_iterator& x)
      : table(x.table), pos(x.pos) {}
    
    template <typename T, typename I, typename R, typename V>
    __compact_hashtable_iterator(const __compact_hashtable_iterator<T,I,R,V>& x)
      : table(x.table), pos(const_cast<pointer>(x.pos)) {}
    
    __compact_hashtable_iterator(const Table& __table,
				 Iterator __pos,
				 bool forward)
      : table(&__table), pos(const_cast<pointer>(__pos))
    {
      if (forward)
	advance();
    }
    
  public:
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
    
  public:
    template <typename T1, typename I1, typename R1, typename V1,
	      typename T2, typename I2, typename R2, typename V2>
    friend
    bool operator==(const __compact_hashtable_iterator<T1,I1,R1,V1>& x,
		    const __compact_hashtable_iterator<T2,I2,R2,V2>& y);

    template <typename T1, typename I1, typename R1, typename V1,
	      typename T2, typename I2, typename R2, typename V2>
    friend
    bool operator!=(const __compact_hashtable_iterator<T1,I1,R1,V1>& x,
		    const __compact_hashtable_iterator<T2,I2,R2,V2>& y);

  private:
    void advance()
    {
      for (/**/; pos != table->__bucket.end() && (table->is_deleted(*this) || table->is_empty(*this)); ++ pos);
    }
    
  private:
    const Table* table;
    pointer pos;
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
  
  template <typename Key, typename Value, typename Empty, typename Deleted,
	    typename ExtractKey, typename Hash, typename Pred, typename Alloc>
  class compact_hashtable : public Empty,
			    public Deleted,
			    public ExtractKey,
                            public Hash,
                            public Pred
  {
    template <typename T, typename I, typename R, typename V>
    friend struct __compact_hashtable_iterator;

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
    typedef compact_hashtable<Key, Value, Empty, Deleted, ExtractKey, Hash, Pred, Alloc> self_type;
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
    static const size_type __cache_size_power2 = (utils::bithack::static_next_largest_power2<__cache_size>::result >> 2);
    static const size_type __cache_minimum = (__cache_size_power2 == 0 ? size_type(1) : __cache_size_power2);
    static const size_type __cache_linear = (__cache_minimum << 2);
    
  public:
    compact_hashtable(size_type hint=0,
		      const Hash& __hash = Hash(),
		      const Pred& __pred = Pred())
      : Hash(__hash),
	Pred(__pred),
	__bucket(),
	__size_element(0),
	__size_deleted(0) {  }
    compact_hashtable(const compact_hashtable& x)
      : __bucket(),
	__size_element(0),
	__size_deleted(0) 
    {
      assign(x);
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

      static_cast<Empty&>(*this) = static_cast<const Empty&>(x);
      static_cast<Deleted&>(*this) = static_cast<const Deleted&>(x);
      extract_key() = x.extract_key();
      hash() = x.hash();
      pred() = x.pred();
      
      __bucket.assign(x.__bucket);
      
      __size_element = x.__size_element;
      __size_deleted = x.__size_deleted;
    }

    void swap(compact_hashtable& x)
    {
      std::swap(static_cast<Empty&>(*this), static_cast<Empty&>(x));
      std::swap(static_cast<Deleted&>(*this), static_cast<Deleted&>(x));
      std::swap(extract_key(), x.extract_key());
      std::swap(hash(), x.hash());
      std::swap(pred(), x.pred());
      
      __bucket.swap(x.__bucket);
      
      std::swap(__size_element, x.__size_element);
      std::swap(__size_deleted, x.__size_deleted);
    }
    
    void clear()
    {
      utils::destroy_range(__bucket.begin(), __bucket.end());
      std::uninitialized_fill(__bucket.begin(), __bucket.end(), Empty::operator()());
      
      __size_element = 0;
      __size_deleted = 0;
    }
    
    bool empty() const { return size() == 0; }
    size_type size() const { return __size_element - __size_deleted; }
    size_type bucket_count() const { return __bucket.size(); }
    size_type occupied_count() const { return __size_element; }
    
    void resize(size_type __n) { rehash(__n); }
    
    void rehash(size_type minimum_size)
    {
      rehash_bucket(minimum_size);
    }
    
    const_iterator begin() const { return const_iterator(*this, __bucket.begin(), true); }
    iterator begin() { return iterator(*this, __bucket.begin(), true); }
      
    const_iterator end() const { return const_iterator(*this, __bucket.end(), false); }
    iterator end() { return iterator(*this, __bucket.end(), false); }
								
    const_iterator find(const key_type& key) const
    {
      if (empty()) return end();
      
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_linear
						   ? find_linear(key)
						   : find_bucket(key));
      
      if (pos.first == size_type(-1))
	return end();
      else
	return const_iterator(*this, __bucket.begin() + pos.first, false);
    }

    iterator find(const key_type& key)
    {
      if (empty()) return end();
      
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_linear
						   ? find_linear(key)
						   : find_bucket(key));
      
      if (pos.first == size_type(-1))
	return end();
      else
	return iterator(*this, __bucket.begin() + pos.first, false);
    }
    
    size_type erase(const key_type& key)
    {
      iterator iter = find(key);
      
      if (iter != end()) {
	copy_value(*iter, Deleted::operator()());
	++ __size_deleted;
	return 1;
      } else
	return 0;
    }
    
    void erase(iterator iter)
    {
      if (iter == end()) return;
      
      const key_type& key = extract_key()(*iter.pos);
      
      if (! pred()(key, extract_key()(Empty::operator()())) && ! pred()(key, extract_key()(Deleted::operator()()))) {
	copy_value(*iter.pos, Deleted::operator()());
	++ __size_deleted;
      }
    }
    
    void erase(iterator first, iterator last)
    {
      for (/**/; first != last; ++ first) {
	const key_type& key = extract_key()(*first.pos);
	
	if (! pred()(key, extract_key()(Empty::operator()()))
	    && ! pred()(key, extract_key()(Deleted::operator()()))) {
	  copy_value(*first.pos, Deleted::operator()());
	  ++ __size_deleted;
	}
      }
    }

    void erase(const_iterator iter)
    {
      if (iter == end()) return;
      
      const key_type& key = extract_key()(*iter.pos);
      
      if (! pred()(key, extract_key()(Empty::operator()())) && ! pred()(key, extract_key()(Deleted::operator()()))) {
	copy_value(*iter.pos, Deleted::operator()());
	++ __size_deleted;
      }
    }
    
    void erase(const_iterator first, const_iterator last)
    {
      for (/**/; first != last; ++ first) {
	const key_type& key = extract_key()(*first.pos);
	
	if (! pred()(key, extract_key()(Empty::operator()()))
	    && ! pred()(key, extract_key()(Deleted::operator()()))) {
	  copy_value(*first.pos, Deleted::operator()());
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
    
    template <typename DefaultValue>
    value_type& insert_default(const key_type& x)
    {
      rehash(__size_element + 1);
      
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_linear
						   ? find_linear(x)
						   : find_bucket(x));
      
      if (pos.first != size_type(-1))
	return __bucket[pos.first];
      else {
	if (pred()(extract_key()(__bucket[pos.second]), extract_key()(Deleted::operator()())))
	  -- __size_deleted;
	else
	  ++ __size_element;
	
	copy_value(__bucket[pos.second], DefaultValue()(x));
	return __bucket[pos.second];
      }
    }

  private:
    bool is_deleted(const_iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(Deleted::operator()()));
    }
    bool is_deleted(iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(Deleted::operator()()));
    }
      
    bool is_empty(const_iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(Empty::operator()()));
    }
    bool is_empty(iterator& x) const
    {
      return pred()(extract_key()(*x.pos), extract_key()(Empty::operator()()));
    }
    
  private:
    bool rehash_bucket(size_type minimum_size)
    {
      if (minimum_size <= __size_element) return false;

      const size_type capacity = capacity_bucket(minimum_size);
      
      if (capacity <= __bucket.size()) return false;
      
      bucket_type bucket_new(capacity, Empty::operator()());
      __bucket.swap(bucket_new);

      __size_element = 0;
      __size_deleted = 0;
      
      if (capacity <= __cache_linear)
	initialize_linear(bucket_new);
      else
	initialize_bucket(bucket_new);
      
      return true;
    }

    std::pair<size_type, size_type> find_linear(const key_type& key) const
    {
      size_type pos_insert = size_type(-1);
      for (size_type pos_buck = 0; pos_buck != __bucket.size(); ++ pos_buck) {
	const key_type& key_buck = extract_key()(__bucket[pos_buck]);
	
	if (pred()(key_buck, extract_key()(Empty::operator()())))
	  return std::make_pair(size_type(-1), utils::bithack::branch(pos_insert == size_type(-1), pos_buck, pos_insert));
	else if (pred()(key_buck, extract_key()(Deleted::operator()())))
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
      
      for (;;) {
	const key_type& key_buck = extract_key()(__bucket[pos_buck]);
	
	if (pred()(key_buck, extract_key()(Empty::operator()()))) // no searching further
	  return std::make_pair(size_type(-1), utils::bithack::branch(pos_insert == size_type(-1), pos_buck, pos_insert));
	else if (pred()(key_buck, extract_key()(Deleted::operator()()))) // searching...
	  pos_insert = utils::bithack::branch(pos_insert == size_type(-1), pos_buck, pos_insert);
	else if (pred()(key_buck, key))
	  return std::make_pair(pos_buck, size_type(-1));
	
	// quadratic probing...
	++ num_probes;
	pos_buck = (pos_buck + num_probes) & (__bucket.size() - 1);
      }
    }

    static inline
    size_type capacity_power2(size_type n)
    {
      return bithack::branch(bithack::is_power2(n), n, static_cast<size_type>(bithack::next_largest_power2(n)));
    }
    
    static inline
    size_type capacity_bucket(size_type n)
    {
      const size_type size_power2 = capacity_power2(n);
      
      return utils::bithack::branch(size_power2 <= __cache_linear,
				    utils::bithack::max(size_power2, __cache_minimum),
				    capacity_power2(n + (n >> 2)));
    }
    
    void initialize_linear(const bucket_type& bucket)
    {
      // we assume that: we are empty, and the old bucket has no duplicates
      // enough bucket allocated

      if (bucket.empty()) return;

      typename bucket_type::iterator iter = __bucket.begin();
      
      typename bucket_type::const_iterator biter_end = bucket.end();
      for (typename bucket_type::const_iterator biter = bucket.begin(); biter != biter_end; ++ biter) {
	const key_type& key = extract_key()(*biter);
	
	if (pred()(key, extract_key()(Empty::operator()()))
	    || pred()(key, extract_key()(Deleted::operator()()))) continue;
	
	copy_value(*iter, *biter);
	++ iter;
	++ __size_element;
      }
    }
    
    void initialize_bucket(const bucket_type& bucket)
    {
      // we assume that: we are empty, and the old bucket has no duplicates
      // enough bucket allocated

      if (bucket.empty()) return;
      
      typename bucket_type::const_iterator biter_end = bucket.end();
      for (typename bucket_type::const_iterator biter = bucket.begin(); biter != biter_end; ++ biter) {
	const key_type& key = extract_key()(*biter);
	
	if (pred()(key, extract_key()(Empty::operator()()))
	    || pred()(key, extract_key()(Deleted::operator()()))) continue;
	
	size_type num_probes = 0;
	size_type pos_buck = hash()(key) & (__bucket.size() - 1);
	
	for (;;) {
	  const key_type& key_buck = extract_key()(__bucket[pos_buck]);
	  
	  if (pred()(key_buck, extract_key()(Empty::operator()()))) break;
	  
	  ++ num_probes;
	  pos_buck = (pos_buck + num_probes) & (__bucket.size() - 1);
	}
	
	copy_value(__bucket[pos_buck], *biter);
	++ __size_element;
      }
    }
    
    // insert when we already know that the storage is big enough
    std::pair<iterator, bool> insert_noresize(const value_type& x)
    {
      const std::pair<size_type, size_type> pos = (__bucket.size() <= __cache_linear
						   ? find_linear(extract_key()(x))
						   : find_bucket(extract_key()(x)));
      if (pos.first != size_type(-1))
	return std::make_pair(iterator(*this, __bucket.begin() + pos.first, false), false);
      else {
	if (pred()(extract_key()(__bucket[pos.second]), extract_key()(Deleted::operator()())))
	  -- __size_deleted;
	else
	  ++ __size_element;
	
	copy_value(__bucket[pos.second], x);
	return std::make_pair(iterator(*this, __bucket.begin() + pos.second, false), true);
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
    extract_key_type& extract_key() { return static_cast<extract_key_type&>(*this); }
    const extract_key_type& extract_key() const { return static_cast<const extract_key_type&>(*this); }

    hash_type& hash() { return static_cast<hash_type&>(*this); }
    const hash_type& hash() const { return static_cast<const hash_type&>(*this); }

    pred_type& pred() { return static_cast<pred_type&>(*this); }
    const pred_type& pred() const { return static_cast<const pred_type&>(*this); }
    
  private:
    void copy_value(value_type& dest, const value_type& x) 
    {
      copy_value(dest, x, boost::has_trivial_assign<value_type>());
    }
    
    void copy_value(value_type& dest, const value_type& x, boost::true_type)
    {
      std::memcpy(&dest, &x, sizeof(value_type));
    }

    void copy_value(value_type& dest, const value_type& x, boost::false_type)
    {
      utils::destroy_object(&dest);
      utils::construct_object(&dest, x);
    }
      
    
  private:
    bucket_type __bucket;
    size_type   __size_element;
    size_type   __size_deleted;
  };
  
};

namespace std
{
  template <typename K, typename V, typename _E, typename _D, typename E, typename H, typename P, typename A>
  inline
  void swap(utils::compact_hashtable<K,V,_E,_D,E,H,P,A>& x,
	    utils::compact_hashtable<K,V,_E,_D,E,H,P,A>& y)
  {
    x.swap(y);
  }
  
};

#endif
