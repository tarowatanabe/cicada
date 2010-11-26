// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SYMBOL__HPP__
#define __UTILS__SYMBOL__HPP__ 1

//
// a thrad-unsafe implementation...
//
// Can we make it lock-free?
//

#include <stdint.h>

#include <stdexcept>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>

#include <boost/functional/hash/hash.hpp>
#include <boost/type_traits.hpp>

namespace utils
{
  template <typename Tp, typename Alloc, bool Trivial>
  struct __symbol_vector_impl {};
  
  template <typename Tp, typename Alloc>
  struct __symbol_vector_impl<Tp, Alloc, true>
  {
    typedef std::vector<Tp, Alloc> vector_type;
  };
  
  template <typename Tp, typename Alloc>
  struct __symbol_vector_impl<Tp, Alloc, false>
  {
    typedef utils::chunk_vector<Tp, 4096/sizeof(Tp), Alloc> vector_type;
  };

  template <typename Tp, bool Trivial>
  struct __symbol_clear {};
  
  template <typename Tp>
  struct __symbol_clear<Tp, true>
  {
    static inline
    void clear(Tp& value) {}
  };
  
  template <typename Tp>
  struct __symbol_clear<Tp, false>
  {
    static inline
    void clear(Tp& value)
    {
      static const Tp __default_value = Tp();
      value = __default_value;
    }
  };


  template <typename Tp, typename Hash=boost::hash<Tp>, typename Equal=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class __symbol_base : public Hash,
			public Equal
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp    value_type;
    typedef Hash  hash_type;
    typedef Equal equal_type;
    
    typedef uint32_t index_type;
    typedef uint32_t counter_type;

  private:
    struct Data
    {
      value_type   value;
      index_type   next;
      counter_type counter;
      
      Data() : value(), next(0), counter(0) {}
      Data(const value_type& _value,
	   const index_type& _next,
	   const counter_type& _counter)
	: value(_value), next(_next), counter(_counter) {}
    };
    typedef Data data_type;
    
    typedef typename Alloc::template rebind<index_type >::other index_allocator_type;
    typedef typename Alloc::template rebind<data_type >::other  data_allocator_type;
    
    typedef __symbol_vector_impl<data_type, data_allocator_type, boost::has_trivial_copy<data_type>::value > impl_type;
    
    typedef typename impl_type::vector_type               data_set_type;
    typedef std::vector<index_type, index_allocator_type> bin_set_type;
    
  public:
    __symbol_base() : bins(8, 0), garbage(0) {}
    
    const Tp& operator[](size_type pos) const { return data[pos].value; }
    size_type size() const { return data.size(); }
    
    counter_type counter(size_type pos) const
    {
      return data[pos].counter;
    }
    void increment(size_type pos)
    {
      ++ data[pos].counter;
    }
    
    void decrement(size_type pos)
    {
      -- data[pos].counter;
    }
    
    index_type insert(const Tp& x)
    {
      const size_t key = hash()(x);
      index_type index = bins[key & (bins.size() - 1)];
      index_type index_empty = 0;
      
      for (/**/; index && ! equal()(data[index - 1].value, x); index = data[index - 1].next) {
	const index_type mask = index_type(index_empty == 0 && data[index - 1].counter == 0) - 1;
	index_empty = ((~mask) & index) | (mask & index_empty);
      }
      
      if (index) {
	++ data[index - 1].counter;
	return index - 1;
      } else if (index_empty) {
	// not found, but there exists empty bucket on this chain...
	
	data[index_empty - 1].value = x;
	data[index_empty - 1].counter = 1;
	
	return index_empty - 1;
      } else {
	// not found and no empty bucket...
	
	if (data.size() > (bins.size() + (bins.size() >> 1)))
	  resize(2 * bins.size());
	
	const size_type pos = key & (bins.size() - 1);
	
	if (garbage) {
	  index = garbage;
	  garbage = data[index - 1].next;
	  
	  data[index - 1].value = x;
	  data[index - 1].next = bins[pos];
	  data[index - 1].counter = 1;
	  
	  bins[pos] = index;
	  
	  return index - 1;
	} else {
	  index = data.size() + 1;
	  
	  data.push_back(data_type(x, bins[pos], 1));
	  
	  bins[pos] = index;
	  
	  return index - 1;
	}
      }
    }
    
    void resize(size_type size)
    {
      if (! utils::bithack::is_power2(size))
	size = utils::bithack::next_largest_power2(size);

      if (size <= bins.size()) return;
      
      bin_set_type bins_new(size, 0);
      bins_new.swap(bins);
      
      garbage = 0;
      
      const size_type hash_mask = bins.size() - 1;
      
      index_type pos = 0;
      typename data_set_type::iterator iter_begin = data.begin();
      typename data_set_type::iterator iter_end = data.end();
      for (typename data_set_type::iterator iter = iter_begin; iter != iter_end; ++ iter, ++ pos) {
	if (iter->counter) {
	  const size_type bin_pos = hash()(iter->value) & hash_mask;
	  
	  iter->next = bins[bin_pos];
	  bins[bin_pos] = pos + 1;
	} else {
	  __symbol_clear<value_type, boost::has_trivial_destructor<value_type>::value >::clear(iter->value);
	  iter->next = garbage;
	  garbage = pos + 1;
	}
      }
    }
    
  private:
    hash_type& hash() { return static_cast<hash_type&>(*this); }
    const hash_type& hash() const { return static_cast<const hash_type&>(*this); }
    equal_type& equal() { return static_cast<equal_type&>(*this); }
    const equal_type& equal() const { return static_cast<const equal_type&>(*this); }
    
  private:
    bin_set_type  bins;
    data_set_type data;
    
    index_type    garbage;
  };
  
  
  template <typename Tp, typename Hash=boost::hash<Tp>, typename Equal=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class symbol
  {
  public:
    typedef uint32_t id_type;
    
  public:
    symbol() : m_id(id_type(-1)) { allocate(Tp()); }
    symbol(const Tp& x) : m_id(id_type(-1)) { allocate(x); }
    symbol(const symbol& x) : m_id(x.m_id) { increment(); }
    ~symbol() { decrement(); m_id = id_type(-1); }
    
    symbol& operator=(const symbol& x)
    {
      assign(x);
      return *this;
    }
    
    void assign(const symbol& x)
    {
      if (m_id == x.m_id) return;
      
      base().increment(x.m_id);
      base().decrement(m_id);
      m_id = x.m_id;
    }
    
    void swap(symbol& x) 
    { 
      std::swap(m_id, x.m_id);
    }
    
    const id_type& id() const { return m_id; }
    const Tp& value() const { return base()[m_id]; }
    operator const Tp&() const { return value(); }
    
    friend
    std::size_t hash_value(symbol const& x) { return utils::hashmurmur<size_t>()(x.m_id, size_t(0)); }
    
    friend
    bool operator<(const symbol& x, const symbol& y) { return x.m_id < y.m_id; }
    friend
    bool operator<=(const symbol& x, const symbol& y) { return x.m_id <= y.m_id; }
    friend
    bool operator>(const symbol& x, const symbol& y) { return x.m_id > y.m_id; }
    friend
    bool operator>=(const symbol& x, const symbol& y) { return x.m_id >= y.m_id; }
    friend
    bool operator==(const symbol& x, const symbol& y) { return x.m_id == y.m_id; }
    friend
    bool operator!=(const symbol& x, const symbol& y) { return x.m_id != y.m_id; }

    size_t counter() const
    {
      return base().counter(m_id);
    }

    size_t allocated() const 
    {
      return base().allocated();
    }
    
    size_t deallocated() const
    {
      return base().deallocated();
    }
    
  private:
    typedef __symbol_base<Tp, Hash, Equal, Alloc> base_type;
    
    base_type& base() const
    {
      static base_type __base;
      return __base;
    }
    
    void increment()
    {
      base().increment(m_id);
    }
    
    void decrement()
    {
      base().decrement(m_id);
    }
    
    void allocate(const Tp& x)
    {
      m_id = base().insert(x);
    }
    
  private:
    id_type m_id;
  };  
  
};

namespace std
{
  template <typename T, typename H, typename E, typename A>
  inline
  void swap(utils::symbol<T,H,E,A>& x, utils::symbol<T,H,E,A>& y)
  {
    x.swap(y);
  }
};

#endif
