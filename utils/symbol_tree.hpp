// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SYMBOL_TREE__HPP__
#define __UTILS__SYMBOL_TREE__HPP__ 1

#include <stdint.h>

#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>

#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>

#include <boost/functional/hash/hash.hpp>
#include <boost/type_traits.hpp>



namespace utils
{

  template <typename Tp, typename Alloc, bool Trivial>
  struct __symbol_tree_vector_impl {};
  
  template <typename Tp, typename Alloc>
  struct __symbol_tree_vector_impl<Tp, Alloc, true>
  {
    typedef std::vector<Tp, Alloc> vector_type;
  };
  
  template <typename Tp, typename Alloc>
  struct __symbol_tree_vector_impl<Tp, Alloc, false>
  {
    typedef utils::chunk_vector<Tp, 4096/sizeof(Tp), Alloc> vector_type;
  };



  template <typename Tp, bool Trivial>
  struct __symbol_tree_clear {};
  
  template <typename Tp>
  struct __symbol_tree_clear<Tp, true>
  {
    static inline
    void clear(Tp& value) {}
  };
  
  template <typename Tp>
  struct __symbol_tree_clear<Tp, false>
  {
    static inline
    void clear(Tp& value)
    {
      static const Tp __default_value = Tp();
      value = __default_value;
    }
  };


  template <typename Tp, typename Hash=boost::hash<Tp>, typename Equal=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class __symbol_tree_base : public Hash,
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
      index_type   parent;
      counter_type counter;
      
      Data() : value(), next(0), parent(), counter(0) {}
      Data(const value_type& _value,
	   const index_type& _next,
	   const index_type& _parent,
	   const counter_type& _counter)
	: value(_value), next(_next), parent(_parent), counter(_counter) {}
    };
    typedef Data data_type;

    typedef typename Alloc::template rebind<index_type >::other index_allocator_type;
    typedef typename Alloc::template rebind<data_type >::other  data_allocator_type;
    
    typedef __symbol_tree_vector_impl<data_type, data_allocator_type, boost::has_trivial_copy<data_type>::value > impl_type;
    
    typedef typename impl_type::vector_type               data_set_type;
    typedef std::vector<index_type, index_allocator_type> bin_set_type;
    
  public:
    __symbol_tree_base() : bins(8, 0), garbage(0) {}

    const Tp& operator[](size_type pos) const { return data[pos].value; }
    size_type size() { return data.size(); }
    
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
      if (data[pos].counter == 0 && data[pos].parent != index_type(-1))
	-- data[data[pos].parent].counter;
    }
    
    index_type parent(size_type pos) const { return data[pos].parent; }
    
    index_type insert(index_type parent_id, const Tp& x)
    {
      if (parent_id != index_type(-1) && data[parent_id].counter == 0)
	throw std::runtime_error("no parent???");
      
      const size_t key = utils::hashmurmur<size_t>()(hash()(x), parent_id);
      index_type index = bins[key & (bins.size() - 1)];
      index_type index_empty = 0;
      
      for (/**/; index && (data[index - 1].parent != parent_id || ! equal()(data[index - 1].value, x)); index = data[index - 1].next) {
	const index_type mask = index_type(index_empty == 0 && data[index - 1].counter == 0) - 1;
	index_empty = ((~mask) & index) | (mask & index_empty);
      }
      
      if (index) {
	++ data[index - 1].counter;
	if (parent_id != index_type(-1) && data[index - 1].counter == 1)
	  ++ data[parent_id].counter;
	
	return index - 1;
      } else if (index_empty) {
	// not found, but there exists empty bucket on this chain...

	
	data[index_empty - 1].value = x;
	data[index_empty - 1].parent = parent_id;
	data[index_empty - 1].counter = 1;
	
	if (parent_id != index_type(-1))
	  ++ data[parent_id].counter;
	
	return index_empty - 1;
      } else {
	// not found and no empty bucket...
	
	if (data.size() > (bins.size() + (bins.size() >> 1)))
	  resize(2 * bins.size());
	
	const size_t pos = key & (bins.size() - 1);
	
	if (garbage) {
	  index = garbage;
	  garbage = data[index - 1].next;
	  
	  data[index - 1].value = x;
	  data[index - 1].next = bins[pos];
	  data[index - 1].parent = parent_id;
	  
	  data[index - 1].counter = 1;
	  if (parent_id != index_type(-1))
	    ++ data[parent_id].counter;
	  
	  bins[pos] = index;
	  
	  return index - 1;
	  
	} else {
	  index = data.size() + 1;
	  
	  data.push_back(data_type(x, bins[pos], parent_id, 1));
	  
	  if (parent_id != index_type(-1))
	    ++ data[parent_id].counter;
	  
	  bins[pos] = index;
	  
	  return index - 1;
	}
      }
    }
    
    void resize(size_t size)
    {
      if (! utils::bithack::is_power2(size))
	size = utils::bithack::next_largest_power2(size);

      if (size <= bins.size()) return;
      
      bin_set_type bins_new(size, 0);
      bins_new.swap(bins);
      
      const size_t hash_mask = bins.size() - 1;
      
      garbage = 0;
      
      index_type pos = 0;
      typename data_set_type::iterator iter_begin = data.begin();
      typename data_set_type::iterator iter_end = data.end();
      for (typename data_set_type::iterator iter = iter_begin; iter != iter_end; ++ iter, ++ pos) {
	if (iter->counter) {
	  const size_type bin_pos = utils::hashmurmur<size_t>()(hash()(iter->value), iter->parent) & hash_mask;
	  
	  iter->next = bins[bin_pos];
	  bins[bin_pos] = pos + 1;
	} else {
	  __symbol_tree_clear<value_type, boost::has_trivial_destructor<value_type>::value>::clear(iter->value);
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
  class symbol_tree
  {
  public:
    typedef uint32_t id_type;
    
  public:
    symbol_tree() : m_id(id_type(-1)) { allocate(id_type(-1), Tp()); }
    explicit symbol_tree(const symbol_tree& parent, const Tp& x) : m_id(id_type(-1)) { allocate(parent.id(), x); }
    symbol_tree(const symbol_tree& x) : m_id(x.m_id) { increment(); }
    ~symbol_tree() { decrement(); m_id = id_type(-1); }
    
  private:
    explicit symbol_tree(const id_type& x) : m_id(x) { increment(); }

  public:
    
    symbol_tree& operator=(const symbol_tree& x)
    {
      assign(x);
      return *this;
    }
    
    void assign(const symbol_tree& x)
    {
      if (m_id == x.m_id) return;

      base().increment(x.m_id);
      base().decrement(m_id);
      m_id = x.m_id;
    }
    
    void swap(symbol_tree& x) 
    { 
      std::swap(m_id, x.m_id);
    }
    
    const id_type& id() const { return m_id; }
    const Tp& value() const { return base()[m_id]; }
    operator const Tp&() const { return value(); }
    bool is_root() const { return base().parent(m_id) == id_type(-1); }
    symbol_tree parent() const { return symbol_tree(base().parent(m_id)); }
    static const symbol_tree& root()
    {
      static symbol_tree __root;
      return __root;
    }
    
    friend
    std::size_t hash_value(symbol_tree const& x) { return utils::hashmurmur<size_t>()(x.m_id, size_t(0)); }
    
    friend
    bool operator<(const symbol_tree& x, const symbol_tree& y) { return x.m_id < y.m_id; }
    friend
    bool operator<=(const symbol_tree& x, const symbol_tree& y) { return x.m_id <= y.m_id; }
    friend
    bool operator>(const symbol_tree& x, const symbol_tree& y) { return x.m_id > y.m_id; }
    friend
    bool operator>=(const symbol_tree& x, const symbol_tree& y) { return x.m_id >= y.m_id; }
    friend
    bool operator==(const symbol_tree& x, const symbol_tree& y) { return x.m_id == y.m_id; }
    friend
    bool operator!=(const symbol_tree& x, const symbol_tree& y) { return x.m_id != y.m_id; }

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
    typedef __symbol_tree_base<Tp, Hash, Equal, Alloc> base_type;
        
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
    
    void allocate(const id_type& parent, const Tp& x)
    {
      m_id = base().insert(parent, x);
    }
    
  private:
    id_type m_id;
  };
    
};

namespace std
{
  template <typename T, typename H, typename E, typename A>
  inline
  void swap(utils::symbol_tree<T,H,E,A>& x, utils::symbol_tree<T,H,E,A>& y)
  {
    x.swap(y);
  }
};

#endif
