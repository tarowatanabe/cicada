// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// this is basically an indexed-hash-table, but allow "erasing"

#ifndef __UTILS__SYMBOL_HASHTABLE__HPP__
#define __UTILS__SYMBOL_HASHTABLE__HPP__

#include <stdint.h>
#include <vector>
#include <climits>
#include <stdexcept>
#include <algorithm>

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/type_traits.hpp>

#include <utils/bithack.hpp>
#include <utils/memory.hpp>
#include <utils/chunk_vector.hpp>


namespace utils
{

  class __symbol_hashtable_power2
  {
  public:
    typedef uint32_t index_type;
    typedef size_t   size_type;
    
  public:
    
    inline size_type power_of_two(size_type n)
    {
      return bithack::branch(bithack::is_power2(n), n, static_cast<size_type>(bithack::next_largest_power2(n)));
    }
  };
  
  // simple bin classes with power2 size assignment
  // we assume bins are integer types (unt32_t)
  // no fill when resizing...

  template <typename __Index, typename __Alloc>
  struct __symbol_hashtable_bin_impl
  {
    typedef typename __Alloc::template rebind<__Index >::other allocator_type;
  };
  
  template <typename __Index, typename __Alloc>
  class __symbol_hashtable_bin : public __symbol_hashtable_bin_impl<__Index, __Alloc>::allocator_type
  {
  public:
    typedef __Index           index_type;
    typedef size_t            size_type;
    
    typedef       index_type* pointer;
    typedef       index_type* iterator;
    typedef const index_type* const_iterator;
    
    typedef       index_type& reference;
    typedef const index_type& const_reference;

    typedef typename __symbol_hashtable_bin_impl<__Index, __Alloc>::allocator_type allocator_type;
    
  public:
    __symbol_hashtable_bin(size_type x=0) : base(0), last(0) { allocate(x); }
    __symbol_hashtable_bin(const __symbol_hashtable_bin& x) : base(0), last(0) { assign(x); }
    ~__symbol_hashtable_bin() { deallocate(); }
    
    __symbol_hashtable_bin& operator=(const __symbol_hashtable_bin& x)
    {
      assign(x);
      return *this;
    }
    
    inline bool empty() const { return base == last; }
    inline size_t size() const { return last - base; }
    
    inline       iterator begin() { return base; }
    inline const_iterator begin() const { return base; }
    inline       iterator end() { return last; }
    inline const_iterator end() const { return last; }
    
    inline       reference operator[](size_type x) { return *(base + x); }
    inline const_reference operator[](size_type x) const { return *(base + x); }

    void clear() 
    {
      deallocate();
    }
    
    void assign(const __symbol_hashtable_bin& x)
    {
      if (x.size() != size())
	resize(x.size());
      
      std::copy(x.begin(), x.end(), base);
    }
    
    void resize(size_type x)
    {
      if (size() == x) return;
      allocate(x);
    }

    void swap(__symbol_hashtable_bin& x)
    {
      using namespace std;
      
      swap(base, x.base);
      swap(last, x.last);
      swap(allocator(), x.allocator());
    }
    
  private:
    void deallocate()
    {
      if (! base) return;
      
      allocator().deallocate(base, size());
      
      base = 0;
      last = 0;
    }
    
    void allocate(size_type x)
    {
      if (base)
	deallocate();
      
      if (x > 0) {
	base = allocator().allocate(x);
	last = base + x;
      }
    }
    
    // allocator...
    inline allocator_type& allocator()
    {
      return static_cast<allocator_type&>(*this);
    }
    
  private:
    iterator base;
    iterator last;
  };
  
  
  // we will differentiate whether we have trivial destructor or not...
  // we support push-back, (like, std::vector)

  
  template <typename Tp, bool Trivial>
  struct __symbol_hashtable_value_clear {};
  
  template <typename Tp>
  struct __symbol_hashtable_value_clear<Tp, true>
  {
    static inline
    void clear(Tp& value) {}
  };
  
  template <typename Tp>
  struct __symbol_hashtable_value_clear<Tp, false>
  {
    static inline
    void clear(Tp& value)
    {
      value = Tp();
    }
  };


  template <typename Value, typename Index, typename Alloc, bool Trivial>
  struct __symbol_hashtable_node_impl {};
  
  template <typename Value, typename Index, typename Alloc>
  struct __symbol_hashtable_node_impl<Value, Index, Alloc, true>
  {
    typedef Value value_type;
    typedef Index index_type;
    
    struct value_index_type
    {
      value_type value;
      index_type next;
      bool       erased;
      
      value_index_type(const value_type& x) : value(x), next(), erased(false) {}
      value_index_type() : value(), next(), erased(false) {}
      
      value_type&       data()       { return value; }
      const value_type& data() const { return value; }
    };
    
    typedef typename Alloc::template rebind<value_index_type >::other allocator_type;
    typedef std::vector<value_index_type, allocator_type> vector_type;
  };
  
  template <typename Value, typename Index, typename Alloc>
  struct __symbol_hashtable_node_impl<Value, Index, Alloc, false>
  {
    typedef Value value_type;
    typedef Index index_type;
    
    struct value_index_type
    {
      value_type value;
      index_type next;
      bool       erased;
      
      value_index_type(const value_type& x) : value(x), next(), erased(false) {}
      value_index_type() : value(), next(), erased(false) {}
      
      value_type&       data()       { return value; }
      const value_type& data() const { return value; }
    };
    
    typedef typename Alloc::template rebind<value_index_type >::other allocator_type;
    typedef utils::chunk_vector<value_index_type, 4096 / sizeof(value_index_type), allocator_type> vector_type;
  };
  
  
  // we will simply use vector...
  // no explicit handling for allocated values...
  template <typename Value, typename Index, typename Alloc>
  class __symbol_hashtable_node : public __symbol_hashtable_node_impl<Value,Index,Alloc,boost::has_trivial_destructor<Value>::value>::vector_type
  {
  public:
    typedef Value value_type;
    typedef Index index_type;
    
    typedef typename __symbol_hashtable_node_impl<Value,Index,Alloc,boost::has_trivial_destructor<Value>::value>::vector_type vector_type;
    typedef typename vector_type::value_type value_index_type;
    typedef typename vector_type::size_type size_type;
    
    typedef typename vector_type::reference       reference;
    typedef typename vector_type::const_reference const_reference;
    typedef typename vector_type::iterator       iterator;
    typedef typename vector_type::const_iterator const_iterator;
    
  public:
    size_type size() const { return nodes().size(); }
    bool empty() const { return nodes().empty(); }
    
    iterator begin() { return nodes().begin(); }
    iterator end() { return nodes().end(); }
    
    const_iterator begin() const { return nodes().begin(); }
    const_iterator end() const { return nodes().end(); }

    reference front() { return nodes().front(); }
    reference back() { return nodes().back(); }
    
    const_reference front() const { return nodes().front(); }
    const_reference back() const { return nodes().back(); }
    
    reference       operator[](size_type x) { return nodes()[x]; }
    const_reference operator[](size_type x) const { return nodes()[x]; }
    
    void push_back(const value_type& x) { nodes().push_back(value_index_type(x)); }
    void clear() { nodes().clear(); }

    void swap(__symbol_hashtable_node& x) { nodes().swap(x.nodes()); }
    
  private:
    inline       vector_type& nodes() { return static_cast<vector_type&>(*this); }
    inline const vector_type& nodes() const { return static_cast<const vector_type&>(*this); }
  };
    
  
  template <typename Key, typename Value, typename ExtractKey, typename Hash, typename Equal, typename Alloc>
  class symbol_hashtable : public __symbol_hashtable_power2,
			    public ExtractKey,
			    public Hash,
			    public Equal
  {
  public:
    typedef Key        key_type;
    typedef Value      value_type;
    typedef ExtractKey extract_key_type;
    typedef Hash       hash_type;
    typedef Equal      equal_type;
    typedef Alloc      allocator_type;
    
    typedef size_t     size_type;
    
  private:
    typedef __symbol_hashtable_bin<index_type,
				    allocator_type> bin_set_type;
    typedef __symbol_hashtable_node<value_type,
				     index_type,
				     allocator_type> node_set_type;
    
    
  public:
    typedef       value_type& reference;
    typedef const value_type& const_reference;
    
    struct const_iterator : public node_set_type::const_iterator
    {
    private:
      typedef typename node_set_type::const_iterator iter_type;

    public:
      const_iterator() : iter_type(), end() {}
      const_iterator(iter_type __iter, iter_type __end) : iter_type(__iter), end(__end) { __move(); }
      
      const value_type& operator*() const { return iter_type::operator*().data(); }
      const value_type* operator->() const { return &(iter_type::operator*().data()); }
      
      const_iterator& operator++()
      {
	iter_type::operator++();
	__move();
	return *this;
      }
      
      const_iterator operator++(int)
      {
	const_iterator __tmp(*this);
	iter_type::operator++();
	__move();
	return __tmp;
      }
      
    private:
      void __move()
      {
	while (iter_type(*this) != end && iter_type::operator*().erased)
	  iter_type::operator++();
      }
      
      iter_type end;
    };
    typedef const_iterator iterator;
    
  private:
    extract_key_type& extract_key() { return static_cast<extract_key_type&>(*this); }
    const extract_key_type& extract_key() const { return static_cast<const extract_key_type&>(*this); }
    hash_type& hash() { return static_cast<hash_type&>(*this); }
    const hash_type& hash() const { return static_cast<const hash_type&>(*this); }
    equal_type& equal() { return static_cast<equal_type&>(*this); }
    const equal_type& equal() const { return static_cast<const equal_type&>(*this); }
    
  private:
    bin_set_type  bins;
    node_set_type nodes;
    index_type    trash;
    
  public:
    symbol_hashtable(const size_type size=8,const hash_type& __hash=hash_type(), const equal_type& __equal=equal_type())
      : Hash(__hash), Equal(__equal), bins(power_of_two(size)), nodes(), trash(0) { using namespace std; fill(bins.begin(), bins.end(), 0); }

    symbol_hashtable(const symbol_hashtable& x)
      : bins(), nodes(), trash(0) { assign(x); }
    
    ~symbol_hashtable() { clear(); }
    
    symbol_hashtable& operator=(const symbol_hashtable& x)
    {
      assign(x);
      return *this;
    }
    
  public:

    void assign(const symbol_hashtable& x)
    {
      clear();
      
      extract_key() = x.extract_key();
      hash() = x.hash();
      equal() = x.equal();
      
      bins = x.bins;
      nodes = x.nodes;
      trash = x.trash;
    }
    
    void swap(symbol_hashtable& x)
    {
      std::swap(extract_key(), x.extract_key());
      std::swap(hash(), x.hash());
      std::swap(equal(), x.equal());
      
      bins.swap(x.bins);
      nodes.swap(x.nodes);
      std::swap(trash, x.trash);
    }
    
    size_type size() const { return nodes.size(); }
    bool      empty() const { return nodes.empty(); }
    
    inline const_reference operator[](const index_type x) const { return nodes[x].data(); }
    inline       reference operator[](const index_type x)       { return nodes[x].data(); }
    
    const_iterator begin() const { return const_iterator(nodes.begin(), nodes.end()); }
    const_iterator end() const { return const_iterator(nodes.end(), nodes.end()); }
    
    void clear() 
    { 
      std::fill(bins.begin(), bins.end(), 0);
      nodes.clear();
      trash = 0;
    }
    
    void resize(size_type size)
    {
      using namespace std;
      
      size = power_of_two(size);

      if (size <= bins.size()) return;
      
      bins.resize(size);
      fill(bins.begin(), bins.end(), 0);

      trash = 0;
      
      const size_type hash_mask = bins.size() - 1;
      
      index_type pos = 0;
      typename node_set_type::iterator iter_end = nodes.end();
      for (typename node_set_type::iterator iter = nodes.begin(); iter != iter_end; ++ iter, ++ pos) {
	if (! iter->erased) {
	  const index_type key = hash()(extract_key()(iter->data())) & hash_mask;
	  iter->next = bins[key];
	  bins[key] = pos + 1;
	} else {
	  __symbol_hashtable_value_clear<value_type, boost::has_trivial_destructor<value_type>::value >::clear(iter->value);
	  iter->next = trash;
	  trash = pos + 1;
	}
      }
    }
    
    const_iterator find(const key_type& x) const
    {
      const size_type key = hash()(x);
      
      index_type index = bins[key & (bins.size() - 1)];
      for (/**/; index && (nodes[index - 1].erased || ! equal()(extract_key()(nodes[index - 1].data()), x)); index = nodes[index - 1].next) {}
      return const_iterator(index ? nodes.begin() + index - 1 : nodes.end(), nodes.end());
    }
    
    std::pair<iterator, bool> insert(const value_type& x)
    {
      const size_type key = hash()(extract_key()(x));
      index_type index = bins[key & (bins.size() - 1)];
      index_type index_empty = 0;
      
      for (/**/; index && ! equal()(extract_key()(nodes[index - 1].data()), extract_key()(x)); index = nodes[index - 1].next)
	if (! index_empty && nodes[index - 1].erased)
	  index_empty = index;
      
      if (index) {
	const bool erased = nodes[index - 1].erased;
	nodes[index - 1].erased = false;
	return std::make_pair(iterator(nodes.begin() + index - 1, nodes.end()), erased);
      } else if (index_empty) {
	
	nodes[index_empty - 1].value = x;
	nodes[index_empty - 1].erased = false;
	
	return std::make_pair(iterator(nodes.begin() + index_empty - 1, nodes.end()), true);
      } else {
	// not found...
      
	// perform resizing if we store many nodes
	if (nodes.size() > bins.size())
	  resize(2 * bins.size());
      
	// insert at new bins position...
	const index_type pos = key & (bins.size() - 1);
	
	if (trash) {
	  index = trash;
	  trash = nodes[index - 1].next;
	  
	  nodes[index - 1].value = x;
	  nodes[index - 1].next = bins[pos];
	  nodes[index - 1].erased = false;
	  
	  bins[pos]= index;
	  
	  return std::make_pair(iterator(nodes.begin() + index - 1, nodes.end()), true);
	} else {
	  nodes.push_back(x);
	  nodes.back().next = bins[pos];
	  
	  bins[pos] = nodes.size();
	  
	  return std::make_pair(iterator(nodes.end() - 1, nodes.end()), true);
	}
      }
    }

    void erase(const index_type x)
    {
      nodes[x].erased = true;
    }
    
    void erase(const key_type& x) const
    {
      const size_type key = hash()(x);
      index_type index = bins[key & (bins.size() - 1)];
      for (/**/; index && (nodes[index - 1].erased || ! equal()(extract_key()(nodes[index - 1].data()), x)); index = nodes[index - 1].next) {}
      
      if (index)
	nodes[index - 1].erased = true;
    }
    
  };
};

namespace std
{
  
  template <typename Key, typename Value, typename ExtractKey, typename Hash, typename Equal, typename Alloc >
  void swap(utils::symbol_hashtable<Key, Value, ExtractKey, Hash, Equal, Alloc>& x,
	    utils::symbol_hashtable<Key, Value, ExtractKey, Hash, Equal, Alloc>& y)
  {
    x.swap(y);
  }
  
}

#endif 
