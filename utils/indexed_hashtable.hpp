// -*- mode: c++ -*-

#ifndef __UTILS__INDEXED_HASHTABLE__HPP__
#define __UTILS__INDEXED_HASHTABLE__HPP__

#include <stdint.h>
#include <vector>
#include <climits>
#include <stdexcept>
#include <algorithm>

#include <boost/numeric/conversion/bounds.hpp>

#include <utils/bithack.hpp>
#include <utils/memory.hpp>
#include <utils/chunk_vector.hpp>


namespace utils
{

  class __indexed_hashtable_power2
  {
  public:
    typedef uint32_t index_type;
    typedef size_t   size_type;
        
  public:
    
    inline size_type power_of_two(size_type n)
    {
      return bithack::is_power2(n) ? n : bithack::next_largest_power2(n);
    }
  };
  
  
  // simple bin classes with power2 size assignment
  // we assume bins are integer types (unt32_t)
  // no fill when resizing...

  template <typename __Index, typename __Alloc>
  struct __indexed_hashtable_bin_impl
  {
    typedef typename __Alloc::template rebind<__Index >::other allocator_type;
  };
  
  template <typename __Index, typename __Alloc>
  class __indexed_hashtable_bin : public __indexed_hashtable_bin_impl<__Index, __Alloc>::allocator_type
  {
  public:
    typedef __Index           index_type;
    typedef size_t            size_type;
    
    typedef       index_type* pointer;
    typedef       index_type* iterator;
    typedef const index_type* const_iterator;
    
    typedef       index_type& reference;
    typedef const index_type& const_reference;

    typedef typename __indexed_hashtable_bin_impl<__Index, __Alloc>::allocator_type allocator_type;
    
  public:
    __indexed_hashtable_bin(size_type x=0) : base(0), last(0) { allocate(x); }
    __indexed_hashtable_bin(const __indexed_hashtable_bin& x) : base(0), last(0) { assign(x); }
    ~__indexed_hashtable_bin() { deallocate(); }
    
    __indexed_hashtable_bin& operator=(const __indexed_hashtable_bin& x)
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
    
    void assign(const __indexed_hashtable_bin& x)
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

    void swap(__indexed_hashtable_bin& x)
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

  template <typename Value, typename Index, typename Alloc, bool Trivial>
  struct __indexed_hashtable_node_impl {};
  
  template <typename Value, typename Index, typename Alloc>
  struct __indexed_hashtable_node_impl<Value, Index, Alloc, true>
  {
    typedef Value value_type;
    typedef Index index_type;
    
    struct value_index_type
    {
      value_type value;
      index_type next;
      
      value_index_type(const value_type& x) : value(x), next() {}
      value_index_type() : value(), next() {}
      
      value_type&       data()       { return value; }
      const value_type& data() const { return value; }
    };
    
    typedef typename Alloc::template rebind<value_index_type >::other allocator_type;
    typedef std::vector<value_index_type, allocator_type> vector_type;
  };
  
  template <typename Value, typename Index, typename Alloc>
  struct __indexed_hashtable_node_impl<Value, Index, Alloc, false>
  {
    typedef Value value_type;
    typedef Index index_type;
    
    struct value_index_type
    {
      value_type value;
      index_type next;
      
      value_index_type(const value_type& x) : value(x), next() {}
      value_index_type() : value(), next() {}
      
      value_type&       data()       { return value; }
      const value_type& data() const { return value; }
    };
    
    typedef typename Alloc::template rebind<value_index_type >::other allocator_type;
    typedef utils::chunk_vector<value_index_type, 4096 / sizeof(value_index_type), allocator_type> vector_type;
  };
  
  
  // we will simply use vector...
  // no explicit handling for allocated values...
  template <typename Value, typename Index, typename Alloc>
  class __indexed_hashtable_node : public __indexed_hashtable_node_impl<Value,Index,Alloc,boost::has_trivial_destructor<Value>::value>::vector_type
  {
  public:
    typedef Value value_type;
    typedef Index index_type;
    
    typedef typename __indexed_hashtable_node_impl<Value,Index,Alloc,boost::has_trivial_destructor<Value>::value>::vector_type vector_type;
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

    void swap(__indexed_hashtable_node& x) { nodes().swap(x.nodes()); }
    
  private:
    inline       vector_type& nodes() { return static_cast<vector_type&>(*this); }
    inline const vector_type& nodes() const { return static_cast<const vector_type&>(*this); }
  };
    
  
  template <typename Key, typename Value, typename ExtractKey, typename Hash, typename Equal, typename Alloc>
  class indexed_hashtable : public __indexed_hashtable_power2,
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
    typedef __indexed_hashtable_bin<index_type,
				    allocator_type> bin_set_type;
    typedef __indexed_hashtable_node<value_type,
				     index_type,
				     allocator_type> node_set_type;
    
    
  public:
    typedef       value_type& reference;
    typedef const value_type& const_reference;
    
    struct const_iterator : public node_set_type::const_iterator
    {
      const_iterator(const typename node_set_type::const_iterator& x) : node_set_type::const_iterator(x) {}
      
      const value_type& operator*() const { return node_set_type::const_iterator::operator*().data(); }
      const value_type* operator->() const { return &(node_set_type::const_iterator::operator*().data()); }
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
    bin_set_type   bins;
    node_set_type nodes;
    
  public:
    indexed_hashtable(const size_type size=8,const hash_type& __hash=hash_type(), const equal_type& __equal=equal_type())
      : Hash(__hash), Equal(__equal), bins(power_of_two(size)), nodes() { using namespace std; fill(bins.begin(), bins.end(), 0); }

    indexed_hashtable(const indexed_hashtable& x)
      : bins(), nodes() { assign(x); }
    
    ~indexed_hashtable() { clear(); }
    
    indexed_hashtable& operator=(const indexed_hashtable& x)
    {
      assign(x);
      return *this;
    }
    
  public:

    void assign(const indexed_hashtable& x)
    {
      clear();
      
      extract_key() = x.extract_key();
      hash() = x.hash();
      equal() = x.equal();
      
      bins = x.bins;
      nodes = x.nodes;
    }
    
    void swap(indexed_hashtable& x)
    {
      std::swap(extract_key(), x.extract_key());
      std::swap(hash(), x.hash());
      std::swap(equal(), x.equal());
      
      bins.swap(x.bins);
      nodes.swap(x.nodes);
    }
    
    size_type size() const { return nodes.size(); }
    bool      empty() const { return nodes.empty(); }
    
    inline const_reference operator[](const index_type x) const { return nodes[x].data(); }
    inline       reference operator[](const index_type x)       { return nodes[x].data(); }
    
    const_iterator begin() const { return const_iterator(nodes.begin()); }
    const_iterator end() const { return const_iterator(nodes.end()); }
    
    void clear() 
    { 
      std::fill(bins.begin(), bins.end(), 0);
      nodes.clear();
    }
    
    void resize(size_type size)
    {
      using namespace std;
      
      size = power_of_two(size);

      if (size <= bins.size()) return;
      
      bins.resize(size);
      fill(bins.begin(), bins.end(), 0);
      
      const size_type hash_mask = bins.size() - 1;
      
      index_type pos = 0;
      typename node_set_type::iterator iter_end = nodes.end();
      for (typename node_set_type::iterator iter = nodes.begin(); iter != iter_end; ++ iter, ++ pos) {
	const index_type key = hash()(extract_key()(iter->data())) & hash_mask;
	iter->next = bins[key];
	bins[key] = pos + 1;
      }
    }
    
    const_iterator find(const key_type& x) const
    {
      const size_type key = hash()(x);
      
      index_type index = bins[key & (bins.size() - 1)];
      for (/**/; index && ! equal()(extract_key()(nodes[index - 1].data()), x); index = nodes[index - 1].next);
      return (index ? begin() + index - 1 : end());
    }
    
    std::pair<iterator, bool> insert(const value_type& x)
    {
      const size_type key = hash()(extract_key()(x));
      index_type index = bins[key & (bins.size() - 1)];
      for (/**/; index && ! equal()(extract_key()(nodes[index - 1].data()), extract_key()(x)); index = nodes[index - 1].next);
      if (index) return std::make_pair(begin() + index - 1, false);
      
      // not found...
      
      // perform resizing if we store many nodes
      if (nodes.size() > bins.size())
	resize(2 * bins.size());
      
      // insert at new bins position...
      const index_type pos = key & (bins.size() - 1);
      
      nodes.push_back(x);
      nodes.back().next = bins[pos];
      bins[pos] = nodes.size();
      
      return std::make_pair(end() - 1, true);
    }

  };
};

namespace std
{
  
  template <typename Key, typename Value, typename ExtractKey, typename Hash, typename Equal, typename Alloc >
  void swap(utils::indexed_hashtable<Key, Value, ExtractKey, Hash, Equal, Alloc>& x,
	    utils::indexed_hashtable<Key, Value, ExtractKey, Hash, Equal, Alloc>& y)
  {
    x.swap(y);
  }
  
}

#endif 
