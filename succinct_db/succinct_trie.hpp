// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __SUCCINCT_DB__SUCCINCT_TRIE__HPP__
#define __SUCCINCT_DB__SUCCINCT_TRIE__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <vector>
#include <deque>
#include <stdexcept>
#include <memory>
#include <utility>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#include <utils/map_file.hpp>
#include <utils/repository.hpp>
#include <utils/succinct_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/filesystem.hpp>
#include <utils/packed_vector.hpp>
#include <utils/packed_device.hpp>
#include <utils/atomicop.hpp>
#include <utils/array_power2.hpp>

// double-array like interface for succinct-trie structure
// key should be integral type. unsigned char (or uint8_t ) is recommented.
// any value for data...

namespace succinctdb
{
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  struct __succinct_trie_cursor
  {
    typedef Key      key_type;
    typedef Data     data_type;
    typedef Data     mapped_type;
    typedef Impl     impl_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef __succinct_trie_cursor<Key,Data,Impl,Alloc> self_type;

    __succinct_trie_cursor()
      : impl(0), node_pos(0) {}
    __succinct_trie_cursor(const size_type& _node_pos, const impl_type* _impl)
      : impl(_impl), node_pos(_node_pos) {}

    size_type node() const { return node_pos; }
    bool exists() const { return impl->exists(node_pos); }
    data_type data() const { return impl->data(node_pos); }
    size_type offset() const { return impl->offset(node_pos); }
    data_type operator*() const { return data(); }
    
    self_type operator++(int)
    {
      self_type tmp = *this;
      increment();
      return tmp;
    }
    
    self_type& operator++()
    {
      increment();
      return *this;
    }
    void increment()
    {
      if (! impl->is_next_sibling(node_pos) || impl->key(node_pos) != impl->key(node_pos + 1)) {
	impl = 0;
	node_pos = 0;
      } else
	++ node_pos;
    }
    
    const impl_type* impl;
    size_type node_pos;
  };
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  inline
  bool operator==(const __succinct_trie_cursor<Key,Data,Impl,Alloc>& x,
		  const __succinct_trie_cursor<Key,Data,Impl,Alloc>& y)
  {
    return x.impl == y.impl && x.node_pos == y.node_pos;
  }
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  inline
  bool operator!=(const __succinct_trie_cursor<Key,Data,Impl,Alloc>& x,
		  const __succinct_trie_cursor<Key,Data,Impl,Alloc>& y)
  {
    return x.impl != y.impl || x.node_pos != y.node_pos;
  }


  template <typename Key, typename Alloc>
  struct __succinct_trie_iterator_node
  {
    typedef Key       key_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef std::vector<key_type, Alloc > buffer_type;
    
    buffer_type buffer;
    size_type   pos;
    
    __succinct_trie_iterator_node() : buffer(), pos() {}
    __succinct_trie_iterator_node(const __succinct_trie_iterator_node& x) : buffer(x.buffer), pos(x.pos) {}
    __succinct_trie_iterator_node(const size_type& _pos)
      : buffer(), pos(_pos) {}
  };
  
  template <typename Key, typename Alloc>
  inline
  bool operator==(const __succinct_trie_iterator_node<Key,Alloc>& x, const __succinct_trie_iterator_node<Key,Alloc>& y)
  {
    return x.pos == y.pos && x.buffer == y.buffer;
  }
  
  template <typename Key, typename Alloc>
  inline
  bool operator!=(const __succinct_trie_iterator_node<Key,Alloc>& x, const __succinct_trie_iterator_node<Key,Alloc>& y)
  {
    return ! (x == y);
  }
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  class __succinct_trie_iterator
  {
  public:
    typedef Key      key_type;
    typedef Data     data_type;
    typedef Data     mapped_type;
    typedef Impl     impl_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef std::input_iterator_tag iterator_category;

    typedef __succinct_trie_cursor<Key,Data,Impl,Alloc> cursor;
    typedef __succinct_trie_cursor<Key,Data,Impl,Alloc> const_cursor;

  private:
    typedef __succinct_trie_iterator<Key,Data,Impl,Alloc> self_type;

    typedef typename Alloc::template rebind<key_type>::other key_alloc_type;
    typedef __succinct_trie_iterator_node<key_type, key_alloc_type>  node_type;
    
    typedef typename Alloc::template rebind<node_type>::other node_alloc_type;
    typedef std::deque<node_type, node_alloc_type > node_set_type;
    
  public:
    typedef typename node_type::buffer_type key_buffer_type;
    
  public:
    __succinct_trie_iterator() : impl(), nodes() {}
    __succinct_trie_iterator(const __succinct_trie_iterator& x) : impl(x.impl), nodes(x.nodes) {}
    __succinct_trie_iterator(const impl_type& _impl) : impl(&_impl), nodes() {}
    __succinct_trie_iterator(const size_type pos, const impl_type& _impl) : impl(&_impl), nodes(1, node_type(pos))
    {
      if (! exists())
	increment();
    }
    
  public:
    size_type node() const { return  nodes.back().pos; }
    
    cursor begin() const { return cursor(nodes.back().pos, impl); }
    cursor end() const { return cursor(); }
    
    bool exists() const { return impl->exists(nodes.back().pos); }
    data_type data() const { return impl->data(nodes.back().pos); }
    size_type offset() const { return impl->offset(nodes.back().pos); }
    data_type operator*() const { return data(); }
    
    const key_buffer_type& key() const { return nodes.back().buffer; }
    
    self_type operator++(int)
    {
      self_type tmp = *this;
      increment();
      return tmp;
    }
    
    self_type& operator++()
    {
      increment();
      return *this;
    }

    template <typename _Key, typename _Data, typename _Impl, typename _Alloc>
    friend
    bool operator==(const __succinct_trie_iterator<_Key,_Data,_Impl,_Alloc>& x,
		    const __succinct_trie_iterator<_Key,_Data,_Impl,_Alloc>& y);
    template <typename _Key, typename _Data, typename _Impl, typename _Alloc>
    friend
    bool operator!=(const __succinct_trie_iterator<_Key,_Data,_Impl,_Alloc>& x,
		    const __succinct_trie_iterator<_Key,_Data,_Impl,_Alloc>& y);
    
  private:
    void increment()
    {
      move();
      while (! nodes.empty() && ! impl->exists(nodes.back().pos))
	move();
    }
    
    void move()
    {
      node_type node = nodes.back();
      nodes.pop_back();
      node.buffer.push_back(key_type());
      
      // we push nodes in reverse order...
      std::pair<size_type, size_type> range = impl->range(node.pos);

      for (size_type node_pos = range.second; node_pos != range.first; -- node_pos) {
	const key_type key = impl->key(node_pos - 1);
	
	// we will exclude duplicates...
	// do we record the duplica?
	if (node_pos == range.second || key != nodes.back().buffer.back())
	  nodes.push_back(node);
	
	nodes.back().pos = node_pos - 1;
	nodes.back().buffer.back() = key;
      }
    }
    
  private:
    const impl_type* impl;
    node_set_type nodes;
  };

  template <typename Key, typename Data, typename Impl, typename Alloc>
  inline
  bool operator==(const __succinct_trie_iterator<Key,Data,Impl,Alloc>& x,
		  const __succinct_trie_iterator<Key,Data,Impl,Alloc>& y)
  {
    return x.impl == y.impl && x.nodes == y.nodes;
  }
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  inline
  bool operator!=(const __succinct_trie_iterator<Key,Data,Impl,Alloc>& x,
		  const __succinct_trie_iterator<Key,Data,Impl,Alloc>& y)
  {
    return x.impl != y.impl || x.nodes != y.nodes;
  }
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  struct __succinct_trie_reverse_iterator
  {
  public:
    typedef Key      key_type;
    typedef Data     data_type;
    typedef Data     mapped_type;
    typedef Impl     impl_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef std::input_iterator_tag iterator_category;
    
    typedef __succinct_trie_cursor<Key,Data,Impl,Alloc> cursor;
    typedef __succinct_trie_cursor<Key,Data,Impl,Alloc> const_cursor;

  private:
    typedef __succinct_trie_reverse_iterator<Key,Data,Impl,Alloc> self_type;

  public:
    __succinct_trie_reverse_iterator() : impl(), node_pos(0) {}
    __succinct_trie_reverse_iterator(const __succinct_trie_reverse_iterator& x) : impl(x.impl), node_pos(x.node_pos) {}
    __succinct_trie_reverse_iterator(const impl_type& _impl) : impl(&_impl), node_pos(0) {}
    __succinct_trie_reverse_iterator(const size_type _node_pos, const impl_type& _impl)
      : impl(&_impl), node_pos(_node_pos) {}
    
  public:
    cursor begin() const { return cursor(node_pos, impl); }
    cursor end() const { return cursor(); }
    
    size_type node() const { return node_pos; }
    bool exists() const { return impl->exists(node_pos); }
    data_type data() const { return impl->data(node_pos); }
    size_type offset() const { return impl->offset(node_pos); }
    data_type operator*() const { return data(); }
    
    key_type key() const
    {
      return impl->key(node_pos);
    }
    
    self_type& operator++()
    {
      increment();
      return *this;
    }
    
    self_type operator++(int)
    {
      self_type tmp = *this;
      increment();
      return *this;
    }
    
    template <typename _Key, typename _Data, typename _Impl, typename _Alloc>
    friend
    bool operator==(const __succinct_trie_reverse_iterator<_Key,_Data,_Impl,_Alloc>& x,
		    const __succinct_trie_reverse_iterator<_Key,_Data,_Impl,_Alloc>& y);
    template <typename _Key, typename _Data, typename _Impl, typename _Alloc>
    friend
    bool operator!=(const __succinct_trie_reverse_iterator<_Key,_Data,_Impl,_Alloc>& x,
		    const __succinct_trie_reverse_iterator<_Key,_Data,_Impl,_Alloc>& y);
    
  private:
    void increment()
    {
      node_pos = impl->parent(node_pos);
    }
    
  private:
    const impl_type* impl;
    size_type        node_pos;
  };

  template <typename Key, typename Data, typename Impl, typename Alloc>
  inline
  bool operator==(const __succinct_trie_reverse_iterator<Key,Data,Impl,Alloc>& x,
		  const __succinct_trie_reverse_iterator<Key,Data,Impl,Alloc>& y)
  {
    return x.impl == y.impl && x.node_pos == y.node_pos;
  }
  
  template <typename Key, typename Data, typename Impl, typename Alloc>
  inline
  bool operator!=(const __succinct_trie_reverse_iterator<Key,Data,Impl,Alloc>& x,
		  const __succinct_trie_reverse_iterator<Key,Data,Impl,Alloc>& y)
  {
    return x.impl != y.impl || x.node_pos != y.node_pos;
  }


  template <typename Key, typename Data, typename Alloc>
  class __succinct_trie_base
  {
  public:
    typedef Key      key_type;
    typedef Data     data_type;
    typedef Data     mapped_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef boost::filesystem::path path_type;
    
  public:
    static inline const size_type out_of_range()
    {
      return size_type(-1);
    }
    
    template <typename Index>
    key_type __key(const Index& index, size_type pos) const { return index[pos]; }
    
    template <typename Index, typename Positions>
    std::pair<size_type, size_type> __range(const Index& index, const Positions& positions, size_type node_pos) const
    {
      const typename Positions::size_type first = positions.select(node_pos + 1, false);
      const typename Positions::size_type last = positions.select(node_pos + 2, false);
      
      if (first == typename Positions::size_type(-1) || last == typename Positions::size_type(-1))
	return std::make_pair(index.size(), index.size());
      else
	return std::make_pair((first + 1) - (node_pos + 1), (last + 1) - (node_pos + 2));
    }
    
    template <typename IndexMap, typename Mapped>
    data_type __data(const IndexMap& index_map, const Mapped& mapped, size_type node_pos) const
    {
      return mapped[index_map.rank(node_pos, true) - 1];
    }

    template <typename IndexMap>
    size_type __offset(const IndexMap& index_map, size_type node_pos) const
    {
      return index_map.rank(node_pos, true) - 1;
    }
    
    template <typename IndexMap>
    bool __exists(const IndexMap& index_map, size_type node_pos) const
    {
      return (node_pos < index_map.size() && index_map[node_pos]);
    }
    
    template <typename Index, typename Positions>
    bool __has_children(const Index& index, const Positions& positions, size_type node_pos) const
    {
      const std::pair<size_type, size_type> result = __range(index, positions, node_pos);
      return result.first < result.second;
    }
    
    template <typename Positions>
    size_type __parent(const Positions& positions, size_type node_pos) const
    {
      // rank0(select1(node_pos + 1)) - 1
      // y = select1(x)  <=> x = rank1(y), (y + 1) - x = rank0(y)
      const size_type pos_one = positions.select(node_pos + 1, true);
      return pos_one - node_pos - 1;
    }
    
    template <typename Positions>
    bool __is_next_sibling(const Positions& positions, size_type node_pos) const
    {
      const size_type pos_one = positions.select(node_pos + 1, true);
      return positions.test(pos_one + 1);
    }
    
    template <typename Index, typename Positions>
    size_type __traverse(const Index& index,
			 const Positions& positions,
			 const key_type* key_buf,
			 size_type& node_pos,
			 size_type& key_pos,
			 size_type key_len) const
    {
      if (key_pos == key_len) return node_pos;
      
      for (/**/; key_pos < key_len; ++ key_pos) {
	const std::pair<size_type, size_type> node_range = __range(index, positions, node_pos);
	
	if (node_range.first == node_range.second) return out_of_range();
	
	const size_type pos = __lower_bound(index, node_range.first, node_range.second, key_buf[key_pos]);
	if (pos == node_range.second || key_buf[key_pos] < index[pos]) return out_of_range();
	
	node_pos = pos;
      }
      
      return node_pos;
    }
    
    template <typename Index>
    size_type __lower_bound(const Index& index, size_type first, size_type last, const key_type& key) const
    {
      size_type length = last - first;
      
      if (length <= 128) {
	typename Index::const_iterator iter = index.begin() + first;
	for (/**/; first != last && (*iter) < key; ++ first, ++ iter);
	return first;
      } else {
	typename Index::const_iterator __first = index.begin() + first;
	typename Index::const_iterator __middle;

	while (length > 0) {
	  const size_type half  = length >> 1;
	  const size_type middle = first + half;
	  __middle = __first + half;
	  
	  if (*__middle < key) {
	    first = middle + 1;
	    __first = __middle + 1;
	    length = length - half - 1;
	  } else 
	    length = half;
	}
	return first;
      }
    }

  };
  
  
  template <typename Key, typename Alloc, size_t KeySize>
  struct __succinct_trie_mapped_index_impl
  {
    typedef utils::packed_vector_mapped<Key, Alloc> mapped_set_type;
  };
  
  template <typename Key, typename Alloc>
  struct __succinct_trie_mapped_index_impl<Key, Alloc, 1>
  {
    typedef utils::map_file<Key, Alloc> mapped_set_type;
  };
  
  template <typename Key, typename Alloc>
  struct __succinct_trie_mapped_index
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;
    typedef typename __succinct_trie_mapped_index_impl<Key, Alloc, sizeof(Key) >::mapped_set_type impl_type;

    typedef typename impl_type::const_iterator const_iterator;
    
    __succinct_trie_mapped_index() {}
    __succinct_trie_mapped_index(const path_type& path) : __impl(path) {}
    
    Key operator[](size_type pos) const { return __impl[pos]; }
    
    const_iterator begin() const { return __impl.begin(); }
    const_iterator end() const { return __impl.end(); }
    
    void open(const path_type& path) { __impl.open(path); }
    void clear() { __impl.clear(); }
    void close() { __impl.close(); }
    
    bool empty() const { return __impl.empty(); }
    size_type size() const { return __impl.size(); }
    path_type path() const { return __impl.path(); }
    
    uint64_t size_bytes() const { return __impl.size_bytes(); }
    uint64_t size_compressed() const { return __impl.size_compressed(); }
    uint64_t size_cache() const { return __impl.size_cache(); }
    
    impl_type __impl;
  };
  
  template <typename Data, typename Alloc, bool Integral>
  struct __succinct_trie_mapped_data_impl {};

  template <typename Data, typename Alloc>
  struct __succinct_trie_mapped_data_impl<Data,Alloc,false>
  {
    typedef utils::map_file<Data, Alloc> mapped_set_type;
  };
  
  template <typename Data, typename Alloc>
  struct __succinct_trie_mapped_data_impl<Data,Alloc,true>
  {
    typedef utils::packed_vector_mapped<Data, Alloc> mapped_set_type;
  };
  
  
  template <typename Data, typename Alloc>
  struct __succinct_trie_mapped_data
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;
    typedef typename __succinct_trie_mapped_data_impl<Data, Alloc, boost::is_integral<Data>::value >::mapped_set_type impl_type;
    
    __succinct_trie_mapped_data() {}
    __succinct_trie_mapped_data(const path_type& path) : __impl(path) {}

    Data operator[](size_type pos) const { return __impl[pos]; }
    
    void open(const path_type& path) { __impl.open(path); }
    void clear() { __impl.clear(); }
    void close() { __impl.close(); }
    
    bool empty() const { return __impl.empty(); }
    size_type size() const { return __impl.size(); }
    path_type path() const { return __impl.path(); }
    
    uint64_t size_bytes() const { return __impl.size_bytes(); }
    uint64_t size_compressed() const { return __impl.size_compressed(); }
    uint64_t size_cache() const { return __impl.size_cache(); }
    
    impl_type __impl;
  };
  
  
  template <typename Key, typename Data, typename Alloc=std::allocator<std::pair<Key, Data> > >
  class succinct_trie_mapped : public __succinct_trie_base<Key,Data,Alloc>
  {
  public:
    typedef Key      key_type;
    typedef Data     data_type;
    typedef Data     mapped_type;
    typedef Alloc    allocator_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef boost::filesystem::path path_type;
    
  private:
    typedef succinct_trie_mapped<Key,Data,Alloc> self_type;
    typedef __succinct_trie_base<Key,Data,Alloc> base_type;

  public:
    typedef __succinct_trie_cursor<Key,Data,self_type,Alloc> cursor;
    typedef __succinct_trie_cursor<Key,Data,self_type,Alloc> const_cursor;

    typedef __succinct_trie_iterator<Key,Data,self_type,Alloc> iterator;
    typedef __succinct_trie_iterator<Key,Data,self_type,Alloc> const_iterator;
    
    typedef __succinct_trie_reverse_iterator<Key,Data,self_type,Alloc> reverse_iterator;
    typedef __succinct_trie_reverse_iterator<Key,Data,self_type,Alloc> const_reverse_iterator;
    
    using base_type::out_of_range;
    
  private:
    typedef typename Alloc::template rebind<uint32_t>::other  bit_alloc_type;
    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<data_type>::other data_alloc_type;
    
    typedef utils::succinct_vector_mapped<bit_alloc_type> position_set_type;
    typedef utils::succinct_vector_mapped<bit_alloc_type> index_map_type;
    //typedef utils::map_file<key_type, key_alloc_type>     index_set_type;
    typedef __succinct_trie_mapped_index<key_type, key_alloc_type> index_set_type;
    //typedef utils::map_file<data_type, data_alloc_type>   mapped_set_type;
    typedef __succinct_trie_mapped_data<data_type, data_alloc_type> mapped_set_type;

  public:
    succinct_trie_mapped() {}
    succinct_trie_mapped(const path_type& path) { open(path); }
    
  public:
    const_iterator begin(const size_type node_pos) const { return iterator(node_pos, *this); }
    const_iterator end(const size_type node_pos) const { return iterator(*this); }
    const_iterator begin() const { return iterator(0, *this); }
    const_iterator end() const { return iterator(*this); }
    
    const_reverse_iterator rbegin(const size_type node_pos) const { return reverse_iterator(node_pos, *this); }
    const_reverse_iterator rend(const size_type node_pos) const { return reverse_iterator(node_pos, *this); }
    const_reverse_iterator rbegin() const { return reverse_iterator(0, *this); }
    const_reverse_iterator rend() const { return reverse_iterator(0, *this); }
    
    const_cursor cbegin(const size_type node_pos) const { return cursor(node_pos, this); }
    const_cursor cend(const size_type node_pos) const { return cursor(); }
    const_cursor cbegin() const { return cursor(0, this); }
    const_cursor cend() const { return cursor(); }
    
  public:
    void read(const path_type& path) { open(path); }
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      close();

      repository_type rep(path, repository_type::read);
      positions.open(rep.path("positions"));
      index_map.open(rep.path("index-map"));
      index.open(rep.path("index"));
      mapped.open(rep.path("mapped"));
    }

    void write(const path_type& file) const
    {
      if (path() == file) return;
      
      // remove first...
      if (boost::filesystem::exists(file) && ! boost::filesystem::is_directory(file))
	boost::filesystem::remove_all(file);
      
      // create directory
      if (! boost::filesystem::exists(file))
	boost::filesystem::create_directories(file);
      
      // remove all the files...
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(file); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
      
      // copy all...
      for (boost::filesystem::directory_iterator iter(path()); iter != iter_end; ++ iter)
	utils::filesystem::copy_files(*iter, file);
    }
    
    bool empty() const { return mapped.empty(); }
    size_type size() const { return mapped.size(); }
    size_type index_size() const { return index.size(); }
    path_type path() const { return mapped.path().parent_path(); }
    
    uint64_t size_bytes() const
    { 
      return positions.size_bytes() + index_map.size_bytes() + index.size_bytes() + mapped.size_bytes();
    }
    uint64_t size_compressed() const
    {
      return positions.size_compressed() + index_map.size_compressed() + index.size_compressed() + mapped.size_compressed();
    }
    uint64_t size_cache() const
    {
      return positions.size_cache() + index_map.size_cache() + index.size_cache() + mapped.size_cache();
    }
    
    void close() { clear(); }
    void clear()
    {
      positions.clear();
      index_map.clear();
      index.clear();
      mapped.clear();
    }
    
    
    std::pair<size_type, size_type> range(const size_type node_pos) const
    {
      return base_type::__range(index, positions, node_pos);
    }
    
    size_type parent(size_type node_pos) const
    {
      return base_type::__parent(positions, node_pos);
    }
    bool is_next_sibling(size_type node_pos) const
    {
      return base_type::__is_next_sibling(positions, node_pos);
    }
    
    key_type key(size_type pos) const
    { 
      return base_type::__key(index, pos);
    }
    
    data_type data(size_type node_pos) const
    {
      return base_type::__data(index_map, mapped, node_pos);
    }
    
    size_type offset(size_type node_pos) const
    {
      return base_type::__offset(index_map, node_pos);
    }
    
    data_type operator[](size_type node_pos) const
    {
      return data(node_pos);
    }
    
    bool exists(size_type node_pos) const
    {
      return base_type::__exists(index_map, node_pos);
    }

    bool has_children(size_type node_pos) const
    {
      return base_type::__has_children(index, positions, node_pos);
    }
    
    bool is_valid(size_type node_pos) const
    {
      return node_pos != out_of_range();
    }
    
    size_type find(const key_type* key_buf, size_type key_size, size_type node_pos) const
    {
      size_type key_pos = 0;
      return traverse(key_buf, node_pos, key_pos, key_size);
    }

    size_type find(const key_type* key_buf, size_type key_size) const
    {
      size_type node_pos = 0;
      size_type key_pos = 0;
      return traverse(key_buf, node_pos, key_pos, key_size);
    }
    
    size_type traverse(const key_type* key_buf, size_type& node_pos, size_type& key_pos, size_type key_len) const
    {
      return base_type::__traverse(index, positions, key_buf, node_pos, key_pos, key_len);
    }
    
  private:
    position_set_type positions;
    index_map_type    index_map;
    index_set_type    index;
    mapped_set_type   mapped;
  };

  template <typename Key, typename Alloc, size_t KeySize>
  struct __succinct_trie_index_impl
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;

    typedef utils::packed_vector<Key, Alloc> mapped_set_type;
    typedef utils::packed_sink<Key, Alloc>   mapped_sink_type;

    typedef typename mapped_set_type::const_iterator const_iterator;
    
    __succinct_trie_index_impl() : __mapped() {}
    
    bool empty() const { return __mapped.empty(); }
    size_type size() const { return __mapped.size(); }
    
    uint64_t size_bytes() const { return __mapped.size_bytes(); }
    uint64_t size_compressed() const { return __mapped.size_compressed(); }
    uint64_t size_cache() const { return __mapped.size_cache(); }
    
    void clear() { __mapped.clear(); }
    
    void write(const path_type& path) const { __mapped.write(path); }
    
    void push_back(const Key& x) { __mapped.push_back(x); }
    void build() { __mapped.build(); }
    
    Key operator[](size_type pos) const { return __mapped[pos]; }
    
    const_iterator begin() const { return __mapped.begin(); }
    const_iterator end() const { return __mapped.end(); }
    
    mapped_set_type __mapped;
  };

  template <typename Key, typename Alloc>
  struct __succinct_trie_index_impl<Key,Alloc,1>
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;

    typedef std::vector<Key, Alloc>       mapped_set_type;
    typedef boost::iostreams::file_sink   mapped_sink_type;

    typedef typename mapped_set_type::const_iterator const_iterator;
    
    __succinct_trie_index_impl() : __mapped() {}
    
    bool empty() const { return __mapped.empty(); }
    size_type size() const { return __mapped.size(); }
    
    uint64_t size_bytes() const { return __mapped.size() * sizeof(Key); }
    uint64_t size_compressed() const { return __mapped.size() * sizeof(Key); }
    uint64_t size_cache() const { return 0; }
    
    void clear() { __mapped.clear(); }
    
    void write(const path_type& path)  const
    {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::file_sink(path.file_string()), 1024 * 1024);
      
      const int64_t file_size = sizeof(Key) * __mapped.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os.write(((char*) &(*__mapped.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("succinct trie: write()");
    }
    
    void push_back(const Key& x) { __mapped.push_back(x); }
    void build() {}
    
    Key operator[](size_type pos) const { return __mapped[pos]; }
    
    const_iterator begin() const { return __mapped.begin(); }
    const_iterator end() const { return __mapped.end(); }
    
    mapped_set_type __mapped;
  };
  

  template <typename Key, typename Alloc>
  struct __succinct_trie_index
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;
    typedef __succinct_trie_index_impl<Key, Alloc, sizeof(Key) >  impl_type;
    typedef typename impl_type::mapped_sink_type  sink_type;

    typedef typename impl_type::const_iterator const_iterator;
    
    __succinct_trie_index() : __impl() {}
    
    Key operator[](size_type pos) const { return __impl[pos]; }
    const_iterator begin() const { return __impl.begin(); }
    const_iterator end() const { return __impl.end(); }

    void push_back(const Key& x) { __impl.push_back(x); }
    void clear() { __impl.clear(); }
    void build() { __impl.build(); }
    
    bool empty() const { return __impl.empty(); }
    size_type size() const { return __impl.size(); }
    
    uint64_t size_bytes() const { return __impl.size_bytes(); }
    uint64_t size_compressed() const { return __impl.size_compressed(); }
    uint64_t size_cache() const { return __impl.size_cache(); }
    
    void write(const path_type& path) const { __impl.write(path); }
    
    impl_type __impl;
  };
  
  template <typename Data, typename Alloc, bool Integral>
  struct __succinct_trie_data_impl {};

  template <typename Data, typename Alloc>
  struct __succinct_trie_data_impl<Data,Alloc,false>
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;

    typedef std::vector<Data, Alloc>      mapped_set_type;
    typedef boost::iostreams::file_sink   mapped_sink_type;
    
    __succinct_trie_data_impl() : __mapped() {}
    
    bool empty() const { return __mapped.empty(); }
    size_type size() const { return __mapped.size(); }
    
    uint64_t size_bytes() const { return __mapped.size() * sizeof(Data); }
    uint64_t size_compressed() const { return __mapped.size() * sizeof(Data); }
    uint64_t size_cache() const { return 0; }
    
    void clear() { __mapped.clear(); }
    
    void write(const path_type& path)  const
    {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::file_sink(path.file_string()), 1024 * 1024);
      
      const int64_t file_size = sizeof(Data) * __mapped.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os.write(((char*) &(*__mapped.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("succinct trie: write()");
    }
    
    void push_back(const Data& x) { __mapped.push_back(x); }
    void build() {}
    
    Data operator[](size_type pos) const { return __mapped[pos]; }

    mapped_set_type __mapped;
  };
  
  template <typename Data, typename Alloc>
  struct __succinct_trie_data_impl<Data,Alloc,true>
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;

    typedef utils::packed_vector<Data, Alloc> mapped_set_type;
    typedef utils::packed_sink<Data, Alloc>   mapped_sink_type;
    
    __succinct_trie_data_impl() : __mapped() {}
    
    bool empty() const { return __mapped.empty(); }
    size_type size() const { return __mapped.size(); }
    
    uint64_t size_bytes() const { return __mapped.size_bytes(); }
    uint64_t size_compressed() const { return __mapped.size_compressed(); }
    uint64_t size_cache() const { return __mapped.size_cache(); }
    
    void clear() { __mapped.clear(); }
    
    void write(const path_type& path) const { __mapped.write(path); }
    
    void push_back(const Data& x) { __mapped.push_back(x); }
    void build() { __mapped.build(); }
    
    Data operator[](size_type pos) const { return __mapped[pos]; }
    
    mapped_set_type __mapped;
  };
  
  
  template <typename Data, typename Alloc>
  struct __succinct_trie_data
  {
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef boost::filesystem::path path_type;
    typedef __succinct_trie_data_impl<Data, Alloc, boost::is_integral<Data>::value >  impl_type;
    typedef typename impl_type::mapped_sink_type  sink_type;
    
    __succinct_trie_data() : __impl() {}

    Data operator[](size_type pos) const { return __impl[pos]; }
    void push_back(const Data& x) { __impl.push_back(x); }
    void clear() { __impl.clear(); }
    void build() { __impl.build(); }
    
    bool empty() const { return __impl.empty(); }
    size_type size() const { return __impl.size(); }
    
    uint64_t size_bytes() const { return __impl.size_bytes(); }
    uint64_t size_compressed() const { return __impl.size_compressed(); }
    uint64_t size_cache() const { return __impl.size_cache(); }
    
    void write(const path_type& path) const { __impl.write(path); }
    
    impl_type __impl;
  };


  template <typename Key, typename Data, typename Alloc=std::allocator<std::pair<Key, Data> > >
  class succinct_trie : public __succinct_trie_base<Key,Data,Alloc>
  {
  public:
    typedef Key      key_type;
    typedef Data     data_type;
    typedef Data     mapped_type;
    typedef Alloc    allocator_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef boost::filesystem::path path_type;
    
  private:
    typedef succinct_trie<Key,Data,Alloc> self_type;
    typedef __succinct_trie_base<Key,Data,Alloc> base_type;

  public:
    typedef __succinct_trie_cursor<Key,Data,self_type,Alloc> cursor;
    typedef __succinct_trie_cursor<Key,Data,self_type,Alloc> const_cursor;

    typedef __succinct_trie_iterator<Key,Data,self_type,Alloc> iterator;
    typedef __succinct_trie_iterator<Key,Data,self_type,Alloc> const_iterator;
    
    typedef __succinct_trie_reverse_iterator<Key,Data,self_type,Alloc> reverse_iterator;
    typedef __succinct_trie_reverse_iterator<Key,Data,self_type,Alloc> const_reverse_iterator;

    using base_type::out_of_range;
        
  private:
    typedef typename Alloc::template rebind<uint32_t>::other  bit_alloc_type;
    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<data_type>::other data_alloc_type;
    
    typedef utils::succinct_vector<bit_alloc_type>  position_set_type;
    typedef utils::succinct_vector<bit_alloc_type>  index_map_type;
    //typedef std::vector<key_type, key_alloc_type>   index_set_type;
    typedef __succinct_trie_index<key_type, key_alloc_type> index_set_type;
    //typedef std::vector<data_type, data_alloc_type> mapped_set_type;
    typedef __succinct_trie_data<data_type, data_alloc_type> mapped_set_type;
    
  public:
    const_iterator begin(const size_type node_pos) const { return iterator(node_pos, *this); }
    const_iterator end(const size_type node_pos) const { return iterator(*this); }
    const_iterator begin() const { return iterator(0, *this); }
    const_iterator end() const { return iterator(*this); }
    
    const_reverse_iterator rbegin(const size_type node_pos) const { return reverse_iterator(node_pos, *this); }
    const_reverse_iterator rend(const size_type node_pos) const { return reverse_iterator(node_pos, *this); }
    const_reverse_iterator rbegin() const { return reverse_iterator(0, *this); }
    const_reverse_iterator rend() const { return reverse_iterator(0, *this); }
    
    const_cursor cbegin(const size_type node_pos) const { return cursor(node_pos, this); }
    const_cursor cend(const size_type node_pos) const { return cursor(); }
    const_cursor cbegin() const { return cursor(0, this); }
    const_cursor cend() const { return cursor(); }
    
  public:
    bool empty() const { return mapped.empty(); }
    size_type size() const { return mapped.size(); }
    size_type index_size() const { return index.size(); }
    
    uint64_t size_bytes() const
    { 
      return positions.size_bytes() + index_map.size_bytes() + index.size() * sizeof(key_type) + mapped.size_bytes();
    }
    uint64_t size_compressed() const
    {
      return positions.size_compressed() + index_map.size_compressed() + index.size() * sizeof(key_type) + mapped.size_compressed();
    }
    uint64_t size_cache() const
    {
      return positions.size_cache() + index_map.size_cache() + mapped.size_cache();
    }

    void close() { clear(); }
    void clear()
    {
      positions.clear();
      index_map.clear();
      index.clear();
      mapped.clear();
    }
      
    std::pair<size_type, size_type> range(const size_type node_pos) const
    {
      return base_type::__range(index, positions, node_pos);
    }
    
    size_type parent(size_type node_pos) const
    {
      return base_type::__parent(positions, node_pos);
    }
    bool is_next_sibling(size_type node_pos) const
    {
      return base_type::__is_next_sibling(positions, node_pos);
    }
    
    key_type key(size_type pos) const
    { 
      return base_type::__key(index, pos);
    }
    
    data_type data(size_type node_pos) const
    {
      return base_type::__data(index_map, mapped, node_pos);
    }
    
    size_type offset(size_type node_pos) const
    {
      return base_type::__offset(index_map, node_pos);
    }
    
    data_type operator[](size_type node_pos) const
    {
      return data(node_pos);
    }
    
    bool exists(size_type node_pos) const
    {
      return base_type::__exists(index_map, node_pos);
    }

    bool has_children(size_type node_pos) const
    {
      return base_type::__has_children(index, positions, node_pos);
    }
    
    bool is_valid(size_type node_pos) const
    {
      return node_pos != out_of_range();
    }

    size_type find(const key_type* key_buf, size_type key_size, size_type node_pos) const
    {
      size_type key_pos = 0;
      return traverse(key_buf, node_pos, key_pos, key_size);
    }

    size_type find(const key_type* key_buf, size_type key_size) const
    {
      size_type node_pos = 0;
      size_type key_pos = 0;
      return traverse(key_buf, node_pos, key_pos, key_size);
    }

    size_type traverse(const key_type* key_buf, size_type& node_pos, size_type& key_pos, size_type key_len) const
    {
      return base_type::__traverse(index, positions, key_buf, node_pos, key_pos, key_len);
    }
    
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      rep["type"] = "succinct-trie";
      
      positions.write(rep.path("positions"));
      index_map.write(rep.path("index-map"));
      index.write(rep.path("index"));
      mapped.write(rep.path("mapped"));
    }
    
    // build!
  private:
    typedef typename Alloc::template rebind<size_type>::other size_alloc_type;
    typedef std::vector<size_type, size_alloc_type> terminal_set_type;

    struct Node
    {
      size_type depth;
      size_type first;
      size_type last;
      terminal_set_type terminals;
      
      key_type key;
      
      Node() : depth(), first(), last(), terminals(), key() {}
      Node(const size_type& _depth, const size_type _first, const size_type _last)
	: depth(_depth), first(_first), last(_last), terminals(), key() {}
    };
    typedef Node node_type;
    
    typedef typename Alloc::template rebind<node_type>::other node_alloc_type;
    typedef std::vector<node_type, node_alloc_type> node_set_type;

    typedef typename Alloc::template rebind<node_set_type>::other node_set_alloc_type;
    typedef std::deque<node_set_type, node_set_alloc_type> queue_type;
    
  private:
    template <typename _Tp, typename Index>
    struct __push_back_vector
    {
      __push_back_vector(Index& __index) : index(__index) {}
      
      void push_back(const _Tp& x)
      {
	index.push_back(x);
      }
      Index& index;
    };
    
    template <typename _Tp, typename Stream>
    struct __push_back_stream
    {
      __push_back_stream(Stream& __stream) : stream(__stream) {}
      
      void push_back(const _Tp& x)
      {
	if (! stream.write((char*) &x, sizeof(_Tp)))
	  throw std::runtime_error("sucinct trie: build()");
      }
      
      Stream& stream;
    };
    
  public:
    template <typename Iterator, typename ExtractKey, typename ExtractData>
    void build(const path_type& path, Iterator first, Iterator last, ExtractKey extract_key, ExtractData extract_data)
    {
      typedef utils::repository repository_type;
      
      // open and close..
      {
	repository_type rep(path, repository_type::write);
	rep["type"] = "succinct-trie";
      }
      
      repository_type rep(path, repository_type::read);
      
      {
	boost::iostreams::filtering_ostream os_index;
	boost::iostreams::filtering_ostream os_mapped;
	os_index.push(typename index_set_type::sink_type(rep.path("index").file_string()), 1024 * 1024);
	os_mapped.push(typename mapped_set_type::sink_type(rep.path("mapped").file_string()), 1024 * 1024);
	
	__push_back_stream<key_type, boost::iostreams::filtering_ostream> __index(os_index);
	__push_back_stream<data_type, boost::iostreams::filtering_ostream> __mapped(os_mapped);
	
	try {
	  __build(first, last, extract_key, extract_data, __index, __mapped);
	}
	catch (const std::exception& err) {
	  throw err;
	}
	catch (...) {
	  throw std::runtime_error("succinct_trie::build() memory error?");
	}
	os_index.pop();
	os_mapped.pop();
      }
      
      positions.write(rep.path("positions"));
      index_map.write(rep.path("index-map"));

      clear();
    }
    
    template <typename Iterator, typename ExtractKey>
    void build(const path_type& path, Iterator first, Iterator last, ExtractKey extract_key)
    {
      typedef utils::repository repository_type;
      
      // open and close...
      {
	repository_type rep(path, repository_type::write);
	rep["type"] = "succinct-trie";
      }
      repository_type rep(path, repository_type::read);
      
      {
	boost::iostreams::filtering_ostream os_index;
	boost::iostreams::filtering_ostream os_mapped;
	os_index.push(typename index_set_type::sink_type(rep.path("index").file_string()), 1024 * 1024);
	os_mapped.push(typename mapped_set_type::sink_type(rep.path("mapped").file_string()), 1024 * 1024);
	
	__push_back_stream<key_type, boost::iostreams::filtering_ostream> __index(os_index);
	__push_back_stream<data_type, boost::iostreams::filtering_ostream> __mapped(os_mapped);
	
	try {
	  __build(first, last, extract_key, __index, __mapped);
	}
	catch (const std::exception& err) {
	  throw err;
	}
	catch (...) {
	  throw std::runtime_error("succinct_trie::build() memory error?");
	}

	os_index.pop();
	os_mapped.pop();
      }
      
      positions.write(rep.path("positions"));
      index_map.write(rep.path("index-map"));

      clear();
    }

    
    template <typename Iterator, typename ExtractKey, typename ExtractData>
    void build(Iterator first, Iterator last, ExtractKey extract_key, ExtractData extract_data)
    {
      __push_back_vector<key_type, index_set_type> __index(index);
      __push_back_vector<data_type, mapped_set_type> __mapped(mapped);
      
      __build(first, last, extract_key, extract_data, __index, __mapped);
      
      index.build();
      mapped.build();
    }
    
    template <typename Iterator, typename ExtractKey>
    void build(Iterator first, Iterator last, ExtractKey extract_key)
    {
      __push_back_vector<key_type, index_set_type> __index(index);
      __push_back_vector<data_type, mapped_set_type> __mapped(mapped);
      
      __build(first, last, extract_key, __index, __mapped);
      
      index.build();
      mapped.build();
    }
    
    
  private:
    template <typename Iterator, typename ExtractKey, typename ExtractData, typename Index, typename Mapped>
    void __build(Iterator first, Iterator last, ExtractKey extract_key, ExtractData extract_data, Index& __index, Mapped& __mapped)
    {
      clear();
      
      const size_type data_size = last - first;
      
      queue_type queue;
      size_type index_size = 0;
      
      queue.push_front(node_set_type());
      children(node_type(0, 0, data_size), queue.front(), first, last, extract_key);
      
      __index.push_back(key_type());
      ++ index_size;
      
      // initial 01
      positions.set(0, true);
      positions.set(1, false);
      
      while (! queue.empty()) {
	const node_set_type& nodes = queue.back();
	
	size_type children_size = 0;
	
	typename node_set_type::const_iterator niter_begin = nodes.begin();
	typename node_set_type::const_iterator niter_end = nodes.end();
	for (typename node_set_type::const_iterator niter = niter_begin; niter != niter_end; ++ niter) {
	  // dump key and index...
	  
	  if (niter->terminals.empty()) {
	    __index.push_back(niter->key);
	    ++ index_size;
	    
	    queue.push_front(node_set_type());
	    children(*niter, queue.front(), first, last, extract_key);
	    
	    ++ children_size;
	  } else {
	    typename terminal_set_type::const_iterator titer_end = niter->terminals.end();
	    for (typename terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer) {
	      __index.push_back(niter->key);
	      ++ index_size;
	      __mapped.push_back(extract_data(*(first + *titer)));
	      index_map.set(index_size - 1, true);
	    }
	    
	    queue.push_front(node_set_type());
	    children(*niter, queue.front(), first, last, extract_key);
	    
	    for (int i = 0; i < int(niter->terminals.size()) - 1; ++ i)
	      queue.push_front(node_set_type());
	    
	    children_size += niter->terminals.size();
	  }
	}
	
	queue.pop_back();
	
	for (size_t i = 0; i != children_size; ++ i)
	  positions.set(positions.size(), true);
	positions.set(positions.size(), false);
      }
      
      index_map.build();
      positions.build();
    }
    
    template <typename Iterator, typename ExtractKey, typename Index, typename Mapped>
    void __build(Iterator first, Iterator last, ExtractKey extract_key, Index& __index, Mapped& __mapped)
    {
      clear();
      
      const size_type data_size = last - first;
      
      queue_type queue;
      size_type index_size = 0;
      
      queue.push_front(node_set_type());
      children(node_type(0, 0, data_size), queue.front(), first, last, extract_key);
      
      __index.push_back(key_type());
      ++ index_size;
      
      // initial 01
      positions.set(0, true);
      positions.set(1, false);
      
      while (! queue.empty()) {
	const node_set_type& nodes = queue.back();
	
	size_type children_size = 0;
	
	typename node_set_type::const_iterator niter_begin = nodes.begin();
	typename node_set_type::const_iterator niter_end = nodes.end();
	for (typename node_set_type::const_iterator niter = niter_begin; niter != niter_end; ++ niter) {
	  // dump key and index...
	  
	  if (niter->terminals.empty()) {
	    __index.push_back(niter->key);
	    ++ index_size;
	    
	    queue.push_front(node_set_type());
	    children(*niter, queue.front(), first, last, extract_key);
	    
	    ++ children_size;
	  } else {
	    typename terminal_set_type::const_iterator titer_end = niter->terminals.end();
	    for (typename terminal_set_type::const_iterator titer = niter->terminals.begin(); titer != titer_end; ++ titer) {
	      __index.push_back(niter->key);
	      ++ index_size;
	      __mapped.push_back(*titer);
	      index_map.set(index_size - 1, true);
	    }
	    
	    
	    queue.push_front(node_set_type());
	    children(*niter, queue.front(), first, last, extract_key);
	    
	    for (int i = 0; i < int(niter->terminals.size()) - 1; ++ i)
	      queue.push_front(node_set_type());
	    
	    children_size += niter->terminals.size();
	  }
	}
	
	queue.pop_back();
	
	for (size_t i = 0; i != children_size; ++ i)
	  positions.set(positions.size(), true);
	positions.set(positions.size(), false);
      }
      
      index_map.build();
      positions.build();
    }


  private:
    template <typename Iterator, typename ExtractKey>
    void children(const node_type& parent,
		  node_set_type&   nodes,
		  Iterator first,
		  Iterator last,
		  ExtractKey extract_key)
    {
      Iterator iter = first + parent.first;
      for (size_type pos = parent.first; pos != parent.last && iter != last; ++ pos, ++ iter)
	if (parent.depth < extract_key(*iter).size()) {
	  
	  const key_type& key = extract_key(*iter)[parent.depth];
	  
	  if (nodes.empty() || key != nodes.back().key) {
	    nodes.push_back(node_type());
	    nodes.back().depth = parent.depth + 1;
	    nodes.back().key = key;
	    nodes.back().first = pos;
	    nodes.back().terminals.clear();
	  }
	  if (parent.depth + 1 == extract_key(*iter).size())
	    nodes.back().terminals.push_back(pos);
	  nodes.back().last = pos + 1;
	}
      }
    
  private:
    position_set_type positions;
    index_map_type    index_map;
    index_set_type    index;
    mapped_set_type   mapped;
  };
  
};

#endif
