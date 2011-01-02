// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__CHUNK_VECTOR__HPP__
#define __UTILS__CHUNK_VECTOR__HPP__ 1

//
// chunk vector inspired by std::deque
//
// However, chunk-vector supports only push_back, but power-of-two optimized for
// faster access...
//

// TODO: revise this implementation, especially, interaction with iterator and uninitialized_fill or copy...
// remarks: What to do about the iterator...? use of map_pointer + offset?

#include <memory>
#include <vector>
#include <iterator>

// use of static-allocator...
#include <utils/static_allocator.hpp>
#include <utils/bithack.hpp>
#include <utils/memory.hpp>

namespace utils
{
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  struct __chunk_vector_iterator
  {
    typedef __chunk_vector_iterator<_Tp, _Tp&, _Tp*, _Size>                   iterator;
    typedef __chunk_vector_iterator<_Tp, const _Tp&, const _Tp*, _Size> const_iterator;
    
    typedef std::random_access_iterator_tag   iterator_category;
    
    typedef _Tp                           value_type;
    typedef _Ptr                          pointer;
    typedef _Ref                          reference;
    typedef size_t                        size_type;
    typedef ptrdiff_t                     difference_type;
    
    typedef _Tp**                         map_pointer;
    typedef __chunk_vector_iterator       self_type;
    
    static const size_type __chunk_size = _Size;
    static const size_type __chunk_mask = __chunk_size - 1;
    static const size_type __chunk_shift = utils::bithack::static_bit_count<__chunk_mask>::result;
    
    _Tp*        __curr;
    _Tp*        __first;
    map_pointer __node;
    
    __chunk_vector_iterator(const iterator& x)
      : __curr(x.__curr), __first(x.__first), __node(x.__node) {}
	
    __chunk_vector_iterator(_Tp* __x, map_pointer __y)
      : __curr(__x), __first(*__y), __node(__y) { }
    __chunk_vector_iterator()
      :  __curr(0), __first(0), __node(0) {}
    
  public:
    reference operator*() const { return *__curr; }
    pointer   operator->() const { return __curr; }
    
    // operator++
    self_type& operator++()
    {
      ++ __curr;
      if (__curr == __first + __chunk_size) {
	set_node(__node + 1);
	__curr = __first;
      }
      return *this;
    }
    
    self_type operator++(int)
    {
      self_type __tmp = *this;
      ++ *this;
      return __tmp;
    }
    
    // operator--
    self_type& operator--()
    {
      if (__curr == __first) {
	set_node(__node - 1);
	__curr = __first + __chunk_size;
      }
      -- __curr;
      return *this;
    }
    
    self_type operator--(int)
    {
      self_type __tmp = *this;
      -- *this;
      return __tmp;
    }

    // operator+=
    self_type& operator+=(difference_type __n)
    {
      const difference_type __offset = __n + (__curr - __first);
      if (__offset >= 0 && __offset < difference_type(__chunk_size))
	__curr += __n;
      else {	
	if (__offset > 0) {
	  // power-of-two optimization...
	  set_node(__node + (__offset >> __chunk_shift));
	  __curr = __first + (__offset & __chunk_mask);
	} else {
	  const difference_type __node_offset = - difference_type((- __offset - 1) / __chunk_size) - 1;
	  set_node(__node + __node_offset);
	  __curr = __first + (__offset - __node_offset * difference_type(__chunk_size));
	}
      }
      return *this;
    }
    
    self_type operator+(difference_type __n) const
    {
      self_type __tmp = *this;
      return __tmp += __n;
    }
    
    // operator-=
    self_type& operator-=(difference_type __n)
    {
      return *this += -__n;
    }
    
    self_type operator-(difference_type __n) const
    {
      self_type __tmp = *this;
      return __tmp -= __n;
    }
    
    reference operator[](difference_type __n) const
    {
      return *(*this + __n);
    }
    
    void set_node(map_pointer __new_node)
    {
      __node = __new_node;
      __first = *__new_node;
    }
  };

  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline bool
  operator==(const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& y)
  {
    return x.__curr == y.__curr;
  }
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline bool
  operator==(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    return x.__curr == y.__curr;
  }
  
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline bool
  operator!=(const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& y)
  {
    return ! (x.__curr == y.__curr);
  }
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline bool
  operator!=(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    return ! (x.__curr == y.__curr);
  }
  
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline bool
  operator<(const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& x,
	    const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& y)
  {
    return (x.__node < y.__node || (! (y.__node < x.__node) && x.__curr < y.__curr));
  }
  
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline bool
  operator<(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	    const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    return (x.__node < y.__node || (! (y.__node < x.__node) && x.__curr < y.__curr));
  }
  
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline bool
  operator>(const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& x,
	    const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& y)
  {
    return y < x;
  }
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline bool
  operator>(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	    const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    return y < x;
  }
  
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline bool
  operator<=(const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& y)
  {
    return ! (y < x);
  }
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline bool
  operator<=(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    return ! (y < x);
  }
  
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline bool
  operator>=(const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& y)
  {
    return ! (x < y);
  }
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline bool
  operator>=(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	     const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    return ! (x < y);
  }
  
  
  template <typename _Tp, typename _RefL, typename _PtrL, typename _RefR, typename _PtrR, size_t _Size>
  inline typename __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>::difference_type
  operator-(const __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size>& x,
	    const __chunk_vector_iterator<_Tp,_RefR,_PtrR,_Size>& y)
  {
    typedef __chunk_vector_iterator<_Tp,_RefL,_PtrL,_Size> iterator;

    return typename iterator::difference_type
      (((x.__node - y.__node - 1) * iterator::__chunk_size)
       + (x.__curr - x.__first)
       + ((y.__first + iterator::__chunk_size) - y.__curr));
  }
  
  template <typename _Tp, typename _Ref, typename _Ptr, size_t _Size>
  inline __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>
  operator+(ptrdiff_t __n, const __chunk_vector_iterator<_Tp,_Ref,_Ptr,_Size>& __x)
  {
    return __x + __n;
  }
  
  
  
  template <typename _Tp, size_t _ChunkSize=128, typename _Alloc=std::allocator<_Tp> >
  struct __chunk_vector_alloc_impl
  {
    typedef _Tp* pointer;
    
    static const bool   __is_power2 = utils::bithack::static_is_power2<_ChunkSize>::result;
    static const size_t __next_power2 = size_t(utils::bithack::static_next_largest_power2<_ChunkSize>::result);
    static const size_t chunk_size = (_ChunkSize == 0 ? size_t(1) : (__is_power2 ? _ChunkSize : __next_power2));
    
    typedef utils::static_allocator<_Tp, chunk_size, _Alloc>   node_alloc_type;
    // use raw allocator for debug purpose...
    // typedef _Alloc   node_alloc_type;
    typedef typename _Alloc::template rebind<_Tp*>::other map_alloc_type;
  };
  

  template <typename _Tp, size_t _ChunkSize, typename _Alloc>
  struct __chunk_vector_impl : public __chunk_vector_alloc_impl<_Tp, _ChunkSize, _Alloc>::node_alloc_type,
			       public __chunk_vector_alloc_impl<_Tp, _ChunkSize, _Alloc>::map_alloc_type
  {
    // we will concentrate on memory management, not construction/destruction

  public:
    
    typedef __chunk_vector_alloc_impl<_Tp, _ChunkSize, _Alloc> alloc_impl_type;

    typedef _Tp        value_type;
    typedef size_t     size_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    
    typedef _Tp**      map_pointer;
    
    typedef typename alloc_impl_type::node_alloc_type node_alloc_type;
    typedef typename alloc_impl_type::map_alloc_type  map_alloc_type;
    
    static const size_t chunk_size = alloc_impl_type::chunk_size;
    static const size_t chunk_mask = chunk_size - 1;
    static const size_t chunk_shift = utils::bithack::static_bit_count<chunk_mask>::result;
    
    typedef __chunk_vector_iterator<_Tp, _Tp&, _Tp*, chunk_size>             iterator;
    typedef __chunk_vector_iterator<_Tp, const _Tp&, const _Tp*, chunk_size> const_iterator;
    
    node_alloc_type& node_allocator() { return static_cast<node_alloc_type&>(*this); }
    map_alloc_type& map_allocator() { return static_cast<map_alloc_type&>(*this); }
    
    __chunk_vector_impl() : __map(0), __node_size(0) {}
    __chunk_vector_impl(size_type x) : __map(0), __node_size(0) { __initialize_map(x); }
    ~__chunk_vector_impl()
    {
      if (__map) {
	const size_type __map_size = (__node_size >> chunk_shift) + 1;
	__destroy_nodes(__map, __map + __map_size);
	map_allocator().deallocate(__map, __map_size);
      }
    }
    
  protected:
    
    void __initialize_map(size_type x)
    {
      const size_type __map_size = (x >> chunk_shift) + 1;
      
      __node_size = x;
      __map = map_allocator().allocate(__map_size);
      
      try {
	__create_nodes(__map, __map + __map_size);
      }
      catch (...) {
	map_allocator().deallocate(__map, __map_size);
	__map = 0;
	__node_size = 0;
	throw;
      }
    }
    
    void __reallocate_map(size_type x)
    {
      const size_type __map_size_curr = (__node_size >> chunk_shift) + 1;
      const size_type __map_size_new = (x >> chunk_shift) + 1;
      
      if (__map_size_new > __map_size_curr) {
	map_pointer __map_new  = map_allocator().allocate(__map_size_new);
	try {
	  __create_nodes(__map_new + __map_size_curr, __map_new + __map_size_new);
	}
	catch (...) {
	  map_allocator().deallocate(__map_new, __map_size_new);
	  throw;
	}
	std::copy(__map, __map + __map_size_curr, __map_new);
	
	map_allocator().deallocate(__map, __map_size_curr);
	__map = __map_new;
      } else if (__map_size_new < __map_size_curr) {
	map_pointer __map_new = map_allocator().allocate(__map_size_new);
	std::copy(__map, __map + __map_size_new, __map_new);
	
	__destroy_nodes(__map + __map_size_new, __map + __map_size_curr);
	
	map_allocator().deallocate(__map, __map_size_curr);
	__map = __map_new;
      }
      
      __node_size = x;
    }
    
    void __create_nodes(map_pointer first, map_pointer last)
    {
      map_pointer curr;
      try {
	for (curr = first; curr < last; ++ curr)
	  *curr = node_allocator().allocate(chunk_size);
      }
      catch (...) {
	__destroy_nodes(first, curr);
	throw;
      }
    }
    
    void __destroy_nodes(map_pointer first, map_pointer last)
    {
      for (/**/; first < last; ++ first)
	node_allocator().deallocate(*first, chunk_size);
    }
    
    
  protected:
    // actual storage... very minimum...
    map_pointer __map;
    size_type   __node_size;
  };
  
  template <typename _Tp, size_t _ChunkSize=128, typename _Alloc=std::allocator<_Tp> >
  class chunk_vector : protected __chunk_vector_impl<_Tp, _ChunkSize, _Alloc>
  {
  protected:
    typedef __chunk_vector_impl<_Tp, _ChunkSize, _Alloc> impl_type;
    
  public:
    typedef _Tp                                   value_type;
    typedef typename impl_type::pointer           pointer;
    typedef typename impl_type::const_pointer     const_pointer;
    typedef typename impl_type::reference         reference;
    typedef typename impl_type::const_reference   const_reference;
    typedef typename impl_type::iterator          iterator;
    typedef typename impl_type::const_iterator    const_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef size_t                                size_type;
    typedef ptrdiff_t                             difference_type;
    
  private:
    using impl_type::__map;
    using impl_type::__node_size;
    
    using impl_type::__initialize_map;
    using impl_type::__reallocate_map;
    
  public:
    chunk_vector() : impl_type(0) {}
    chunk_vector(size_type x) : impl_type(x) { __fill_initialize(value_type()); }
    chunk_vector(size_type x, const value_type& __value) : impl_type(x) { __fill_initialize(__value); }
    chunk_vector(const chunk_vector& x) : impl_type(x.size()) { __copy_initialize(x); }
    
    // note that we do not assign anything...
    template <typename Iterator>
    chunk_vector(Iterator first, Iterator last) : impl_type()
    {
      typedef typename boost::is_integral<Iterator>::type __integral;
      __initialize_dispatch(first, last, __integral());
    }
    
    ~chunk_vector() throw() { __destroy_range(begin(), end()); }
    
    chunk_vector& operator=(const chunk_vector& x)
    {
      __copy_assign(x);
      return *this;
    }
    
  public:
    void assign(const chunk_vector& x)
    {
      __copy_assign(x);
    }
    
    void assign(size_type x, const value_type& __value)
    {
      __fill_assign(x, __value);
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      typedef typename boost::is_integral<Iterator>::type __integral;
      __assign_dispatch(first, last, __integral());
    }
    
  public:

    inline bool empty() const { return __node_size == 0; }
    inline size_type size() const { return __node_size; }
    
    inline const_iterator begin() const { return const_iterator(*__map, __map); }
    inline       iterator begin()       { return iterator(*__map, __map); }
    inline const_iterator end() const { return __iterator(__node_size); }
    inline       iterator end()       { return __iterator(__node_size); }
    
    inline const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    inline       reverse_iterator rbegin()       { return reverse_iterator(end()); }
    inline const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
    inline       reverse_iterator rend()       { return reverse_iterator(begin()); }
    
    
    inline const_reference front() const { return *(*__map); }
    inline       reference front()       { return *(*__map); }
    
    inline const_reference back() const { return operator[](__node_size - 1); }
    inline       reference back()       { return operator[](__node_size - 1); }
    
    inline const_reference operator[](size_type __x) const 
    {
      return *(__map[__x >> impl_type::chunk_shift] + (__x & impl_type::chunk_mask));
    }
    inline       reference operator[](size_type __x)
    {
      return *(__map[__x >> impl_type::chunk_shift] + (__x & impl_type::chunk_mask));
    }
    
  public:
    
    void push_back(const _Tp& __x)
    {
      const size_type __map_back = __node_size >> impl_type::chunk_shift;
      const size_type __node_back = __node_size & impl_type::chunk_mask;
      
      __reallocate_map(__node_size + 1);
      utils::construct_object(*(__map + __map_back) + __node_back, __x);
    }

    void pop_back()
    {
      const size_type __map_back = (__node_size - 1) >> impl_type::chunk_shift;
      const size_type __node_back = (__node_size - 1) & impl_type::chunk_mask;
      
      utils::destroy_object(*(__map + __map_back) + __node_back);
      __reallocate_map(__node_size - 1);
    }
    
    
    void clear()
    {
      __destroy_range(begin(), end());
      __reallocate_map(0);
    }
    
    
    
    
    void resize(size_type x) { resize(x, value_type()); }
    void resize(size_type x, const value_type& __value)
    {
      if (x < __node_size) {
	__destroy_range(__iterator(x), __iterator(__node_size));
	__reallocate_map(x);
      } else if (x > __node_size) {
	const size_type __size_old = __node_size;
	__reallocate_map(x);
	__uninitialized_fill(__iterator(__size_old), end(), __value);
      }
    }
    
    void swap(chunk_vector& x)
    {
      std::swap(__map, x.__map);
      std::swap(__node_size, x.__node_size);
    }
    
  private:
    
    const_iterator __iterator(size_type x) const
    {
      const size_type __map_pos = x >> impl_type::chunk_shift;
      const size_type __node_pos = x & impl_type::chunk_mask;
      
      return const_iterator(*(__map + __map_pos) + __node_pos, __map + __map_pos);
    }
    iterator __iterator(size_type x)
    {
      const size_type __map_pos = x >> impl_type::chunk_shift;
      const size_type __node_pos = x & impl_type::chunk_mask;
      
      return iterator(*(__map + __map_pos) + __node_pos, __map + __map_pos);
    }
    
    // assignment...
    template <typename Integer>
    void __assign_dispatch(Integer n, Integer x, boost::true_type)
    {
      __fill_assign(n, x);
    }

    template <typename _InputIterator>
    void __assign_dispatch(_InputIterator first, _InputIterator last, boost::false_type)
    {
      typedef typename std::iterator_traits<_InputIterator>::iterator_category __category;
      __range_assign(first, last, __category());
    }
    
    void __fill_assign(size_type x, const value_type& __value)
    {
      __fill(begin(), __iterator(std::min(x, __node_size)), __value);
      resize(x, __value);
    }
    
    void __copy_assign(const chunk_vector& x)
    {
      const size_type copy_size = std::min(size(), x.size());
      
      __copy(x.begin(), x.__iterator(copy_size), begin());
      
      if (x.size() < __node_size) {
	__destroy_range(__iterator(copy_size), end());
	__reallocate_map(x.size());
      } else if (x.size() > __node_size) {
	__reallocate_map(x.size());
	__uninitialized_copy(x.__iterator(copy_size), x.end(), __iterator(copy_size));
      }
    }
    
    template <typename _InputIterator>
    void __range_assign(_InputIterator first, _InputIterator last, std::input_iterator_tag)
    {
      clear();
      
      // push-back...
      try {
	for (/**/; first != last; ++ first)
	  push_back(*first);
      }
      catch (...) {
	clear();
	throw;
      }
    }
    
    template <typename _ForwardIterator>
    void __range_assign(_ForwardIterator first, _ForwardIterator last, std::forward_iterator_tag)
    {
      const size_type __n = std::distance(first, last);
      const size_type copy_size = std::min(__node_size, __n);
      
      // trivial copy...
      // we will introducce specialized copy...
      _ForwardIterator mid = first;
      std::advance(mid, copy_size);
      __copy(first, mid, begin());
      first = mid;
      
      if (__n < __node_size) {
	__destroy_range(__iterator(__n), __iterator(__node_size));
	__reallocate_map(__n);
      } else if (__n > __node_size) {	
	const size_type __size_old = __node_size;
	__reallocate_map(__n);
	__uninitialized_copy(first, last, __iterator(__size_old));
      }
    }

    // initialize...
    
    template <typename Integer>
    void __initialize_dispatch(Integer n, Integer x, boost::true_type)
    {
      __initialize_map(n);
      __fill_initialize(x);
    }
    template <typename _InputIterator>
    void __initialize_dispatch(_InputIterator first, _InputIterator last, boost::false_type)
    {
      typedef typename std::iterator_traits<_InputIterator>::iterator_category __category;
      __range_initialize(first, last, __category());
    }
    
    template <typename _InputIterator>
    void __range_initialize(_InputIterator first, _InputIterator last, std::input_iterator_tag)
    {
      __initialize_map(0);
      // push-back...
      try {
	for (/**/; first != last; ++ first)
	  push_back(*first);
      }
      catch (...) {
	clear();
	throw;
      }
    }
    
    template <typename _ForwardIterator>
    void __range_initialize(_ForwardIterator first, _ForwardIterator last, std::forward_iterator_tag)
    {
      // TODO: exception safety...
      const size_type __n = std::distance(first, last);
      __initialize_map(__n);
      
      iterator __curr = begin();
      iterator __last = end();

      try {
	for (/**/; __curr.__node < __last.__node; ++ __curr.__node) {
	  _ForwardIterator mid = first;
	  std::advance(mid, impl_type::chunk_size);
	  std::uninitialized_copy(first, mid, *(__curr.__node));
	  first = mid;
	}
	std::uninitialized_copy(first, last, *(__last.__node));
      }
      catch (...) {
	__destroy_range(begin(), iterator(*(__curr.__node), __curr.__node));
	throw;
      }
    }
    
    void __copy_initialize(const chunk_vector& x)
    {
      __uninitialized_copy(x.begin(), x.end(), begin());
    }
    
    void __fill_initialize(const value_type& __value)
    {
      __uninitialized_fill(begin(), end(), __value);
    }
    
    // auxiliary functions...
    // we assume forward-iterator...
    // we also assume that necessary memory region is already allocated...
    template <typename Iterator>
    void __copy(Iterator first, Iterator last, iterator iter)
    {
      size_type __n = std::distance(first, last);
      while (__n) {
	const difference_type dist = std::min(__n, (iter.__first + impl_type::chunk_size) - iter.__curr);
	Iterator mid = first;
	std::advance(mid, dist);
	
	// actual copy
	std::copy(first, mid, iter.__curr);
	
	// advance
	first = mid;
	__n -= dist;
	iter += dist;
      }
    }
    
    template <typename Iterator>
    void __uninitialized_copy(Iterator first, Iterator last, iterator iter)
    {
      size_type __n = std::distance(first, last);
      while (__n) {
	const difference_type dist = std::min(__n, (iter.__first + impl_type::chunk_size) - iter.__curr);
	Iterator mid = first;
	std::advance(mid, dist);
	
	// actual copy
	std::uninitialized_copy(first, mid, iter.__curr);
	
	// advance
	first = mid;
	__n -= dist;
	iter += dist;
      }
    }
    
    void __copy(const_iterator first, const_iterator last, iterator iter)
    {
      if (first.__node == last.__node) {
	std::copy(first.__curr, last.__curr, iter.__curr);
      } else {
	std::copy(first.__curr, first.__first + impl_type::chunk_size, iter.__curr);
	
	++ first.__node;
	++ iter.__node;
	
	for (/**/; first.__node < last.__node; ++ first.__node, ++ iter.__node)
	  std::copy(*(first.__node), *(first.__node) + impl_type::chunk_size, *(iter.__node));
	std::copy(last.__first, last.__curr, *(iter.__node));
      }
    }
    
    void __uninitialized_copy(const_iterator first, const_iterator last, iterator iter)
    {
      if (first.__node == last.__node) {
	std::uninitialized_copy(first.__curr, last.__curr, iter.__curr);
      } else {
	std::uninitialized_copy(first.__curr, first.__first + impl_type::chunk_size, iter.__curr);
	
	++ first.__node;
	++ iter.__node;
	
	for (/**/; first.__node < last.__node; ++ first.__node, ++ iter.__node)
	  std::uninitialized_copy(*(first.__node), *(first.__node) + impl_type::chunk_size, *(iter.__node));
	std::uninitialized_copy(last.__first, last.__curr, *(iter.__node));
      }
    }
    
    void __fill(iterator first, iterator last, const value_type& __value)
    {
      // TODO: support for safer exception...
      if (first.__node == last.__node)  {
	std::fill(first.__curr, last.__curr, __value);
      } else {
	std::fill(first.__curr, first.__first + impl_type::chunk_size, __value);
	
	++ first.__node;
	
	for (/**/; first.__node < last.__node; ++ first.__node)
	  std::fill(*(first.__node), *(first.__node) + impl_type::chunk_size, __value);
	
	std::fill(last.__first, last.__curr, __value);
      }
    }
    
    void __uninitialized_fill(iterator first, iterator last, const value_type& __value)
    {
      // TODO: support for safer exception...
      
      if (first.__node == last.__node) {
	std::uninitialized_fill(first.__curr, last.__curr, __value);
      } else {
	std::uninitialized_fill(first.__curr, first.__first + impl_type::chunk_size, __value);
	
	++ first.__node;
	
	for (/**/; first.__node < last.__node; ++ first.__node)
	  std::uninitialized_fill(*(first.__node), *(first.__node) + impl_type::chunk_size, __value);
	std::uninitialized_fill(last.__first, last.__curr, __value);
      }
    }

    
    void __destroy_range(iterator first, iterator last)
    {
      __destroy_range_aux(first, last, boost::has_trivial_destructor<value_type>());
    }
    
    void __destroy_range_aux(iterator first, iterator last, const boost::true_type&)
    {
      
    }
    
    void __destroy_range_aux(iterator first, iterator last, const boost::false_type&)
    {
      if (first.__node == last.__node)
	utils::destroy_range(first.__curr, last.__curr);
      else {
	utils::destroy_range(first.__curr, first.__first + impl_type::chunk_size);
	
	++ first.__node;
	
	for (/**/; first.__node < last.__node; ++ first.__node)
	  utils::destroy_range(*(first.__node), *(first.__node) + impl_type::chunk_size);
	utils::destroy_range(last.__first, last.__curr);
      }
    }
  };
};

namespace std
{
  // specialization for std::swap...
  
  template <typename _Tp, size_t _ChunkSize, typename _Alloc >
  inline
  void swap(utils::chunk_vector<_Tp, _ChunkSize, _Alloc>& x,
	    utils::chunk_vector<_Tp, _ChunkSize, _Alloc>& y) 
  {
    x.swap(y);
  }
  
};

#endif
