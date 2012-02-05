// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SMALL_VECTOR__HPP__
#define __UTILS__SMALL_VECTOR__HPP__ 1

#include <memory>
#include <iterator>
#include <vector>
#include <stdexcept>

#include <boost/type_traits.hpp>

#include <utils/memory.hpp>
#include <utils/bithack.hpp>

namespace utils {
  
  
  template <typename Tp, typename _Alloc>
  struct __small_vector_base : public _Alloc
  {
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;
    
    typedef       Tp  value_type;
    typedef       Tp* iterator;
    typedef const Tp* const_iterator;
    typedef       Tp& reference;
    typedef const Tp& const_reference;
    typedef       Tp* pointer;
    
    typedef _Alloc allocator_type;
    
    static const size_type small_threshold = sizeof(pointer) / sizeof(Tp);
    
    __small_vector_base() : __base(0), __size(0) {}
    __small_vector_base(size_type __n) : __base(0), __size(0) { initialize_vector(__n); }
    ~__small_vector_base()
    {
      if (__size > small_threshold && __base)
	allocator().deallocate(__base, capacity(__size));
    }
    
    void swap(__small_vector_base& x)
    {
      std::swap(__base, x.__base);
      std::swap(__size, x.__size);
    }
    
    inline const_iterator begin() const { return (__size > small_threshold ? __base : reinterpret_cast<const_iterator>(&__base)); }
    inline       iterator begin() { return (__size > small_threshold ? __base : reinterpret_cast<iterator>(&__base)); }
    inline const_iterator end() const { return begin() + __size; }
    inline       iterator end()       { return begin() + __size; }

    static inline size_type capacity(const size_type& size)
    {
      const size_t power2 = bithack::branch(bithack::is_power2(size),
					    size,
					    static_cast<size_type>(bithack::next_largest_power2(size)));
      
      const size_t size_alloc        = size * sizeof(value_type);
      const size_t size_power2_alloc = power2 * sizeof(value_type);
      const size_t size_256          = size_t(256) / sizeof(value_type);
      
      return bithack::branch(size_alloc > 256 || size == 0, size, bithack::branch(size_power2_alloc > 256, size_256, power2));
    }
    
    bool empty() const { return ! __size; }
    size_type size() const { return __size; }
    size_type& size() { return __size; }
    size_type capacity() const { return capacity(__size); }
    
    void initialize_vector(size_type __n)
    {
      __size = __n;
      if (__size > small_threshold)
	__base = allocator().allocate(capacity(__size));
    }
    
    const allocator_type& allocator() const { return static_cast<allocator_type&>(*this); }
    allocator_type& allocator() { return static_cast<allocator_type&>(*this); }
    
  private:
    pointer   __base;
    size_type __size;
  };
  
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class small_vector
  {
  private:
    typedef __small_vector_base<Tp, Alloc> base_type;

  public:
    typedef Tp                                  value_type;
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;
    
    typedef typename base_type::reference       reference;
    typedef typename base_type::const_reference const_reference;

    typedef typename base_type::iterator       iterator;
    typedef typename base_type::const_iterator const_iterator;

    typedef       Tp* pointer;
    typedef const Tp* const_pointer;
    
    typedef typename std::reverse_iterator<iterator>       reverse_iterator;
    typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;
    
  public:
    small_vector() : __base(0) {}
    small_vector(size_type __n) : __base(__n) { std::uninitialized_fill(__base.begin(), __base.end(), value_type()); }
    small_vector(const size_type __n, const Tp& __value) : __base(__n) { std::uninitialized_fill(__base.begin(), __base.end(), __value); }
    
    template <typename Iterator>
    small_vector(Iterator first, Iterator last) : __base() {
      typedef typename boost::is_integral<Iterator>::type __integral;
      __initialize_dispatch(first, last, __integral());
    }
    small_vector(const small_vector& x) : __base(x.size()) { std::uninitialized_copy(x.begin(), x.end(), __base.begin()); }
    
    ~small_vector() { utils::destroy_range(__base.begin(), __base.end()); }
    
    small_vector& operator=(const small_vector& x)
    {
      __range_assign(x.begin(), x.end());
      return *this;
    }
    
  public:
    void assign(size_type __n, const Tp& __value) { __fill_assign(__n, __value); }
    void assign(const small_vector& x) { __range_assign(x.begin(), x.end()); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    { 
      typedef typename boost::is_integral<Iterator>::type __integral;
      __assign_dispatch(first, last, __integral());
    }
    
  public:
    bool empty() const { return __base.empty(); }
    size_type size() const { return __base.size(); }
    
    reference       operator[](size_type __n)       { return *(begin() + __n); }
    const_reference operator[](size_type __n) const { return *(begin() + __n); }
    
    reference front() { return *(begin()); }
    const_reference front() const { return *(begin()); }
    reference back() { return *(end() - 1); }
    const_reference back() const { return *(end() - 1); }
    
    const_iterator begin() const { return __base.begin(); }
    iterator begin() { return __base.begin(); }
    
    const_iterator end() const { return __base.end(); }
    iterator end() { return __base.end(); }
    
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
    
  public:
    void swap(small_vector& x) { __base.swap(x.__base); }
    
    void clear() { resize(0); }
    void resize(size_type __n) { resize(__n, value_type()); }
    void resize(size_type __n, const Tp& __value) {
      const size_type __size = size();
      
      if (__n < __size) {
	if (base_type::capacity(__n) == __base.capacity()) {
	  utils::destroy_range(__base.begin() + __n, __base.begin() + __size);
	  __base.size() = __n;
	} else {
	  base_type __base_new(__n);
	  std::uninitialized_copy(__base.begin(), __base.begin() + __n, __base_new.begin());
	  __base.swap(__base_new);
	  utils::destroy_range(__base_new.begin(), __base_new.end());
	}
      } else if (__n > __size) {
	if (base_type::capacity(__n) == __base.capacity()) {
	  std::uninitialized_fill(__base.begin() + __size, __base.begin() + __n, __value);
	  __base.size() = __n;
	} else {
	  base_type __base_new(__n);
	  std::uninitialized_copy(__base.begin(), __base.end(), __base_new.begin());
	  std::uninitialized_fill(__base_new.begin() + __size, __base_new.end(), __value);
	  __base.swap(__base_new);
	  utils::destroy_range(__base_new.begin(), __base_new.end());
	}
      }
    }
    
    void pop_back()
    {
      // nochecking...
      resize(size() - 1);
    }

    void push_back(const Tp& __value)
    {
      const size_type __size = size();
      
      if (base_type::capacity(__size + 1) == __base.capacity()) {
	utils::construct_object(__base.begin() + __size, __value);
	++ __base.size();
      } else {
	base_type __base_new(__size + 1);
	std::uninitialized_copy(__base.begin(), __base.end(), __base_new.begin());
	utils::construct_object(__base_new.begin() + __size, __value);
	
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
    }
    
    iterator insert(iterator position, const Tp& x)
    {
      const size_type __n = position - begin();
      const size_type __size = size();

      if (base_type::capacity(__size + 1) == __base.capacity()) {
	if (__n == __size) {
	  utils::construct_object(position, x);
	  ++ __base.size();
	} else {
	  utils::construct_object(end(), back());
	  ++ __base.size();
	  std::copy_backward(position, end() - 2, end() - 1);
	  *position = x;
	}
      } else {
	base_type __base_new(__size + 1);
	std::uninitialized_copy(begin(), position, __base_new.begin());
	utils::construct_object(__base_new.begin() + __n, x);
	std::uninitialized_copy(position, end(), __base_new.begin() + __n + 1);
	
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
      
      return begin() + __n;
    }

    template <typename Iterator>
    void insert(iterator position, Iterator first, Iterator last)
    {
      typedef typename std::iterator_traits<Iterator>::iterator_category __category;
      __range_insert(position, first, last, __category());
    }

    template <typename Iterator>
    void __range_insert(iterator position, Iterator first, Iterator last, std::forward_iterator_tag)
    {
      const size_type __n = position - begin();
      const size_type __d = std::distance(first, last);
      const size_type __size = size();
      
      if (__n == __size && base_type::capacity(__size + __d) == __base.capacity()) {
	std::uninitialized_copy(first, last, __base.end());
	__base.size() += __d;
      } else {
	base_type __base_new(__size + __d);
	std::uninitialized_copy(begin(), position, __base_new.begin());
	std::uninitialized_copy(first, last, __base_new.begin() + __n);
	std::uninitialized_copy(position, end(), __base_new.begin() + __n + __d);
	
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
    }

    template <typename Iterator>
    void __range_insert(iterator position, Iterator first, Iterator last, std::input_iterator_tag)
    {
      std::vector<Tp, Alloc> vec(first, last);
      __range_insert(position, vec.begin(), vec.end(), std::forward_iterator_tag());
    }

    iterator erase(iterator position)
    {
      const size_type __n = position - begin();
      const size_type __size = size();

      if (__base.capacity() == base_type::capacity(__size - 1)) {
	if (__n + 1 != __size)
	  std::copy(position + 1, end(), position);
	utils::destroy_object(end() - 1);
	-- __base.size();
      } else {
	base_type __base_new(__size - 1);
	std::uninitialized_copy(begin(), position, __base_new.begin());
	std::uninitialized_copy(position + 1, end(), __base_new.begin() + __n);
	
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
      
      return begin() + __n;
    }

    iterator erase(iterator first, iterator last)
    {
      const size_type __n = first - begin();
      const size_type __d = last - first;
      const size_type __size = size();

      if (__base.capacity() == base_type::capacity(__size - __d)) {
	if (last != end())
	  std::copy(last, end(), first);
	utils::destroy_range(end() - __d, end());
	__base.size() -= __d;
      } else {
	base_type __base_new(__size - __d);
	std::uninitialized_copy(begin(), first, __base_new.begin());
	std::uninitialized_copy(last, end(), __base_new.begin() + __n);
	
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
      
      return begin() + __n;
    }
    
  private:

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
    
    void __fill_assign(size_type __n, const Tp& __value)
    {
      const size_type __size = size();
      
      if (base_type::capacity(__n) == base_type::capacity(__size)) {
	if (__n == __size)
	  std::fill(begin(), end(), __value);
	else if (__n < __size) { 
	  std::fill(begin(), begin() + __n, __value);
	  utils::destroy_range(begin() + __n, end());
	} else {
	  std::fill(begin(), end(), __value);
	  std::uninitialized_fill(end(), begin() + __n, __value);
	}
	__base.size() = __n;
      } else {
	base_type __base_new(__n);
	std::uninitialized_fill(__base_new.begin(), __base_new.end(), __value);
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
    }
    
    template <typename ForwardIterator>
    void __range_assign(ForwardIterator first, ForwardIterator last, std::forward_iterator_tag)
    {
      const size_type __n = std::distance(first, last);
      const size_type __size = size();
      
      if (base_type::capacity(__n) == base_type::capacity(__size)) {
	if (__n == __size)
	  std::copy(first, last, begin());
	else if (__n < __size) {
	  std::copy(first, last, begin());
	  utils::destroy_range(begin() + __n, end());
	} else {
	  utils::destroy_range(begin(), end());
	  std::uninitialized_copy(first, last, begin());
	}
	__base.size() = __n;
      } else {
	base_type __base_new(__n);
	std::uninitialized_copy(first, last, __base_new.begin());
	__base.swap(__base_new);
	utils::destroy_range(__base_new.begin(), __base_new.end());
      }
    }

    template <typename InputIterator>
    void __range_assign(InputIterator first, InputIterator last, std::input_iterator_tag)
    {
      // we will first allocate at temporary buffer...
      std::vector<Tp, Alloc> vec(first, last);
      __range_assign(vec.begin(), vec.end(), std::forward_iterator_tag());
    }
    
    template <typename Iterator>
    void __range_assign(Iterator first, Iterator last)
    {
      typedef typename std::iterator_traits<Iterator>::iterator_category __category;
      __range_assign(first, last, __category());
    }
    
    template <typename _InputIterator>
    void __range_initialize(_InputIterator first, _InputIterator last, std::input_iterator_tag)
    {
      std::vector<Tp, Alloc> vec(first, last);
      __range_initialize(vec.begin(), vec.end(), std::forward_iterator_tag());
    }
    
    template <typename _ForwardIterator>
    void __range_initialize(_ForwardIterator first, _ForwardIterator last, std::forward_iterator_tag)
    {
      const size_type __n = std::distance(first, last);
      __base.initialize_vector(__n);
      std::uninitialized_copy(first, last, __base.begin());
    }

    template <typename Integer>
    void __initialize_dispatch(Integer n, Integer x, boost::true_type)
    {
      __base.initialize_vector(n);
      std::uninitialized_fill(__base.begin(), __base.end(), x);
    }
    
    template <typename _InputIterator>
    void __initialize_dispatch(_InputIterator first, _InputIterator last, boost::false_type)
    {
      typedef typename std::iterator_traits<_InputIterator>::iterator_category __category;
      __range_initialize(first, last, __category());
    }
    
  private:
    base_type __base;
  };
  
  template <typename Tp, typename Alloc>
  inline
  bool operator==(const small_vector<Tp, Alloc>& x,
		  const small_vector<Tp, Alloc>& y)
  {
    return x.size() == y.size() && std::equal(x.begin(), x.end(), y.begin());
  }
  template <typename Tp, typename Alloc>
  inline
  bool operator!=(const small_vector<Tp, Alloc>& x,
		  const small_vector<Tp, Alloc>& y)
  {
    return !(x == y);
  }
  template <typename Tp, typename Alloc>
  inline
  bool operator<(const small_vector<Tp, Alloc>& x,
		 const small_vector<Tp, Alloc>& y)
  {
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
  }
  template <typename Tp, typename Alloc>
  inline
  bool operator>(const small_vector<Tp, Alloc>& x,
		 const small_vector<Tp, Alloc>& y)
  {
    return y < x;
  }
  template <typename Tp, typename Alloc>
  inline
  bool operator<=(const small_vector<Tp, Alloc>& x,
		  const small_vector<Tp, Alloc>& y)
  {
    return ! (y < x);
  }
  template <typename Tp, typename Alloc>
  inline
  bool operator>=(const small_vector<Tp, Alloc>& x,
		  const small_vector<Tp, Alloc>& y)
  {
    return ! (x < y);
  }
};

namespace std {
  template<typename _Tp, typename _Alloc>
  inline void
  swap(utils::small_vector<_Tp,_Alloc>& x,
       utils::small_vector<_Tp,_Alloc>& y)
  {
    x.swap(y);
  }
};

#endif
