// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// byte aligned vector
//

#ifndef __UTILS__BYTE_ALIGNED_VECTOR__HPP__
#define __UTILS__BYTE_ALIGNED_VECTOR__HPP__ 1

#include <vector>
#include <algorithm>

#include <utils/simple_vector.hpp>
#include <utils/memory.hpp>
#include <utils/byte_aligned_code.hpp>

#include <boost/numeric/conversion/bounds.hpp>

namespace utils
{

  template <typename Tp, size_t Size>
  struct __byte_aligned_vector_value_impl {};
  
  template <typename Tp>
  struct __byte_aligned_vector_value_impl<Tp, 1>
  {
    typedef int8_t signed_type;
    typedef uint8_t unsigned_type;
  };
  
  template <typename Tp>
  struct __byte_aligned_vector_value_impl<Tp, 2>
  {
    typedef int16_t signed_type;
    typedef uint16_t unsigned_type;
  };
  
  template <typename Tp>
  struct __byte_aligned_vector_value_impl<Tp, 4>
  {
    typedef int32_t signed_type;
    typedef uint32_t unsigned_type;
  };
  
  template <typename Tp>
  struct __byte_aligned_vector_value_impl<Tp, 8>
  {
    typedef int64_t signed_type;
    typedef uint64_t unsigned_type;
  };

  template <typename Tp>
  struct __byte_aligned_vector_value
  {
    typedef __byte_aligned_vector_value_impl<Tp, sizeof(Tp)> __impl_type;
    
    typedef typename __impl_type::signed_type   signed_type;
    typedef typename __impl_type::unsigned_type unsigned_type;
  };


  template <typename Tp>
  struct __byte_aligned_vector_iterator
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp  value_type;
    typedef Tp* pointer;
    typedef const Tp& reference;
    
    typedef std::input_iterator_tag   iterator_category;
    
    typedef const char* __pointer_type;
    typedef __byte_aligned_vector_iterator<Tp> self_type;


    typedef typename __byte_aligned_vector_value<Tp>::signed_type signed_value_type;
    typedef typename __byte_aligned_vector_value<Tp>::unsigned_type unsigned_value_type;
    
  public:
    __byte_aligned_vector_iterator() : first(0), last(0) {}
    __byte_aligned_vector_iterator(__pointer_type __first,
					 __pointer_type __last)
      : first(__first), last(__last) 
    {
      if (first == last) {
	first = 0;
	last = 0;
      } else {
	first += utils::byte_aligned_decode(__min, first);
	first += utils::byte_aligned_decode((signed_value_type&) __value, first);
	((signed_value_type&) __value) += __min;
      }
    }
    
    const value_type& operator*() const { return __value; }
    const value_type* operator->() const { return &__value; }
    
    self_type& operator++() 
    {
      if (first == last) {
	first = 0;
	last = 0;
      } else {
	first += utils::byte_aligned_decode((signed_value_type&) __value, first);
	((signed_value_type&) __value) += __min;
      }
      
      return *this;
    }
    
    self_type operator++(int)
    {
      self_type tmp = *this;
      ++ *this;
      return tmp;
    }



  public:
    __pointer_type first;
    __pointer_type last;

    value_type        __value;
    signed_value_type __min;
  };
  
  template <typename _Tp>
  inline
  bool operator==(const __byte_aligned_vector_iterator<_Tp>& x,
		  const __byte_aligned_vector_iterator<_Tp>& y)
  {
    return x.first == y.first && x.last == y.last;
  }
  
  template <typename _Tp>
  inline
  bool operator!=(const __byte_aligned_vector_iterator<_Tp>& x,
		  const __byte_aligned_vector_iterator<_Tp>& y)
  {
    return x.first != y.first || x.last != y.last;
  }
  
  
  
  template <typename _Tp>
  inline
  bool operator<(const __byte_aligned_vector_iterator<_Tp>& x,
		 const __byte_aligned_vector_iterator<_Tp>& y)
  {
    return (x.first < y.first && (!(y.first < x.first) && x.last < y.last));
  }
  
  template <typename _Tp>
  inline
  bool operator>(const __byte_aligned_vector_iterator<_Tp>& x,
		 const __byte_aligned_vector_iterator<_Tp>& y)
  {
    return y < x;
  }
  
  template <typename _Tp>
  inline
  bool operator<=(const __byte_aligned_vector_iterator<_Tp>& x,
		  const __byte_aligned_vector_iterator<_Tp>& y)
  {
    return ! (y < x);
  }
  
  template <typename _Tp>
  inline
  bool operator>=(const __byte_aligned_vector_iterator<_Tp>& x,
		  const __byte_aligned_vector_iterator<_Tp>& y)
  {
    return ! (x < y);
  }

  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class byte_aligned_vector
  {
  private:
    typedef char byte_type;
    typedef typename Alloc::template rebind<byte_type>::other byte_alloc_type;
    
    typedef utils::simple_vector<byte_type, byte_alloc_type> base_type;
    
  public:
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;

    typedef Tp value_type;

    typedef __byte_aligned_vector_iterator<Tp> iterator;
    typedef __byte_aligned_vector_iterator<Tp> const_iterator;
    
  public:
    byte_aligned_vector() : base() {}
    template <typename Iterator>
    byte_aligned_vector(Iterator first, Iterator last) : base() { assign(first, last); }
    
    const_iterator begin() const { return const_iterator(&(*base.begin()), &(*base.end()));}
    const_iterator end() const { return const_iterator(); }
    
    bool empty() const { return base.empty(); }
    size_type compressed_size() const { return base.size(); }
    
    void clear() { base.clear(); }
    void swap(byte_aligned_vector& x) { base.swap(x.base); }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      typedef typename __byte_aligned_vector_value<Tp>::signed_type   signed_value_type;
      typedef typename __byte_aligned_vector_value<Tp>::unsigned_type unsigned_value_type;
      
      typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;
      
      base.clear();
      
      if (first == last) return;

      size_type size_vector = 0;
      signed_value_type value_min = boost::numeric::bounds<signed_value_type>::highest();
      for (Iterator iter = first; iter != last; ++ iter, ++ size_vector)
	value_min = std::min(value_min, (const signed_value_type&) *iter);
      
      // # of diff-elements + lower_bound for diff value
      buffer_type buffer(size_vector * 16 + 16);
      
      typename buffer_type::iterator biter = buffer.begin();
      biter += utils::byte_aligned_encode(value_min, &(*biter));
      
      for (Iterator iter = first; iter != last; ++ iter)
	biter += utils::byte_aligned_encode((const signed_value_type&) *iter - value_min, &(*biter));
      
      base.assign(buffer.begin(), biter);
    }

  public:
    template <typename T, typename A>
    friend
    bool operator==(const byte_aligned_vector<T,A>& x, const byte_aligned_vector<T,A>& y) { return x.base == y.base; }
    
    template <typename T, typename A>
    friend
    bool operator!=(const byte_aligned_vector<T,A>& x, const byte_aligned_vector<T,A>& y) { return x.base != y.base; }
    
    template <typename T, typename A>
    friend
    bool operator<(const byte_aligned_vector<T,A>& x, const byte_aligned_vector<T,A>& y) { return x.base < y.base; }
    
    template <typename T, typename A>
    friend
    bool operator>(const byte_aligned_vector<T,A>& x, const byte_aligned_vector<T,A>& y) { return x.base > y.base; }
    
    template <typename T, typename A>
    friend
    bool operator<=(const byte_aligned_vector<T,A>& x, const byte_aligned_vector<T,A>& y) { return x.base <= y.base; }
    
    template <typename T, typename A>
    friend
    bool operator>=(const byte_aligned_vector<T,A>& x, const byte_aligned_vector<T,A>& y) { return x.base >= y.base; }    
    
  private:
    base_type base;
  };
  
};

namespace std
{

  template <typename Tp, typename Alloc>
  inline
  void swap(utils::byte_aligned_vector<Tp,Alloc>& x,
	    utils::byte_aligned_vector<Tp,Alloc>& y)
  {
    x.swap(y);
  }

};

#endif
