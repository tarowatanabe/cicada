// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// byte aligned "delta" vector
//

#ifndef __UTILS__BYTE_ALIGNED_PAIR_DELTA_VECTOR__HPP__
#define __UTILS__BYTE_ALIGNED_PAIR_DELTA_VECTOR__HPP__ 1

#include <vector>
#include <algorithm>

#include <utils/simple_vector.hpp>
#include <utils/memory.hpp>
#include <utils/byte_aligned_code.hpp>

#include <boost/numeric/conversion/bounds.hpp>

namespace utils
{

  template <typename Tp, size_t Size>
  struct __byte_aligned_pair_delta_vector_value_impl {};
  
  template <typename Tp>
  struct __byte_aligned_pair_delta_vector_value_impl<Tp, 1>
  {
    typedef int8_t signed_type;
    typedef uint8_t unsigned_type;
  };
  
  template <typename Tp>
  struct __byte_aligned_pair_delta_vector_value_impl<Tp, 2>
  {
    typedef int16_t signed_type;
    typedef uint16_t unsigned_type;
  };
  
  template <typename Tp>
  struct __byte_aligned_pair_delta_vector_value_impl<Tp, 4>
  {
    typedef int32_t signed_type;
    typedef uint32_t unsigned_type;
  };
  
  template <typename Tp>
  struct __byte_aligned_pair_delta_vector_value_impl<Tp, 8>
  {
    typedef int64_t signed_type;
    typedef uint64_t unsigned_type;
  };

  template <typename Tp>
  struct __byte_aligned_pair_delta_vector_value
  {
    typedef __byte_aligned_pair_delta_vector_value_impl<Tp, sizeof(Tp)> __impl_type;
    
    typedef typename __impl_type::signed_type   signed_type;
    typedef typename __impl_type::unsigned_type unsigned_type;
  };


  template <typename Tp1, typename Tp2>
  struct __byte_aligned_pair_delta_vector_iterator
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp1 first_type;
    typedef Tp2 second_type;
    
    typedef std::pair<Tp1, Tp2> value_type;
    typedef value_type* pointer;
    typedef const value_type& reference;
    
    typedef std::input_iterator_tag   iterator_category;
    
    typedef const char* __pointer_type;
    typedef __byte_aligned_pair_delta_vector_iterator<Tp1, Tp2> self_type;


    typedef typename __byte_aligned_pair_delta_vector_value<Tp1>::signed_type   signed_first_type;
    typedef typename __byte_aligned_pair_delta_vector_value<Tp1>::unsigned_type unsigned_first_type;

    typedef typename __byte_aligned_pair_delta_vector_value<Tp2>::signed_type   signed_second_type;
    typedef typename __byte_aligned_pair_delta_vector_value<Tp2>::unsigned_type unsigned_second_type;
    
  public:
    __byte_aligned_pair_delta_vector_iterator() : first(0), last(0) {}
    __byte_aligned_pair_delta_vector_iterator(__pointer_type __first,
					      __pointer_type __last)
      : first(__first), last(__last) 
    {
      if (first == last) {
	first = 0;
	last = 0;
      } else {
	first += utils::byte_aligned_decode(__min_first, first);
	first += utils::byte_aligned_decode(__min_second, first);
	
	first += utils::byte_aligned_decode((signed_first_type&) __value.first, first);
	first += utils::byte_aligned_decode((signed_second_type&) __value.second, first);
	
	((signed_second_type&) __value.second) += __min_second;
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
	signed_first_type __diff;
	first += utils::byte_aligned_decode(__diff, first);
	first += utils::byte_aligned_decode((signed_second_type&) __value.second, first);
	
	((signed_first_type&) __value.first)   += __diff + __min_first;
	((signed_second_type&) __value.second) += __min_second;
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
    
    value_type         __value;
    signed_first_type  __min_first;
    signed_second_type __min_second;
  };
  
  template <typename _Tp1, typename _Tp2>
  inline
  bool operator==(const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& x,
		  const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& y)
  {
    return x.first == y.first && x.last == y.last;
  }
  
  template <typename _Tp1, typename _Tp2>
  inline
  bool operator!=(const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& x,
		  const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& y)
  {
    return x.first != y.first || x.last != y.last;
  }
  
  
  
  template <typename _Tp1, typename _Tp2>
  inline
  bool operator<(const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& x,
		 const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& y)
  {
    return (x.first < y.first && (!(y.first < x.first) && x.last < y.last));
  }
  
  template <typename _Tp1, typename _Tp2>
  inline
  bool operator>(const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& x,
		 const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& y)
  {
    return y < x;
  }
  
  template <typename _Tp1, typename _Tp2>
  inline
  bool operator<=(const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& x,
		  const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& y)
  {
    return ! (y < x);
  }
  
  template <typename _Tp1, typename _Tp2>
  inline
  bool operator>=(const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& x,
		  const __byte_aligned_pair_delta_vector_iterator<_Tp1,_Tp2>& y)
  {
    return ! (x < y);
  }

  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class byte_aligned_pair_delta_vector
  {
  private:
    typedef char byte_type;
    typedef typename Alloc::template rebind<byte_type>::other byte_alloc_type;
    
    typedef utils::simple_vector<byte_type, byte_alloc_type> base_type;
    
  public:
    typedef typename base_type::size_type       size_type;
    typedef typename base_type::difference_type difference_type;

    typedef Tp value_type;

    typedef typename value_type::first_type  first_type;
    typedef typename value_type::second_type second_type;

    typedef __byte_aligned_pair_delta_vector_iterator<first_type, second_type> iterator;
    typedef __byte_aligned_pair_delta_vector_iterator<first_type, second_type> const_iterator;
    
  public:
    byte_aligned_pair_delta_vector() : base() {}
    template <typename Iterator>
    byte_aligned_pair_delta_vector(Iterator first, Iterator last) : base() { assign(first, last); }
    
    const_iterator begin() const { return const_iterator(&(*base.begin()), &(*base.end()));}
    const_iterator end() const { return const_iterator(); }
    
    bool empty() const { return base.empty(); }
    size_type compressed_size() const { return base.size(); }
    
    void clear() { base.clear(); }
    void swap(byte_aligned_pair_delta_vector& x) { base.swap(x.base); }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      typedef typename __byte_aligned_pair_delta_vector_value<first_type>::signed_type   signed_first_type;
      typedef typename __byte_aligned_pair_delta_vector_value<first_type>::unsigned_type unsigned_first_type;

      typedef typename __byte_aligned_pair_delta_vector_value<second_type>::signed_type   signed_second_type;
      typedef typename __byte_aligned_pair_delta_vector_value<second_type>::unsigned_type unsigned_second_type;
      
      typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;
      
      base.clear();
      
      if (first == last) return;

      size_type size_vector = 0;
      
      const signed_first_type  first_init  = (const signed_first_type&) (*first).first;
      const signed_second_type second_init = (const signed_second_type&) (*first).second;

      signed_first_type  first_min  = boost::numeric::bounds<signed_first_type>::highest();
      signed_second_type second_min = second_init;
      
      {
	signed_first_type first_prev = first_init;
	Iterator iter = first;
	++ iter;
	for (/**/; iter != last; ++ iter, size_vector += 2) {
	  const signed_first_type& value_first = (const signed_first_type&) (*iter).first;
	  first_min = std::min(first_min, signed_first_type(value_first - first_prev));
	  second_min = std::min(second_min, (const signed_second_type&) (*iter).second);
	  first_prev = value_first;
	}
      }
      
      
      // # of diff-elements + initial value + lower_bound for diff value
      buffer_type buffer(size_vector * 16 + 16 * 2 + 16 * 2);
      
      typename buffer_type::iterator biter = buffer.begin();
      biter += utils::byte_aligned_encode(first_min, &(*biter));
      biter += utils::byte_aligned_encode(second_min, &(*biter));
      
      biter += utils::byte_aligned_encode(first_init, &(*biter));
      biter += utils::byte_aligned_encode(second_init - second_min, &(*biter));
      
      {
	signed_first_type first_prev = first_init;
	Iterator iter = first;
	++ iter;
	for (/**/; iter != last; ++ iter) {
	  const signed_first_type& value_first = (const signed_first_type&) (*iter).first;
	  biter += utils::byte_aligned_encode((value_first - first_prev) - first_min, &(*biter));
	  biter += utils::byte_aligned_encode((const signed_second_type&) (*iter).second - second_min, &(*biter));
	  first_prev = value_first;
	}
      }
      
      base.assign(buffer.begin(), biter);
    }

  public:
    template <typename T, typename A>
    friend
    bool operator==(const byte_aligned_pair_delta_vector<T,A>& x, const byte_aligned_pair_delta_vector<T,A>& y) { return x.base == y.base; }
    
    template <typename T, typename A>
    friend
    bool operator!=(const byte_aligned_pair_delta_vector<T,A>& x, const byte_aligned_pair_delta_vector<T,A>& y) { return x.base != y.base; }
    
    template <typename T, typename A>
    friend
    bool operator<(const byte_aligned_pair_delta_vector<T,A>& x, const byte_aligned_pair_delta_vector<T,A>& y) { return x.base < y.base; }
    
    template <typename T, typename A>
    friend
    bool operator>(const byte_aligned_pair_delta_vector<T,A>& x, const byte_aligned_pair_delta_vector<T,A>& y) { return x.base > y.base; }
    
    template <typename T, typename A>
    friend
    bool operator<=(const byte_aligned_pair_delta_vector<T,A>& x, const byte_aligned_pair_delta_vector<T,A>& y) { return x.base <= y.base; }
    
    template <typename T, typename A>
    friend
    bool operator>=(const byte_aligned_pair_delta_vector<T,A>& x, const byte_aligned_pair_delta_vector<T,A>& y) { return x.base >= y.base; }    
    
  private:
    base_type base;
  };
  
};

namespace std
{

  template <typename Tp, typename Alloc>
  inline
  void swap(utils::byte_aligned_pair_delta_vector<Tp,Alloc>& x,
	    utils::byte_aligned_pair_delta_vector<Tp,Alloc>& y)
  {
    x.swap(y);
  }

};

#endif
