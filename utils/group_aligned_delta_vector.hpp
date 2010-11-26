// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// TODO...!

#ifndef __UTILS__GROUP_ALIGNED_DELTA_VECTOR__HPP__
#define __UTILS__GROUP_ALIGNED_DELTA_VECTOR__HPP__ 1

#include <vector>

#include <utils/group_aligned_code.hpp>
#include <utils/simple_vector.hpp>

namespace utils
{
  
  template <typename Tp>
  struct __group_aligned_delta_iterator
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Tp  value_type;
    typedef Tp* pointer;
    typedef const Tp& reference;
    
    typedef std::input_iterator_tag   iterator_category;
    
    typedef char* __pointer_type;
    typedef __group_aligned_delta_vector_iterator<Tp> self_type;
    
    
  public:
    __group_aligned_delta_vector_iterator() : value(), first(0), last(0), pos(0) {}
    __group_aligned_delta_vector_iterator(__pointer_type __first,
					  __pointer_type __last)
      : value(), first(__first), last(__last), pos(0)
    {
      if (first == last) {
	first = 0;
	last = 0;
      } else
	first += utils::group_aligned_decode(value, first);
    }
    
    const value_type& operator*() const { return value; }
    const value_type* operator->() const { return &value; }
    
    self_type& operator++() 
    {
      if (! first) return *this;
      if (first == last) {
	first = 0;
	last = 0;
	return *this;
      }
      
      value_type delta;
      first += utils::group_aligned_decode(delta, first);
      value += delta;
      
      return *this;
    }
    
    self_type operator++(int)
    {
      self_type tmp = *this;
      ++ *this;
      return tmp;
    }



  private:    
    value_type     value;
    __pointer_type first;
    __pointer_type last;
    size_type      pos;
  };

  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class group_aligned_delta_vector
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp        value_type;
    
  private:
    typedef char byte_type:
    typedef typename Alloc::template rebind<byte_type>::other byte_alloc_type;
    typedef utils::simple_vector<byte_type, byte_alloc_type> base_type;

  public:
    typedef __group_aligned_delta_iterator<value_type> iterator;
    typedef __group_aligned_delta_iterator<value_type> const_iterator;
    
  public:
    group_aligned_deleta_vector() {}
    template <typename Iterator>
    group_aligned_delta_vector(Iterator first, Iterator last) { assign(first, last); }
    
    
    const_iterator begin() const { return const_iterator(&(*base.begin()), &(*base.end())); }
    cosnt_iterator end() const { return const_iterator(); }
    
    void clear() { base.clear(); }
    bool empty() const { return base.empty(); }
    
    void swap(group_aligned_delta_vector& x) { base.swap(x.base); }

    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;

      clear();
      
      if (first == last) return;
      
      buffer_type buffer((last - first) * sizeof(value_type) * 2);
      buffer_type::iterator biter = buffer.begin();
      buffer_type::iterator hiter = buffer.begin();
      
      size_type pos = 0;
      const size_type offset = utils::group_aligned_encode(*first, &(*hiter), pos);
      biter = hiter + offset;
      hiter += offset & (- size_type((pos & 0x03) == 0x03));
      ++ pos;
      
      Iterator prev = first;
      ++ first;
      for (/**/; first != last; ++ first, ++ pos) {
	const size_type offset = utils::group_aligned_encode(*first - *prev, &(*hiter), pos);
	biter = hiter + offset;
	hiter += offset & (- size_type((pos & 0x03) == 0x03));
	
	prev = first;
      }
      
      base.resize(biter - buffer.begin());
      std::copy(buffer.begin(), biter, base.begin());
    }
    
  private:
    base_type base;
  };
  
};

namespace std
{
  
  template <typename Tp, typename Alloc>
  inline
  void swap(utils::group_aligned_delta_vector<Tp,Alloc>& x,
	    utils::group_aligned_delta_vector<Tp,Alloc>& y)
  {
    x.swap(y);
  }

};

#endif
