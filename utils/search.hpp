// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SEARCH__HPP__
#define __UTILS__SEARCH__HPP__ 1

#include <algorithm>

#include <utils/bithack.hpp>

namespace utils
{
  template <typename Iterator, typename Tp>
  Iterator linear_search(Iterator first, Iterator last, const Tp& val)
  {
    for (/**/; first != last && *first < val; ++ first) {}
    
    return (first != last && !(val < *first) ? first : last);
  }

  template <typename Iterator, typename Tp, typename Compare>
  Iterator linear_search(Iterator first, Iterator last, const Tp& val, Compare comp)
  {
    for (/**/; first != last && comp(*first, val); ++ first) {}
    
    return (first != last && ! comp(val, *first) ? first : last);
  }

  template <typename Iterator, typename Tp>
  Iterator binary_search(Iterator first, Iterator last, const Tp& val)
  {
    Iterator iter = std::lower_bound(first, last, val);
    
    return (iter != last && !(val < *iter) ? iter : last);
  }
  
  template <typename Iterator, typename Tp, typename Compare>
  Iterator binary_search(Iterator first, Iterator last, const Tp& val, Compare comp)
  {
    Iterator iter = std::lower_bound(first, last, val, comp);
    
    return (iter != last && ! comp(val, *iter) ? iter : last);
  }
  
  template <size_t Bytes>
  struct __interpolation_pivot
  {
    size_t operator()(uint64_t offset, uint64_t range, uint64_t width) const
    {
      return (offset * width) / (range + 1);
    }
  };
  
  template <> struct __interpolation_pivot<8>
  {
    size_t operator()(uint64_t offset, uint64_t range, uint64_t width) const
    {
      const uint64_t ret = static_cast<uint64_t>(float(offset) / float(range) * float(width));
      
      return utils::bithack::branch(ret < width, ret, width - 1);
    }
  };
  
  template <typename Iterator, typename Tp>
  Iterator interpolation_search(Iterator first, Iterator last, const Tp& val)
  {
    typedef typename std::iterator_traits<Iterator>::value_type      value_type;
    typedef typename std::iterator_traits<Iterator>::difference_type difference_type;
    
    if (first == last) return first;
    
    value_type front = *first;
    
    // val <= front
    if (! (front < val))
      return (val < front ? last : first);
    
    // we will work with the range of [first, last], not a conventional [first, last)
    Iterator __last = last;
    -- last;
    value_type back = *last;
    
    // back <= val
    if (! (val < back))
      return (back < val ? __last : last);
    
    // search the range [first + 1, last - 1]
    
    value_type pivot;
    Iterator middle;
    difference_type diff;
    difference_type length = std::distance(first, last);
    __interpolation_pivot<sizeof(value_type)> interpolation;
    
    while (length > 1) {
      diff = 1 + interpolation(val - front, back - front, length - 1);
      middle = first;
      std::advance(middle, diff);
      pivot = *middle;
      
      if (pivot < val) {
	first = middle;
	front = pivot;
	length -= diff;
      } else if (val < pivot) {
	last = middle;
	back = pivot;
	length = diff;
      } else
	return middle;
    }
    
    return __last;
  }


};


#endif
