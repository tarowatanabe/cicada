// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MEMORY__HPP__
#define __UTILS__MEMORY__HPP__ 1

#include <iterator>

#include <boost/type_traits.hpp>


namespace utils
{
  
  // construct object

  template <typename Tp>
  inline void __construct_object(Tp* pointer, const Tp& value, boost::true_type)
  {
    *pointer = value;
  }
  
  template <typename Tp>
  inline void __construct_object(Tp* pointer, const Tp& value, boost::false_type)
  {
    ::new(static_cast<void*>(pointer)) Tp(value);
  }
  
  template <typename Tp>
  inline void construct_object(Tp* pointer, const Tp& value)
  {
    __construct_object(pointer, value, boost::has_trivial_constructor<Tp>());
  }
  
  template <typename Tp>
  inline void construct_object(Tp* pointer)
  {
    construct_object(pointer, Tp());
  }

  // destroy object

  template <typename Tp>
  inline void __destroy_object(Tp* pointer, boost::true_type)
  {
    
  }

  template <typename Tp>
  inline void __destroy_object(Tp* pointer, boost::false_type)
  {
    pointer->~Tp();
  }
  
  
  template <typename Tp>
  inline void destroy_object(Tp* pointer)
  {
    __destroy_object(pointer, boost::has_trivial_destructor<Tp>());
  }

  // destroy range

  template <typename Iterator>
  inline
  void __destroy_range(Iterator first, Iterator last, boost::true_type)
  {
    
  }
  
  template <typename Iterator>
  inline
  void __destroy_range(Iterator first, Iterator last, boost::false_type)
  {
    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    
    for (/**/; first != last; ++ first)
      utils::destroy_object(&(*first));
  }
  
  
  template <typename Iterator>
  inline void destroy_range(Iterator first, Iterator last)
  {
    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    
    __destroy_range(first, last, boost::has_trivial_destructor<value_type>());
  }
};

#endif
