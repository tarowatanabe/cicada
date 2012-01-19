// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS_ISTREAM_LINE_ITERATOR__HPP__
#define __UTILS_ISTREAM_LINE_ITERATOR__HPP__ 1

#include <cstddef>
#include <string>
#include <iostream>
#include <iterator>

namespace utils
{
  
  struct istream_line_iterator
  {
  public:
    typedef std::input_iterator_tag iterator_category;
    
    typedef std::istream stream_type;
    typedef std::string  value_type;
    
    typedef std::string& reference;
    typedef std::string* pointer;

    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

  public:
    istream_line_iterator() : __line(), __stream() {}
    istream_line_iterator(stream_type& stream) : __line() , __stream(&stream) { __read(); }

  public:
    const value_type& operator*()  const { return __line; }
    const value_type* operator->() const { return &__line; }

    istream_line_iterator& operator++()
    {
      __read();
      return *this;
    }
    
    istream_line_iterator operator++(int)
    {
      istream_line_iterator __tmp = *this;
      __read();
      return __tmp;
    }
    
  public:
    void __read()
    {
      if (! __stream) return;
      
      if (! std::getline(*__stream, __line)) {
	__line.clear();
	__stream = 0;
      }
    }


  public:    
    value_type   __line;
    stream_type* __stream;
  };

  inline
  bool operator==(const istream_line_iterator& x,
		  const istream_line_iterator& y)
  {
    return x.__stream == y.__stream;
  }
  
  inline
  bool operator!=(const istream_line_iterator& x,
		  const istream_line_iterator& y)
  {
    return x.__stream != y.__stream;
  }
};

#endif
