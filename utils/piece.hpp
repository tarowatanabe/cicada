// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__PIECE__HPP__
#define __UTILS__PIECE__HPP__ 1

// piece, an implementation of string-like object
// inspired by google's string-piece
//

#include <string>
#include <cstring>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <utils/hashmurmur.hpp>

namespace utils
{
  class piece
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef const char  value_type;
    typedef const char* pointer;
    typedef const char* iterator;
    typedef const char* const_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef const char& reference;
    typedef const char& const_reference;

    static const size_type npos;
    
  public:
    piece() : first_(0), last_(0) {}
    piece(const char* str) : first_(str), last_(str == 0 ? 0 : str + ::strlen(str)) {}
    piece(const std::string& str) : first_(str.c_str()), last_(str.c_str() + str.size()) {}
    piece(const char* offset, const size_type len) : first_(offset), last_(offset + len) {}
    piece(const char* __first, const char* __last) : first_(__first), last_(__last) {}
    piece(std::string::const_iterator __first, std::string::const_iterator __last)
      : first_(&(*__first)), last_(&(*__last)) {}

    pointer data() const { return first_; }
    pointer c_str() const { return first_; }
    
    bool empty() const { return first_ == last_; }
    size_type size() const { return last_ - first_; }
    size_type length() const { return last_ - first_; }
    
    void clear() { first_ = 0; last_ = 0; }

    void assign(const char* str)
    {
      first_ = str;
      last_  = (str == 0 ? 0 : str + ::strlen(str));
    }
    void assign(const std::string& str)
    {
      first_ = str.c_str();
      last_  = str.c_str() + str.size();
    }
    
    void assign(const char* offset, const size_type len)
    {
      first_ = offset;
      last_  = offset + len;
    }
    
    void assign(const char* __first, const char* __last)
    {
      first_ = __first;
      last_  = __last;
    }
    
    void assign(std::string::const_iterator __first, std::string::const_iterator __last)
    {
      first_ = &(*__first);
      last_  = &(*__last);
    }
    
    reference operator[](size_type pos) const { return first_[pos]; }
    
    operator std::string() const { return std::string(first_, last_); }
    
    const_iterator begin() const { return first_; }
    const_iterator end() const { return last_; }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(last_); }
    const_reverse_iterator rend() const { return const_reverse_iterator(first_); }
    
    size_type max_size() const { return last_ - first_; }
    size_type capacity() const { return last_ - first_; }
    
    piece substr(size_type pos = 0, size_type n = npos) const
    {
      const_iterator first = first_ + pos;
      const_iterator last  = (n == npos ? last_ : std::min(first_ + pos + n, last_));
      
      if (first > last_)
	throw std::out_of_range("piece::substr");
      
      return piece(first, last);
    }
    
  private:
    const_iterator first_;
    const_iterator last_;
  };
  
  inline
  size_t hash_value(piece const& x)
  {
    return utils::hashmurmur<size_t>()(x.begin(), x.end(), 0);
  }

  inline
  bool operator==(const piece& x, const piece& y)
  {
    return x.size() == y.size() && std::equal(x.begin(), x.end(), y.begin());
  }
  
  inline
  bool operator!=(const piece& x, const piece& y)
  {
    return ! (x == y);
  }
  
  inline
  bool operator<(const piece& x, const piece& y)
  {
    return std::lexicographical_compare(x.begin(), x.end(), y.begin(), y.end());
  }
  
  inline
  bool operator>(const piece& x, const piece& y)
  {
    return y < x;
  }
  
  inline
  bool operator<=(const piece& x, const piece& y)
  {
    return ! (y < x);
  }
  
  inline
  bool operator>=(const piece& x, const piece& y)
  {
    return ! (x < y);
  }
  
  inline
  std::ostream& operator<<(std::ostream& os, const piece& x)
  {
    return os.write(x.c_str(), x.size());
  }
  
};

#endif
