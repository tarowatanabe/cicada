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
#include <cctype>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>

namespace utils
{
  template <typename _Traits>
  class basic_piece
  {
  public:
    typedef _Traits   traits_type;
    
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

    static const size_type npos()
    {
      return size_type(-1);
    }
    
  public:
    basic_piece() : first_(0), last_(0) {}
    
    // construct from basic_piece
    template <typename _T>
    basic_piece(const basic_piece<_T>& str)
      : first_(str.begin()), last_(str.end()) {}
    template <typename _T>
    basic_piece(const basic_piece<_T>& str, size_type pos, size_type n = npos())
      : first_(str.begin() + pos), last_(n == npos() ? str.end() : std::min(str.begin() + pos + n, str.end()))
    {
      if (first_ > str.end())
	throw std::out_of_range("basic_piece::basic_piece");
    }
    
    // construc from string
    basic_piece(const std::string& str)
      : first_(str.c_str()), last_(str.c_str() + str.size()) {}
    basic_piece(const std::string& str, size_type pos, size_type n = npos())
      : first_(str.c_str() + pos), last_(n == npos() ? str.c_str() + str.size() : std::min(str.c_str() + pos + n, str.c_str() + str.size()))
    {
      if (first_ > str.c_str() + str.size())
	throw std::out_of_range("basic_piece::basic_piece");
    }
    
    // consruct from bare string
    basic_piece(const char* str) : first_(str), last_(str == 0 ? 0 : str + ::strlen(str)) {}
    basic_piece(const char* offset, const size_type len) : first_(offset), last_(offset + len) {}
    
    // construct by iterator
    basic_piece(const char* __first, const char* __last) : first_(__first), last_(__last) {}
    basic_piece(std::string::const_iterator __first, std::string::const_iterator __last)
      : first_(&(*__first)), last_(&(*__last)) {}
    basic_piece(std::vector<char>::const_iterator __first, std::vector<char>::const_iterator __last)
      : first_(&(*__first)), last_(&(*__last)) {}
    
    pointer data() const { return first_; }
    pointer c_str() const { return first_; }
    
    bool empty() const { return first_ == last_; }
    size_type size() const { return last_ - first_; }
    size_type length() const { return last_ - first_; }
    
    void clear() { first_ = 0; last_ = 0; }
    
    // assign from basic_piece
    template <typename _T>
    void assign(const basic_piece<_T>& str)
    {
      first_ = str.begin();
      last_  = str.end();
    }
    template <typename _T>
    void assign(const basic_piece<_T>& str, size_type pos, size_type n = npos())
    {
      first_ = str.begin() + pos;
      last_  = (n == npos() ? str.end() : std::min(str.begin() + pos + n, str.end()));

      if (first_ > str.end())
	throw std::out_of_range("basic_piece::assign");
    }
    
    // assign from string
    void assign(const std::string& str)
    {
      first_ = str.c_str();
      last_  = str.c_str() + str.size();
    }
    void assign(const std::string& str, size_type pos, size_type n = npos())
    {
      first_ = str.c_str() + pos;
      last_  = (n == npos() ? str.c_str() + str.size() : std::min(str.c_str() + pos + n, str.c_str() + str.size()));

      if (first_ > str.c_str() + str.size())
	throw std::out_of_range("basic_piece::assign");
    }
    
    
    // assignment by bare string
    void assign(const char* str)
    {
      first_ = str;
      last_  = (str == 0 ? 0 : str + ::strlen(str));
    }
    void assign(const char* offset, const size_type len)
    {
      first_ = offset;
      last_  = offset + len;
    }
    
    // assignment by iterators
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
    
    basic_piece substr(size_type pos = 0, size_type n = npos()) const
    {
      const_iterator first = first_ + pos;
      const_iterator last  = (n == npos() ? last_ : std::min(first_ + pos + n, last_));
      
      if (first > last_)
	throw std::out_of_range("basic_piece::substr");
      
      return basic_piece(first, last);
    }

    size_type find(const char* __s, size_type __pos, size_type __n) const
    {
      const size_type __size = size();
      const char* __data = first_;
      for (/**/; __pos + __n <= __size; ++__pos)
        if (traits_type::compare(__data + __pos, __s, __n) == 0)
          return __pos;
      return npos();
    }
    
    size_type find(char __c, size_type __pos = 0) const
    {
      size_type __ret = npos();
      const size_type __size = size();
      if (__pos < __size) {
	const char* __data = first_;
	const size_type __n = __size - __pos;
	const char* __p = traits_type::find(__data + __pos, __n, __c);
	if (__p)
	  __ret = __p - __data;
      }
      return __ret;
    }

    int compare(const basic_piece& x) const
    {
      const difference_type __x = size();
      const difference_type __y = x.size();
      
      const int __result = traits_type::compare(data(), x.data(), utils::bithack::min(__x, __y));
      if (__result)
	return __result;
      else
	return __x - __y;
    }
    
  private:
    const_iterator first_;
    const_iterator last_;
  };
  
  template <typename _T>
  inline
  size_t hash_value(basic_piece<_T> const& x)
  {
    return utils::hashmurmur<size_t>()(x.begin(), x.end(), 0);
  }

  template <typename _T>
  inline
  bool operator==(const basic_piece<_T>& x, const basic_piece<_T>& y)
  {
    return x.compare(y) == 0;
  }
    
  template <typename _T>
  inline
  bool operator!=(const basic_piece<_T>& x, const basic_piece<_T>& y)
  {
    return x.compare(y) != 0;
  }
  
  template <typename _T>
  inline
  bool operator<(const basic_piece<_T>& x, const basic_piece<_T>& y)
  {
    return x.compare(y) < 0;
  }
  
  template <typename _T>
  inline
  bool operator>(const basic_piece<_T>& x, const basic_piece<_T>& y)
  {
    return x.compare(y) > 0;
  }
  
  
  template <typename _T>
  inline
  bool operator<=(const basic_piece<_T>& x, const basic_piece<_T>& y)
  {
    return x.compare(y) <= 0;
  }
  
  template <typename _T>
  inline
  bool operator>=(const basic_piece<_T>& x, const basic_piece<_T>& y)
  {
    return x.compare(y) >= 0;
  }
  
  
  template <typename _T>
  inline
  bool operator==(const basic_piece<_T>& x, const char* y)
  {
    return x == basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator!=(const basic_piece<_T>& x, const char* y)
  {
    return x != basic_piece<_T>(y);
  }

  template <typename _T>
  inline
  bool operator<(const basic_piece<_T>& x, const char* y)
  {
    return x < basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator>(const basic_piece<_T>& x, const char* y)
  {
    return x > basic_piece<_T>(y);
  }

  template <typename _T>
  inline
  bool operator<=(const basic_piece<_T>& x, const char* y)
  {
    return x <= basic_piece<_T>(y);
  }

  template <typename _T>
  inline
  bool operator>=(const basic_piece<_T>& x, const char* y)
  {
    return x >= basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator==(const basic_piece<_T>& x, const std::string& y)
  {
    return x == basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator!=(const basic_piece<_T>& x, const std::string& y)
  {
    return x != basic_piece<_T>(y);
  }

  template <typename _T>
  inline
  bool operator<(const basic_piece<_T>& x, const std::string& y)
  {
    return x < basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator>(const basic_piece<_T>& x, const std::string& y)
  {
    return x > basic_piece<_T>(y);
  }

  template <typename _T>
  inline
  bool operator<=(const basic_piece<_T>& x, const std::string& y)
  {
    return x <= basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator>=(const basic_piece<_T>& x, const std::string& y)
  {
    return x >= basic_piece<_T>(y);
  }
  
  template <typename _T>
  inline
  bool operator==(const char* x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) == y;
  }
  
  template <typename _T>
  inline
  bool operator!=(const char* x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) != y;
  }
  
  template <typename _T>
  inline
  bool operator<(const char* x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) < y;
  }
  
  template <typename _T>
  inline
  bool operator>(const char* x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) > y;
  }
  
  template <typename _T>
  inline
  bool operator<=(const char* x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) <= y;
  }
  
  template <typename _T>
  inline
  bool operator>=(const char* x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) >= y;
  }

  template <typename _T>
  inline
  bool operator==(const std::string& x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) == y;
  }
  
  template <typename _T>
  inline
  bool operator!=(const std::string& x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) != y;
  }
  
  template <typename _T>
  inline
  bool operator<(const std::string& x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) < y;
  }
  
  template <typename _T>
  inline
  bool operator>(const std::string& x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) > y;
  }
  
  template <typename _T>
  inline
  bool operator<=(const std::string& x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) <= y;
  }
  
  template <typename _T>
  inline
  bool operator>=(const std::string& x, const basic_piece<_T>& y)
  {
    return basic_piece<_T>(x) >= y;
  }

  template <typename _T>
  inline
  std::string operator+(const std::string& x, const basic_piece<_T>& y)
  {
    return x + static_cast<std::string>(y);
  }
  
  template <typename _T>
  inline
  std::string operator+(const char* x, const basic_piece<_T>& y)
  {
    return x + static_cast<std::string>(y);
  }

  template <typename _T>
  inline
  std::string operator+(const basic_piece<_T>& x, const std::string& y)
  {
    return static_cast<std::string>(x) + y;
  }
  
  template <typename _T>
  inline
  std::string operator+(const basic_piece<_T>& x, const char* y)
  {
    return static_cast<std::string>(x) + y;
  }
  
  template <typename _T>
  inline
  std::ostream& operator<<(std::ostream& os, const basic_piece<_T>& x)
  {
    return os.write(x.c_str(), x.size());
  }
  
  
  struct __piece_ichar_traits
  {
    
    static const char* find(const char* __s, size_t __n, const char& __a)
    {
      const char __u = ::toupper(__a);
      for (std::size_t __i = 0; __i < __n; ++__i)
        if (::toupper(__s[__i]) == __u)
          return __s + __i;
      return 0;
    }

    static int compare(const char* x, const char* y, size_t n)
    {
      return ::strncasecmp(x, y, n);
    }

  };

  typedef basic_piece<std::char_traits<char> > piece;
  typedef basic_piece<__piece_ichar_traits > ipiece;
};

#endif
