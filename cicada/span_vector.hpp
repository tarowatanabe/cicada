// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SPAN_VECTOR__HPP__
#define __CICADA__SPAN_VECTOR__HPP__ 1

#include <stdint.h>

#include <stdexcept>
#include <iostream>
#include <vector>
#include <utility>
#include <string>

#include <cicada/symbol.hpp>

#include <utils/hashmurmur.hpp>

namespace cicada
{
  class SpanVector
  {
  public:
    typedef int32_t   index_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Symbol label_type;
    
    struct Span
    {
      Span() : first(0), last(0), label() {}
      template <typename Integral>
      Span(const std::pair<Integral, Integral>& x)
	: first(x.first), last(x.second), label() {}
      template <typename Integral>
      Span(const std::pair<Integral, Integral>& x, const label_type& __label)
	: first(x.first), last(x.second), label(__label) {}
      
      Span(index_type _first, index_type _last)
	: first(_first), last(_last), label() {}
      Span(index_type _first, index_type _last, const label_type& __label)
	: first(_first), last(_last), label(__label) {}
      
      difference_type size() const { return difference_type(last) - difference_type(first); }
      bool empty() const { return first == last; }
      void clear() { first = 0; last = 0; label = label_type(); }

      index_type first;
      index_type last;
      label_type label;
    };

    typedef Span span_type;
    
  private:
    typedef std::vector<span_type, std::allocator<span_type> > spans_type;
    
  public:
    typedef spans_type::value_type      value_type;

    typedef spans_type::const_iterator         const_iterator;
    typedef spans_type::iterator               iterator;
    typedef spans_type::const_reverse_iterator const_reverse_iterator;
    typedef spans_type::reverse_iterator       reverse_iterator;
    typedef spans_type::const_reference        const_reference;
    typedef spans_type::reference              reference;

  public:
    SpanVector() : __spans() {}
    SpanVector(const std::string& x) : __spans() { assign(x); }
    SpanVector(size_type __n) : __spans(__n) {}
    SpanVector(size_type __n, const value_type& __value) : __spans(__n, __value) {}
    template <typename Iterator>
    SpanVector(Iterator first, Iterator last) : __spans(first, last) {}
    
    void assign(size_type __n, const value_type& __value) { __spans.assign(__n, __value); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __spans.assign(first, last); }
    
    void assign(const std::string& line);
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    
    // insert/erase...
    iterator insert(iterator pos, const value_type& value) { return __spans.insert(pos, value); }
    void insert(iterator pos, size_type n, const value_type& value) { __spans.insert(pos, n, value); }
    template <typename Iterator>
    void insert(iterator pos, Iterator first, Iterator last) { __spans.insert(pos, first, last); }
    
    iterator erase(iterator pos) { return __spans.erase(pos); }
    iterator erase(iterator first, iterator last) { return __spans.erase(first, last); }
    
    void push_back(const value_type& value) { __spans.push_back(value); }
    void pop_back() { __spans.pop_back(); }
    
    void swap(SpanVector& __x) { __spans.swap(__x.__spans); }
    
    void clear() { __spans.clear(); }
    void reserve(size_type __n) { __spans.reserve(__n); }
    void resize(size_type __n) { __spans.resize(__n); }
    void resize(size_type __n, const value_type& __value) { __spans.resize(__n, __value); }

    size_type size() const { return __spans.size(); }
    bool empty() const { return __spans.empty(); }
    
    inline const_iterator begin() const { return __spans.begin(); }
    inline       iterator begin()       { return __spans.begin(); }
    
    inline const_iterator end() const { return __spans.end(); }
    inline       iterator end()       { return __spans.end(); }
    
    inline const_reverse_iterator rbegin() const { return __spans.rbegin(); }
    inline       reverse_iterator rbegin()       { return __spans.rbegin(); }
    
    inline const_reverse_iterator rend() const { return __spans.rend(); }
    inline       reverse_iterator rend()       { return __spans.rend(); }
    
    inline const_reference operator[](size_type pos) const { return __spans[pos]; }
    inline       reference operator[](size_type pos)       { return __spans[pos]; }

    inline const_reference front() const { return __spans.front(); }
    inline       reference front()       { return __spans.front(); }
    
    inline const_reference back() const { return __spans.back(); }
    inline       reference back()       { return __spans.back(); }

  public:
    
    friend
    size_t hash_value(SpanVector const& x);
    
    friend
    std::ostream& operator<<(std::ostream& os, const SpanVector& x);
    friend
    std::istream& operator>>(std::istream& is, SpanVector& x);

    friend
    std::ostream& operator<<(std::ostream& os, const span_type& x);
    friend
    std::istream& operator>>(std::istream& is, span_type& x);
    
    friend
    bool operator==(const SpanVector& x, const SpanVector& y);
    friend
    bool operator!=(const SpanVector& x, const SpanVector& y);
    friend
    bool operator<(const SpanVector& x, const SpanVector& y);
    friend
    bool operator>(const SpanVector& x, const SpanVector& y);
    friend
    bool operator<=(const SpanVector& x, const SpanVector& y);
    friend
    bool operator>=(const SpanVector& x, const SpanVector& y);
    
  private:
    spans_type __spans;
  };
  
  
  inline
  size_t hash_value(SpanVector const& x)
  { 
    return utils::hashmurmur<size_t>()(x.__spans.begin(), x.__spans.end(), 0);
  }

  inline
  bool operator==(const SpanVector::Span& x, const SpanVector::Span& y)
  {
    return x.first == y.first && x.last == y.last && x.label == y.label;
  }
  
  inline
  bool operator!=(const SpanVector::Span& x, const SpanVector::Span& y)
  {
    return x.first != y.first || x.last != y.last || x.label != y.label;
  }
  
  inline
  bool operator<(const SpanVector::Span& x, const SpanVector::Span& y)
  {
    return (x.first < y.first
	    || (! (y.first < x.first)
		&& (x.last < y.last
		    || (! (y.last < x.last)
			&& x.label < y.label))));
  }
  
  inline
  bool operator>(const SpanVector::Span& x, const SpanVector::Span& y) { return y < x; }
  
  inline
  bool operator<=(const SpanVector::Span& x, const SpanVector::Span& y) { return ! (y < x); }
  
  inline
  bool operator>=(const SpanVector::Span& x, const SpanVector::Span& y) { return ! (x < y); }

  inline
  bool operator==(const SpanVector& x, const SpanVector& y) { return x.__spans == y.__spans; }
  inline
  bool operator!=(const SpanVector& x, const SpanVector& y) { return x.__spans != y.__spans; }
  inline
  bool operator<(const SpanVector& x, const SpanVector& y) { return x.__spans < y.__spans; }
  inline
  bool operator>(const SpanVector& x, const SpanVector& y) { return x.__spans > y.__spans; }
  inline
  bool operator<=(const SpanVector& x, const SpanVector& y) { return x.__spans <= y.__spans; }
  inline
  bool operator>=(const SpanVector& x, const SpanVector& y) { return x.__spans >= y.__spans; }
  
};

namespace std
{
  inline
  void swap(cicada::SpanVector& x, cicada::SpanVector& y)
  {
    x.swap(y);
  }
};



#endif
