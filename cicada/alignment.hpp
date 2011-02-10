// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__ALIGNMENT__HPP__
#define __CICADA__ALIGNMENT__HPP__ 1

// an alignment implementation.
// The major difference to other implementation, such as those in
// nicttm is that point is sorted by target-side ordering....

#include <stdint.h>

#include <cstdlib>

#include <stdexcept>
#include <iostream>
#include <vector>
#include <utility>
#include <string>

#include <utils/hashmurmur.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  class Alignment
  {
  public:
    typedef int32_t index_type;
    struct Point
    {
      index_type source;
      index_type target;
      
      Point() : source(), target() {}
      template <typename Integral>
      Point(const std::pair<Integral, Integral>& x) : source(x.first), target(x.second) {}
      Point(const utils::piece& x) { assign(x); }
      Point(const index_type& _source, const index_type& _target) : source(_source), target(_target) {}

      template <typename Integral>
      void assign(const std::pair<Integral, Integral>& x)
      {
	source = x.first;
	target = x.second;
      }
      void assign(const Point& x)
      {
	source = x.source;
	target = x.target;
      }
      void assign(const utils::piece& x);

      friend
      size_t hash_value(Point const& x);
      friend
      std::ostream& operator<<(std::ostream& os, const Point& x);
      friend
      std::istream& operator>>(std::istream& is, Point& x);
    };
    typedef Point point_type;
    
  private:
    typedef std::vector<point_type, std::allocator<point_type> > align_type;
    
  public:
    typedef align_type::value_type      value_type;
    typedef align_type::size_type       size_type;
    typedef align_type::difference_type difference_type;

    typedef align_type::const_iterator         const_iterator;
    typedef align_type::iterator               iterator;
    typedef align_type::const_reverse_iterator const_reverse_iterator;
    typedef align_type::reverse_iterator       reverse_iterator;
    typedef align_type::const_reference        const_reference;
    typedef align_type::reference              reference;

  public:
    Alignment() : __align() {}
    Alignment(const utils::piece& x) : __align() { assign(x); }
    Alignment(size_type __n) : __align(__n) {}
    Alignment(size_type __n, const value_type& __value) : __align(__n, __value) {}
    template <typename Iterator>
    Alignment(Iterator first, Iterator last) : __align(first, last) {}
    
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    void assign(const utils::piece& line);
    void assign(size_type __n, const value_type& __value) { __align.assign(__n, __value); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __align.assign(first, last); }
    
    // insert/erase...
    iterator insert(iterator pos, const value_type& value) { return __align.insert(pos, value); }
    void insert(iterator pos, size_type n, const value_type& value) { __align.insert(pos, n, value); }
    template <typename Iterator>
    void insert(iterator pos, Iterator first, Iterator last) { __align.insert(pos, first, last); }
    
    iterator erase(iterator pos) { return __align.erase(pos); }
    iterator erase(iterator first, iterator last) { return __align.erase(first, last); }
    
    void push_back(const value_type& value) { __align.push_back(value); }
    void pop_back() { __align.pop_back(); }
    
    void swap(Alignment& __x) { __align.swap(__x.__align); }
    
    void clear() { __align.clear(); }
    void reserve(size_type __n) { __align.reserve(__n); }
    void resize(size_type __n) { __align.resize(__n); }
    void resize(size_type __n, const value_type& __value) { __align.resize(__n, __value); }

    size_type size() const { return __align.size(); }
    bool empty() const { return __align.empty(); }


    void inverse()
    {
      // compute inversed alignment...
      iterator iter_end = __align.end();
      for (iterator iter = __align.begin(); iter != iter_end; ++ iter)
	*iter = point_type(iter->target, iter->source);
    }
    
    inline const_iterator begin() const { return __align.begin(); }
    inline       iterator begin()       { return __align.begin(); }
    
    inline const_iterator end() const { return __align.end(); }
    inline       iterator end()       { return __align.end(); }
    
    inline const_reverse_iterator rbegin() const { return __align.rbegin(); }
    inline       reverse_iterator rbegin()       { return __align.rbegin(); }
    
    inline const_reverse_iterator rend() const { return __align.rend(); }
    inline       reverse_iterator rend()       { return __align.rend(); }
    
    inline const_reference operator[](size_type pos) const { return __align[pos]; }
    inline       reference operator[](size_type pos)       { return __align[pos]; }

    inline const_reference front() const { return __align.front(); }
    inline       reference front()       { return __align.front(); }
    
    inline const_reference back() const { return __align.back(); }
    inline       reference back()       { return __align.back(); }

  public:
    
    friend
    size_t hash_value(Alignment const& x);
    
    friend
    std::ostream& operator<<(std::ostream& os, const Alignment& x);
    friend
    std::istream& operator>>(std::istream& is, Alignment& x);
    
    friend
    bool operator==(const Alignment& x, const Alignment& y);
    friend
    bool operator!=(const Alignment& x, const Alignment& y);
    friend
    bool operator<(const Alignment& x, const Alignment& y);
    friend
    bool operator>(const Alignment& x, const Alignment& y);
    friend
    bool operator<=(const Alignment& x, const Alignment& y);
    friend
    bool operator>=(const Alignment& x, const Alignment& y);
    
  private:
    align_type __align;
  };
  
  inline
  size_t hash_value(Alignment::point_type const& x)
  {
    return utils::hashmurmur<size_t>()(x, 0);
  }
  
  inline
  size_t hash_value(Alignment const& x)
  { 
    return utils::hashmurmur<size_t>()(x.__align.begin(), x.__align.end(), 0);
  }
  
  inline
  bool operator==(const Alignment::point_type& x, const Alignment::point_type& y) { return x.source == y.source && x.target == y.target; }
  inline
  bool operator!=(const Alignment::point_type& x, const Alignment::point_type& y) { return x.source != y.source || x.target != y.target; }
  inline
  bool operator<(const Alignment::point_type& x, const Alignment::point_type& y) { return (x.source < y.source || (!(y.source < x.source) && x.target < y.target)); }
  inline
  bool operator>(const Alignment::point_type& x, const Alignment::point_type& y) { return y < x; }
  inline
  bool operator<=(const Alignment::point_type& x, const Alignment::point_type& y) { return ! (y < x); }
  inline
  bool operator>=(const Alignment::point_type& x, const Alignment::point_type& y) { return ! (x < y); }
  
  
  inline
  bool operator==(const Alignment& x, const Alignment& y) { return x.__align == y.__align; }
  inline
  bool operator!=(const Alignment& x, const Alignment& y) { return x.__align != y.__align; }
  inline
  bool operator<(const Alignment& x, const Alignment& y) { return x.__align < y.__align; }
  inline
  bool operator>(const Alignment& x, const Alignment& y) { return x.__align > y.__align; }
  inline
  bool operator<=(const Alignment& x, const Alignment& y) { return x.__align <= y.__align; }
  inline
  bool operator>=(const Alignment& x, const Alignment& y) { return x.__align >= y.__align; }
  
};

namespace std
{
  inline
  void swap(cicada::Alignment& x, cicada::Alignment& y)
  {
    x.swap(y);
  }
};
  
#endif
