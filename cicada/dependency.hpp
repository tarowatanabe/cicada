// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DEPENDENCY__HPP__
#define __CICADA__DEPENDENCY__HPP__ 1

// dependency implementation

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
  class Dependency
  {
  public:
    typedef int32_t index_type;
    
  private:
    typedef std::vector<index_type, std::allocator<index_type> > dep_type;

  public:
    typedef dep_type::value_type      value_type;
    typedef dep_type::size_type       size_type;
    typedef dep_type::difference_type difference_type;

    typedef dep_type::const_iterator         const_iterator;
    typedef dep_type::iterator               iterator;
    typedef dep_type::const_reverse_iterator const_reverse_iterator;
    typedef dep_type::reverse_iterator       reverse_iterator;
    typedef dep_type::const_reference        const_reference;
    typedef dep_type::reference              reference;

  public:
    Dependency() : __dep() {}
    Dependency(const utils::piece& x) : __dep() { assign(x); }
    Dependency(size_type __n) : __dep(__n) {}
    Dependency(size_type __n, const value_type& __value) : __dep(__n, __value) {}
    template <typename Iterator>
    Dependency(Iterator first, Iterator last) : __dep(first, last) {}
    
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    void assign(const utils::piece& line);
    void assign(size_type __n, const value_type& __value) { __dep.assign(__n, __value); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __dep.assign(first, last); }
    
    // insert/erase...
    iterator insert(iterator pos, const value_type& value) { return __dep.insert(pos, value); }
    void insert(iterator pos, size_type n, const value_type& value) { __dep.insert(pos, n, value); }
    template <typename Iterator>
    void insert(iterator pos, Iterator first, Iterator last) { __dep.insert(pos, first, last); }
    
    iterator erase(iterator pos) { return __dep.erase(pos); }
    iterator erase(iterator first, iterator last) { return __dep.erase(first, last); }
    
    void push_back(const value_type& value) { __dep.push_back(value); }
    void pop_back() { __dep.pop_back(); }
    
    void swap(Dependency& __x) { __dep.swap(__x.__dep); }
    
    void clear() { __dep.clear(); }
    void reserve(size_type __n) { __dep.reserve(__n); }
    void resize(size_type __n) { __dep.resize(__n); }
    void resize(size_type __n, const value_type& __value) { __dep.resize(__n, __value); }

    size_type size() const { return __dep.size(); }
    bool empty() const { return __dep.empty(); }
    
    inline const_iterator begin() const { return __dep.begin(); }
    inline       iterator begin()       { return __dep.begin(); }
    
    inline const_iterator end() const { return __dep.end(); }
    inline       iterator end()       { return __dep.end(); }
    
    inline const_reverse_iterator rbegin() const { return __dep.rbegin(); }
    inline       reverse_iterator rbegin()       { return __dep.rbegin(); }
    
    inline const_reverse_iterator rend() const { return __dep.rend(); }
    inline       reverse_iterator rend()       { return __dep.rend(); }
    
    inline const_reference operator[](size_type pos) const { return __dep[pos]; }
    inline       reference operator[](size_type pos)       { return __dep[pos]; }

    inline const_reference front() const { return __dep.front(); }
    inline       reference front()       { return __dep.front(); }
    
    inline const_reference back() const { return __dep.back(); }
    inline       reference back()       { return __dep.back(); }

  public:
    
    friend
    size_t hash_value(Dependency const& x);
    
    friend
    std::ostream& operator<<(std::ostream& os, const Dependency& x);
    friend
    std::istream& operator>>(std::istream& is, Dependency& x);
    
    friend
    bool operator==(const Dependency& x, const Dependency& y);
    friend
    bool operator!=(const Dependency& x, const Dependency& y);
    friend
    bool operator<(const Dependency& x, const Dependency& y);
    friend
    bool operator>(const Dependency& x, const Dependency& y);
    friend
    bool operator<=(const Dependency& x, const Dependency& y);
    friend
    bool operator>=(const Dependency& x, const Dependency& y);

    dep_type __dep;
  };

  inline
  size_t hash_value(Dependency const& x)
  { 
    return utils::hashmurmur<size_t>()(x.__dep.begin(), x.__dep.end(), 0);
  }

  inline
  bool operator==(const Dependency& x, const Dependency& y) { return x.__dep == y.__dep; }
  inline
  bool operator!=(const Dependency& x, const Dependency& y) { return x.__dep != y.__dep; }
  inline
  bool operator<(const Dependency& x, const Dependency& y) { return x.__dep < y.__dep; }
  inline
  bool operator>(const Dependency& x, const Dependency& y) { return x.__dep > y.__dep; }
  inline
  bool operator<=(const Dependency& x, const Dependency& y) { return x.__dep <= y.__dep; }
  inline
  bool operator>=(const Dependency& x, const Dependency& y) { return x.__dep >= y.__dep; }

};

namespace std
{
  inline
  void swap(cicada::Dependency& x, cicada::Dependency& y)
  {
    x.swap(y);
  }
};

#endif
