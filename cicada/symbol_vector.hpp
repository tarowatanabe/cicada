// -*- mode: c++ -*-

#ifndef __CICADA__SYMBOL_VECTOR__HPP__
#define __CICADA__SYMBOL_VECTOR__HPP__ 1

// reference counted non-thread-safe phrase implementation
// Use with care if you mix with threaded environment... (locking etc.)

#include <iostream>

#include <cicada/symbol.hpp>

#include <utils/simple_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/symbol_sequence.hpp>

namespace cicada
{
  class SymbolVector
  {
  public:
    typedef cicada::Symbol  symbol_type;
    
  private:
    typedef utils::simple_vector<symbol_type, std::allocator<symbol_type> > symbol_vector_impl_type;

    typedef symbol_vector_impl_type impl_type;
    
  public:
    typedef symbol_vector_impl_type::size_type              size_type;
    typedef symbol_vector_impl_type::difference_type        difference_type;
    
    typedef symbol_vector_impl_type::value_type             value_type;
    
    typedef symbol_vector_impl_type::const_iterator         const_iterator;
    typedef symbol_vector_impl_type::iterator               iterator;
    typedef symbol_vector_impl_type::const_reverse_iterator const_reverse_iterator;
    typedef symbol_vector_impl_type::reverse_iterator       reverse_iterator;
    typedef symbol_vector_impl_type::const_reference        const_reference;
    typedef symbol_vector_impl_type::reference              reference;

  public:
    SymbolVector() : __impl() {}
    SymbolVector(size_type size, const symbol_type& word) : __impl(size, word) {}
    template <typename Iterator>
    SymbolVector(Iterator first, Iterator last) : __impl(first, last) {}

  public:
    void clear()
    {
      __impl.clear();
    }

    void assign(const SymbolVector& x)
    {
      __impl = x.__impl;
    }
    
    void assign(size_type size, const symbol_type& word)
    {
      __impl.assign(size, word);
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      __impl.assign(first, last);
    }

    void swap(SymbolVector& x) { __impl.swap(x.__impl); }
    
  public:
    size_type arity() const
    {
      size_type count = 0;
      const_iterator iter_end = end();
      for (const_iterator iter = begin(); iter != iter_end; ++ iter)
	count += iter->is_non_terminal();
      return count;
    }
    
    const_iterator begin() const { return __impl.begin(); }
    iterator begin() { return __impl.begin(); }
    const_iterator end() const { return __impl.end(); }
    iterator end() { return __impl.end(); }
    
    const_reverse_iterator rbegin() const { return __impl.rbegin(); }
    reverse_iterator rbegin() { return __impl.rbegin(); }
    const_reverse_iterator rend() const { return __impl.rend(); }
    reverse_iterator rend() { return __impl.rend(); }
    
    const_reference front() const { return __impl.front(); }
    reference front() { return __impl.front(); }
    const_reference back() const { return __impl.back(); }
    reference back() { return __impl.back(); }
    const_reference operator[](size_type pos) const { return __impl.operator[](pos); }
    reference operator[](size_type pos) { return __impl.operator[](pos); }
    
    size_type size() const { return __impl.size(); }
    bool empty() const { return __impl.empty(); }
    
    
    template <typename ResultIterator>
    void terminals(ResultIterator riter) const
    {
      const_iterator first = begin();
      const_iterator last = end();
      const_iterator iter = first;
      for (/**/; iter != last && iter->is_terminal(); ++ iter);
      *riter = std::make_pair(first, iter);
      ++ riter;
      first = iter;
      
      while (first != last) {
	++ first;
	
	const_iterator iter = first;
	for (/**/; iter != last && iter->is_terminal(); ++ iter);
	*riter = std::make_pair(first, iter);
	++ riter;
	first = iter;
      }
    }


  public:
    friend
    size_t hash_value(SymbolVector const& x);
    
    friend
    std::ostream& operator<<(std::ostream& os, const SymbolVector& x);
    friend
    std::istream& operator>>(std::istream& is, SymbolVector& x);
    
    friend
    bool operator==(const SymbolVector& x, const SymbolVector& y);
    friend
    bool operator!=(const SymbolVector& x, const SymbolVector& y);
    friend
    bool operator<(const SymbolVector& x, const SymbolVector& y);
    friend
    bool operator>(const SymbolVector& x, const SymbolVector& y);
    friend
    bool operator<=(const SymbolVector& x, const SymbolVector& y);
    friend
    bool operator>=(const SymbolVector& x, const SymbolVector& y);

    
    
  private:
    impl_type __impl;
  };
  
  inline
  size_t hash_value(SymbolVector const& x)
  {
    return utils::hashmurmur<size_t>()(x.begin(), x.end(), 0);
  }
    
  inline
  bool operator==(const SymbolVector& x, const SymbolVector& y) { return x.__impl == y.__impl; }
  inline
  bool operator!=(const SymbolVector& x, const SymbolVector& y) { return x.__impl != y.__impl; }
  inline
  bool operator<(const SymbolVector& x, const SymbolVector& y) { return x.__impl < y.__impl; }
  inline
  bool operator>(const SymbolVector& x, const SymbolVector& y) { return x.__impl > y.__impl; }
  inline
  bool operator<=(const SymbolVector& x, const SymbolVector& y) { return x.__impl <= y.__impl; }
  inline
  bool operator>=(const SymbolVector& x, const SymbolVector& y) { return x.__impl >= y.__impl; }

  
  
};

namespace std
{
  inline
  void swap(cicada::SymbolVector& x, cicada::SymbolVector& y)
  {
    x.swap(y);
  }
};

#endif
