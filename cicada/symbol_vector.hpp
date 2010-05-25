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
    struct symbol_vector_impl_hash_type
    {
      utils::hashmurmur<size_t> __hasher;
      
      template <typename Iterator>
      size_t operator()(Iterator first, Iterator last) const
      {
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	
	return __hash_dispatch(first, last, value_type());
      }
      
      size_t operator()(const symbol_vector_impl_type& x) const
      {
	size_t seed = 0;
	symbol_vector_impl_type::const_iterator iter_end = x.end();
	for (symbol_vector_impl_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	  seed = __hasher(*iter, seed);
	return seed;
      }
      
    private:
      template <typename Iterator, typename _Word>
      size_t __hash_dispatch(Iterator first, Iterator last, _Word) const
      {
	size_t seed = 0;
	for (/**/; first != last; ++ first)
	  seed = __hasher(symbol_type(*first), seed);
	return seed;
      }

      template <typename Iterator>
      size_t __hash_dispatch(Iterator first, Iterator last, symbol_type) const
      {
	size_t seed = 0;
	for (/**/; first != last; ++ first)
	  seed = __hasher(*first, seed);
	return seed;
      }
    };
    
    struct symbol_vector_impl_equal_type
    {

      template <typename Iterator>
      bool operator()(Iterator first, Iterator last, const symbol_vector_impl_type& x) const
      {
	return x.size() == std::distance(first, last) && std::equal(x.begin(), x.end(), first);
      }

      template <typename Iterator>
      bool operator()(const symbol_vector_impl_type& x, Iterator first, Iterator last) const
      {
	return x.size() == std::distance(first, last) && std::equal(x.begin(), x.end(), first);
      }

      bool operator()(const symbol_vector_impl_type& x, const symbol_vector_impl_type& y) const
      {
	return x == y;
      }
    };

    typedef std::allocator<symbol_vector_impl_type> symbol_vector_impl_alloc_type;
    
    typedef utils::symbol_sequence<symbol_vector_impl_type, symbol_vector_impl_hash_type, symbol_vector_impl_equal_type, symbol_vector_impl_alloc_type> impl_type;
    
  public:
    typedef symbol_vector_impl_type::size_type              size_type;
    typedef symbol_vector_impl_type::difference_type        difference_type;
    
    typedef symbol_vector_impl_type::value_type             value_type;
    
    typedef symbol_vector_impl_type::const_iterator         const_iterator;
    typedef symbol_vector_impl_type::const_reverse_iterator const_reverse_iterator;
    typedef symbol_vector_impl_type::const_reference        const_reference;
    
    typedef impl_type::id_type id_type;

  public:
    SymbolVector() : __impl(__default()) {}
    SymbolVector(size_type size, const symbol_type& word) : __impl(symbol_vector_impl_type(size, word)) {}
    template <typename Iterator>
    SymbolVector(Iterator first, Iterator last) : __impl(__default())
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      __assign_dispatch(first, last, value_type());
    }

  public:
    void clear()
    {
      __impl = __default();
    }

    void assign(const SymbolVector& x)
    {
      __impl = x.__impl;
    }
    
    void assign(size_type size, const symbol_type& word)
    {
      __impl = impl_type(symbol_vector_impl_type(size, word));
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      __assign_dispatch(first, last, value_type());
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

    id_type id() const { return __impl.id(); }
    
    const_iterator begin() const { return symbol_vector().begin(); }
    const_iterator end() const { return symbol_vector().end(); }

    const_reverse_iterator rbegin() const { return symbol_vector().rbegin(); }
    const_reverse_iterator rend() const { return symbol_vector().rend(); }
    
    const_reference front() const { return symbol_vector().front(); }
    const_reference back() const { return symbol_vector().back(); }
    const_reference operator[](size_type pos) const { return symbol_vector().operator[](pos); }
    
    size_type size() const { return symbol_vector().size(); }
    bool empty() const { return symbol_vector().empty(); }

    

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
    const symbol_vector_impl_type& symbol_vector() const { return __impl.value(); }

  private:
    template <typename Iterator, typename _Word>
    void __assign_dispatch(Iterator first, Iterator last, _Word)
    {
      __impl = impl_type(symbol_vector_impl_type(first, last));
    }
    
    template <typename Iterator>
    void __assign_dispatch(Iterator first, Iterator last, symbol_type)
    {
      __impl = impl_type(first, last);
    }
    
    
  private:
    static const impl_type& __default()
    {
      static const impl_type __impl = impl_type();
      return __impl;
    }
    
  private:
    impl_type __impl;
  };
  
  inline
  size_t hash_value(SymbolVector const& x) { return hash_value(x.__impl); }
    
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
