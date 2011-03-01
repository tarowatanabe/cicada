// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SYMBOL__HPP__
#define __CICADA__SYMBOL__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include <utils/indexed_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/spinlock.hpp>
#include <utils/piece.hpp>
#include <utils/chunk_vector.hpp>

namespace cicada
{
  
  struct SymbolImpl;

  class Symbol
  {
  private:
    friend struct SymbolImpl;

  public:
    typedef std::string  symbol_type;
    typedef utils::piece piece_type;
    typedef uint32_t     id_type;
    
    typedef symbol_type::size_type              size_type;
    typedef symbol_type::difference_type        difference_type;
    
    typedef symbol_type::value_type             value_type;
    typedef symbol_type::const_iterator         const_iterator;
    typedef symbol_type::const_reverse_iterator const_reverse_iterator;
    typedef symbol_type::const_reference        const_reference;

    typedef boost::filesystem::path path_type;
    
  private:
    typedef utils::spinlock             mutex_type;
    typedef mutex_type::scoped_lock     lock_type;
    typedef mutex_type::scoped_try_lock trylock_type;
    
  public:
    Symbol() : __id(__allocate_empty()) { }
    Symbol(const piece_type& x) : __id(__allocate(x)) { }
    Symbol(const symbol_type& x) : __id(__allocate(x)) { }
    Symbol(const char* x) : __id(__allocate(x)) { }
    Symbol(const id_type& x) : __id(x) {}
    template <typename Iterator>
    Symbol(Iterator first, Iterator last) : __id(__allocate(piece_type(first, last))) { }
    
    void assign(const piece_type& x) { __id = __allocate(x); }
    void assign(const symbol_type& x) { __id = __allocate(x); }
    void assign(const char* x) { __id = __allocate(x); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __id = __allocate(piece_type(first, last)); }
    
  public:
    void swap(Symbol& x) { std::swap(__id, x.__id); }
    
    const id_type   id() const { return __id; }
    operator const symbol_type&() const { return symbol(); }
    operator utils::piece() const { return symbol(); }
    
    const symbol_type& symbol() const
    {
      symbol_map_type& maps = __symbol_maps();
      
      if (__id >= maps.size())
	maps.resize(__id + 1, 0);
      if (! maps[__id]) {
	lock_type lock(__mutex);
	maps[__id] = &(__symbols()[__id]);
      }
      
      return *maps[__id];
    }
    
    const_iterator begin() const { return symbol().begin(); }
    const_iterator end() const { return symbol().end(); }
    
    const_reverse_iterator rbegin() const { return symbol().rbegin(); }
    const_reverse_iterator rend() const { return symbol().rend(); }
    
    const_reference operator[](size_type x) const { return symbol()[x]; }
    
    size_type size() const { return symbol().size(); }
    bool empty() const { return symbol().empty(); }
    
    
    bool is_terminal() const { return ! is_non_terminal(); }
    bool is_non_terminal() const;
    
    
    piece_type non_terminal_strip() const;
    int non_terminal_index() const;
    Symbol non_terminal() const;
    Symbol non_terminal(const int index) const;
    
  public:
    // boost hash
    friend
    size_t  hash_value(Symbol const& x);
    
    // iostreams
    friend
    std::ostream& operator<<(std::ostream& os, const Symbol& x);
    friend
    std::istream& operator>>(std::istream& is, Symbol& x);
    
    // comparison...
    friend
    bool operator==(const Symbol& x, const Symbol& y);
    friend
    bool operator!=(const Symbol& x, const Symbol& y);
    friend
    bool operator<(const Symbol& x, const Symbol& y);
    friend
    bool operator>(const Symbol& x, const Symbol& y);
    friend
    bool operator<=(const Symbol& x, const Symbol& y);
    friend
    bool operator>=(const Symbol& x, const Symbol& y);
    
    
  private:
    typedef utils::indexed_set<piece_type, boost::hash<piece_type>, std::equal_to<piece_type>, std::allocator<piece_type> > symbol_index_type;
    typedef utils::chunk_vector<symbol_type, 4096 / sizeof(symbol_type), std::allocator<symbol_type> > symbol_set_type;
    typedef std::vector<const symbol_type*, std::allocator<const symbol_type*> > symbol_map_type;
    
  public:
    static bool exists(const piece_type& x)
    {
      lock_type lock(__mutex);
      
      const symbol_index_type& index = __index();
      
      return index.find(x) != index.end();
    }
    static size_t allocated()
    {
      lock_type lock(__mutex);
      
      return __symbols().size();
    }
    static void write(const path_type& path);
    
  private:
    static mutex_type    __mutex;
    
    static symbol_map_type& __symbol_maps();
    
    static symbol_set_type& __symbols()
    {
      static symbol_set_type syms;
      return syms;
    }
    
    static symbol_index_type& __index()
    {
      static symbol_index_type index;
      return index;
    }
    
    static const id_type& __allocate_empty()
    {
      static const id_type __id = __allocate("");
      return __id;
    }
    
    static id_type __allocate(const piece_type& x)
    {
      lock_type lock(__mutex);

      symbol_index_type& index = __index();
      
      std::pair<symbol_index_type::iterator, bool> result = index.insert(x);
      
      if (result.second) {
	symbol_set_type& symbols = __symbols();
	symbols.push_back(x);
	const_cast<piece_type&>(*result.first) = symbols.back();
      }
      
      return result.first - index.begin();
    }
    
  private:
    id_type __id;
  };
  
  inline
  size_t hash_value(Symbol const& x)
  {
    return utils::hashmurmur<size_t>()(x.__id);
  }
  
  inline
  std::ostream& operator<<(std::ostream& os, const Symbol& x)
  {
    os << x.symbol();
    return os;
  }
  
  inline
  std::istream& operator>>(std::istream& is, Symbol& x)
  {
    std::string symbol;
    is >> symbol;
    x.assign(symbol);
    return is;
  }
  
  inline
   bool operator==(const Symbol& x, const Symbol& y)
  {
    return x.__id == y.__id;
  }
  inline
  bool operator!=(const Symbol& x, const Symbol& y)
  {
    return x.__id != y.__id;
  }
  inline
  bool operator<(const Symbol& x, const Symbol& y)
  {
    return x.__id < y.__id;
  }
  inline
  bool operator>(const Symbol& x, const Symbol& y)
  {
    return x.__id > y.__id;
  }
  inline
  bool operator<=(const Symbol& x, const Symbol& y)
  {
    return x.__id <= y.__id;
  }
  inline
  bool operator>=(const Symbol& x, const Symbol& y)
  {
    return x.__id >= y.__id;
  }

};

namespace std
{
  inline
  void swap(cicada::Symbol& x, cicada::Symbol& y)
  {
    x.swap(y);
  }
};

#endif

