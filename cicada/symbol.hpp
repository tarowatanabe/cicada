// -*- mode: c++ -*-

#ifndef __CICADA__SYMBOL__HPP__
#define __CICADA__SYMBOL__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>

#include <utils/indexed_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/spinlock.hpp>

namespace cicada
{
  class Symbol
  {
  public:
    typedef std::string symbol_type;
    typedef uint32_t    id_type;
    
    typedef symbol_type::size_type              size_type;
    typedef symbol_type::difference_type        difference_type;
    
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
    Symbol(const symbol_type& x) : __id(__allocate(x)) { }
    Symbol(const char* x) : __id(__allocate(x)) { }
    Symbol(const id_type& x) : __id(x) {}
    template <typename Iterator>
    Symbol(Iterator first, Iterator last) : __id(__allocate(symbol_type(first, last))) { }
    
    void assign(const symbol_type& x) { __id = __allocate(x); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __id = __allocate(symbol_type(first, last)); }
    
  public:
    void swap(Symbol& x) { std::swap(__id, x.__id); }
    
    const id_type   id() const { return __id; }
    operator const symbol_type&() const { return symbol(); }
    
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
    bool is_non_terminal() const
    {
      const size_type scan_pos = (__id << 1);
      const size_type flag_pos = (__id << 1) + 1;

      non_terminal_map_type& maps = __non_terminal_maps();
      if (flag_pos >= maps.size())
	maps.resize(flag_pos + 1, false);
      
      if (! maps[scan_pos]) {
	const symbol_type& word = symbol();
	const size_type size = word.size();
	
	// at least we have [<char>]
	maps[scan_pos] = true;
	maps[flag_pos] = (size >= 3 && word[0] == '[' && word[size - 1] == ']');
      }
      
      return maps[flag_pos];
    }
    
    int non_terminal_index() const
    {
      if (! is_non_terminal())
	return 0;
      
      index_map_type& maps = __index_maps();
      if (__id >= maps.size())
	maps.resize(__id + 1, -1);

      if (maps[__id] < 0) {
	const symbol_type& word = symbol();
	const size_type size = word.size();
	
	// at leas we have [<char>,<digit>]
	maps[__id] = 0;
	if (size >= 5 && word[0] == '[' && word[size - 1] == ']') {
	  // first-index and # of char...
	  symbol_type::size_type pos = word.find(',');
	  maps[__id] = (pos != symbol_type::npos ? atoi(word.c_str() + pos + 1) : 0);
	}
      }
      
      return maps[__id];
    }
    
    Symbol non_terminal() const
    {
      if (! is_non_terminal())
	return *this;
      
      non_terminal_symbol_map_type& maps = __non_terminal_symbol_maps();
      if (__id >= maps.size())
	maps.resize(__id + 1, id_type(-1));
      
      if (maps[__id] == id_type(-1)) {
	const symbol_type& word = symbol();
	const size_type size = word.size();
	
	symbol_type::size_type pos = word.find(',');
	
	maps[__id] = __id;
	if (pos != symbol_type::npos)
	  maps[__id] = Symbol(word.substr(0, pos) + ']').id();
      }
      
      return maps[__id];
    }

    Symbol non_terminal(const int index) const
    {
      if (! is_non_terminal() || index < 0)
	return *this;
      
      if (index == 0)
	return non_terminal();
      
      const symbol_type& word = symbol();
      const symbol_type::size_type pos = word.find(',');
      
      if (pos == symbol_type::npos)
	return Symbol(word.substr(0, word.size() - 1) + ',' + boost::lexical_cast<std::string>(index) + ']');
      else
	return Symbol(word.substr(0, pos) + ',' + boost::lexical_cast<std::string>(index) + ']');
    }
    
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
    struct hasher
    {
      size_t operator()(const symbol_type& x) const
      {
	return __hasher(x.begin(), x.end(), size_t(0));
      }
      utils::hashmurmur<size_t> __hasher;
    };
    typedef utils::indexed_set<symbol_type, hasher, std::equal_to<symbol_type>, std::allocator<symbol_type> > symbol_set_type;

    typedef std::vector<const symbol_type*, std::allocator<const symbol_type*> > symbol_map_type;
    typedef std::vector<int, std::allocator<int> >   index_map_type;
    typedef std::vector<bool, std::allocator<bool> > non_terminal_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > non_terminal_symbol_map_type;

    
  public:
    static bool exists(const symbol_type& x)
    {
      lock_type lock(__mutex);
      return __symbols().find(x) == __symbols().end();
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
    static index_map_type& __index_maps();
    static non_terminal_map_type& __non_terminal_maps();
    static non_terminal_symbol_map_type& __non_terminal_symbol_maps();
    
    static symbol_set_type& __symbols()
    {
      static symbol_set_type symbols;
      return symbols;
    }
    
    static const id_type& __allocate_empty()
    {
      static const id_type __id = __allocate(symbol_type());
      return __id;
    }
    
    static id_type __allocate(const symbol_type& x)
    {
      lock_type lock(__mutex);
      symbol_set_type::iterator siter = __symbols().insert(x).first;
      return siter - __symbols().begin();
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

#endif

