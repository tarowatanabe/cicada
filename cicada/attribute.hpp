// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__ATTRIBUTE__HPP__
#define __CICADA__ATTRIBUTE__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/thread.hpp>

#include <utils/indexed_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/spinlock.hpp>
#include <utils/piece.hpp>
#include <utils/chunk_vector.hpp>

namespace cicada
{
  struct AttributeImpl;

  class Attribute
  {
  private:
    friend struct AttributeImpl;
    
  public:
    typedef std::string  attribute_type;
    typedef utils::piece piece_type;
    typedef uint32_t     id_type;
    
    typedef attribute_type::size_type              size_type;
    typedef attribute_type::difference_type        difference_type;
    
    typedef attribute_type::const_iterator         const_iterator;
    typedef attribute_type::const_reverse_iterator const_reverse_iterator;
    typedef attribute_type::const_reference        const_reference;

  private:
    typedef utils::spinlock              mutex_type;
    typedef utils::spinlock::scoped_lock lock_type;
    
  public:
    Attribute() : __id(__allocate_empty()) { }
    Attribute(const utils::piece& x) : __id(__allocate(x)) { }
    Attribute(const attribute_type& x) : __id(__allocate(x)) { }
    Attribute(const char* x) : __id(__allocate(x)) { }
    Attribute(const id_type& x) : __id(x) {}
    template <typename Iterator>
    Attribute(Iterator first, Iterator last) : __id(__allocate(piece_type(first, last))) { }
    
    void assign(const piece_type& x) { __id = __allocate(x); }
    void assign(const attribute_type& x) { __id = __allocate(x); }
    void assign(const char* x) { __id = __allocate(x); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __id = __allocate(piece_type(first, last)); }
    
  public:
    void swap(Attribute& x) { std::swap(__id, x.__id); }
    
    const id_type   id() const { return __id; }
    operator const attribute_type&() const { return attribute(); }
    operator utils::piece() const { return attribute(); }
    
    const attribute_type& attribute() const
    {
      attribute_map_type& maps = __attribute_maps();
      
      if (__id >= maps.size())
	maps.resize(__id + 1, 0);
      if (! maps[__id]) {
	lock_type lock(__mutex_data);
	maps[__id] = &(__attributes()[__id]);
      }
      
      return *maps[__id];
    }
    
    const_iterator begin() const { return attribute().begin(); }
    const_iterator end() const { return attribute().end(); }
    
    const_reverse_iterator rbegin() const { return attribute().rbegin(); }
    const_reverse_iterator rend() const { return attribute().rend(); }
    
    const_reference operator[](size_type x) const { return attribute()[x]; }
    
    size_type size() const { return attribute().size(); }
    bool empty() const { return attribute().empty(); }
        
  public:
    // boost hash
    friend
    size_t  hash_value(Attribute const& x);
    
    // iostreams
    friend
    std::ostream& operator<<(std::ostream& os, const Attribute& x);
    friend
    std::istream& operator>>(std::istream& is, Attribute& x);
    
    // comparison...
    friend
    bool operator==(const Attribute& x, const Attribute& y);
    friend
    bool operator!=(const Attribute& x, const Attribute& y);
    friend
    bool operator<(const Attribute& x, const Attribute& y);
    friend
    bool operator>(const Attribute& x, const Attribute& y);
    friend
    bool operator<=(const Attribute& x, const Attribute& y);
    friend
    bool operator>=(const Attribute& x, const Attribute& y);

  private:
    typedef utils::indexed_set<piece_type, boost::hash<piece_type>, std::equal_to<piece_type>, std::allocator<piece_type> > attribute_index_type;
    typedef utils::chunk_vector<attribute_type, 4096 / sizeof(attribute_type), std::allocator<attribute_type> > attribute_set_type;
    typedef std::vector<const attribute_type*, std::allocator<const attribute_type*> > attribute_map_type;
    
  public:
    static bool exists(const piece_type& x)
    {
      lock_type lock(__mutex_index);
      
      const attribute_index_type& index = __index();
      
      return index.find(x) != index.end();
    }
    
    static size_t allocated()
    {
      lock_type lock(__mutex_data);
      return __attributes().size();
    }
    
  private:
    static mutex_type    __mutex_index;
    static mutex_type    __mutex_data;
    
    static attribute_map_type& __attribute_maps();
    
    static attribute_set_type& __attributes()
    {
      static attribute_set_type feats;
      return feats;
    }
    
    static attribute_index_type& __index()
    {
      static attribute_index_type index;
      return index;
    }
    
    static const id_type& __allocate_empty()
    {
      static const id_type __id = __allocate("");
      return __id;
    }
    
    static id_type __allocate(const piece_type& x)
    {
      lock_type lock(__mutex_index);
      
      attribute_index_type& index = __index();
      
      std::pair<attribute_index_type::iterator, bool> result = index.insert(x);
      
      if (result.second) {
	lock_type lock(__mutex_data);

	attribute_set_type& attributes = __attributes();
	attributes.push_back(attribute_type(x.begin(), x.end()));
	const_cast<piece_type&>(*result.first) = attributes.back();	
      }
      
      return result.first - index.begin();
    }
    
  private:
    id_type __id;
  };
  
  inline
  size_t hash_value(Attribute const& x)
  {
    return utils::hashmurmur<size_t>()(x.__id);
  }
  
  inline
  std::ostream& operator<<(std::ostream& os, const Attribute& x)
  {
    os << x.attribute();
    return os;
  }
  
  inline
  std::istream& operator>>(std::istream& is, Attribute& x)
  {
    std::string attribute;
    is >> attribute;
    x.assign(attribute);
    return is;
  }
  
  inline
   bool operator==(const Attribute& x, const Attribute& y)
  {
    return x.__id == y.__id;
  }
  inline
  bool operator!=(const Attribute& x, const Attribute& y)
  {
    return x.__id != y.__id;
  }
  inline
  bool operator<(const Attribute& x, const Attribute& y)
  {
    return x.__id < y.__id;
  }
  inline
  bool operator>(const Attribute& x, const Attribute& y)
  {
    return x.__id > y.__id;
  }
  inline
  bool operator<=(const Attribute& x, const Attribute& y)
  {
    return x.__id <= y.__id;
  }
  inline
  bool operator>=(const Attribute& x, const Attribute& y)
  {
    return x.__id >= y.__id;
  }

};

namespace std
{
  inline
  void swap(cicada::Attribute& x, cicada::Attribute& y)
  {
    x.swap(y);
  }
};

#endif

