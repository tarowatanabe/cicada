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

namespace cicada
{
  struct AttributeImpl;

  class Attribute
  {
  private:
    friend struct AttributeImpl;
    
  public:
    typedef std::string attribute_type;
    typedef uint32_t    id_type;
    
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
    Attribute(Iterator first, Iterator last) : __id(__allocate(attribute_type(first, last))) { }
    
    void assign(const attribute_type& x) { __id = __allocate(x); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __id = __allocate(attribute_type(first, last)); }
    
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
	lock_type lock(__mutex);
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
    struct hasher : public utils::hashmurmur<size_t>
    {
      typedef utils::hashmurmur<size_t> hasher_type;

      size_t operator()(const attribute_type& x) const
      {
	return hasher_type::operator()(x.begin(), x.end(), size_t(0));
      }
    };

  private:
    typedef utils::indexed_set<attribute_type, hasher, std::equal_to<attribute_type>, std::allocator<attribute_type> > attribute_set_type;

    typedef std::vector<const attribute_type*, std::allocator<const attribute_type*> > attribute_map_type;
    
  public:
    static bool exists(const attribute_type& x)
    {
      lock_type lock(__mutex);
      return __attributes().find(x) != __attributes().end();
    }
    static size_t allocated()
    {
      lock_type lock(__mutex);
      return __attributes().size();
    }
    
  private:
    static mutex_type    __mutex;
    
    static attribute_map_type& __attribute_maps();
    
    static attribute_set_type& __attributes()
    {
      static attribute_set_type attrs;
      return attrs;
    }
    
    static const id_type& __allocate_empty()
    {
      static const id_type __id = __allocate(attribute_type());
      return __id;
    }
    
    static id_type __allocate(const attribute_type& x)
    {
      lock_type lock(__mutex);
      attribute_set_type::iterator siter = __attributes().insert(x).first;
      return siter - __attributes().begin();
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

