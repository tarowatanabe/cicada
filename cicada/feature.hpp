// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__HPP__
#define __CICADA__FEATURE__HPP__ 1

#include <stdint.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <cicada/symbol.hpp>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include <utils/indexed_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/spinlock.hpp>


namespace cicada
{
  class Feature
  {
  public:
    typedef cicada::Symbol symbol_type;
    
  public:
    typedef std::string feature_type;
    typedef uint32_t    id_type;
    
    typedef feature_type::size_type              size_type;
    typedef feature_type::difference_type        difference_type;
    
    typedef feature_type::const_iterator         const_iterator;
    typedef feature_type::const_reverse_iterator const_reverse_iterator;
    typedef feature_type::const_reference        const_reference;

    typedef boost::filesystem::path path_type;
    
  private:
    typedef utils::spinlock              mutex_type;
    typedef utils::spinlock::scoped_lock lock_type;
    
  public:
    Feature() : __id(__allocate_empty()) { }
    Feature(const symbol_type& x) : __id(__allocate(static_cast<const std::string&>(x))) {}
    Feature(const feature_type& x) : __id(__allocate(x)) { }
    Feature(const id_type& x) : __id(x) {}
    template <typename Iterator>
    Feature(Iterator first, Iterator last) : __id(__allocate(feature_type(first, last))) { }
    
    void assign(const feature_type& x) { __id = __allocate(x); }
    template <typename Iterator>
    void assign(Iterator first, Iterator last) { __id = __allocate(feature_type(first, last)); }
    
  public:
    void swap(Feature& x) { std::swap(__id, x.__id); }
    
    const id_type   id() const { return __id; }
    operator const feature_type&() const { return feature(); }
    
    const feature_type& feature() const
    {
      feature_map_type& maps = __feature_maps();
      
      if (__id >= maps.size())
	maps.resize(__id + 1, 0);
      if (! maps[__id]) {
	lock_type lock(__mutex);
	maps[__id] = &(__features()[__id]);
      }
      
      return *maps[__id];
    }
    
    const_iterator begin() const { return feature().begin(); }
    const_iterator end() const { return feature().end(); }
    
    const_reverse_iterator rbegin() const { return feature().rbegin(); }
    const_reverse_iterator rend() const { return feature().rend(); }
    
    const_reference operator[](size_type x) const { return feature()[x]; }
    
    size_type size() const { return feature().size(); }
    bool empty() const { return feature().empty(); }
        
  public:
    // boost hash
    friend
    size_t  hash_value(Feature const& x);
    
    // iostreams
    friend
    std::ostream& operator<<(std::ostream& os, const Feature& x);
    friend
    std::istream& operator>>(std::istream& is, Feature& x);
    
    // comparison...
    friend
    bool operator==(const Feature& x, const Feature& y);
    friend
    bool operator!=(const Feature& x, const Feature& y);
    friend
    bool operator<(const Feature& x, const Feature& y);
    friend
    bool operator>(const Feature& x, const Feature& y);
    friend
    bool operator<=(const Feature& x, const Feature& y);
    friend
    bool operator>=(const Feature& x, const Feature& y);
    
  private:
    struct hasher
    {
      size_t operator()(const feature_type& x) const
      {
	return __hasher(x.begin(), x.end(), size_t(0));
      }
      utils::hashmurmur<size_t> __hasher;
    };
    typedef utils::indexed_set<feature_type, hasher, std::equal_to<feature_type>, std::allocator<feature_type> > feature_set_type;

    typedef std::vector<const feature_type*, std::allocator<const feature_type*> > feature_map_type;
    
  public:
    static bool exists(const feature_type& x)
    {
      lock_type lock(__mutex);
      return __features().find(x) == __features().end();
    }
    static size_t allocated()
    {
      lock_type lock(__mutex);
      return __features().size();
    }
    
  private:
    static mutex_type    __mutex;
    
    static feature_map_type& __feature_maps();
    
    static feature_set_type& __features()
    {
      static feature_set_type features;
      return features;
    }
    
    static const id_type& __allocate_empty()
    {
      static const id_type __id = __allocate(feature_type());
      return __id;
    }
    
    static id_type __allocate(const feature_type& x)
    {
      lock_type lock(__mutex);
      feature_set_type::iterator siter = __features().insert(x).first;
      return siter - __features().begin();
    }
    
  private:
    id_type __id;
  };
  
  inline
  size_t hash_value(Feature const& x)
  {
    return utils::hashmurmur<size_t>()(x.__id);
  }
  
  inline
  std::ostream& operator<<(std::ostream& os, const Feature& x)
  {
    os << x.feature();
    return os;
  }
  
  inline
  std::istream& operator>>(std::istream& is, Feature& x)
  {
    std::string feature;
    is >> feature;
    x.assign(feature);
    return is;
  }
  
  inline
   bool operator==(const Feature& x, const Feature& y)
  {
    return x.__id == y.__id;
  }
  inline
  bool operator!=(const Feature& x, const Feature& y)
  {
    return x.__id != y.__id;
  }
  inline
  bool operator<(const Feature& x, const Feature& y)
  {
    return x.__id < y.__id;
  }
  inline
  bool operator>(const Feature& x, const Feature& y)
  {
    return x.__id > y.__id;
  }
  inline
  bool operator<=(const Feature& x, const Feature& y)
  {
    return x.__id <= y.__id;
  }
  inline
  bool operator>=(const Feature& x, const Feature& y)
  {
    return x.__id >= y.__id;
  }

};

#endif

