// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__ATTRIBUTE_VECTOR__HPP__
#define __CICADA__ATTRIBUTE_VECTOR__HPP__ 1

#include <map>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <cicada/attribute.hpp>

#include <boost/variant.hpp>

namespace cicada
{
  
  class AttributeVector
  {
  public:
    typedef cicada::Attribute attribute_type;
    typedef cicada::Attribute key_type;

    typedef int64_t     int_type;
    typedef double      float_type;
    typedef std::string string_type;
    
    typedef boost::variant<int_type, float_type, string_type> mapped_type;
    typedef mapped_type                                       data_type;
    
    typedef std::pair<const attribute_type, data_type> value_type;
    
  public:
    typedef std::map<key_type, data_type, std::less<key_type>, std::allocator<value_type> > attribute_vector_type;
    typedef AttributeVector self_type;
    
  public:
    typedef attribute_vector_type::size_type       size_type;
    typedef attribute_vector_type::difference_type difference_type;
    
    typedef attribute_vector_type::const_iterator  const_iterator;
    typedef attribute_vector_type::iterator        iterator;
    
    typedef attribute_vector_type::const_reverse_iterator  const_reverse_iterator;
    typedef attribute_vector_type::reverse_iterator        reverse_iterator;
    
    typedef attribute_vector_type::const_reference const_reference;
    typedef attribute_vector_type::reference       reference;
    
  public:
    AttributeVector() {}
    AttributeVector(const AttributeVector& x) : __values(x.__values) {}
    template <typename Iterator>
    AttributeVector(Iterator first, Iterator last) : __values(first, last) { }
    AttributeVector(const std::string& x) : __values() { assign(x); }

  public:
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    void assign(const std::string& x);

    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      __values.clear();
      __values.insert(first, last);
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      __values.insert(first, last);
    }
    
    std::pair<iterator, bool> insert(const value_type& x)
    {
      return __values.insert(x);
    }

    size_type size() const { return __values.size(); }
    bool empty() const { return __values.empty(); }
    
    void reserve(size_type x) { }
    
    void clear() { __values.clear(); }
    
    data_type& operator[](const key_type& x)
    {
      return __values[x];
    }
    
    const_iterator find(const key_type& x) const
    {
      return __values.find(x);
    }
    
    iterator find(const key_type& x)
    {
      return __values.find(x);
    }
    
    void erase(const key_type& x)
    {
      __values.erase(x);
    }

    void erase(iterator x)
    {
      __values.erase(x);
    }
    
    const_iterator begin() const { return __values.begin(); }
    iterator begin() { return __values.begin(); }
    const_iterator end() const { return __values.end(); }
    iterator end() { return __values.end(); }

    const_reverse_iterator rbegin() const { return __values.rbegin(); }
    reverse_iterator rbegin() { return __values.rbegin(); }
    const_reverse_iterator rend() const { return __values.rend(); }
    reverse_iterator rend() { return __values.rend(); }
    
    const_reference front() const { return *__values.begin(); }
    reference front() { return *__values.begin(); }
    
    const_reference back() const { return *(-- __values.end());}
    reference back() { return *(-- __values.end());}
    
    void swap(AttributeVector& x) { __values.swap(x.__values); }

  private:
    template <typename Tp, typename Iterator>
    struct __find_visitor : public boost::static_visitor<Iterator>
    {
      __find_visitor(Iterator __iter, Iterator __iter_end) : iter(__iter), iter_end(__iter_end) {}
      
      Iterator iter;
      Iterator iter_end;
      
      Iterator operator()(Tp& i) const
      {
	return iter;
      }
      
      template <typename T>
      Iterator operator()(T& x) const
      {
	return iter_end;
      }
    };

  public:
    const_iterator find_int(const key_type& x) const
    {
      const_iterator iter = __values.find(x);
      if (iter == __values.end())
	return iter;
      else
	return boost::apply_visitor(__find_visitor<int_type, const_iterator>(iter, __values.end()), iter->second);
    }
    
    iterator find_int(const key_type& x)
    {
      iterator iter = __values.find(x);
      if (iter == __values.end())
	return iter;
      else
	return boost::apply_visitor(__find_visitor<int_type, iterator>(iter, __values.end()), iter->second);
    }
    

    const_iterator find_float(const key_type& x) const
    {
      const_iterator iter = __values.find(x);
      
      if (iter == __values.end())
	return iter;
      else
	return boost::apply_visitor(__find_visitor<float_type, const_iterator>(iter, __values.end()), iter->second);
    }

    iterator find_float(const key_type& x)
    {
      iterator iter = __values.find(x);
      
      if (iter == __values.end())
	return iter;
      else
	return boost::apply_visitor(__find_visitor<float_type, iterator>(iter, __values.end()), iter->second);
    }

    const_iterator find_string(const key_type& x) const
    {
      const_iterator iter = __values.find(x);
      
      if (iter == __values.end())
	return iter;
      else
	return boost::apply_visitor(__find_visitor<string_type, const_iterator>(iter, __values.end()), iter->second);
    }

    iterator find_string(const key_type& x)
    {
      iterator iter = __values.find(x);
      
      if (iter == __values.end())
	return iter;
      else
	return boost::apply_visitor(__find_visitor<string_type, iterator>(iter, __values.end()), iter->second);
    }
    
  public:
    // comparison
    friend
    bool operator==(const AttributeVector& x, const AttributeVector& y)
    {
      return x.__values == y.__values;
    }

    friend
    bool operator!=(const AttributeVector& x, const AttributeVector& y)
    {
      return x.__values != y.__values;
    }

    friend
    bool operator<(const AttributeVector& x, const AttributeVector& y)
    {
      return x.__values < y.__values;
    }

    friend
    bool operator<=(const AttributeVector& x, const AttributeVector& y)
    {
      return x.__values <= y.__values;
    }

    friend
    bool operator>(const AttributeVector& x, const AttributeVector& y)
    {
      return x.__values > y.__values;
    }
    
    friend
    bool operator>=(const AttributeVector& x, const AttributeVector& y)
    {
      return x.__values >= y.__values;
    }


  public:
    friend
    std::ostream& operator<<(std::ostream& os, const AttributeVector& x);
    
    friend
    std::istream& operator>>(std::istream& is, AttributeVector& x);
    
  private:
    attribute_vector_type __values;
  };
  
};

namespace std
{
  inline
  void swap(cicada::AttributeVector& x, cicada::AttributeVector& y)
  {
    x.swap(y);
  }
};

#endif

