// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARAMETER__HPP__
#define __CICADA__PARAMETER__HPP__ 1

#include <string>
#include <vector>
#include <iostream>

#include <utils/piece.hpp>

namespace cicada
{
  struct Parameter
  {
  public:
    typedef std::string attribute_type;
    typedef std::string key_type;
    typedef std::string data_type;
    typedef std::string mapped_type;
    typedef std::pair<std::string, std::string> value_type;

  private:
    typedef std::vector<value_type, std::allocator<value_type> > value_set_type;      

  public:
    typedef value_set_type::size_type       size_type;
    typedef value_set_type::difference_type difference_type;
    
    typedef value_set_type::const_iterator       iterator;
    typedef value_set_type::const_iterator const_iterator;
      
  public:
    Parameter(const std::string& parameter)
      : __attr(), __values()  { parse(parameter); }
    Parameter() : __attr(), __values() {}
    
    operator attribute_type() const { return __attr; }
    
    const attribute_type& name() const { return __attr; }
    attribute_type& name() { return __attr; }
    
    void push_back(const value_type& x) { __values.push_back(x); }
    
    const_iterator begin() const { return __values.begin(); }
    const_iterator end() const { return __values.end(); }
    
    bool empty() const { return __values.empty(); }
    size_type size() const { return __values.size(); }
    
    const_iterator find(const utils::piece& key) const
    {
      for (const_iterator iter = begin(); iter != end(); ++ iter)
	if (utils::piece(iter->first) == key)
	  return iter;
      return end();
    }

    void erase(const utils::piece& key) 
    {
      while (! __values.empty()) {
	bool found = false;
	for (value_set_type::iterator iter = __values.begin(); iter != __values.end(); ++ iter)
	  if (utils::piece(iter->first) == key) {
	    __values.erase(iter);
	    found = true;
	    break;
	  }
	
	if (! found) break;
      }
    }
    
  public:
    friend
    std::ostream& operator<<(std::ostream& os, const Parameter& x);
      
  private:
    void parse(const std::string& parameter);
      
  private:
    attribute_type __attr;
    value_set_type __values;
  };
  
};


#endif
