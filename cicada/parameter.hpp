// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARAMETER__HPP__
#define __CICADA__PARAMETER__HPP__ 1

#include <string>
#include <vector>

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
    typedef value_set_type::const_iterator       iterator;
    typedef value_set_type::const_iterator const_iterator;
      
  public:
    Parameter(const std::string& parameter)
      : __attr(), __values()  { parse(parameter); }
      
    operator attribute_type() const { return __attr; }

    const attribute_type& name() const { return __attr; }
    const_iterator begin() const { return __values.begin(); }
    const_iterator end() const { return __values.end(); }
    bool empty() const { return __values.empty(); }

    const_iterator find(const std::string& key) const
    {
      for (const_iterator iter = begin(); iter != end(); ++ iter)
	if (iter->first == key)
	  return iter;
      return end();
    }
      
  private:
    void parse(const std::string& parameter);
      
  private:
    attribute_type __attr;
    value_set_type __values;
  };
  
};


#endif
