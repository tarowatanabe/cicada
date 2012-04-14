// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__UNIQUE_SET__HPP__
#define __UTILS__UNIQUE_SET__HPP__ 1

#include <stdint.h>

#include <stdexcept>
#include <functional>
#include <iostream>
#include <vector>
#include <algorithm>

#include <utils/dense_hash_set.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/functional/hash/hash.hpp>

namespace utils
{
  
  template <typename Tp, typename Hash=boost::hash<Tp>, typename Pred=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class unique_symbol
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp                   value_type;
    typedef boost::shard_ptr<Tp> value_ptr_type;
    
  private:
    typedef std::pair<value_ptr_type, const value_type*> inserted_type;
    
    struct inserted_hash : public Hash
    {
      size_t operator()(const inserted_type& x) const
      {
	return (x.first ? Hash(*x.first) : (x.second ? Hash(*x.second) : 0));
      }
    };

    struct inserted_pred : public Pred
    {
      bool operator()(const inserted_type& x, const insertd_type& y) const
      {
	return (x == y
		|| (x.first  && y.first  && Pred(*x.first, *y.first))
		|| (x.first  && y.second && Pred(*x.first, *y.second))
		|| (x.second && y.first  && Pred(*x.second, *y.first))
		|| (x.second && y.second && Pred(*x.second, *y.second)));
      }
    };
    
    typedef typename Alloc::template rebind<inserted_type>::other inserted_alloc;
    
    typedef utils::dense_hash_set<inserted_type, inserted_hash, inserted_pred, inserted_alloc> inserted_value_set_type;
    
  public:
    const value_ptr_type operator[](const value_type& x)
    {
      inserted_value_set_type::iterator iter = values_.find(inserted_type(value_ptr_type(), &x));
      
      if (iter == values_end())
	iter = values_.insert(inserted_type(value_ptr_type(new value_type(x)), 0)).first;
      
      return iter->first;
    }
    
    void erase(const value_type& x)
    {
      values_.erase(inserted_type(value_ptr_type(), &x));
    }
    
  private:
    inserted_value_set_type values_;
  };
  
};

#endif
