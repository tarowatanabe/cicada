// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// this is not debugged...

#ifndef __UTILS__UNIQUE_SET__HPP__
#define __UTILS__UNIQUE_SET__HPP__ 1

#include <stdint.h>

#include <stdexcept>
#include <iterator>
#include <algorithm>

#include <utils/compact_set.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/functional/hash/hash.hpp>

namespace utils
{
  
  template <typename Tp, typename Hash=boost::hash<Tp>, typename Pred=std::equal_to<Tp>, typename Alloc=std::allocator<Tp> >
  class unique_set
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Tp                    value_type;
    typedef boost::shared_ptr<Tp> value_ptr_type;
    
  private:
    typedef std::pair<value_ptr_type, const value_type*> inserted_type;
    
    struct inserted_hash : public Hash
    {
      size_t operator()(const inserted_type& x) const
      {
	return (x.first ? Hash::operator()(*x.first) : (x.second ? Hash::operator()(*x.second) : size_t(0)));
      }
    };

    struct inserted_pred : public Pred
    {
      bool operator()(const inserted_type& x, const inserted_type& y) const
      {
	return (x == y
		|| (x.first  && y.first  && Pred::operator()(*x.first, *y.first))
		|| (x.first  && y.second && Pred::operator()(*x.first, *y.second))
		|| (x.second && y.first  && Pred::operator()(*x.second, *y.first))
		|| (x.second && y.second && Pred::operator()(*x.second, *y.second)));
      }
    };
    
    typedef typename Alloc::template rebind<inserted_type>::other inserted_alloc;

    struct inserted_unassigned
    {
      inserted_type operator()() const
      {
	return inserted_type(value_ptr_type(), 0);
      }
    };
    
    typedef utils::compact_set<inserted_type,
			       inserted_unassigned, inserted_unassigned, 
			       inserted_hash, inserted_pred, inserted_alloc> inserted_value_set_type;

  public:
    class const_iterator : public inserted_value_set_type::const_iterator
    {
    public:
      typedef std::forward_iterator_tag iterator_category;  // very little defined!
      typedef Tp        value_type;
      typedef const Tp& reference;
      typedef const Tp* pointer;
      typedef typename inserted_value_set_type::difference_type difference_type;
      typedef typename inserted_value_set_type::size_type size_type;

    private:
      typedef typename inserted_value_set_type::const_iterator iterator_type;
      
    public:
      const_iterator(const iterator_type& x)  : iterator_type(x) {}
    public:
      reference operator*() const { return *iterator_type::operator*().first; }
      pointer   operator->() const { return iterator_type::operator*().first.get(); }
    };

    typedef const_iterator iterator;

  public:
    unique_set() : values_() {  }
    
  public:
    const value_ptr_type operator[](const value_type& x)
    {
      typename inserted_value_set_type::iterator iter = values_.find(inserted_type(value_ptr_type(), &x));
      
      if (iter == values_.end())
	iter = values_.insert(inserted_type(value_ptr_type(new value_type(x)), 0)).first;
      
      return iter->first;
    }

    const value_ptr_type operator[](const value_ptr_type& x)
    {
      return values_.insert(inserted_type(x, 0)).first->first;
    }
    
#if 0
    void erase(const value_type& x)
    {
      values_.erase(inserted_type(value_ptr_type(), &x));
    }
    
    void erase(const value_ptr_type& x)
    {
      values_.erase(inserted_type(x, 0));
    }
#endif

    void clear() { values_.clear(); }

    iterator find(const value_type& x) const { return values_.find(inserted_type(value_ptr_type(), &x)); }

    iterator begin() const { return values_.begin(); }
    iterator end() const { return values_.end(); }
    
    size_type size() const { return values_.size(); }
    bool empty() const { return values_.empty(); }
    
  private:
    inserted_value_set_type values_;
  };
  
};

#endif
