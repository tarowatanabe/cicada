// -*- mode: c++ -*-

#ifndef __UTILS__ARC_LIST__HPP__
#define __UTILS__ARC_LIST__HPP__ 1

#include <utility>
#include <list>
#include <algorithm>
#include <memory>

#include <utils/static_allocator.hpp>
#include <utils/atomicop.hpp>
#include <utils/memory.hpp>

//#include <boost/intrusive/list.hpp>

namespace utils
{
  // key-data cache...
  // find will return iteraot, bool pair.
  // the second argument, bool indicate whether the iterator is cached or not...
  // since we do not perform any "erase" all the iterator will be valid, as long as list-object exists...
  
  template <typename _Key,
	    typename _Data,
	    size_t _CacheSize,
	    typename _KeyEqual=std::equal_to<_Key>,
	    typename _Alloc=std::allocator<std::pair<_Key, _Data> > >
  class arc_list : public _KeyEqual
  {
  public:
    typedef size_t                         size_type;
    typedef ptrdiff_t                      difference_type;
    typedef _Key                           key_type;
    typedef _Data                          data_type;
    
    typedef _KeyEqual                      key_equal;
    
    typedef std::pair<key_type, data_type> value_type;
    
  private:
    
    typedef typename _Alloc::template rebind<value_type>::other value_alloc_type;
    typedef std::list<value_type, utils::static_allocator<value_type, 1, value_alloc_type> > list_type;
    
  public:
    typedef typename list_type::iterator iterator;
    
  public:
    static const size_type cache_size = _CacheSize;

  public:

    arc_list()
      : list_t1(), list_t2(), list_b1(), list_b2(),
	list_tmp(),
	target_t1(0),
	size_t1(0), size_t2(0), size_b1(0), size_b2(0) {}
    
  public:

    size_type size() const { return size_t1 + size_t2; }
    size_type capacity() const { return size_t1 + size_t2 + size_b1 + size_b2; }
    
    void clear()
    {
      list_t1.clear();
      list_t2.clear();
      
      list_b1.clear();
      list_b2.clear();
      
      list_tmp.clear();
      
      target_t1 = 0;
      
      size_t1 = 0;
      size_t2 = 0;
      
      size_b1 = 0;
      size_b2 = 0;
    }
    
    void swap(arc_list& x)
    {
      using namespace std;
      
      list_t1.swap(x.list_t1);
      list_t2.swap(x.list_t2);
      
      list_b1.swap(x.list_b1);
      list_b2.swap(x.list_b2);
      
      list_tmp.swap(x.list_tmp);
      
      swap(target_t1, x.target_t1);
      
      swap(size_t1, x.size_t1);
      swap(size_t2, x.size_t2);
      swap(size_b1, x.size_b1);
      swap(size_b2, x.size_b2);
    }

    std::pair<iterator, bool> find(const key_type& key)
    {
      typename list_type::iterator iter = find_list(list_t1.begin(), list_t1.end(), key);
      if (iter != list_t1.end()) {
	// found in recency list...
	
	-- size_t1;
	++ size_t2;
	list_t2.splice(list_t2.begin(), list_t1, iter);
	return std::make_pair(iterator(list_t2.begin()), true);
      }
      
      iter = find_list(list_t2.begin(), list_t2.end(), key);
      if (iter != list_t2.end()) {
	
	// found in frequency list...
	list_t2.splice(list_t2.begin(), list_t2, iter);
	return std::make_pair(iterator(list_t2.begin()), true);
      }
      
      iter = find_list(list_b1.begin(), list_b1.end(), key);
      if (iter != list_b1.end()) {
	// favor recency...
	
	// increment...
	target_t1 = std::min(difference_type(target_t1) + std::max(difference_type(size_b2 / size_b1), difference_type(1)), difference_type(cache_size));
	
	// remove from b1
	list_tmp.splice(list_tmp.begin(), list_b1, iter);
	-- size_b1;
	
	// adaptive replacement!
	replace();
	
	// add...
	list_t2.splice(list_t2.begin(), list_tmp, list_tmp.begin());
	++ size_t2;
	
	return std::make_pair(iterator(list_t2.begin()), true);
      }

      iter = find_list(list_b2.begin(), list_b2.end(), key);
      if (iter != list_b2.end()) {
	// favor frequency...
	
	// decrement...
	target_t1 = std::max(difference_type(target_t1) - std::max(difference_type(size_b1 / size_b2), difference_type(1)), difference_type(0));
	
	// remove from b2
	list_tmp.splice(list_tmp.begin(), list_b2, iter);
	-- size_b2;
	
	// adaptive replacement!
	replace();
	
	// add...
	list_t2.splice(list_t2.begin(), list_tmp, list_tmp.begin());
	++ size_t2;
	
	return std::make_pair(iterator(list_t2.begin()), true);
      }
      
      // not found in any list!
      if (size_t1 + size_b1 == cache_size) {
	// t1 + b1 is full
	if (size_t1 < cache_size) {
	  // still room in t1...
	  typename list_type::iterator iter = list_b1.end();
	  -- iter;
	  
	  iter->first = key;
	  
	  // remove from b1
	  list_tmp.splice(list_tmp.begin(), list_b1, iter);
	  -- size_b1;
	  
	  // replace!
	  replace();
	  
	  list_t1.splice(list_t1.begin(), list_tmp, list_tmp.begin());
	  ++ size_t1;
	  
	} else {
	  // recycle from the last of list_t1
	  typename list_type::iterator iter = list_t1.end();
	  -- iter;
	  
	  iter->first = key;
	  
	  list_t1.splice(list_t1.begin(), list_t1, iter);
	}
      } else {
	if (size_t1 + size_t2 + size_b1 + size_b2 >= cache_size) {
	  
	  // cache full...
	  if (size_t1 + size_t2 + size_b1 + size_b2 == cache_size * 2) {
	    // directory is full...
	    // extract from B2....
	    
	    typename list_type::iterator iter = list_b2.end();
	    -- iter;
	    
	    iter->first = key;
	    
	    list_tmp.splice(list_tmp.begin(), list_b2, iter);
	    -- size_b2;
	    
	    // replace!
	    replace();
	    
	    list_t1.splice(list_t1.begin(), list_tmp, list_tmp.begin());
	    ++ size_t1;
	    
	  } else {
	    // replace!
	    replace();
	    
	    // simply allocate...
	    list_t1.insert(list_t1.begin(), value_type(key, data_type()));
	    ++ size_t1;
	  }
	  
	} else {
	  // simply allocate...
	  list_t1.insert(list_t1.begin(), value_type(key, data_type()));
	  ++ size_t1;
	}
      }      
      
      return std::make_pair(iterator(list_t1.begin()), false);
    }
    
  private:
    

    void replace()
    {
      // is this correct?
      if (size_t1 >= std::max(size_type(1), target_t1)) {
	// t1 is too big...
	typename list_type::iterator iter = list_t1.end();
	-- iter;
	list_b1.splice(list_b1.begin(), list_t1, iter);
	
	-- size_t1;
	++ size_b1;
      } else {
	
	// t2 is too big...
	typename list_type::iterator iter = list_t2.end();
	-- iter;
	list_b2.splice(list_b2.begin(), list_t2, iter);
	
	-- size_t2;
	++ size_b2;
      }
    }
    
  private:
    
    typename list_type::iterator find_list(typename list_type::iterator first, typename list_type::iterator last, const key_type& key)
    {
      while (first != last && ! equal_to()(first->first, key))
	++ first;
      return first;
    }
    
    key_equal& equal_to() { return static_cast<key_equal&>(*this); }
    
  private:
    list_type list_t1;
    list_type list_t2;
    list_type list_b1;
    list_type list_b2;
    list_type list_tmp;
    
    size_type target_t1;
    
    size_type size_t1;
    size_type size_t2;
    size_type size_b1;
    size_type size_b2;
  };
};

namespace std
{
  
  template <typename _Key, typename _Data, size_t _CacheSize, typename _KeyEqual, typename _Alloc>
  void swap(utils::arc_list<_Key,_Data,_CacheSize,_KeyEqual,_Alloc>& x,
	    utils::arc_list<_Key,_Data,_CacheSize,_KeyEqual,_Alloc>& y) 
  {
    x.swap(y);
  }
};

#endif
