// -*- mode: c++ -*-

// lock-free implementation of LIFO stack


#ifndef __UTILS__LOCKFREE_STACK__HPP__
#define __UTILS__LOCKFREE_STACK__HPP__ 1

#include <vector>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace utils
{
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class lockfree_stack
  {
  private:
    typedef std::vector<Tp, Alloc> buffer_type;
    
  public:
    typedef Tp     value_type;
    typedef size_t size_type;
    
  private:
    typedef int32_t border_type;
        
  public:
    lockfree_stack(size_type size)
      : buffer(size),
	border_push(0),
	border_pop(0),
    {}
    
  public:
    bool push(const value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __push(x);
      else {
	while (! __push(x))
	  boost::thread::yield();
	return true;
      }
    }
    
    bool push_swap(value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __push_swap(x);
      else {
	while (! __push_swap(x))
	  boost::thread::yield();
	return true;
      }
    }
    
    bool pop(value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __pop(x);
      else {
	while (! __pop(x))
	  boost::thread::yield();
	return true;
      }
    }
    
    bool pop_swap(value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __pop_swap(x);
      else {
	while (! __pop_swap(x))
	  boost::thread::yield();
	return true;
      }
    }
    
  private:
    bool __push(const value_type& x)
    {
      border_type prev_size;
      
      atomicop::memory_barrier();
      prev_size = border_push;
      
      while (prev_size < buffer.size()) {
	
	if (atomicop::compare_and_swap(border_push, prev_size, prev_size + 1)) {
	  buffer[prev_size] = x;
	  atomicop::fetch_and_add(border_pop, 1);
	  return true;
	}
	
	atomicop::memory_barrier();
	prev_size = border_push;
      }
      
      return false;
    }
    
    bool __push_swap(value_type& x)
    {
      border_type prev_size;
      
      atomicop::memory_barrier();
      prev_size = border_push;
      
      while (prev_size < buffer.size()) {
	
	if (atomicop::compare_and_swap(border_push, prev_size, prev_size + 1)) {
	  using namespace boost;
	  using namespace std;
	  swap(buffer[prev_size], x);
	  atomicop::fetch_and_add(border_pop, 1);
	  return true;
	}
	
	atomicop::memory_barrier();
	prev_size = border_push;
      }
      
      return false;
    }
    
    bool __pop(value_type& x)
    {
      border_type prev_size;
      
      atomicop::memory_barrier();
      prev_size = border_pop;
      
      while (prev_size > 0) {
	
	if (atomicop::compare_and_swap(border_pop, prev_size, prev_size - 1)) {
	  x = buffer[prev_size - 1];
	  atomicop::fetch_and_add(border_push, -1);
	  return true;
	}
	
	atomicop::memory_barrier();
	prev_size = border_pop;
      }
      return false;
    }
    
    bool __pop_swap(value_type& x)
    {
      border_type prev_size;
      
      atomicop::memory_barrier();
      prev_size = border_pop;
      
      while (prev_size > 0) {
	
	if (atomicop::compare_and_swap(border_pop, prev_size, prev_size - 1)) {
	  using namespace boost;
	  using namespace std;
	  swap(x, buffer[prev_size - 1]);
	  atomicop::fetch_and_add(border_push, -1);
	  return true;
	}
	
	atomicop::memory_barrier();
	prev_size = border_pop;
      }
      return false;
    }
    
  private:
    buffer_type buffer;
    volatile border_type border_push;
    volatile border_type border_pop;
  };
};

#endif
