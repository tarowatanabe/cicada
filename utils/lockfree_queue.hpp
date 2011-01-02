// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// lock-free implementation of FIFO queue

#ifndef __UTILS__LOCKFREE_QUEUE__HPP__
#define __UTILS__LOCKFREE_QUEUE__HPP__ 1

#include <stdint.h>

#include <vector>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>
#include <utils/spinlock.hpp>

namespace utils
{
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class lockfree_queue
  {
  public:
    typedef Tp     value_type;
    typedef size_t size_type;
    
  private:
    typedef std::vector<value_type, Alloc> buffer_type;
    
    typedef utils::spinlock            spinlock_type;
    typedef spinlock_type::scoped_lock lock_type;
    
    
    struct border_type
    {
      struct border_pair_type
      {
	int32_t front;
	int32_t back;
      };
      
      union {
	border_pair_type border;
	int64_t value;
      } data;
      
      border_type() {}
      
      inline volatile int64_t& value() volatile { return data.value; }
      inline const    int64_t& value() const    { return data.value; }
      inline          int64_t& value()          { return data.value; }
      
      inline volatile int32_t& front() volatile { return data.border.front; }
      inline const    int32_t& front() const    { return data.border.front; }
      inline          int32_t& front()          { return data.border.front; }
      
      inline volatile int32_t& back() volatile { return data.border.back; }
      inline const    int32_t& back() const    { return data.border.back; }
      inline          int32_t& back()          { return data.border.back; }
    };
    
  public:
    lockfree_queue(size_type size)
      : buffer(size),
	border_push(),
	border_pop()
    {
      border_push.value() = 0;
      border_pop.value() = 0;
    }

  private:
    struct assign_equal
    {
      void operator()(value_type& x, const value_type& y) const
      {
	x = y;
      }
    };
    
    struct assign_swap
    {
      void operator()(value_type& x, const value_type& y) const
      {
	using namespace boost;
	using namespace std;
	swap(x, const_cast<value_type&>(y));
      }
    };
    
  public:
    bool push(const value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __push(x, assign_equal());
      else {
	for (;;) {
	  for (int i = 0; i < 50; ++ i) {
	    if (__push(x, assign_equal()))
	      return true;
	    else
	      boost::thread::yield();
	  }
	  __sleep_long();
	}
	
	return true;
      }
    }
    
    bool push_swap(value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __push(x, assign_swap());
      else {
	for (;;) {
	  for (int i = 0; i < 50; ++ i) {
	    if (__push(x, assign_swap()))
	      return true;
	    else
	      boost::thread::yield();
	  }
	  __sleep_long();
	}
	return true;
      }
    }
    
    bool pop(value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __pop(x, assign_equal());
      else {
	for (;;) {
	  for (int i = 0; i < 50; ++ i) {
	    if (__pop(x, assign_equal()))
	      return true;
	    else
	      boost::thread::yield();
	  }
	  __sleep_long();
	}
	
	return true;
      }
    }
    
    bool pop_swap(value_type& x, const bool no_wait=false)
    {
      if (no_wait)
	return __pop(x, assign_swap());
      else {
	for (;;) {
	  for (int i = 0; i < 50; ++ i) {
	    if (__pop(x, assign_swap()))
	      return true;
	    else
	      boost::thread::yield();
	  }
	  __sleep_long();
	}
	
	return true;
      }
    }
    
  private:
    
    void __sleep_long()
    {
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    }

    template <typename Assigner>
    bool __push(const value_type& x, Assigner __assign)
    {
      border_type border_prev;
      border_type border_post;
      
      border_prev.value() = atomicop::fetch_and_add(border_push.value(), int64_t(0));
      
      while (border_prev.back() - border_prev.front() < buffer.size()) {
	
	border_post.front() = border_prev.front();
	border_post.back()  = border_prev.back() + 1;
	
	{
	  lock_type lock(spinlock_push);
	  if (atomicop::compare_and_swap(border_push.value(), border_prev.value(), border_post.value())) {
	    atomicop::memory_barrier();
	    __assign(buffer[border_prev.back() % buffer.size()], x);
	    
	    atomicop::fetch_and_add(border_pop.back(), 1);
	    return true;
	  }
	}
	
	border_prev.value() = atomicop::fetch_and_add(border_push.value(), int64_t(0));
      }
      return false;
      
    }
    
    template <typename Assigner>
    bool __pop(value_type& x, Assigner __assign)
    {
      border_type border_prev;
      border_type border_post;
      
      border_prev.value() = atomicop::fetch_and_add(border_pop.value(), int64_t(0));
      
      while (border_prev.back() - border_prev.front() > 0) {
	border_post.front() = border_prev.front() + 1;
	border_post.back()  = border_prev.back();
	
	{
	  lock_type lock(spinlock_pop);
	  if (atomicop::compare_and_swap(border_pop.value(), border_prev.value(), border_post.value())) {
	    atomicop::memory_barrier();
	    __assign(x, buffer[border_prev.front() % buffer.size()]);
	    
	    atomicop::fetch_and_add(border_push.front(), 1);
	    return true;
	  }
	}
	
	border_prev.value() = atomicop::fetch_and_add(border_pop.value(), int64_t(0));
      }
      return false;
    }
    

    
  private:
    buffer_type buffer;
    volatile border_type border_push;
    volatile border_type border_pop;

    spinlock_type spinlock_push;
    spinlock_type spinlock_pop;
  };
};

#endif
