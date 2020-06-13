// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// an implementation from: http://locklessinc.com/articles/locks/

#ifndef __UTILS__RWSPINLOCK__HPP__
#define __UTILS__RWSPINLOCK__HPP__ 1

#if 0
#define atomic_xadd(P, V) __sync_fetch_and_add((P), (V))
#define cmpxchg(P, O, N) __sync_val_compare_and_swap((P), (O), (N))
#define atomic_inc(P) __sync_add_and_fetch((P), 1)
#define atomic_dec(P) __sync_add_and_fetch((P), -1) 
#define atomic_add(P, V) __sync_add_and_fetch((P), (V))
#define atomic_set_bit(P, V) __sync_or_and_fetch((P), 1<<(V))
#define atomic_clear_bit(P, V) __sync_and_and_fetch((P), ~(1<<(V)))
#endif

#include <stdexcept>

#include <boost/core/noncopyable.hpp>
#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace utils
{
  class rwspinlock : private boost::noncopyable
  {
  private:
    enum {
      RW_WAIT_BIT  = 0,
      RW_WRITE_BIT = 1,
      RW_READ_BIT  = 2,
    };
    
    enum {
      RW_WAIT  = 1,
      RW_WRITE = 2,
      RW_READ  = 4,
    };
  public:
    
    struct scoped_writer_lock
    {
      scoped_writer_lock(rwspinlock& x) : lock(x) { lock.lock_writer(); }
      ~scoped_writer_lock() { lock.unlock_writer(); }
      
    private:
      rwspinlock& lock;
    };
    
    struct scoped_reader_lock
    {
      scoped_reader_lock(rwspinlock& x) : lock(x) { lock.lock_reader(); }
      ~scoped_reader_lock() { lock.unlock_reader(); }
      
    private:
      rwspinlock& lock;
    };
    
    rwspinlock() : lock_(0) {}
    
  public:
    void lock_writer()
    {
      for (;;) {
	uint32_t state = lock_;
	
	if (state < RW_WRITE) {
	  if (utils::atomicop::compare_and_swap(lock_, state, uint32_t(RW_WRITE))) return;
	  
	  state = lock_;
	}
	
	if (! (state & RW_WAIT))
	  __sync_or_and_fetch(&lock_, uint32_t(1) << RW_WAIT_BIT);
	
	while (lock_ > RW_WAIT)
	  boost::thread::yield();
      }
    }

    bool trylock_writer()
    {
      uint32_t state = lock_;
      
      return (state < RW_WRITE && utils::atomicop::compare_and_swap(&lock_, state, uint32_t(RW_WRITE)));
    }
    
    void unlock_writer()
    {
      utils::atomicop::fetch_and_add(lock_, uint32_t(- RW_WRITE));
    }
    
    void lock_reader()
    {
      for (;;) {
	while (lock_ & (RW_WAIT | RW_WRITE))
	  boost::thread::yield();
	
	if (! (utils::atomicop::fetch_and_add(lock_, uint32_t(RW_READ)) & (RW_WAIT | RW_WRITE))) return;
	
	utils::atomicop::fetch_and_add(lock_, uint32_t(- RW_READ));
      }
    }

    bool trylock_reader()
    {
      uint32_t state = utils::atomicop::fetch_and_add(lock_, uint32_t(RW_READ));
      
      if (!(state & (RW_WAIT | RW_WRITE))) return true;
      
      utils::atomicop::fetch_and_add(lock_, uint32_t(- RW_READ));
      
      return false;
    }

    void unlock_reader()
    {
      utils::atomicop::fetch_and_add(lock_, uint32_t(- RW_READ));
    }
    
  private:
    volatile uint32_t lock_;
  };
};

#endif
