// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// an implementation from: http://locklessinc.com/articles/locks/

#ifndef __UTILS__RWTICKET2__HPP__
#define __UTILS__RWTICKET2__HPP__ 1

#include <stdexcept>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

#if 0
namespace utils
{
  class rwticket2 : private boost::noncopyable
  {
  public:
    
    struct scoped_writer_lock
    {
      scoped_writer_lock(rwticket2& x) : lock(x) { lock.lock_writer(); }
      ~scoped_writer_lock() { lock.unlock_writer(); }
      
    private:
      rwticket2& lock;
    };
    
    struct scoped_reader_lock
    {
      scoped_reader_lock(rwticket2& x) : lock(x) { lock.lock_reader(); }
      ~scoped_reader_lock() { lock.unlock_reader(); }
      
    private:
      rwticket2& lock;
    };
    
    rwticket2()
      : lock_(0), pending_(0) {}
    
  public:
    void lock_writer()
    {
      // write lock
      while (! utils::atomicop::compare_and_swap(lock_, 0, 1))
	boost::thread::yield();
      
      while (pending_)
	boost::thread::yield();
    }

    bool trylock_writer()
    {
      if (pending_) return false;
      
      if (! utils::atomicop::compare_and_swap(lock_, 0, 1)) return false;
      
      if (pending_) {
	// unlock!
	__sync_lock_test_and_set(&lock_, 0);
	return false;
      }
      
      return true;
    }
    
    void unlock_writer()
    {
      // test and set...
      __sync_lock_test_and_set(&lock_, 0);
    }
    
    void lock_reader()
    {
      for (;;) {
	utils::atomicop::fetch_and_add(pending_, 1);
	
	if (! lock_) return;
	
	// release
	utils::atomicop::fetch_and_add(pending_, -1);
	
	while (lock_)
	  boost::thread::yield();
      }
    }

    bool trylock_reader()
    {
      utils::atomicop::fetch_and_add(pending_, 1);
      
      if (! lock_) return true;
      
      utils::atomicop::fetch_and_add(pending_, -1);
      
      return false;
    }
    
    void unlock_reader()
    {
      utils::atomicop::fetch_and_add(pending_, -1);
    }
    
  private:
    volatile int lock_;
    volatile int pending_;
  };
};
#endif

namespace utils
{
  class rwticket2 : private boost::noncopyable
  {
  private:
    typedef union 
    {
      uint64_t u;
      uint32_t us;
      struct {
	uint16_t write;
	uint16_t read;
	uint16_t users;
      } s;
    } ticket_type;
    
  public:
    
    struct scoped_writer_lock
    {
      scoped_writer_lock(rwticket2& x) : lock(x) { lock.lock_writer(); }
      ~scoped_writer_lock() { lock.unlock_writer(); }
      
    private:
      rwticket2& lock;
    };
    
    struct scoped_reader_lock
    {
      scoped_reader_lock(rwticket2& x) : lock(x) { lock.lock_reader(); }
      ~scoped_reader_lock() { lock.unlock_reader(); }
      
    private:
      rwticket2& lock;
    };
    
    rwticket2() 
    {
      ticket_.u = 0;
    }
    
  public:
    void lock_writer()
    {
      const uint64_t me = __sync_fetch_and_add(&ticket_.u, uint64_t(1) << 32);
      const uint16_t val = me >> 32;
      
      for (int loop = 0; val != ticket_.s.write; ++ loop) {
	if ((loop & 0x63) == 0x63) {
	  struct timespec tm;
	  tm.tv_sec = 0;
	  tm.tv_nsec = 2000001;
	  nanosleep(&tm, NULL);
	} else
	  boost::thread::yield();
      }
    }

    bool trylock_writer()
    {
      const uint64_t me = ticket_.s.users;
      const uint16_t menew = me + 1;
      const uint64_t read = uint64_t(ticket_.s.read) << 16;
      const uint64_t cmp    = (me << 16) + read + me;
      const uint64_t cmpnew = (uint64_t(menew) << 16) + read + me;
      
      return utils::atomicop::compare_and_swap(ticket_.u, cmp, cmpnew);
    }
    
    void unlock_writer()
    {
      ticket_type t = ticket_;
      
      utils::atomicop::memory_barrier();
      
      ++ t.s.write;
      ++ t.s.read;
      
      ticket_.us = t.us;
    }
    
    void lock_reader()
    {
      const uint64_t me = __sync_fetch_and_add(&ticket_.u, uint64_t(1) << 32);
      const uint16_t val = me >> 32;

      for (int loop = 0; val != ticket_.s.read; ++ loop) {
	if ((loop & 0x63) == 0x63) {
	  struct timespec tm;
	  tm.tv_sec = 0;
	  tm.tv_nsec = 2000001;
	  nanosleep(&tm, NULL);
	} else
	  boost::thread::yield();
      }
      
      __sync_add_and_fetch(&ticket_.s.read, uint16_t(1));
    }

    bool trylock_reader()
    {
      const uint64_t me    = ticket_.s.users;
      const uint16_t menew = me + 1;
      const uint64_t write = ticket_.s.write;
      const uint64_t cmp    = (me << 32) + (me << 16) + write;
      const uint64_t cmpnew = (uint64_t(menew) << 32) + (menew << 16) + write;
      
      return utils::atomicop::compare_and_swap(ticket_.u, cmp, cmpnew);
    }

    void unlock_reader()
    {
      __sync_add_and_fetch(&ticket_.s.write, uint16_t(1));
    }
    
  private:
    ticket_type ticket_;
  };
};

#endif
