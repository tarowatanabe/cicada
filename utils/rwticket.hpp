// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// an implementation from: http://locklessinc.com/articles/locks/

#ifndef __UTILS__RWTICKET__HPP__
#define __UTILS__RWTICKET__HPP__ 1

#include <stdexcept>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace utils
{
  class rwticket : private boost::noncopyable
  {
  private:
    struct mutex_type
    {
      typedef union {
	uint32_t u;
	struct
	{
	  uint16_t ticket;
	  uint16_t users;
	} s;
      } ticket_type;
      
      mutex_type() { ticket_.u = 0; }
      
      void lock()
      {
	const uint16_t me = __sync_fetch_and_add(&ticket_.s.users, uint16_t(1));
	
	while (ticket_.s.ticket != me)
	  boost::thread::yield();
      }

      void unlock()
      {
	utils::atomicop::memory_barrier();
      
	++ ticket_.s.ticket;
      }

      bool try_lock()
      {
	uint16_t me = ticket_.s.users;
	uint16_t menew = me + 1;
	uint32_t cmp    = (uint32_t(me) << 16) + me;
	uint32_t cmpnew = (uint32_t(menew) << 16) + me;
	
	return utils::atomicop::compare_and_swap(ticket_.u, cmp, cmpnew);
      }
      
      bool locked()
      {
	ticket_type u = ticket_;
	
	utils::atomicop::memory_barrier();
	
	return u.s.ticket != u.s.users;
      }
      
      ticket_type ticket_;
    };
    

  public:
    
    struct scoped_writer_lock
    {
      scoped_writer_lock(rwticket& x) : lock(x) { lock.lock_writer(); }
      ~scoped_writer_lock() { lock.unlock_writer(); }
      
    private:
      rwticket& lock;
    };
    
    struct scoped_reader_lock
    {
      scoped_reader_lock(rwticket& x) : lock(x) { lock.lock_reader(); }
      ~scoped_reader_lock() { lock.unlock_reader(); }
      
    private:
      rwticket& lock;
    };
    
    rwticket()
      : mutex_(), pending_(0) { }
    
  public:
    void lock_writer()
    {
      mutex_.lock();
      
      while (pending_)
	boost::thread::yield();
    }

    bool trylock_writer()
    {
      if (pending_) return false;
      
      if (! mutex_.try_lock()) return false;
      
      if (pending_) {
	mutex_.unlock();
	
	return false;
      }
      
      return true;
    }
    
    void unlock_writer()
    {
      mutex_.unlock();
    }
    
    void lock_reader()
    {
      for (;;) {
	utils::atomicop::fetch_and_add(pending_, 1);
	
	if (! mutex_.locked()) return;
	
	// release
	utils::atomicop::fetch_and_add(pending_, -1);
	
	while (mutex_.locked())
	  boost::thread::yield();
      }
    }

    bool trylock_reader()
    {
      utils::atomicop::fetch_and_add(pending_, 1);
      
      if (! mutex_.locked()) return true;
      
      utils::atomicop::fetch_and_add(pending_, -1);
      
      return false;
    }
    
    void unlock_reader()
    {
      utils::atomicop::fetch_and_add(pending_, -1);
    }
    
  private:
    mutex_type mutex_;
    int        pending_;
  };
};

#if 0
namespace utils
{
  class rwticket : private boost::noncopyable
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
      scoped_writer_lock(rwticket& x) : lock(x) { lock.lock_writer(); }
      ~scoped_writer_lock() { lock.unlock_writer(); }
      
    private:
      rwticket& lock;
    };
    
    struct scoped_reader_lock
    {
      scoped_reader_lock(rwticket& x) : lock(x) { lock.lock_reader(); }
      ~scoped_reader_lock() { lock.unlock_reader(); }
      
    private:
      rwticket& lock;
    };
    
    rwticket() 
    {
      ticket_.u = 0;
    }
    
  public:
    void lock_writer()
    {
      const uint64_t me = __sync_fetch_and_add(&ticket_.u, uint64_t(1) << 32);
      const uint16_t val = me >> 32;

      while (val != ticket_.s.write)
	boost::thread::yield();
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

      while (val != ticket_.s.read)
	boost::thread::yield();
      
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

#endif
