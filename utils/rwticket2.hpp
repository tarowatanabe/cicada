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
#include <utils/spinlock.hpp>

namespace utils
{
  class rwticket2 : private boost::noncopyable
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

#endif
