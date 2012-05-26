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

      for (;;) {
	if (val == ticket_.s.write) break;
	
	boost::thread::yield();
	utils::atomicop::memory_barrier();
      }
    }

    bool trylock_writer()
    {
      const uint64_t me = ticket_.s.users;
      const uint16_t menew = me + 1;
      const uint64_t read = ticket_.s.read << 16;
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

      for (;;) {
	if (val == ticket_.s.read) break;
	
	boost::thread::yield();
	utils::atomicop::memory_barrier();
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
