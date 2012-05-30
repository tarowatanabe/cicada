// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// an implementation from: http://www.cs.unc.edu/~anderson/papers/rtsj10-for-web.pdf
// B. Brandenburg and J. Anderson,
// " Spin-Based Reader-Writer Synchronization for Multiprocessor Real-Time Systems",
// Real-Time Systems, special issue on selected papers from the 21st Euromicro Conference on Real-Time Systems,
//  Volume 46, Issue 1, pp. 25-87, 2010.


#ifndef __UTILS__RWPHASE__HPP__
#define __UTILS__RWPHASE__HPP__ 1

#include <stdexcept>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace utils
{
  class rwphase : private boost::noncopyable
  {
  private:
    typedef union 
    {
      uint64_t u;
      struct {
	uint16_t rin;
	uint16_t rout;
	uint16_t win;
	uint16_t wout;
      } s;
    } ticket_type;

  public:
    
    struct scoped_writer_lock
    {
      scoped_writer_lock(rwphase& x) : lock(x) { lock.lock_writer(); }
      ~scoped_writer_lock() { lock.unlock_writer(); }
      
    private:
      rwphase& lock;
    };
    
    struct scoped_reader_lock
    {
      scoped_reader_lock(rwphase& x) : lock(x) { lock.lock_reader(); }
      ~scoped_reader_lock() { lock.unlock_reader(); }
      
    private:
      rwphase& lock;
    };
    
    rwphase() { ticket_.u = 0; }
    
  public:
    void lock_writer()
    {
      uint16_t ticket = __sync_fetch_and_add(&ticket_.s.win, uint16_t(1));
      
      while (ticket != ticket_.s.wout)
	boost::thread::yield();
      
      const uint16_t w = 0x02 | (ticket & 0x01);
      
      ticket = __sync_fetch_and_add(&ticket_.s.rin, w);
      
      while (ticket != ticket_.s.rout)
	boost::thread::yield();
    }
    
    void unlock_writer()
    {
      ticket_.s.rin &= uint16_t(0xfffc);
      
      __sync_fetch_and_add(&ticket_.s.wout, uint16_t(1));
    }
    
    void lock_reader()
    {
      const uint16_t w = __sync_fetch_and_add(&ticket_.s.rin, uint16_t(0x04)) & 0x03;
      
      while ((w != 0) && (w == (ticket_.s.rin & 0x03)))
	boost::thread::yield();
    }
    
    void unlock_reader()
    {
      __sync_fetch_and_add(&ticket_.s.rout, uint16_t(0x04));
    }
    
  private:
    ticket_type ticket_;
  };
};


#endif
