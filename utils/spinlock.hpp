// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SPINLOCK__HPP__
#define __UTILS__SPINLOCK__HPP__ 1

#include <stdexcept>
#include <cassert>
#include <errno.h>
#include <pthread.h>

#include <boost/thread/detail/config.hpp>
#include <boost/utility.hpp>

#include <boost/thread.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/exceptions.hpp>

#include <utils/config.hpp>
#include <utils/atomicop.hpp>

#ifdef HAVE_LIBKERN_OSATOMIC_H
  #include <libkern/OSAtomic.h>
#endif

#if defined(__GNUC__) && ( (__GNUC__ > 4) || ((__GNUC__ >= 4) && (__GNUC_MINOR__ >= 1)) )
#define __UTILS_SPINLOCK_GCC_CAS__ 1
#endif

namespace utils
{
  class spinlock : private boost::noncopyable
  {
  public:
    typedef boost::unique_lock<spinlock>       scoped_lock;
    typedef boost::detail::try_lock_wrapper<spinlock> scoped_try_lock;
    
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
    typedef union {
      uint32_t u;
      struct
      {
	uint16_t ticket;
	uint16_t users;
      } s;
    } ticket_type;
#endif
    
    spinlock()
    {
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
      m_spinlock.u = 0;
#elif defined(HAVE_OSSPINLOCK)
      m_spinlock = OS_SPINLOCK_INIT;
#elif defined(HAVE_PTHREAD_SPINLOCK)
      const int res = pthread_spin_init(&m_spinlock, PTHREAD_PROCESS_SHARED);
      if (res != 0)
	throw boost::thread_resource_error();
#else
      const int res = pthread_mutex_init(&m_spinlock, 0);
      if (res != 0)
	throw boost::thread_resource_error();
#endif
      
    }
    ~spinlock() { 
#if defined(__UTILS_SPINLOCK_GCC_CAS__)

#elif defined(HAVE_OSSPINLOCK)
      // do nothing...
#elif defined(HAVE_PTHREAD_SPINLOCK)
      const int res = pthread_spin_destroy(&m_spinlock);
      assert(res == 0);
#else
      const int res = pthread_mutex_destroy(&m_spinlock);
      assert(res == 0);
#endif
    }
    
  public:
    bool try_lock()
    {
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
      uint16_t me = m_spinlock.s.users;
      uint16_t menew = me + 1;
      uint32_t cmp    = (uint32_t(me) << 16) + me;
      uint32_t cmpnew = (uint32_t(menew) << 16) + me;
      
      return utils::atomicop::compare_and_swap(m_spinlock.u, cmp, cmpnew);
      
#elif defined(HAVE_OSSPINLOCK)
      return OSSpinLockTry(&m_spinlock);
#elif defined(HAVE_PTHREAD_SPINLOCK)
      return ! pthread_spin_trylock(&m_spinlock);
#else
      return ! pthread_mutex_trylock(&m_spinlock);
#endif
    }
    void lock()
    {
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
      const uint16_t me = __sync_fetch_and_add(&m_spinlock.s.users, uint16_t(1));
      
      while (m_spinlock.s.ticket != me)
	boost::thread::yield();
#elif defined(HAVE_OSSPINLOCK)
      OSSpinLockLock(&m_spinlock);
#elif defined(HAVE_PTHREAD_SPINLOCK)
      const int res = pthread_spin_lock(&m_spinlock);
      if (res == EDEADLK)
	throw boost::lock_error();
      assert(res == 0);
#else
      const int res = pthread_mutex_lock(&m_spinlock); 
      if (res == EDEADLK)
	throw boost::lock_error();
      assert(res == 0);
#endif
    }
    
    void unlock()
    {
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
      utils::atomicop::memory_barrier();
      
      ++ m_spinlock.s.ticket;
#elif defined(HAVE_OSSPINLOCK)
      OSSpinLockUnlock(&m_spinlock);
#elif defined(HAVE_PTHREAD_SPINLOCK)
      const int res = pthread_spin_unlock(&m_spinlock);
      if (res == EPERM)
	throw boost::lock_error();
      assert(res == 0);
#else
      const int res = pthread_mutex_unlock(&m_spinlock);
      if (res == EPERM)
	throw boost::lock_error();
      assert(res == 0);
#endif
    }
    
  private:
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
    ticket_type m_spinlock;
#elif defined(HAVE_OSSPINLOCK)
    OSSpinLock m_spinlock;
#elif defined(HAVE_PTHREAD_SPINLOCK)
    pthread_spinlock_t m_spinlock;
#else
    pthread_mutex_t m_spinlock;
#endif
  };
};

#endif
