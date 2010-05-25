// -*- mode: c++ -*-

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
    
    
    spinlock() : m_spinlock()
    {
#if defined(__UTILS_SPINLOCK_GCC_CAS__)
      m_spinlock = 0;
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
      return __sync_bool_compare_and_swap(&m_spinlock, 0, 1);
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
      while (! __sync_bool_compare_and_swap(&m_spinlock, 0, 1))
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
      __sync_lock_test_and_set(&m_spinlock, 0);
      //__sync_synchronize();
      //m_spinlock = 0;
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
    volatile int m_spinlock;
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
