// -*- mode: c++ -*-

#ifndef __UTILS__THREAD_SPECIFIC_PTR__HPP__
#define __UTILS__THREAD_SPECIFIC_PTR__HPP__ 1

#include <errno.h>
#include <pthread.h>

#include <boost/thread.hpp>
#include <utils/config.hpp>

namespace utils
{
  
  template <typename Tp>
  class thread_specific_ptr
  {
  public:
#ifdef HAVE_TLS
    thread_specific_ptr() : key(), initialized(false), instance(0) { initialize(); }
#else
    thread_specific_ptr() : key(), initialized(false) { initialize(); }
#endif
  
  public:
    Tp& operator*() const { return *get(); }
    Tp* operator->() const { return get(); }
  
    Tp* get() const 
    {
#ifdef HAVE_TLS
      return instance;
#else
      return static_cast<Tp*>(pthread_getspecific(key));
#endif
    }
  
    void reset(Tp* new_value=0)
    {
#ifdef HAVE_TLS
      std::auto_ptr<Tp> value(get());
      pthread_setspecific(key, new_value);
      instance = new_value;
#else
      std::auto_ptr<Tp> value(get());
      pthread_setspecific(key, new_value);
#endif
    }

  private:
    void initialize()
    {
      std::cerr << "initializing" << std::endl;

      volatile bool tmp = initialized;
      __sync_synchronize();
      if (! initialized) {
	boost::mutex::scoped_lock lock(mutex);
      
	if (! initialized) {
	  std::cerr << "creating" << std::endl;
	  pthread_key_create(&key, delete_value);
	
	  tmp = true;
	  __sync_synchronize();
	  initialized = tmp;
	}
      }
    }

  private:
    thread_specific_ptr(const thread_specific_ptr& x) {}
    thread_specific_ptr& operator=(const thread_specific_ptr& x) {}
  
  private:
    static void delete_value(void* data)
    {
      if (data)
	delete static_cast<Tp*>(data);
    }
  
  private:
    boost::mutex mutex;
    pthread_key_t key;
    bool initialized;
  
#ifdef HAVE_TLS
    __thread Tp* instance;
#endif
  };
  
};

#endif
