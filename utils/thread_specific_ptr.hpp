// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__THREAD_SPECIFIC_PTR__HPP__
#define __UTILS__THREAD_SPECIFIC_PTR__HPP__ 1

#include <errno.h>
#include <pthread.h>

#include <boost/thread.hpp>

#include <utils/atomicop.hpp>

namespace utils
{
  
  template <typename Tp>
  class thread_specific_ptr
  {
  public:
    thread_specific_ptr() : key(), initialized(false) { initialize(); }
  
  public:
    Tp& operator*() const { return *get(); }
    Tp* operator->() const { return get(); }
  
    Tp* get() const 
    {
      return static_cast<Tp*>(pthread_getspecific(key));
    }
  
    void reset(Tp* new_value=0)
    {
      std::auto_ptr<Tp> value(get());
      pthread_setspecific(key, new_value);
    }

  private:
    void initialize()
    {
      volatile bool tmp = initialized;
      utils::atomicop::memory_barrier();
      if (! initialized) {
	boost::mutex::scoped_lock lock(mutex);
      
	if (! initialized) {
	  pthread_key_create(&key, delete_value);
	
	  tmp = true;
	  utils::atomicop::memory_barrier();
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
  };
  
};

#endif
