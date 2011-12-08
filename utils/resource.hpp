// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__RESOURCE__HPP__
#define __UTILS__RESOURCE__HPP__ 1

#include <cstdio>

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <utils/config.hpp>

#ifdef HAVE_CLOCK_GETTIME
#include <pthread.h>
#endif

#ifdef HAVE_THREAD_INFO
#include <mach/mach_init.h>
#include <mach/thread_info.h>
#include <mach/thread_act.h>
#endif

namespace utils
{
  struct resource
  {
  public:
    resource()
    {
#if defined RUSAGE_THREAD
      struct rusage  ruse;
      struct rusage  ruse_thread;
      struct timeval utime;
    
      gettimeofday(&utime, NULL);
      getrusage(RUSAGE_SELF, &ruse);
      getrusage(RUSAGE_THREAD, &ruse_thread);
      
      __cpu_time = (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
		    + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec));
      __user_time = double(utime.tv_sec) + 1e-6 * utime.tv_usec;
      __thread_time = (double(ruse_thread.ru_utime.tv_sec + ruse_thread.ru_stime.tv_sec)
		       + 1e-6 * (ruse_thread.ru_utime.tv_usec + ruse_thread.ru_stime.tv_usec));
#elif defined HAVE_CLOCK_GETTIME
      struct rusage   ruse;
      struct timeval  utime;
      struct timespec tspec;
      
      gettimeofday(&utime, NULL);
      getrusage(RUSAGE_SELF, &ruse);
#if defined CLOCK_THREAD_CPUTIME_ID
      ::clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tspec);
#else
      // get the current thread
      pthread_t pth=pthread_self();
      // get the clock_id associated to the current thread
      clockid_t clock_id;
      pthread_getcpuclockid(pth, &clock_id);
      // get the timespec associated to the thread clock
      ::clock_gettime(clock_id, &tspec);
#endif
      
      __cpu_time = (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
		    + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec));
      __user_time = double(utime.tv_sec) + 1e-6 * utime.tv_usec;
      __thread_time = double(tspec.tv_sec) + 1e-9 * tspec.tv_nsec;
#elif defined HAVE_THREAD_INFO
      struct timeval utime;
      struct rusage  ruse;
      struct thread_basic_info th_info;
      mach_msg_type_number_t th_info_count = THREAD_BASIC_INFO_COUNT;
      
      gettimeofday(&utime, NULL);
      getrusage(RUSAGE_SELF, &ruse);
      thread_info(mach_thread_self(), THREAD_BASIC_INFO, (thread_info_t)&th_info, &th_info_count);
      
      __cpu_time = (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
		    + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec));
      __user_time = double(utime.tv_sec) + 1e-6 * utime.tv_usec;;
      __thread_time = (double(th_info.user_time.seconds + th_info.system_time.seconds)
		       + 1e-6 * (th_info.user_time.microseconds + th_info.system_time.microseconds));
#else
      struct rusage  ruse;
      struct timeval utime;
      
      gettimeofday(&utime, NULL);
      getrusage(RUSAGE_SELF, &ruse);
      
      __cpu_time = (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
		    + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec));
      __user_time = double(utime.tv_sec) + 1e-6 * utime.tv_usec;
      __thread_time = __cpu_time;
#endif
    }
    
  public:
    double cpu_time() const { return __cpu_time; }
    double user_time() const { return __user_time; }
    double thread_time() const { return __thread_time; }
    
  public:
    resource& operator+=(const resource& x)
    {
      __cpu_time    += x.__cpu_time;
      __user_time   += x.__user_time;
      __thread_time += x.__thread_time;
      return *this;
    }
    
    resource& operator-=(const resource& x)
    {
      __cpu_time    -= x.__cpu_time;
      __user_time   -= x.__user_time;
      __thread_time -= x.__thread_time;
      return *this;
    }
    
  private:
    double __cpu_time;
    double __user_time;
    double __thread_time;
  };
  
  inline
  resource operator+(const resource& x, const resource& y)
  {
    resource ret = x;
    ret += y;
    return ret;
  }
  
  inline
  resource operator-(const resource& x, const resource& y)
  {
    resource ret = x;
    ret -= y;
    return ret;
  }
};

#endif
