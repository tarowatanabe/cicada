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

#if defined HAVE_CLOCK_GETTIME

#include <pthread.h>

namespace utils
{
  class resource
  {
  public:
    resource()
    {
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
    }
    
  public:
    double cpu_time() const { return (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
				      + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec)); }
    double user_time() const { return double(utime.tv_sec) + 1e-6 * utime.tv_usec; }
    double thread_time() const { return double(tspec.tv_sec) + 1e-9 * tspec.tv_nsec; }
    
  private:
    struct rusage   ruse;
    struct timeval  utime;
    struct timespec tspec;
  };
};

#elif defined HAVE_THREAD_INFO

#include <mach/mach_init.h>
#include <mach/thread_info.h>
#include <mach/thread_act.h>

namespace utils
{
  class resource
  {
  public:
    resource()
    {
      gettimeofday(&utime, NULL);
      getrusage(RUSAGE_SELF, &ruse);

      struct thread_basic_info th_info;
      mach_msg_type_number_t th_info_count = THREAD_BASIC_INFO_COUNT;
      
      if (thread_info(mach_thread_self(), THREAD_BASIC_INFO, (thread_info_t)&th_info, &th_info_count) == KERN_SUCCESS)
	__thread_time = (double(th_info.user_time.seconds + th_info.system_time.seconds)
			 + 1e-6 * (th_info.user_time.microseconds + th_info.system_time.microseconds));
      else
	__thread_time = 0.0;
    }
    
  public:
    double cpu_time() const { return (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
				      + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec)); }
    double user_time() const { return double(utime.tv_sec) + 1e-6 * utime.tv_usec; }
    double thread_time() const { return __thread_time; }
    
  private:
    struct rusage  ruse;
    struct timeval utime;
    double         __thread_time;
  };
};
#else
namespace utils
{
  class resource
  {
  public:
    resource()
    {
      gettimeofday(&utime, NULL);
      getrusage(RUSAGE_SELF, &ruse);
    }
    
  public:
    double cpu_time() const { return (double(ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec)
				      + 1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec)); }
    double user_time() const { return double(utime.tv_sec) + 1e-6 * utime.tv_usec; }
    double thread_time() const { return cpu_time(); }
    
  private:
    struct rusage  ruse;
    struct timeval utime;
  };
};
#endif

#endif
