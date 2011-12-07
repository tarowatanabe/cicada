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
    
  private:
    struct rusage  ruse;
    struct timeval utime;
  };
  
};

#endif
