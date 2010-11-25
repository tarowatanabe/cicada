// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__BOUNDED_QUEUE__HPP__
#define __UTILS__BOUNDED_QUEUE__HPP__ 1

#include <vector>
#include <algorithm>

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

namespace utils
{
  template <typename Tp, typename Alloc=std::allocator<Tp> >
  class bounded_queue : private boost::noncopyable
  {
  private:
    typedef std::vector<Tp, Alloc>    buffer_type;
    typedef boost::mutex              mutex_type;
    typedef boost::condition          condition_type;
    typedef boost::mutex::scoped_lock lock_type;
    
  private:
    buffer_type buffer;
    int begin;
    int end;
    int buffered;
    
    mutex_type mutex;
    condition_type not_full, not_empty;
    
  public:
    bounded_queue(const int n)
      : buffer(n), begin(), end(), buffered() {}
    
  public:
    size_t size() { lock_type lock(mutex); return buffered; }
    bool empty() { lock_type lock(mutex); return buffered == 0; }
    bool full() { lock_type lock(mutex); return buffered == buffer.size(); }
    
  public:    
    bool push_swap(Tp& x, const bool no_wait=false)
    {
      lock_type lock(mutex);

      if (no_wait) {
	if (buffered == buffer.size())
	  return false;
      } else {
	while (buffered == buffer.size())
	  not_full.wait(lock);
      }
      
      using namespace boost;
      using namespace std;
      
      swap(buffer[end], x);
      end = (end + 1) % buffer.size();
      ++ buffered;
      
      not_empty.notify_one();
      return true;
    }

    bool push(const Tp& x, const bool no_wait=false)
    {
      lock_type lock(mutex);
      
      if (no_wait) {
	if (buffered == buffer.size())
	  return false;
      } else {
	while (buffered == buffer.size())
	  not_full.wait(lock);
      }
      
      buffer[end] = x;
      end = (end + 1) % buffer.size();
      ++ buffered;
      
      not_empty.notify_one();
      return true;
    }
    
    bool pop_swap(Tp& x, const bool no_wait=false)
    {
      lock_type lock(mutex);
      
      if (no_wait) {
	if (buffered == 0)
	  return false;
      } else {
	while (buffered == 0)
	  not_empty.wait(lock);
      }
      
      using namespace boost;
      using namespace std;
      
      swap(x, buffer[begin]);
      buffer[begin] = Tp();
      begin = (begin + 1) % buffer.size();
      -- buffered;
      
      not_full.notify_one();
      return true;
    }
    
    bool pop(Tp& x, const bool no_wait=false)
    {
      lock_type lock(mutex);
      
      if (no_wait) {
	if (buffered == 0)
	  return false;
      } else {
	while (buffered == 0)
	  not_empty.wait(lock);
      }
      
      x = buffer[begin];
      buffer[begin] = Tp();
      begin = (begin + 1) % buffer.size();
      -- buffered;
      
      not_full.notify_one();
      return true;
    }
  };
};

#endif
