

#include <pthread.h>

#include <iostream>
#include <utility>
#include <memory>

#include <boost/thread.hpp>

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

template <typename Tp>
class singleton
{
public:
  singleton() : instance(0) {}
  
  Tp* get() const
  {
    // double lock
    Tp* tmp = instance;
    //utils::atomicop::memory_barrier();
    __sync_synchronize();
    if (! instance) {
      boost::mutex::scoped_lock lock(const_cast<boost::mutex&>(mutex));
      if (! instance) {
	tmp = new Tp();
	//utils::atomicop::memory_barrier();
	__sync_synchronize();
	const_cast<Tp*>(instance) = tmp;
      }
    }
    return instance;
  }
  
  Tp& operator*() const { return *get(); }
  Tp* operator->() const { return get(); }

  
private:
  singleton(const singleton& ){}
  singleton& operator=(const singleton& ) {}
  
private:
  Tp*          instance;
  boost::mutex mutex;
};

struct Task
{
  
  Task(int __i) : i(__i) {}
  
  void operator()()
  {
    for (int j = 0; j < 1024 * 1024; ++ j) {
      static thread_specific_ptr<int> value;

      if (! value.get())
	value.reset(new int(i));
      
      if (*value != i)
	std::cout << i << " " << *value << std::endl;
    }
  }

  
  
  int i;
};


int main(int argc, char** argv)
{
  
  boost::thread_group workers;
  for (int i = 0; i < 24; ++ i)
    workers.add_thread(new boost::thread(Task(i)));
  
  workers.join_all();

}
