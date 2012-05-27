#include <iostream>

#include "spinlock.hpp"



struct writer_type
{
  writer_type(utils::spinlock& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    for (int i = 0; i != 1024; ++ i) {
      utils::spinlock::scoped_lock lock(mutex);
      
      if (value == 1024) break;
      
      ++ value;
      
      std::cerr << "writer: " << value << std::endl;

      boost::thread::yield();
    }
  }
  
  utils::spinlock& mutex;
  int& value;
};

struct reader_type
{
  reader_type(utils::spinlock& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    size_t accessed = 0;
    while (value != 1024) {
      utils::spinlock::scoped_lock lock(mutex);
      
      boost::thread::yield();

      ++ accessed;
    }
    
    utils::spinlock::scoped_lock lock(mutex);
    std::cerr << "accessed: " << accessed << std::endl;
  }
  
  utils::spinlock& mutex;
  int& value;
};

int main(int argc, char** argv)
{
  utils::spinlock test;
  
  if (! test.try_lock())
    std::cerr << "try read failed?" << std::endl;
  
  if (test.try_lock())
    std::cerr << "try write succeed?" << std::endl;
  
  if (! test.try_lock())
    std::cerr << "try read failed?" << std::endl;
  
  test.unlock();
  
  test.lock();
  
  if (test.try_lock())
    std::cerr << "try read succeed?" << std::endl;

  if (test.try_lock())
    std::cerr << "try writer succeed?" << std::endl;
  
  test.unlock();

  boost::thread_group workers;
  
  utils::spinlock mutex;
  int value = 0;
  
  for (int i = 0; i != 4; ++ i)
    workers.add_thread(new boost::thread(reader_type(mutex, value)));
  
  workers.add_thread(new boost::thread(writer_type(mutex, value)));
  workers.add_thread(new boost::thread(writer_type(mutex, value)));

  workers.join_all();
}
