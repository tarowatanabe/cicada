#include <iostream>

#include "rwspinlock.hpp"



struct writer_type
{
  writer_type(utils::rwspinlock& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    for (int i = 0; i != 1024; ++ i) {
      utils::rwspinlock::scoped_writer_lock lock(mutex);
      
      if (value == 1024) break;
      
      ++ value;
      
      std::cerr << "writer: " << value << std::endl;

      boost::thread::yield();
    }
  }
  
  utils::rwspinlock& mutex;
  int& value;
};

struct reader_type
{
  reader_type(utils::rwspinlock& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    size_t accessed = 0;
    while (value != 1024) {
      utils::rwspinlock::scoped_reader_lock lock(mutex);
      
      boost::thread::yield();

      ++ accessed;
    }
    
    utils::rwspinlock::scoped_writer_lock lock(mutex);
    std::cerr << "accessed: " << accessed << std::endl;
  }
  
  utils::rwspinlock& mutex;
  int& value;
};

int main(int argc, char** argv)
{
  utils::rwspinlock test;
  
  if (! test.trylock_reader())
    std::cerr << "try read failed?" << std::endl;

  if (test.trylock_writer())
    std::cerr << "try write succeed?" << std::endl;

  if (! test.trylock_reader())
    std::cerr << "try read failed?" << std::endl;

  test.unlock_reader();
  test.unlock_reader();
  
  test.lock_writer();
  
  if (test.trylock_reader())
    std::cerr << "try read succeed?" << std::endl;

  if (test.trylock_writer())
    std::cerr << "try writer succeed?" << std::endl;
  
  test.unlock_writer();

  boost::thread_group workers;
  
  utils::rwspinlock mutex;
  int value = 0;
  
  for (int i = 0; i != 4; ++ i)
    workers.add_thread(new boost::thread(reader_type(mutex, value)));
  
  workers.add_thread(new boost::thread(writer_type(mutex, value)));
  workers.add_thread(new boost::thread(writer_type(mutex, value)));

  workers.join_all();
}
