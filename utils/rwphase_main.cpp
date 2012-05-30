#include <iostream>

#include "rwphase.hpp"



struct writer_type
{
  writer_type(utils::rwphase& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    for (int i = 0; i != 1024; ++ i) {
      utils::rwphase::scoped_writer_lock lock(mutex);
      
      if (value == 1024) break;
      
      ++ value;
      
      std::cerr << "writer: " << value << std::endl;

      boost::thread::yield();
    }
  }
  
  utils::rwphase& mutex;
  int& value;
};

struct reader_type
{
  reader_type(utils::rwphase& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    size_t accessed = 0;
    while (value != 1024) {
      utils::rwphase::scoped_reader_lock lock(mutex);
      
      boost::thread::yield();

      ++ accessed;
    }
    
    utils::rwphase::scoped_writer_lock lock(mutex);
    std::cerr << "accessed: " << accessed << std::endl;
  }
  
  utils::rwphase& mutex;
  int& value;
};

int main(int argc, char** argv)
{
  utils::rwphase test;
  

  test.lock_reader();
  test.lock_reader();
  test.unlock_reader();
  test.unlock_reader();
  
  test.lock_writer();
  test.unlock_writer();

  boost::thread_group workers;
  
  utils::rwphase mutex;
  int value = 0;
  
  for (int i = 0; i != 4; ++ i)
    workers.add_thread(new boost::thread(reader_type(mutex, value)));
  
  workers.add_thread(new boost::thread(writer_type(mutex, value)));
  workers.add_thread(new boost::thread(writer_type(mutex, value)));

  workers.join_all();
}
