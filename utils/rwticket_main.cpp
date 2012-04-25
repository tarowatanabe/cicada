#include <iostream>

#include "rwticket.hpp"



struct writer_type
{
  writer_type(utils::rwticket& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    for (int i = 0; i != 1024; ++ i) {
      utils::rwticket::scoped_writer_lock lock(mutex);
      
      ++ value;
      
      std::cerr << "writer: " << value << std::endl;
    }
  }
  
  utils::rwticket& mutex;
  int& value;
};

struct reader_type
{
  reader_type(utils::rwticket& __mutex,
	      int& __value) : mutex(__mutex), value(__value) {}
  
  void operator()()
  {
    while (value != 1024) {
      utils::rwticket::scoped_reader_lock lock(mutex);
      
      std::cerr << "value = " << value << std::endl;

      boost::thread::yield();
    }
  }
  
  utils::rwticket& mutex;
  int& value;
};

int main(int argc, char** argv)
{
  boost::thread_group workers;
  
  utils::rwticket mutex;
  int value = 0;
  
  for (int i = 0; i != 3; ++ i)
    workers.add_thread(new boost::thread(reader_type(mutex, value)));
  
  workers.add_thread(new boost::thread(writer_type(mutex, value)));

  workers.join_all();
}
