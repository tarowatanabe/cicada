//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <memory>
#include <utils/map_file_allocator.hpp>
#include <vector>
#include <iostream>

#include <boost/thread.hpp>

typedef std::vector<int, utils::map_file_allocator<int> > vector_type;

struct Task
{
  vector_type& array;
  
  Task(vector_type& _array) : array(_array) {}
  
  void operator()() throw()
  {
    for (int i = 0; i < 1024 * 1024 * 128; ++ i)
      array.push_back(i);
    for (int i = 0; i < 1024 * 1024 * 128; ++ i)
      if (array[i] != i)
	std::cerr << "DIFFER?" << std::endl;
  }
};

int main(int argc, char** argv)
{
  vector_type array;
  vector_type array2;
  
  for (int i = 0; i < 4; ++ i) {
    
    std::cout << "iteration: " << i << std::endl;

    std::auto_ptr<boost::thread> thread( new boost::thread(Task(array)));
    
    Task task(array2);
    task();
    thread->join();
    
    for (int i = 0; i < array.size(); ++ i)
      if (array2[i] != array[i])
	std::cerr << "differ! i = " << array[i] << " " << array2[i] << std::endl;
    
    array.clear();
    array2.clear();
    
    vector_type(array).swap(array);
    vector_type(array2).swap(array2);
  }
  
  std::cerr << "finished!" << std::endl;
}
