//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>

#include "resource.hpp"

#include "vertical_coded_vector.hpp"
#include "vertical_coded_device.hpp"

#include <boost/thread.hpp>

struct Task
{
  const utils::vertical_coded_vector<int>& coded;
  utils::vertical_coded_vector_mapped<int>& coded_mapped;
  utils::vertical_coded_vector_mapped<int>& coded_mapped_stream;
  
  Task(const utils::vertical_coded_vector<int>& __coded,
       utils::vertical_coded_vector_mapped<int>& __coded_mapped,
       utils::vertical_coded_vector_mapped<int>& __coded_mapped_stream)
    : coded(__coded), 
      coded_mapped(__coded_mapped),
      coded_mapped_stream(__coded_mapped_stream) {}

  void operator()()
  {
    utils::resource start;

    for (int j = 0; j < coded.size(); ++ j) {
      const int __coded = coded[j];
      const int __mapped = coded_mapped[j];
      const int __stream = coded_mapped_stream[j];
      
      if (__coded != __mapped)
	std::cerr << "DIFFER... coded: " << __coded << " mapped: " << __mapped << std::endl;
      if (__mapped != __stream)
	std::cerr << "DIFFER... mapped: " << __mapped << " stream: " << __stream << std::endl;
    }

    utils::resource end;

    std::cerr << "user time: " << (end.user_time() - start.user_time())
	      << " cpu time: " << (end.cpu_time() - start.cpu_time())
	      << std::endl;
  }
  
};

int main(int argc, char** argv)
{
  srandom(time(0) * getpid());
  
  for (int i = 0; i < 4; ++ i) {

    std::cerr << "iteration: " << i << std::endl;
    
    std::vector<int> integers;
    for (int j = 0; j < 1024 * 64; ++ j) {
      integers.push_back(random() % (1024 * 128));
      int val = random() % (1024 * 128);
      integers.push_back(val);
      integers.push_back(val);
    }
    std::sort(integers.begin(), integers.end());

    {
      utils::vertical_coded_vector<int> coded_empty;
      coded_empty.write("tmptmptmp.vc");
      utils::vertical_coded_vector_mapped<int> __coded_mapped("tmptmptmp.vc");
    }
    
    utils::vertical_coded_vector<int> coded;
    
    
    boost::iostreams::filtering_ostream os;
    os.push(utils::vertical_coded_sink<int>("tmptmptmp.vc.iostream"));
    
    for (int j = 0; j < integers.size(); ++ j) {
      coded.push_back(integers[j]);
      os.write((char*) &integers[j], sizeof(int));
    }
    coded.build();
    coded.write("tmptmptmp.vc");
    os.pop();
    
    const utils::vertical_coded_vector<int>& __coded = coded;
    utils::vertical_coded_vector_mapped<int> __coded_mapped("tmptmptmp.vc");
    utils::vertical_coded_vector_mapped<int> __coded_mapped_stream("tmptmptmp.vc.iostream");
    for (int j = 0; j < integers.size(); ++ j) {
      
      
      if (__coded[j] != integers[j])
	std::cerr << "DIFFER... coded: " <<  __coded[j] << " orig: " << integers[j] << std::endl;
      if (__coded_mapped[j] != integers[j])
	std::cerr << "DIFFER... mapped: " << __coded_mapped[j] << " orig: " << integers[j] << std::endl;
      if (__coded_mapped_stream[j] != integers[j])
	std::cerr << "DIFFER... stream: " << __coded_mapped_stream[j] << " orig: " << integers[j] << std::endl;
    }
    
    std::cout << "START THREAD" << std::endl;

    utils::resource start;

    std::vector<boost::thread*> threads(6);
  
    for (int i = 0; i < threads.size(); ++ i) {
      threads[i] = new boost::thread(Task(__coded, __coded_mapped, __coded_mapped_stream));
    }
    for (int i = 0; i < threads.size(); ++ i) {
      threads[i]->join();
      delete threads[i];
    }

    utils::resource end;

    std::cerr << "user time: " << (end.user_time() - start.user_time())
	      << " cpu time: " << (end.cpu_time() - start.cpu_time())
	      << std::endl;
  }

  
  
}
