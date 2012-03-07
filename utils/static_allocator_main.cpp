//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//


// TLS support is disabled for MacOSX... (impossible???)
#include <vector>
#include <iostream>

#include <utils/static_allocator.hpp>

#include <boost/thread.hpp>

template <typename Tp, size_t Size>
struct worker
{
  
  void operator()() 
  {
    typedef utils::static_allocator<Tp, Size> alloc_type;
    typedef std::vector<Tp*> vec_type;
    
    alloc_type alloc;
    vec_type vec;
    for (int i = 0; i < 1024 * 16; ++ i)
      vec.push_back(alloc.allocate(Size));
    
    for (typename vec_type::iterator iter = vec.begin(); iter != vec.end(); ++ iter)
      alloc.deallocate(*iter, Size);				      
  }
};

int main(int argc, char** argv)
{
  std::auto_ptr<boost::thread> thread1(new boost::thread(worker<int, 512>()));
  std::auto_ptr<boost::thread> thread2(new boost::thread(worker<int, 512>()));
  
  thread1->join();
  thread2->join();
}
