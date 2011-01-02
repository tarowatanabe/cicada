//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>
#include <set>

#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

#include "utils/lockfree_list_queue.hpp"
#include "utils/hashmurmur.hpp"

typedef boost::shared_ptr<int> ptr_int_type;
typedef std::set<ptr_int_type> integers_type;
typedef utils::lockfree_list_queue<ptr_int_type> lockfree_queue_type;

typedef std::vector<boost::shared_ptr<lockfree_queue_type> > lockfree_queue_set_type;

struct Mapper
{
  integers_type& integers;
  lockfree_queue_type& queue;
  lockfree_queue_set_type& queue_reduce;
  
  Mapper(integers_type& __integers,
	 lockfree_queue_type& __queue,
	 lockfree_queue_set_type& __queue_reduce)
    : integers(__integers), queue(__queue), queue_reduce(__queue_reduce) {}
  
  void operator()()
  {
    utils::hashmurmur<uint32_t> hasher;

    while (1) {
      ptr_int_type value;
      queue.pop_swap(value);
      
      if (! value) break;
      
      integers.insert(value);
      
      queue_reduce[hasher(*value) % queue_reduce.size()]->push_swap(value);
    }
  }
};

struct Reducer
{
  integers_type& integers;
  lockfree_queue_type& queue;
  
  Reducer(integers_type& __integers,
	  lockfree_queue_type& __queue)
    : integers(__integers), queue(__queue) {}
  
  void operator()()
  {
    while (1) {
      ptr_int_type value;
      queue.pop_swap(value);
      
      if (! value) break;
      
      integers.insert(value);
    }
  }
  
};

int main(int argc, char** argv)
{
  typedef utils::lockfree_list_queue<int> queue_type;
  typedef std::vector<int> vector_type;
  
  srandom(time(0) * getpid());
  
  {
    queue_type queue(0);
    vector_type vector;
    
    std::cout << "START PUSH" << std::endl;
    for (int i = 0; i < 1024; ++ i) {
      vector.push_back(random());
      queue.push(vector.back());
    }
  
    std::cout << "START POP" << std::endl;
    vector_type vector_pop;
    for (;;) {
      int value;
      if (queue.pop(value, true))
	vector_pop.push_back(value);
      else
	break;
    }
    
    if (vector != vector_pop)
      std::cerr << "DIFFER!" << std::endl;
  }
  
  std::cerr << "start threads" << std::endl;
  
  // now run thread version...
  int num_threads = 2;
  if (argc >= 2)
    num_threads = std::max(num_threads, atoi(argv[1]));
  
  srandom(time(0) * getpid());
  
  std::vector<boost::shared_ptr<boost::thread> > threads_mapper(num_threads);
  std::vector<boost::shared_ptr<boost::thread> > threads_reducer(num_threads);
  
  std::vector<integers_type>  integers_mapper(num_threads);
  std::vector<integers_type>  integers_reducer(num_threads);
  
  lockfree_queue_type     queue_mapper(256);
  lockfree_queue_set_type queue_reducer(num_threads);
  for (int i = 0; i < num_threads; ++ i)
    queue_reducer[i].reset(new lockfree_queue_type(256));
  
  for (int i = 0; i < num_threads; ++ i) {
    
    threads_mapper[i].reset(new boost::thread(Mapper(integers_mapper[i],
						     queue_mapper,
						     queue_reducer)));
    
    threads_reducer[i].reset(new boost::thread(Reducer(integers_reducer[i],
						       *queue_reducer[i])));
  }
  
  integers_type inserted;
  for (int i = 0; i < 1024 * 128; ++ i) {
    ptr_int_type value(new int(random() & (1024 * 1024 - 1)));
    
    inserted.insert(value);
    queue_mapper.push(value);
  }
    
  for (int i = 0; i < num_threads; ++ i)
    queue_mapper.push(ptr_int_type());
  
  for (int i = 0; i < num_threads; ++ i)
    threads_mapper[i]->join();
  threads_mapper.clear();
  
  for (int i = 0; i < num_threads; ++ i)
    queue_reducer[i]->push(ptr_int_type());
  
  for (int i = 0; i < num_threads; ++ i)
    threads_reducer[i]->join();
  threads_reducer.clear();
  
  
  integers_type received_mapper;
  integers_type received_reducer;
  for (int i = 0; i < num_threads; ++ i) {
    received_mapper.insert(integers_mapper[i].begin(), integers_mapper[i].end());
    received_reducer.insert(integers_reducer[i].begin(), integers_reducer[i].end());
  }
  
  
  if (received_mapper != inserted)
    std::cerr << "WARNING differently sent.."
	      << " received(mapper): " << received_mapper.size()
	      << " inserted: " << inserted.size()
	      << std::endl;
  
  if (received_reducer != inserted)
    std::cerr << "WARNING differently sent.."
	      << " received(reducder): " << received_reducer.size()
	      << " inserted: " << inserted.size()
	      << std::endl;
}
