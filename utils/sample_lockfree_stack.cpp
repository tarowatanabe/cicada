
#include <iostream>

#include <set>

#include <boost/thread.hpp>

#include <utils/lockfree_stack.hpp>

typedef utils::lockfree_stack<int> stack_type;
typedef std::multiset<int> integers_type;

struct Producer
{
  stack_type&    stack;
  integers_type& integers;
  
  Producer(stack_type& __stack,
	   integers_type& __integers)
    : stack(__stack), integers(__integers) {}

  void operator()()
  {
    for (int i = 0; i < 1024 * 128; ++ i) {
      const int value = random() % (1024 * 256);
      integers.insert(value);
      stack.push(value);
    }
    stack.push(-1);
  }
};

struct Consumer
{
  stack_type&    stack;
  integers_type& integers;
  
  Consumer(stack_type& __stack,
	   integers_type& __integers)
    : stack(__stack), integers(__integers) {}

  void operator()()
  {
    int value;
    while (1) {
      if (! stack.pop(value)) {
	boost::thread::yield();
	continue;
      }
      
      if (value < 0) break;
      
      integers.insert(value);
    }
  }
};

int main(int argc, char** argv)
{
  
  stack_type stack;
  integers_type integers_producer;
  integers_type integers_consumer;

  srandom(time() * getpid());

  std::auto_ptr<boost::thread> consumer(new boost::thread(Consumer(stack, integers_consumer)));
  std::auto_ptr<boost::thread> producer(new boost::thread(Producer(stack, integers_producer)));
  
  if (integers_producer != integers_consumer)
    std::cerr << "different!" << std::endl;
    
}
