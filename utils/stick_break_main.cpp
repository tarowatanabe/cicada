#include "utils/stick_break.hpp"
#include "utils/sampler.hpp"

#include <iostream>

int main(int argc, char** argv)
{
  typedef utils::sampler<boost::mt19937> sampler_type;
  typedef utils::stick_break<> breaker_type;

  sampler_type sampler;
  
  breaker_type breaker(0.5, 10.0);
  
  for (ptrdiff_t i = 0; i != 10; ++ i) {
    breaker.increment(sampler);
    
    std::cerr << "i=" << i << " size: " << breaker.size() << std::endl;
    
    for (size_t stick = 0; stick != breaker.size(); ++ stick)
      std::cerr << "stick=" << stick << " " << breaker[stick] << std::endl;
    std::cerr << "remain=" << breaker[breaker.size()] << std::endl;
    
  }
  
}
