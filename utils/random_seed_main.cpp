#include <cmath>
#include <iostream>

#include <boost/random.hpp>

#include "utils/random_seed.hpp"

int main(int argc, char** argv)
{
  std::cout << utils::random_seed() << std::endl;
  
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  
  boost::random::uniform_real_distribution<double> gen(-1, 1);
  
  for (size_t i = 0; i != 1024; ++ i) {
    const float rand = gen(generator);
    const float rounded1 = roundf(rand * 128) / float(128);
    const float rounded2 = roundf(rounded1 * 128) / float(128);
    const float rounded3 = roundf(rounded2 * 128) / float(128);
    const float rounded4 = roundf(rounded3 * 128) / float(128);
    
    if (rounded1 != rounded2)
      std::cerr << "differ: " << rand << " 1st = " << rounded1 << " 2nd = " << rounded2 << std::endl;
    if (rounded2 != rounded3)
      std::cerr << "differ: " << rand << " 2nd = " << rounded2 << " 3rd = " << rounded3 << std::endl;
    if (rounded4 != rounded4)
      std::cerr << "differ: " << rand << " 3rd = " << rounded3 << " 4th = " << rounded4 << std::endl;
  }
}
