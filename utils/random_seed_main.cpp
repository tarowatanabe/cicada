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
    const double rand = gen(generator);
    const double rounded1 = int(rand * 128) / double(128);
    const double rounded2 = int(rounded1 * 128) / double(128);

    const int value = int(rand * 128);
    
    if (rounded1 != rounded2)
      std::cerr << "differ: " << rand << " 1st = " << rounded1 << " 2nd = " << rounded2 << std::endl;
  }
}
