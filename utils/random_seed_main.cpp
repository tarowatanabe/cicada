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
    const double rounded1 = int(rand * 100) * double(0.01);
    const double rounded2 = int(rounded1 * 100) * double(0.01);
    const double rounded3 = int(rounded2 * 100) * double(0.01);
    const double rounded4 = int(rounded3 * 100) * double(0.01);
    
    if (rounded2 != rounded3)
      std::cerr << "differ: " << rand << " 2nd = " << rounded2 << " 3rd = " << rounded3 << std::endl;
  }
}
