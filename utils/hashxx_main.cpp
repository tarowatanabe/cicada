#include <iostream>
#include <string>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <utils/search.hpp>
#include <utils/random_seed.hpp>
#include <utils/unordered_set.hpp>

#include "hashxx.hpp"

int main(int argc, char** argv)
{
  std::string line;

  utils::hashxx<uint64_t> hasher64;
  utils::hashxx<uint32_t> hasher32;
  
  // test for constants...
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  boost::uniform_int<uint32_t> range(0, 1024 * 1024 * 3);
  boost::variate_generator<boost::mt19937, boost::uniform_int<uint32_t> > gen(generator, range);
  
  // random queries...
  for (size_t i = 0; i != 1024 * 4; ++ i) {
    {
      uint64_t key[4];
      
      for (size_t j = 0; j != 4; ++ j) {
	key[j] = gen();
	key[j] <<= 32;
	key[j] |= gen();
      }
      
      if (hasher32(key) != hasher32(key, key + 4, 0))
	std::cerr << "different 32-bit hash...?" << std::endl;
      if (hasher64(key) != hasher64(key, key + 4, 0))
	std::cerr << "different 64-bit hash...?" << std::endl;
    }
    {
      uint64_t key[3];
      
      for (size_t j = 0; j != 3; ++ j) {
	key[j] = gen();
	key[j] <<= 32;
	key[j] |= gen();
      }
      
      if (hasher32(key) != hasher32(key, key + 3, 0))
	std::cerr << "different 32-bit hash...?" << std::endl;
      if (hasher64(key) != hasher64(key, key + 3, 0))
	std::cerr << "different 64-bit hash...?" << std::endl;      
    }

    {
      uint64_t key;
      
      key = gen();
      key <<= 32;
      key |= gen();
      
      if (hasher32(key) != hasher32(&key, &key + 1, 0))
	std::cerr << "different 32-bit hash...?" << std::endl;
      if (hasher64(key) != hasher64(&key, &key + 1, 0))
	std::cerr << "different 64-bit hash...?" << std::endl;            
    }
  }

  while (std::getline(std::cin, line))
    std::cout << hasher64(line.begin(), line.end(), 0)
	      << ' '
	      << hasher32(line.begin(), line.end(), 0)
	      << std::endl;
}
