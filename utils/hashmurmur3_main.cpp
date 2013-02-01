#include <iostream>
#include <string>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <utils/search.hpp>
#include <utils/random_seed.hpp>
#include <utils/unordered_set.hpp>

#include "hashmurmur3.hpp"

template <size_t Size, typename Gen, typename Hasher32, typename Hasher64>
void test_hash(Gen& gen, const Hasher32& hasher32, const Hasher64& hasher64)
{
  uint8_t key[Size];
  
  for (size_t i = 0; i != 1024 * 4; ++ i) {
    for (size_t j = 0; j != Size; ++ j)
      key[j] = gen();
    
    if (hasher32(key) != hasher32(key, key + Size, 0))
      std::cerr << "different 32-bit hash...?" << std::endl;
    if (hasher64(key) != hasher64(key, key + Size, 0))
      std::cerr << "different 64-bit hash...?" << std::endl;
  }
}

int main(int argc, char** argv)
{
  std::string line;

  utils::hashmurmur3<uint64_t> hasher64;
  utils::hashmurmur3<uint32_t> hasher32;
  
  // test for constants...
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  boost::uniform_int<uint32_t> range(0, 1024 * 1024 * 3);
  boost::variate_generator<boost::mt19937, boost::uniform_int<uint32_t> > gen(generator, range);
  
  // random queries...
  test_hash<1>(gen, hasher32, hasher64);
  test_hash<2>(gen, hasher32, hasher64);
  test_hash<3>(gen, hasher32, hasher64);
  test_hash<4>(gen, hasher32, hasher64);
  test_hash<5>(gen, hasher32, hasher64);
  test_hash<6>(gen, hasher32, hasher64);
  test_hash<7>(gen, hasher32, hasher64);
  test_hash<8>(gen, hasher32, hasher64);
  test_hash<9>(gen, hasher32, hasher64);
  test_hash<10>(gen, hasher32, hasher64);
  test_hash<11>(gen, hasher32, hasher64);
  test_hash<12>(gen, hasher32, hasher64);
  test_hash<13>(gen, hasher32, hasher64);
  test_hash<14>(gen, hasher32, hasher64);
  test_hash<15>(gen, hasher32, hasher64);
  test_hash<16>(gen, hasher32, hasher64);
  test_hash<17>(gen, hasher32, hasher64);
  test_hash<18>(gen, hasher32, hasher64);
  test_hash<19>(gen, hasher32, hasher64);
  test_hash<20>(gen, hasher32, hasher64);
  test_hash<21>(gen, hasher32, hasher64);
  test_hash<22>(gen, hasher32, hasher64);
  test_hash<23>(gen, hasher32, hasher64);
  test_hash<24>(gen, hasher32, hasher64);
  test_hash<25>(gen, hasher32, hasher64);
  test_hash<26>(gen, hasher32, hasher64);
  test_hash<27>(gen, hasher32, hasher64);
  test_hash<28>(gen, hasher32, hasher64);
  test_hash<29>(gen, hasher32, hasher64);
  test_hash<30>(gen, hasher32, hasher64);
  test_hash<31>(gen, hasher32, hasher64);
  test_hash<32>(gen, hasher32, hasher64);
  test_hash<33>(gen, hasher32, hasher64);
  test_hash<34>(gen, hasher32, hasher64);
  test_hash<35>(gen, hasher32, hasher64);
  test_hash<36>(gen, hasher32, hasher64);
  test_hash<37>(gen, hasher32, hasher64);
  test_hash<38>(gen, hasher32, hasher64);
  test_hash<39>(gen, hasher32, hasher64);

  // read and compute
  while (std::getline(std::cin, line))
    std::cout << hasher64(line.begin(), line.end(), 0)
	      << ' '
	      << hasher32(line.begin(), line.end(), 0)
	      << std::endl;
}
