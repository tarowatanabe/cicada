#include <iostream>
#include <string>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <utils/search.hpp>
#include <utils/random_seed.hpp>
#include <utils/unordered_set.hpp>
#include <utils/resource.hpp>

#include "hashxx.hpp"
#include "hashmurmur.hpp"
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

template <size_t Num>
struct test_hash_loop
{
  template <typename Gen, typename Hasher32, typename Hasher64>
  static inline
  void test(Gen& gen, const Hasher32& hasher32, const Hasher64& hasher64)
  {
    test_hash<Num>(gen, hasher32, hasher64);
    
    test_hash_loop<Num-1>::test(gen, hasher32, hasher64);
  }
};

template <>
struct test_hash_loop<0>
{
  template <typename Gen, typename Hasher32, typename Hasher64>
  static inline
  void test(Gen& gen, const Hasher32& hasher32, const Hasher64& hasher64)
  {
  }
};

int main(int argc, char** argv)
{
  std::string line;

  
  // test for constants...
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  boost::uniform_int<uint32_t> range(0, 255);
  boost::variate_generator<boost::mt19937, boost::uniform_int<uint32_t> > gen(generator, range);
  
  {
    utils::resource start;
    test_hash_loop<256>::test(gen, utils::hashxx<uint32_t>(), utils::hashxx<uint64_t>());
    utils::resource end;
    
    std::cout << "hashxx cpu time: " << (end.cpu_time() - start.cpu_time())
	      << " user time: " << (end.user_time() - start.user_time())
	      << std::endl;
  }

  {
    utils::resource start;
    test_hash_loop<256>::test(gen, utils::hashmurmur<uint32_t>(), utils::hashmurmur<uint64_t>());
    utils::resource end;
    
    std::cout << "hashmurmur cpu time: " << (end.cpu_time() - start.cpu_time())
	      << " user time: " << (end.user_time() - start.user_time())
	      << std::endl;
  }

  {
    utils::resource start;
    test_hash_loop<256>::test(gen, utils::hashmurmur3<uint32_t>(), utils::hashmurmur3<uint64_t>());
    utils::resource end;
    
    std::cout << "hashmurmur3 cpu time: " << (end.cpu_time() - start.cpu_time())
	      << " user time: " << (end.user_time() - start.user_time())
	      << std::endl;
  }

}
