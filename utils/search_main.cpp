#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

#include <utils/search.hpp>
#include <utils/random_seed.hpp>
#include <utils/unordered_set.hpp>

template <typename Tp>
void random_test(Tp upper, size_t entries, size_t queries)
{
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  boost::uniform_int<Tp> range(0, upper);
  boost::variate_generator<boost::mt19937, boost::uniform_int<Tp> > gen(generator, range);
  
  typedef std::vector<Tp> bucket_type;
  
  typename utils::unordered_set<Tp>::type uniques;
  
  for (size_t i = 0; i != entries; ++ i)
    uniques.insert(gen());
  
  bucket_type buckets(uniques.begin(), uniques.end());
  std::sort(buckets.begin(), buckets.end());
  
  // random queries...
  for (size_t i = 0; i != queries; ++ i) {
    const Tp key = gen();
    
    typename bucket_type::iterator biter = utils::binary_search(buckets.begin(), buckets.end(), key);
    typename bucket_type::iterator iiter = utils::interpolation_search(buckets.begin(), buckets.end(), key);

    if (biter != iiter) {
      std::cerr << "different result...?" << std::endl;

      if (biter == buckets.end())
	std::cerr << "binary search failed" << std::endl;
      
      if (iiter == buckets.end())
	std::cerr << "interpolation search failed" << std::endl;
      
      if (biter != buckets.end() && iiter != buckets.end()) {
	std::cerr << "binary pos: " << (biter - buckets.begin()) << std::endl
		  << "interpolation pos: " << (iiter - buckets.begin()) << std::endl;
	
      }
    }
  }
  
  // check if always succeed...!
  for (size_t i = 0; i != buckets.size(); ++ i) {
    typename bucket_type::iterator biter = utils::binary_search(buckets.begin(), buckets.end(), buckets[i]);
    typename bucket_type::iterator iiter = utils::interpolation_search(buckets.begin(), buckets.end(), buckets[i]);

    if (biter != iiter) {
      std::cerr << "different result...?" << std::endl;

      if (biter == buckets.end())
	std::cerr << "binary search failed" << std::endl;
      
      if (iiter == buckets.end())
	std::cerr << "interpolation search failed" << std::endl;
      
      if (biter != buckets.end() && iiter != buckets.end()) {
	std::cerr << "binary pos: " << (biter - buckets.begin()) << std::endl
		  << "interpolation pos: " << (iiter - buckets.begin()) << std::endl;
	
      }
    }
  }

}

int main(int argc, char** argv)
{
  std::cerr << "100, 100, 200" << std::endl;
  random_test<int>(100, 100, 200);

  std::cerr << "1000, 1000, 2000" << std::endl;
  random_test<int>(1000, 1000, 2000);
  
  std::cerr << "32000, 1000, 2000" << std::endl;
  random_test<int>(32000, 1000, 2000);

  std::cerr << "32000000, 10000, 2000" << std::endl;
  random_test<int>(32000000, 10000, 2000);
}
