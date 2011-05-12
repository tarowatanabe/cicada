
#include <string>
#include <iostream>

#include <utils/dense_set.hpp>
#include <utils/dense_map.hpp>
#include <utils/resource.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>

#include <boost/random.hpp>

#include <map>

template <typename Map, typename Generator>
void process_insert(const char* name, Map& map, Generator& generator)
{
  utils::resource start;

  boost::uniform_int<int> gen;
  for (int i = 0; i != 1024 * 4; ++ i)
    map.insert(std::make_pair(gen(generator), gen(generator)));

  utils::resource end;

  std::cout << "insert:" << name
	    << " cpu time: " << (end.cpu_time() - start.cpu_time())
	    << " user time: " << (end.user_time() - start.user_time())
	    << std::endl;
}

template <typename Map, typename Generator>
void process_erase(const char* name, Map& map, Generator& generator)
{
  utils::resource start;
  
  boost::uniform_int<int> gen(0, 1);
  while (! map.empty()) {
    typename Map::iterator iter = map.begin();
    if (gen(generator))
      ++ iter;
    if (iter != map.end())
      map.erase(iter);
  }
  
  utils::resource end;

  std::cout << "erase:" << name
	    << " cpu time: " << (end.cpu_time() - start.cpu_time())
	    << " user time: " << (end.user_time() - start.user_time())
	    << std::endl;
}



int main(int argc, char** argv)
{
  boost::mt19937 generator;
  generator.seed(time(0) * getpid());

  {
    
#ifdef HAVE_TR1_UNORDERED_MAP
    const char* hashname = "unordered-map";
    std::tr1::unordered_map<int, int, utils::hashmurmur<size_t> > hashmap;
#else
    const char* hashname = "hash-map";
    sgi::hash_map<int, int, utils::hashmurmur<size_t> > hashmap;
#endif
    const char* densename = "dense-map";
    utils::dense_map<int, int, utils::hashmurmur<size_t> > densemap;
    
    process_insert(hashname, hashmap, generator);
    process_insert(densename, densemap, generator);

    process_erase(hashname, hashmap, generator);
    process_erase(densename, densemap, generator);
    
    process_insert(hashname, hashmap, generator);
    process_insert(densename, densemap, generator);
  }
  
  
  
  utils::dense_map<std::string, int> densemap;
  std::map<std::string, int> stdmap;
  
  std::string word;
  while (std::cin >> word) {
    ++ densemap[word];
    ++ stdmap[word];
  }
  
  std::map<std::string, int>::const_iterator iter_end = stdmap.end();
  for (std::map<std::string, int>::const_iterator iter = stdmap.begin(); iter != iter_end; ++ iter) {
    
    utils::dense_map<std::string, int>::iterator siter = densemap.find(iter->first);
    if (siter == densemap.end())
      std::cerr << "no key? " << iter->first << std::endl;
    else if (iter->second != siter->second)
      std::cerr << "different data ? " << iter->second << ' ' << siter->second << std::endl;
  }
  
  // erasing
  for (std::map<std::string, int>::const_iterator iter = stdmap.begin(); iter != iter_end; ++ iter) {
    
    if (densemap.find(iter->first) == densemap.end())
      std::cerr << "no key? " << iter->first << std::endl;
    else {
      densemap.erase(iter->first);
      
      if (densemap.find(iter->first) != densemap.end())
	std::cerr << "not erased! " << iter->first << std::endl;
    }
  }
  
  if (! densemap.empty())
    std::cerr << "after erasing, not empty!" << std::endl;

  densemap.insert(stdmap.begin(), stdmap.end());
  for (std::map<std::string, int>::const_iterator iter = stdmap.begin(); iter != iter_end; ++ iter) {
    
    utils::dense_map<std::string, int>::iterator siter = densemap.find(iter->first);
    if (siter == densemap.end())
      std::cerr << "no key? " << iter->first << std::endl;
    else if (iter->second != siter->second)
      std::cerr << "different data ? " << iter->second << ' ' << siter->second << std::endl;
  }

  stdmap.clear();
  {
    for (utils::dense_map<std::string, int>::const_iterator iter = densemap.begin(); iter != densemap.end(); ++ iter) {
      std::string data = iter->first;
      
    }
  }

}
