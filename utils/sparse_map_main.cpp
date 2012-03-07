
#include <string>
#include <iostream>

#include <utils/sparse_set.hpp>
#include <utils/sparse_map.hpp>

#include <map>

int main(int argc, char** argv)
{
  utils::sparse_map<std::string, int> sparsemap;
  std::map<std::string, int> stdmap;
  
  std::string word;
  while (std::cin >> word) {
    ++ sparsemap[word];
    ++ stdmap[word];
  }
  
  std::map<std::string, int>::const_iterator iter_end = stdmap.end();
  for (std::map<std::string, int>::const_iterator iter = stdmap.begin(); iter != iter_end; ++ iter) {
    
    utils::sparse_map<std::string, int>::iterator siter = sparsemap.find(iter->first);
    if (siter == sparsemap.end())
      std::cerr << "no key? " << iter->first << std::endl;
    else if (iter->second != siter->second)
      std::cerr << "different data ? " << iter->second << ' ' << siter->second << std::endl;
    
  }
  
  // erasing
  for (std::map<std::string, int>::const_iterator iter = stdmap.begin(); iter != iter_end; ++ iter) {
    
    if (sparsemap.find(iter->first) == sparsemap.end())
      std::cerr << "no key? " << iter->first << std::endl;
    else {
      sparsemap.erase(iter->first);
      
      if (sparsemap.find(iter->first) != sparsemap.end())
	std::cerr << "not erased! " << iter->first << std::endl;
    }
  }
  
  if (! sparsemap.empty())
    std::cerr << "after erasing, not empty!" << std::endl;

  sparsemap.insert(stdmap.begin(), stdmap.end());
  for (std::map<std::string, int>::const_iterator iter = stdmap.begin(); iter != iter_end; ++ iter) {
    
    utils::sparse_map<std::string, int>::iterator siter = sparsemap.find(iter->first);
    if (siter == sparsemap.end())
      std::cerr << "no key? " << iter->first << std::endl;
    else if (iter->second != siter->second)
      std::cerr << "different data ? " << iter->second << ' ' << siter->second << std::endl;
  }

  stdmap.clear();
  {
    for (utils::sparse_map<std::string, int>::const_iterator iter = sparsemap.begin(); iter != sparsemap.end(); ++ iter) {
      std::string data = iter->first;
      
    }
  }
  
}
