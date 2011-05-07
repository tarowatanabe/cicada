
#include <string>
#include <iostream>

#include <utils/sparse_set.hpp>
#include <utils/sparse_map.hpp>

#include <utils/dense_set.hpp>
#include <utils/dense_map.hpp>

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
  
  
}
