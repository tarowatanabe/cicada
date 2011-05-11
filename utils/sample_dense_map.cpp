
#include <string>
#include <iostream>

#include <utils/dense_set.hpp>
#include <utils/dense_map.hpp>

#include <map>

int main(int argc, char** argv)
{
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
