#include <iostream>
#include <string>
#include <map>

#include <utils/vector_map.hpp>

int main(int argc, char** argv)
{
  typedef std::map<std::string, int>          map_map_type;
  typedef utils::vector_map<std::string, int> vec_map_type;

  std::cerr << "size: " << sizeof(vec_map_type) << std::endl;

  map_map_type map_map;
  vec_map_type vec_map;

  std::string token;
  while (std::cin >> token) {
    ++ map_map[token];
    ++ vec_map[token];
  }
  
  std::cerr << "map size: " << map_map.size() << std::endl
	    << "vec size: " << vec_map.size() << std::endl;
  
  {
    map_map_type::const_iterator miter = map_map.begin();
    for (vec_map_type::const_iterator viter = vec_map.begin(); viter != vec_map.end(); ++ viter, ++ miter) {
      if (viter->first != miter->first || viter->second != miter->second)
	std::cerr << "differ?"
		  << "\tmap: " << miter->first << ": " << miter->second << std::endl
		  << "\tvec: " << viter->first << ": " << viter->second << std::endl;
    }
  }

  {
    for (map_map_type::const_iterator miter = map_map.begin(); miter != map_map.end(); ++ miter) {
      vec_map_type::iterator viter = vec_map.find(miter->first);
      if (viter == vec_map.end())
	throw std::runtime_error("not found?");
      
      vec_map.erase(viter);
    }
  }

  std::cerr << "erased vec map" << std::endl
	    << "map size: " << map_map.size() << std::endl
	    << "vec size: " << vec_map.size() << std::endl;
  
  vec_map.insert(map_map.begin(), map_map.end());
  
  std::cerr << "map size: " << map_map.size() << std::endl
	    << "vec size: " << vec_map.size() << std::endl;
  
  {
    map_map_type::const_iterator miter = map_map.begin();
    for (vec_map_type::const_iterator viter = vec_map.begin(); viter != vec_map.end(); ++ viter, ++ miter) {
      if (viter->first != miter->first || viter->second != miter->second)
	std::cerr << "differ?"
		  << "\tmap: " << miter->first << ": " << miter->second << std::endl
		  << "\tvec: " << viter->first << ": " << viter->second << std::endl;
    }
  }

  vec_map.clear();
  {
    for (map_map_type::const_iterator miter = map_map.begin(); miter != map_map.end(); ++ miter)
      vec_map.insert(vec_map.end(), *miter);
  }

  std::cerr << "map size: " << map_map.size() << std::endl
	    << "vec size: " << vec_map.size() << std::endl;
  
  {
    map_map_type::const_iterator miter = map_map.begin();
    for (vec_map_type::const_iterator viter = vec_map.begin(); viter != vec_map.end(); ++ viter, ++ miter) {
      if (viter->first != miter->first || viter->second != miter->second)
	std::cerr << "differ?"
		  << "\tmap: " << miter->first << ": " << miter->second << std::endl
		  << "\tvec: " << viter->first << ": " << viter->second << std::endl;
    }
  }

  if (vec_map.size() > 2) {
    std::cerr << "erase range" << std::endl;
    vec_map.erase(vec_map.begin() + 1, vec_map.end() - 1);
    std::cerr << "map size: " << map_map.size() << std::endl
	      << "vec size: " << vec_map.size() << std::endl;

    map_map_type::iterator end = map_map.end();
    -- end;
    
    if (end->first != vec_map.back().first || end->second != vec_map.back().second)
      std::cerr << "differ?" << std::endl
		<< "\tmap: "<< end->first << ": " << end->second << std::endl
		<< "\tvec: "<< vec_map.back().first << ": " << vec_map.back().second << std::endl;
    
    if (map_map.begin()->first != vec_map.front().first || map_map.begin()->second != vec_map.front().second)
      std::cerr << "differ?" << std::endl
		<< "\tmap: "<< map_map.begin()->first << ": " << map_map.begin()->second << std::endl
		<< "\tvec: "<< vec_map.front().first << ": " << vec_map.front().second << std::endl;
  }
}
