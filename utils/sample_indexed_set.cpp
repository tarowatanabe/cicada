//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <string>

#include <set>
#include <map>
#include <sstream>

#include <utils/indexed_set.hpp>


int main (int argc, char** argv)
{
  srandom(time(0));

  // int ...
  {
    typedef utils::indexed_set<int> indexed_set_type;
    typedef std::set<int> set_type;
    
    indexed_set_type indexed_set;
    set_type         set;
    
    std::cout << "sizeof int: " << sizeof(indexed_set) << std::endl;
    
    for (int k = 0; k < 4; ++ k) {
    
      for (int i = 0; i < 1024 * 2; ++ i) {
	if (! indexed_set.insert(i).second)
	  std::cerr << "ALREADY INSERTED?" << std::endl;
	set.insert(i);
      }
    
      std::cout << "1st size:" << " indexed: " << indexed_set.size() << " set: " << set.size() << std::endl;
    
      for (int i = 0; i < 1024 * 2; ++ i) {
	if (indexed_set.insert(i).second)
	  std::cerr << "NEW INSERTION?" << std::endl;
	set.insert(i);
      }
      
      for (int i = 0; i < 1024 * 2; ++ i) {
	if (indexed_set.find(i) == indexed_set.end())
	  std::cerr << "NOT FOUND?" << std::endl;
      
	if (indexed_set[i] != i)
	  std::cerr << "different position?" << std::endl;
      }
    
      std::cout << "2nd size:" << " indexed: " << indexed_set.size() << " set: " << set.size() << std::endl;
    
      for (int j = 0; j < 1024 * 4; ++ j) {
	int i = random() % (1024 * 2);
      
	if (indexed_set.find(i) == indexed_set.end())
	  std::cerr << "NOT INSERTED?" << std::endl;
	else {
	  if (indexed_set.find(i) - indexed_set.begin() != i)
	    std::cerr << "different position?" << std::endl;
	}
      
	if (indexed_set.insert(i).second)
	  std::cerr << "NEW INSERTION?" << std::endl;
      
	if (indexed_set.insert(i).first - indexed_set.begin() != i)
	  std::cerr << "different position?" << std::endl;
      }
    
      for (int j = 0; j < 1024 * 4; ++ j) {
	int i = random() % (1024 * 2) + 1024 * 4;
      
	if (indexed_set.find(i) != indexed_set.end())
	  std::cerr << "found?" << std::endl;
      }
    
      indexed_set.clear();
      set.clear();
    
      std::cout << "CLEAR size:" << " indexed: " << indexed_set.size() << " set: " << set.size() << std::endl;
    }
  }
  
  // string..
  {
    typedef utils::indexed_set<std::string> indexed_set_type;
    typedef std::map<int, std::string> int_string_map_type;
    typedef std::map<std::string, int> string_int_map_type;
    
    indexed_set_type indexed_set;
    string_int_map_type set;
    
    int_string_map_type indexed_inverse;
    int_string_map_type set_inverse;

    std::cout << "sizeof string: " << sizeof(indexed_set) << std::endl;
    
    for (int k = 0; k < 4; ++ k) {
      
      // we will randomly generate...
      for (int i = 0; i < 1024 * 4; ++ i) {
	
	std::ostringstream stream;
	stream << (random() % 1024 * 2);
	const std::string value = stream.str();
	
	if (set.find(value) == set.end()) {
	  const int id = set.size();
	  set.insert(std::make_pair(value, id));
	  set_inverse.insert(std::make_pair(id, value));
	}
      
	std::pair<indexed_set_type::iterator, bool> result = indexed_set.insert(value);
	if (result.second)
	  indexed_inverse.insert(std::make_pair(result.first - indexed_set.begin(), value));
      }
    
      std::cerr << "1st size: indexed: " << indexed_set.size() << " map: " << set.size() << std::endl;
    
      // now verify...
      if (set_inverse != indexed_inverse)
	std::cerr << "differ ?" << std::endl;
      
      indexed_set.clear();
      indexed_inverse.clear();
      
      set.clear();
      set_inverse.clear();
    
      std::cerr << "clear size: indexed: " << indexed_set.size() << " map: " << set.size() << std::endl;
    }
    
  }
}
