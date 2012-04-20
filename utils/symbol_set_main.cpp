//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <string>

#include <set>
#include <map>
#include <sstream>

#include <utils/symbol_set.hpp>
#include <utils/symbol_map.hpp>


int main (int argc, char** argv)
{
  srandom(time(0));

  {
    typedef std::vector<int> vec_type;
    typedef utils::symbol_map<vec_type, vec_type> symbol_set_type;
    typedef std::vector<symbol_set_type> symbol_map_type;

    symbol_map_type symbol_map;
    symbol_map.clear();
    symbol_map.resize(500);
  }
  
  // int ...
  {
    typedef utils::symbol_set<int> symbol_set_type;
    typedef std::set<int> set_type;
    
    symbol_set_type symbol_set;
    set_type         set;
    
    std::cout << "sizeof int: " << sizeof(symbol_set) << std::endl;
    
    for (int k = 0; k < 4; ++ k) {
    
      for (int i = 0; i < 1024 * 2; ++ i) {
	if (! symbol_set.insert(i).second)
	  std::cerr << "ALREADY INSERTED?" << std::endl;
	set.insert(i);
      }
    
      std::cout << "1st size:" << " symbol: " << symbol_set.size() << " set: " << set.size() << std::endl;
    
      for (int i = 0; i < 1024 * 2; ++ i) {
	if (symbol_set.insert(i).second)
	  std::cerr << "NEW INSERTION?" << std::endl;
	set.insert(i);
      }
      
      for (int i = 0; i < 1024 * 2; ++ i) {
	if (symbol_set.find(i) == symbol_set.end())
	  std::cerr << "NOT FOUND?" << std::endl;
      
	if (symbol_set[i] != i)
	  std::cerr << "different position?" << std::endl;
      }
    
      std::cout << "2nd size:" << " symbol: " << symbol_set.size() << " set: " << set.size() << std::endl;
    
      for (int j = 0; j < 1024 * 4; ++ j) {
	int i = random() % (1024 * 2);
      
	if (symbol_set.find(i) == symbol_set.end())
	  std::cerr << "NOT INSERTED?" << std::endl;
      
	if (symbol_set.insert(i).second)
	  std::cerr << "NEW INSERTION?" << std::endl;
      }
    
      for (int j = 0; j < 1024 * 4; ++ j) {
	int i = random() % (1024 * 2) + 1024 * 4;
      
	if (symbol_set.find(i) != symbol_set.end())
	  std::cerr << "found?" << std::endl;
      }
      
      std::cout << "try erase: " << std::endl;
      
      for (int j = 0; j < 1024 * 2; ++ j) {
	int i = random() % (1024 * 2);
	
	symbol_set.erase(symbol_set_type::index_type(i));
	set.erase(i);
      }
      
      for (int i = 0; i < 1024 * 2; ++ i)  {
	if (set.find(i) == set.end())
	  if (symbol_set.find(i) != symbol_set.end())
	    std::cerr << "erased but exists?" << std::endl;
	
	if (symbol_set.find(i) == symbol_set.end())
	  if (set.find(i) != set.end())
	    std::cerr << "erased but exists?" << std::endl;
      }
      
      std::cout << "3rd size:" << " symbol: " << symbol_set.size() << " set: " << set.size() << std::endl;
    
      symbol_set.clear();
      set.clear();
      
      std::cout << "CLEAR size:" << " symbol: " << symbol_set.size() << " set: " << set.size() << std::endl;
    }
  }
  
}
