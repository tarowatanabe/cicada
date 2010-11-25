//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <sstream>
#include <string>

#include <map>

#include "succinct_trie_db.hpp"
#include "utils/byte_aligned_code.hpp"

int main(int argc, char** argv)
{
  typedef succinctdb::succinct_trie_db<char, double> succinct_db_type;
  typedef std::multimap<std::string, double> map_db_type;

  srandom(time(0) * getpid());

  {
    // empty...
    succinct_db_type succinct_db("tmptmp.db.fixed", succinct_db_type::WRITE);
    succinct_db.close();
    succinct_db.open("tmptmp.db.fixed");
    
    std::cerr << "db size: " << succinct_db.size() << std::endl;
  }
  
  for (int iter = 0; iter < 4; ++ iter) {
    
    succinct_db_type succinct_db("tmptmp.db.fixed", succinct_db_type::WRITE);
    map_db_type      map_db;
    
    for (int i = 0; i < 1024 * 16; ++ i) {
      const int value = random();
      const double value_double = double(random()) /  random();
      
      std::ostringstream stream;
      stream << value;
      const std::string value_str = stream.str();
      
      if (! value_str.empty()) {
	map_db.insert(std::make_pair(value_str, value_double));
	const size_t size = succinct_db.insert(value_str.c_str(), value_str.size(), &value_double);
	if (size + 1 != map_db.size())
	  std::cerr << "different size...?" << std::endl;
      }
    }
    succinct_db.close();
    succinct_db.open("tmptmp.db.fixed", succinct_db_type::READ);

    std::cerr << "db size: " << succinct_db.size() << std::endl;
    
    for (map_db_type::const_iterator iter = map_db.begin(); iter != map_db.end(); ++ iter) {
      
      const succinct_db_type::size_type node_pos = succinct_db.find(iter->first.c_str(), iter->first.size());

      if (! succinct_db.is_valid(node_pos))
	std::cerr << "out of range..?" << std::endl;
      
      if (! succinct_db.exists(node_pos))
	std::cerr << "NO KEY FOUND?" << std::endl;
      else {
	
	bool found = false;
	for (succinct_db_type::const_cursor citer = succinct_db.cbegin(node_pos); citer != succinct_db.cend(); ++ citer)
	  if (*citer == iter->second)
	    found = true;
	if (! found)
	  std::cerr << "no data?" << std::endl;
      }
	
    }
  }

  {
    typedef std::vector<char> codes_type;

    typedef succinctdb::succinct_trie_db<char, double> succinct_db_type;
    typedef std::multimap<codes_type, double> map_db_type;
    
    for (int iter = 0; iter < 4; ++ iter) {
    
      succinct_db_type succinct_db("tmptmp.db.fixed", succinct_db_type::WRITE);
      map_db_type      map_db;

      char code[16];
      
      for (int i = 0; i < 1024 * 32; ++ i) {
	
	const int clip1 = random() & 0x03;
	const int clip2 = random() & 0x03;
	
	int value1 = random();
	switch (clip1) {
	case 3: value1 &= 0xffffffff; break;
	case 2: value1 &= 0x00ffffff; break;
	case 1: value1 &= 0x0000ffff; break;
	case 0: value1 &= 0x000000ff; break;
	}
	
	int value2 = random();
	switch (clip2) {
	case 3: value2 &= 0xffffffff; break;
	case 2: value2 &= 0x00ffffff; break;
	case 1: value2 &= 0x0000ffff; break;
	case 0: value2 &= 0x000000ff; break;
	}
	
	const double value_double = double(random()) /  random();
	
	const size_t code_size1 = utils::byte_aligned_encode(value1, code);
	const size_t code_size2 = utils::byte_aligned_encode(value2, code + code_size1);
	
	succinct_db.insert(code, code_size1 + code_size2, &value_double);
	map_db.insert(std::make_pair(codes_type(code, code + code_size1 + code_size2), value_double));
      }
      
      succinct_db.close();
      succinct_db.open("tmptmp.db.fixed", succinct_db_type::READ);
      
      std::cerr << "db size: " << succinct_db.size() << std::endl;
    
      for (map_db_type::const_iterator iter = map_db.begin(); iter != map_db.end(); ++ iter) {
	
	const succinct_db_type::size_type node_pos = succinct_db.find(&(*iter->first.begin()), iter->first.size());
	
	if (! succinct_db.is_valid(node_pos))
	  std::cerr << "out of range..?" << std::endl;
	if (! succinct_db.exists(node_pos))
	  std::cerr << "NO KEY FOUND?" << std::endl;
	else {
	  bool found = false;
	  for (succinct_db_type::const_cursor citer = succinct_db.cbegin(node_pos); citer != succinct_db.cend(); ++ citer)
	    if (*citer == iter->second)
	      found = true;
	  if (! found)
	    std::cerr << "no data?" << std::endl;
	}
      }
    }
    
  }

  std::cerr << "integer" << std::endl;
  
  for (int i = 0; i < 4; ++ i) {
    typedef succinctdb::succinct_trie_db<int, double> succinct_db_type;
    typedef std::multimap<int, double> map_db_type;
    
    succinct_db_type succinct_db("tmptmp.db.fixed", succinct_db_type::WRITE);
    map_db_type      map_db;
    
    for (int i = 0; i < 1024 * 16; ++ i) {
      const int value = random();
      const double value_double = double(random()) /  random();
      
      map_db.insert(std::make_pair(value, value_double));

      const size_t size = succinct_db.insert(&value, 1, &value_double);
      if (size + 1 != map_db.size())
	std::cerr << "differentn size...?" << std::endl;
    }
    
    succinct_db.close();
    succinct_db.open("tmptmp.db.fixed", succinct_db_type::READ);

    std::cerr << "db size: " << succinct_db.size() << std::endl;
    
    for (map_db_type::const_iterator iter = map_db.begin(); iter != map_db.end(); ++ iter) {
      
      const succinct_db_type::size_type node_pos = succinct_db.find(&iter->first, 1);

      if (! succinct_db.is_valid(node_pos))
	std::cerr << "out of range..?" << std::endl;
      
      if (! succinct_db.exists(node_pos))
	std::cerr << "NO KEY FOUND?" << std::endl;
      else {
	
	bool found = false;
	for (succinct_db_type::const_cursor citer = succinct_db.cbegin(node_pos); citer != succinct_db.cend(); ++ citer)
	  if (*citer == iter->second)
	    found = true;
	if (! found)
	  std::cerr << "no data?" << std::endl;
      }
	
    }
  }


}
