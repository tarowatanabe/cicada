
#include <iostream>
#include <sstream>
#include <string>

#include <map>

#include "succinct_trie_database.hpp"

int main(int argc, char** argv)
{
  typedef succinctdb::succinct_trie_database<char, char> succinct_db_type;
  typedef std::multimap<std::string, std::string> map_db_type;

  srandom(time(0) * getpid());

  {
    // empty...
    succinct_db_type succinct_db("tmptmp.database", succinct_db_type::WRITE);
    succinct_db.close();
    succinct_db.open("tmptmp.database");
    
    std::cerr << "db size: " << succinct_db.size() << std::endl;
  }
  
  for (int iter = 0; iter < 4; ++ iter) {
    
    succinct_db_type succinct_db("tmptmp.database", succinct_db_type::WRITE);
    map_db_type      map_db;
    
    for (int i = 0; i < 1024 * 16; ++ i) {
      const int value_key  = random();
      const int value_data = random();
      
      std::ostringstream stream_key;
      std::ostringstream stream_data;
      stream_key  << value_key;
      stream_data << value_data;
      const std::string str_key = stream_key.str();
      const std::string str_data = stream_data.str();
      
      succinct_db.insert(str_key.c_str(), str_key.size(), str_data.c_str(), str_data.size());
      map_db.insert(std::make_pair(str_key, str_data));
    }
    succinct_db.close();
    succinct_db.open("tmptmp.database", succinct_db_type::READ);

    std::cerr << "db size: " << succinct_db.size() << std::endl;
    
    for (map_db_type::const_iterator iter = map_db.begin(); iter != map_db.end(); ++ iter) {
      
      const succinct_db_type::size_type node_pos = succinct_db.find(iter->first.c_str(), iter->first.size());
      
      if (! succinct_db.is_valid(node_pos))
	std::cerr << "out of range..?" << std::endl;
      
      if (! succinct_db.exists(node_pos))
	std::cerr << "NO KEY FOUND?" << std::endl;
      else {
	
	bool found = false;
	
	succinct_db_type::const_iterator siter = succinct_db.begin(node_pos);
	succinct_db_type::const_cursor csiter = siter.begin();

	for (succinct_db_type::const_cursor citer = succinct_db.cbegin(node_pos); citer != succinct_db.cend(); ++ citer, ++ csiter) {
	  
	  if ((*citer).size() == iter->second.size() && std::equal(iter->second.begin(), iter->second.end(), (*citer).begin()))
	    found = true;

	  if (citer.node() != csiter.node())
	    std::cerr << "different cursor?" << std::endl;
	}
	if (! found)
	  std::cerr << "no data?" << std::endl;
      }
      
    }
  }

  std::cerr << "integers" << std::endl;
  
  
  for (int iter = 0; iter < 4; ++ iter) {
    
    typedef succinctdb::succinct_trie_database<int, char> succinct_db_type;
    typedef std::multimap<std::vector<int>, std::string> map_db_type;
    
    succinct_db_type succinct_db("tmptmp.database", succinct_db_type::WRITE);
    map_db_type      map_db;
    
    for (int i = 0; i < 1024 * 16; ++ i) {
      std::vector<int> value_key(2);
      
      value_key[0]  = random() & 0xffff;
      value_key[1]  = random() & 0xffff;
      
      const int value_data = random();
      
      std::ostringstream stream_data;
      stream_data << value_data;
      const std::string str_data = stream_data.str();
      
      succinct_db.insert(&(*value_key.begin()), value_key.size(), str_data.c_str(), str_data.size());
      map_db.insert(std::make_pair(value_key, str_data));
    }
    succinct_db.close();
    succinct_db.open("tmptmp.database", succinct_db_type::READ);

    std::cerr << "db size: " << succinct_db.size() << std::endl;
    
    for (map_db_type::const_iterator iter = map_db.begin(); iter != map_db.end(); ++ iter) {
      
      const succinct_db_type::size_type node_pos = succinct_db.find(&(*iter->first.begin()), iter->first.size());
      
      if (! succinct_db.is_valid(node_pos))
	std::cerr << "out of range..?" << std::endl;
      
      if (! succinct_db.exists(node_pos))
	std::cerr << "NO KEY FOUND?" << std::endl;
      else {
	
	bool found = false;
	
	succinct_db_type::const_iterator siter = succinct_db.begin(node_pos);
	succinct_db_type::const_cursor csiter = siter.begin();

	for (succinct_db_type::const_cursor citer = succinct_db.cbegin(node_pos); citer != succinct_db.cend(); ++ citer, ++ csiter) {
	  
	  if ((*citer).size() == iter->second.size() && std::equal(iter->second.begin(), iter->second.end(), (*citer).begin()))
	    found = true;

	  if (citer.node() != csiter.node())
	    std::cerr << "different cursor?" << std::endl;
	}
	if (! found)
	  std::cerr << "no data?" << std::endl;
      }
      
    }
  }
}
