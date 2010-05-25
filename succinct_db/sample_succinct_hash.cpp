
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include "succinct_hash.hpp"

#include <utils/hashmurmur.hpp>

int main(int argc, char** argv)
{
  std::vector<std::string> strings;
  succinctdb::succinct_hash<char> succinct_hash;

  utils::hashmurmur<uint64_t> hasher;

  std::string word;
  while (std::cin >> word) {
    strings.push_back(word);
    succinct_hash.insert(word.c_str(), word.size(), hasher(word.begin(), word.end(), 0));
  }
  
  for (int i = 0; i < strings.size(); ++ i)
    if (succinct_hash.find(strings[i].c_str(), strings[i].size(), hasher(strings[i].begin(), strings[i].end(), 0)) == uint32_t(-1))
      std::cerr << "NOT FOUND?" << std::endl;

  succinct_hash.write("tmpmtp.hash");
  
  succinctdb::succinct_hash_mapped<char> succinct_hash_mapped("tmpmtp.hash");
  for (int i = 0; i < strings.size(); ++ i)
    if (succinct_hash_mapped.find(strings[i].c_str(), strings[i].size(), hasher(strings[i].begin(), strings[i].end(), 0)) == uint32_t(-1))
      std::cerr << "NOT FOUND? mapped: " << strings[i] << std::endl;

  succinctdb::succinct_hash_stream<char> succinct_hash_stream("tmptmp.hash.stream");

  succinctdb::succinct_hash<char>::const_iterator siter = succinct_hash.begin();
  succinctdb::succinct_hash_mapped<char>::const_iterator miter = succinct_hash_mapped.begin();
  for (/**/; siter != succinct_hash.end(); ++ siter, ++ miter) {
    if (siter.size() != miter.size())
      std::cerr << "size differ!" << std::endl;
    if (! std::equal(siter.begin(), siter.end(), miter.begin()))
      std::cerr << "key differ!" << std::endl;
    
    const std::string key(siter.begin(), siter.end());
    succinct_hash_stream.insert(key.c_str(), key.size(), hasher(key.begin(), key.end(), 0));
  }
  succinct_hash_stream.close();
  
  succinctdb::succinct_hash_mapped<char> succinct_hash_mapped_stream("tmptmp.hash.stream");
  for (int i = 0; i < succinct_hash.size(); ++ i) {
    if (succinct_hash[i].size() != succinct_hash_mapped[i].size())
      std::cerr << "size differ!" << std::endl;
    if (! std::equal(succinct_hash[i].begin(), succinct_hash[i].end(), succinct_hash_mapped[i].begin()))
      std::cerr << "key differ!" << std::endl;
    
    if (succinct_hash[i].size() != succinct_hash_mapped_stream[i].size())
      std::cerr << "size differ!" << std::endl;
    if (! std::equal(succinct_hash[i].begin(), succinct_hash[i].end(), succinct_hash_mapped_stream[i].begin()))
      std::cerr << "key differ!" << std::endl;
  }
  
  for (int i = 0; i < strings.size(); ++ i)
    if (succinct_hash_mapped_stream.find(strings[i].c_str(), strings[i].size(), hasher(strings[i].begin(), strings[i].end(), 0)) == uint32_t(-1))
      std::cerr << "NOT FOUND? mapped: " << strings[i] << std::endl;
}
