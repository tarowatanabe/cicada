
#include <stdint.h>

#include <iostream>
#include <sstream>
#include <algorithm>

#include "succinct_trie.hpp"

struct extract_key
{
  const std::string& operator()(const std::string& x) const { return x; }
};

int main(int argc, char** argv)
{
  typedef succinctdb::succinct_trie<uint8_t, int> succinct_trie_type;
  typedef succinctdb::succinct_trie_mapped<uint8_t, int> succinct_trie_mapped_type;
  typedef succinct_trie_type::size_type size_type;
  
  std::vector<std::string> months;
  months.push_back("January");
  months.push_back("February");
  months.push_back("March");
  months.push_back("Marche");
  months.push_back("April");
  months.push_back("May");
  months.push_back("June");
  months.push_back("July");
  months.push_back("Juck");
  months.push_back("Juck");
  months.push_back("Jucky");
  months.push_back("Jucky");
  months.push_back("August");
  months.push_back("September");
  months.push_back("October");
  months.push_back("November");
  months.push_back("December");
  
  std::sort(months.begin(), months.end());
  
  const uint8_t* data[months.size()];
  size_type length[months.size()];
  
  for (int i = 0; i < months.size(); ++ i) {
    data[i] = (const uint8_t*) months[i].c_str();
    length[i] = months[i].size();
    std::cout << "i = " << i << " " << data[i] << std::endl;
  }
  
  succinct_trie_type succinct_trie;
  
  //succinct_trie.build(months.size(), (const uint8_t**) data, length);
  succinct_trie.build(months.begin(), months.end(), extract_key());
  succinct_trie.write("tmptmp.trie");
  
  {
    succinct_trie_type succinct_trie_stream;
    succinct_trie_stream.build("tmptmp.trie.stream", months.begin(), months.end(), extract_key());
  }
  
  succinct_trie_mapped_type succinct_trie_mapped("tmptmp.trie");
  succinct_trie_mapped_type succinct_trie_mapped_stream("tmptmp.trie.stream");
  

  std::cout << "size: " << succinct_trie.size() << std::endl;
  std::cout << "index size: " << succinct_trie.index_size() << std::endl;
  
  std::cout << "size: " << succinct_trie_mapped.size() << std::endl;
  std::cout << "index size: " << succinct_trie_mapped.index_size() << std::endl;
  
  std::cout << "size: " << succinct_trie_mapped_stream.size() << std::endl;
  std::cout << "index size: " << succinct_trie_mapped_stream.index_size() << std::endl;

  for (int i = 0; i < months.size(); ++ i) {
    size_type node_pos = 0;
    size_type key_pos = 0;
    
    const int result = succinct_trie.traverse((const uint8_t*) data[i], node_pos, key_pos, length[i]);
    if (succinct_trie.exists(result)) {
      if (succinct_trie.is_next_sibling(node_pos))
	std::cerr << "next is sibling!" << std::endl;
      else
	std::cerr << "no sibling" << std::endl;
    }
    
    
    size_type node_pos_mapped = 0;
    size_type key_pos_mapped = 0;
    const int result_mapped = succinct_trie_mapped.traverse((const uint8_t*) data[i], node_pos_mapped, key_pos_mapped, length[i]);
    
    size_type node_pos_mapped_stream = 0;
    size_type key_pos_mapped_stream = 0;
    const int result_mapped_stream = succinct_trie_mapped_stream.traverse((const uint8_t*) data[i], node_pos_mapped_stream, key_pos_mapped_stream, length[i]);
    
    std::cout << "result: " << result << std::endl;
    std::cout << data[i] << std::endl;
    if (succinct_trie.exists(result))
      std::cout << data[succinct_trie.data(result)] << std::endl;
    else
      std::cerr << "no result?" << std::endl;

    if (result != result_mapped)
      std::cerr << "mapped differ?" << std::endl;
    
    if (result != result_mapped_stream)
      std::cerr << "mapped differ?" << std::endl;
  }
  
  {
    // access April iteratively
    size_type node_pos = 0;
    size_type key_pos = 0;
    
    const int result = succinct_trie.traverse((const uint8_t*) "A", node_pos, key_pos, 1);
    std::cout << node_pos << ' ' << result << ' ' << succinct_trie.parent(node_pos) << std::endl;

    key_pos = 0;
    const int result2 = succinct_trie.traverse((const uint8_t*) "p", node_pos, key_pos, 1);
    std::cout << node_pos << ' ' << result2 << ' ' << succinct_trie.parent(node_pos) << std::endl;
    
    key_pos = 0;
    const int result3 = succinct_trie.traverse((const uint8_t*) "r", node_pos, key_pos, 1);
    std::cout << node_pos << ' ' << result3 << ' ' << succinct_trie.parent(node_pos) << std::endl;

    key_pos = 0;
    const int result4 = succinct_trie.traverse((const uint8_t*) "i", node_pos, key_pos, 1);
    std::cout << node_pos << ' ' << result4 << ' ' << succinct_trie.parent(node_pos) << std::endl;
    
    succinct_trie_type::const_reverse_iterator riter_end = succinct_trie.rend();
    for (succinct_trie_type::const_reverse_iterator riter = succinct_trie.rbegin(node_pos); riter != riter_end; ++ riter)
      std::cout << riter.key() << ' ';
    std::cout << std::endl;
    
    key_pos = 0;
    const int result5 = succinct_trie.traverse((const uint8_t*) "l", node_pos, key_pos, 1);
    std::cout << node_pos << ' ' << result5 << ' ' << succinct_trie.parent(node_pos) << std::endl;    
  }
  
  
  std::ostringstream os_memory;
  for (succinct_trie_type::const_iterator iter = succinct_trie.begin(); iter != succinct_trie.end(); ++ iter) {
    std::vector<uint8_t> buffer(iter.key().begin(), iter.key().end());
    
    os_memory << std::string(buffer.begin(), buffer.end()) << " : " << iter.data() << " - " << iter.node() << std::endl;

    for (succinct_trie_type::cursor citer = iter.begin(); citer != iter.end(); ++ citer)
      os_memory << citer.data() << std::endl;
  }
  std::cout << os_memory.str();
  
  std::ostringstream os_mapped;
  {
    std::cout << "MAPPED START" << std::endl;
    
    for (succinct_trie_mapped_type::const_iterator iter = succinct_trie_mapped.begin(); iter != succinct_trie_mapped.end(); ++ iter) {
      std::vector<uint8_t> buffer(iter.key().begin(), iter.key().end());
      os_mapped << std::string(buffer.begin(), buffer.end()) << " : " << iter.data() << " - " << iter.node() << std::endl;
      
      for (succinct_trie_mapped_type::cursor citer = iter.begin(); citer != iter.end(); ++ citer)
	os_mapped << citer.data() << std::endl;
    }
    
    std::cout << "MAPPED END" << std::endl;
  }
  
  if (os_memory.str() != os_mapped.str())
    std::cerr << "different traversal?" << std::endl;
  
  std::ostringstream os_mapped_stream;
  
  {
    std::cout << "MAPPED STREAM START" << std::endl;
    
    for (succinct_trie_mapped_type::const_iterator iter = succinct_trie_mapped_stream.begin(); iter != succinct_trie_mapped_stream.end(); ++ iter) {
      std::vector<uint8_t> buffer(iter.key().begin(), iter.key().end());
      os_mapped_stream << std::string(buffer.begin(), buffer.end()) << " : " << iter.data() << " - " << iter.node() << std::endl;
      
      for (succinct_trie_mapped_type::cursor citer = iter.begin(); citer != iter.end(); ++ citer)
	os_mapped_stream << citer.data() << std::endl;
    }
    
    std::cout << "MAPPED STREAM END" << std::endl;
  }
  
  if (os_memory.str() != os_mapped_stream.str())
    std::cerr << "different traversal?" << std::endl;
  

  {
    size_type node_pos = 0;
    size_type key_pos = 0;
    
    const int result = succinct_trie.traverse((const uint8_t*) "J", node_pos, key_pos, 1);
    std::cout << result << std::endl;
    
    for (succinct_trie_type::const_iterator iter = succinct_trie.begin(node_pos); iter != succinct_trie.end(); ++ iter) {
      std::vector<uint8_t> buffer(iter.key().begin(), iter.key().end());
      std::cout << std::string(buffer.begin(), buffer.end()) << " : " << iter.data() << std::endl;
    }
  }
}
