#include <iostream>
#include <string>
#include <vector>

#include "alloc_vector.hpp"
#include "allocinfo_allocator.hpp"

typedef std::vector<char, utils::allocinfo_allocator<char, std::allocator<char> > > string_type;
typedef utils::alloc_vector<string_type, utils::allocinfo_allocator<string_type, std::allocator<string_type> > > string_set_type;

int main(int argc, char** argv)
{
  srandom(time(0));
  
  {
    string_set_type strings;
  
    {
      std::string word;
      while (std::cin >> word)
	strings[random() % (1024 * 1024)] = string_type(word.begin(), word.end());
    }
  
    size_t exists = 0;
    for (size_t i = 0; i != strings.size(); ++ i)
      exists += strings.exists(i);
  
    std::cout << "size: " << strings.size() << " exist: " << exists << " allocated: "<< utils::allocinfo().allocated() << std::endl;
  
    string_set_type strings1;
    string_set_type strings2;
    string_set_type strings3;
  
    strings1 = strings;
    
    strings.clear();
    std::cout << "size: " << strings1.size() << " allocated: "<< utils::allocinfo().allocated() << std::endl;

    strings2 = strings1;
    
    std::cout << "size: " << strings2.size() << " allocated: "<< utils::allocinfo().allocated() << std::endl;
  
    strings3 = strings;

    std::cout << "size: " << strings3.size() << " allocated: "<< utils::allocinfo().allocated() << std::endl;
  }
  
  std::cout << "allocated: "<< utils::allocinfo().allocated() << std::endl;
}
