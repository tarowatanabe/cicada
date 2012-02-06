//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <vector>
#include <string>

#include "utils/small_vector.hpp"
#include "utils/allocinfo_allocator.hpp"

int main(int argc, char** argv)
{
  srandom(time(0) * getpid());
  
  std::cout << "allocated: " << utils::allocinfo().allocated() << std::endl;

  utils::small_vector<std::string, utils::allocinfo_allocator<std::string, std::allocator<std::string> > > str_vec(12);
  utils::small_vector<std::string, utils::allocinfo_allocator<std::string, std::allocator<std::string> > > str_vec2;

  utils::small_vector<int> int_vec;
  
  std::cout << "allocated: " << utils::allocinfo().allocated() << std::endl;
  std::cout << "sizeof small_vector<int>: " << sizeof(utils::small_vector<int>) << std::endl;
  std::cout << "sizeof small_vector<string>: " << sizeof(utils::small_vector<std::string>) << std::endl;

  std::cout << "begin: " << int_vec.begin() << " begin-1: " << (int_vec.begin() - 1) << " end: " << int_vec.end() << std::endl;
  
  str_vec[0] = "January";
  str_vec[1] = "February";
  str_vec[2] = "March";
  str_vec[3] = "April";
  str_vec[4] = "May";
  str_vec[5] = "June";
  str_vec[6] = "July";
  str_vec[7] = "August";
  str_vec[8] = "September";
  str_vec[9] = "October";
  str_vec[10] = "November";
  str_vec[11] = "December";
  
  for (int i = 0; i < str_vec.size(); ++ i)
    std::cout << "i = " << i << " " << str_vec[i] << std::endl;
  for (utils::small_vector<std::string, utils::allocinfo_allocator<std::string, std::allocator<std::string> > >::const_iterator
	 iter = str_vec.begin(); iter != str_vec.end(); ++ iter)
    std::cout << "i = " << (iter - str_vec.begin()) << " " << *iter << std::endl;
  
  std::cout << "back: " << str_vec.back() << std::endl;

  str_vec2 = str_vec;
  for (int i = 0; i < str_vec2.size(); ++ i) {
    std::cout << "i = " << i << " " << str_vec2[i] << std::endl;
    str_vec2[i] = "";
  }
  str_vec2 = str_vec;
  for (int i = 0; i < str_vec2.size(); ++ i)
    std::cout << "i = " << i << " " << str_vec2[i] << std::endl;

  std::vector<std::string> vec;
  vec.push_back("1");
  vec.push_back("2");
  vec.push_back("3");
  
  str_vec.assign(vec.begin(), vec.end());
  
  for (int i = 0; i < str_vec.size(); ++ i)
    std::cout << "i = " << i << " " << str_vec[i] << std::endl;
  std::cout << "allocated: " << utils::allocinfo().allocated() << std::endl;
  
  str_vec.resize(500);
  std::cout << "resize: " << str_vec.size() << std::endl;
  std::cout << "allocated: " << utils::allocinfo().allocated() << std::endl;

  str_vec.clear();
  str_vec2.clear();
  std::cout << "allocated (clear): " << utils::allocinfo().allocated() << std::endl;

  std::vector<int> intvec;
  for (int i = 0; i < 1024 * 16; ++ i)
    intvec.push_back(i);

  int_vec.assign(intvec.begin(), intvec.end());
  
  std::vector<int>::const_iterator iter = intvec.begin();
  utils::small_vector<int>::const_iterator siter = int_vec.begin();
  for (int i = 0; i < 1024 * 16; ++ i, ++ iter, ++ siter)
    if (*iter != *siter)
      std::cout << "differ: " << i << " " << *iter << " " << *siter << std::endl;
  
  std::cout << "size: " << int_vec.size() << std::endl;
  
  {
    size_t size = int_vec.size();
    for (size_t i = 0; i != size; ++ i) {
      const size_t pos = random() % (size - i);
      intvec.erase(intvec.begin() + pos);
      int_vec.erase(int_vec.begin() + pos);

      if (intvec.size() != int_vec.size())
	std::cout << "size differ: " << intvec.size() << ' ' << int_vec.size() << std::endl;

      if (! std::equal(int_vec.begin(), int_vec.end(), intvec.begin()))
	std::cerr << "content differ" << std::endl;
      
      if (size - i - 1 != int_vec.size())
	std::cout << "differ: " << (size - i - 1) << ' ' << int_vec.size() << std::endl;
    }
  }
  
  {
    std::cout << "size: " << intvec.size() << ' ' << int_vec.size() << std::endl;
    
    const size_t size = 1024 * 16;
    
    for (size_t i = 0; i != size; ++ i) {
      const int value = random();
      const int pos   = random() % (i + 1);
      
      intvec.insert(intvec.begin() + pos, value);
      int_vec.insert(int_vec.begin() + pos, value);
      
      if (intvec.size() != int_vec.size())
	std::cout << "size differ: " << intvec.size() << ' ' << int_vec.size() << std::endl;
      
      if (! std::equal(int_vec.begin(), int_vec.end(), intvec.begin()))
	std::cerr << "content differ" << std::endl;
    }
  }
  
  utils::small_vector<int> intvec2(12, 666);
  utils::small_vector<int> intvec3(12);

  intvec3.clear();
  intvec2.swap(intvec3);
  
  std::cerr << "size: " << intvec2.size() << " " << intvec3.size() << std::endl;

  utils::small_vector<int> intvec4(2, 666);
  utils::small_vector<int> intvec5(2);

  intvec5.clear();
  intvec4.swap(intvec5);
  
  std::cerr << "size: " << intvec4.size() << " " << intvec5.size() << std::endl;

  std::cerr << "pointers: " << intvec4.begin() << " " << &(intvec4)  << std::endl;
  std::cerr << "pointers: " << intvec5.begin() << " " << &(intvec5)  << std::endl;
}
