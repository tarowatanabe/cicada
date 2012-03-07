
#include <utils/filesystem.hpp>

#include <iostream>

int main(int argc, char** argv)
{
  typedef boost::filesystem::path path_type;
  
  if (argc <= 1) return 0;
  
  std::cout << path_type(argv[1]) << std::endl;
  std::cout << path_type(argv[1]).string() << std::endl;
  std::cout << path_type(argv[1]).filename() << std::endl;
  std::cout << path_type(argv[1]).filename().string() << std::endl;
}
