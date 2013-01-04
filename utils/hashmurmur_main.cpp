#include <iostream>
#include <string>

#include "hashmurmur.hpp"

int main(int argc, char** argv)
{
  std::string line;

  utils::hashmurmur<uint32_t> hasher;
  
  while (std::getline(std::cin, line))
    std::cout << hasher(line.begin(), line.end(), 0) << std::endl;
}
