#include <iostream>
#include <string>

#include "hashxx.hpp"

int main(int argc, char** argv)
{
  std::string line;

  utils::hashxx<uint64_t> hasher64;
  utils::hashxx<uint32_t> hasher32;
  
  while (std::getline(std::cin, line))
    std::cout << hasher64(line.begin(), line.end(), 0)
	      << ' '
	      << hasher32(line.begin(), line.end(), 0)
	      << std::endl;
}
