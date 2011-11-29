#include <iostream>

#include "utils/random_seed.hpp"

int main(int argc, char** argv)
{
  std::cout << utils::random_seed() << std::endl;
}
