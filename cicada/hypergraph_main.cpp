#include <iostream>

#include "hypergraph.hpp"

int main(int argc, char** argv)
{
  cicada::HyperGraph graph;
  
  std::cin >> graph;
  std::cout << graph << std::endl;
}
