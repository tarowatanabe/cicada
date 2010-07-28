
#include <iostream>

#include "lattice.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Lattice lattice_type;
  
  lattice_type lattice("((('ein\\'\"en',1.0,1),),(('wettbewerbsbedingten',0.5,2),('wettbewerbs',0.25,1), ('wettbewerb',0.25, 1),),(('bedingten',1.0,1),),(('preissturz',0.5,2), ('preis',0.5,1),),(('sturz',1.0,1),),)");

  std::cout << "lattice size: " << lattice.size() << std::endl;

  std::cout << lattice << std::endl;


  lattice_type input;
  while (std::cin >> input)
    std::cout << input << std::endl;
}
