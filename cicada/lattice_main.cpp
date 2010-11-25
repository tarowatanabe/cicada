//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "lattice.hpp"
#include "remove_epsilon.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Lattice lattice_type;
  
  lattice_type lattice("((('ein\\'\"en',1.0,1),),(('wettbewerbsbedingten',0.5,2),('wettbewerbs',0.25,1), ('wettbewerb',0.25, 1),),(('bedingten',1.0,1),),(('preissturz',0.5,2), ('preis',0.5,1),),(('sturz',1.0,1),),)");

  std::cout << "lattice size: " << lattice.size() << std::endl;

  std::cout << lattice << std::endl;

  lattice_type input;
  lattice_type removed;
  while (std::cin >> input) {
    std::cout << "shortest: " << input.shortest_distance()
	      << " longest: " << input.longest_distance()
	      << " nodes: " << input.size()
	      << std::endl;
    
    std::cout << input << std::endl;
    
    cicada::remove_epsilon(input, removed);

    std::cout << "shortest: " << removed.shortest_distance()
	      << " longest: " << removed.longest_distance()
	      << " nodes: " << removed.size()
	      << std::endl;
    
    std::cout << removed << std::endl;
  }
}
