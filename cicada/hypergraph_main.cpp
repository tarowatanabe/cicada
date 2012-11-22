//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "hypergraph.hpp"

#include <cicada/msgpack/hypergraph.hpp>

#include "msgpack_main_impl.hpp"

int main(int argc, char** argv)
{
  try {
    cicada::HyperGraph graph;
    
    std::cin >> graph;
    std::cout << graph << std::endl;
    
    msgpack_test(graph);
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
}
