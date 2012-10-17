//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "hypergraph.hpp"

#ifdef HAVE_MSGPACK_HPP
#include <msgpack.hpp>
#include <cicada/msgpack/hypergraph.hpp>
#endif

int main(int argc, char** argv)
{
  cicada::HyperGraph graph;
  
  std::cin >> graph;
  std::cout << graph << std::endl;
}
