// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREEBANK__HPP__
#define __CICADA__TREEBANK__HPP__ 1

#include <iostream>
#include <cicada/hypergraph.hpp>

namespace cicada
{
  std::ostream& treebank(std::ostream& os, const HyperGraph& hypergraph);
};

#endif
