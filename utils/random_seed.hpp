// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__RANDOM_SEED__HPP__
#define __UTILS__RANDOM_SEED__HPP__ 1

#include <ctime>
#include <fstream>

namespace utils
{
  inline
  unsigned int random_seed()
  {
    unsigned int seed;
    
    std::ifstream udev("/dev/urandom");
    
    if (udev)
      udev.read((char*) &seed, sizeof(unsigned int));
    
    if (udev.fail() || ! udev) {
      std::ifstream dev("/dev/random");
      
      if (dev)
	dev.read((char*) &seed, sizeof(unsigned int));
      
      if (dev.fail() || ! dev)
	seed = std::time(0);
    }
    return seed;
  }
};

#endif
