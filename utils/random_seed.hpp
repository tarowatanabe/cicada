// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__RANDOM_SEED__HPP__
#define __UTILS__RANDOM_SEED__HPP__ 1

#include <fstream>
#include <ctime>

namespace utils
{
  inline
  unsigned long random_seed()
  {
    unsigned long seed;
    
    std::ifstream udev("/dev/urandom");
    
    if (udev)
      udev.read((char*) &seed, sizeof(unsigned long));
    
    if (udev.fail() || ! udev) {
      std::ifstream dev("/dev/random");
      
      if (dev)
	dev.read((char*) &seed, sizeof(unsigned long));
      
      if (dev.fail() || ! dev)
	seed = std::time(0);
    }
    return seed;
  }
};

#endif
