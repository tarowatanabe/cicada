//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <time.h>
#include <unistd.h>

#include <iostream>
#include <iterator>
#include <time.h>
#include <unistd.h>

#include <iostream>

#include "weight_vector.hpp"

#include "utils/random_seed.hpp"

typedef cicada::WeightVector<double> weights_double_type;
typedef cicada::WeightVector<float>  weights_float_type;

int main(int argc, char** argv)
{
  srandom(utils::random_seed());
  srand(utils::random_seed());

  weights_double_type weights1;
  weights_double_type weights2;
  
  weights1["ngram1"] = 1.0;
  weights1["ngram2"] = 2.0;
  weights1["ngram3"] = 3.0;
  weights1["ngram4"] = 4.0;

  weights2["ngram2"] = 2.0;
  weights2["ngram3"] = 1.0;
  
  std::cout << "weights1" << std::endl
	    << weights1
	    << std::endl;

  std::cout << "weights2" << std::endl
	    << weights2
	    << std::endl;

  weights1 *= 2.5;
  
  std::cout << "weights1" << std::endl
	    << weights1
	    << std::endl;
  
  weights1 += weights2;
  
  std::cout << "weights1" << std::endl
	    << weights1
	    << std::endl;
  
}
