//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <time.h>
#include <unistd.h>

#include <iostream>

#include "feature_vector.hpp"
#include "utils/lexical_cast.hpp"

typedef cicada::FeatureVector<double> feature_set_type;

int main(int argc, char** argv)
{
  feature_set_type features1;
  feature_set_type features2;
  feature_set_type features3;
  
  features1["ngram1"] = 1.0;
  features1["ngram2"] = 2.0;
  features1["ngram3"] = 3.0;

  features2["ngram2"] = 2.0;
  features2["ngram3"] = 1.0;
  
  features3 = features1;
  
  for (int i = 0; i != 8; ++ i) {
    std::string feat = "bad:" + utils::lexical_cast<std::string>(i);
    
    features1[feat] = i;
  }
  
  std::cout << "feature1" << std::endl;
  std::cout << features1;

  std::cout << "feature2" << std::endl;
  std::cout << features2;

  std::cout << "feature1 + feature2" << std::endl;
  std::cout << features1 + features2;
  
  std::cout << "feature1 * feature2" << std::endl;
  std::cout << features1 * features2;
  
  std::cout << "feature1 - feature3" << std::endl;
  std::cout << features1 - features3;

  std::cout << "feature3 - feature1" << std::endl;
  std::cout << features3 - features1;

  features1.erase_prefix(std::string("bad"));
  std::cout << "erased prefix for bad" << std::endl;
  std::cout << features1;

  std::cout << "inserted many" << std::endl;
  for (int i = 0; i != 32; ++ i) {
    std::string feat = "good:" + utils::lexical_cast<std::string>(i);
    
    features1[feat] = i;
  }
  
  std::cout << "feature1" << std::endl;
  std::cout << features1;
  
  std::cout << "feature1 + feature2" << std::endl;
  std::cout << features1 + features2;
  
  std::cout << "feature1 * feature2" << std::endl;
  std::cout << features1 * features2;
  
  std::cout << "feature1 - feature3" << std::endl;
  std::cout << features1 - features3;
  
  std::cout << "feature3 - feature1" << std::endl;
  std::cout << features3 - features1;
}
