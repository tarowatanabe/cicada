//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <time.h>
#include <unistd.h>

#include <iostream>

#include "feature_vector.hpp"
#include "feature_vector_compact.hpp"
#include "dot_product.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"

typedef cicada::FeatureVector<double> feature_set_type;
typedef cicada::FeatureVectorCompact  feature_compact_type;

void check_compact(const feature_set_type& features, const feature_compact_type& feats)
{
  feature_set_type::const_iterator iter = features.begin();
  feature_set_type::const_iterator iter_end = features.end();
  
  feature_compact_type::const_iterator citer = feats.begin();
  feature_compact_type::const_iterator citer_end = feats.end();

  std::cerr << "size: " << features.size() * sizeof(feature_set_type::value_type)
	    << " compressed: " << feats.size_compressed()
	    << std::endl;

  while (iter != iter_end && citer != citer_end) {
    if (iter->first < citer->first) {
      std::cerr << "differ for the original vector!" << std::endl;
      ++ iter;
    } else if (citer->first < iter->first) {
      std::cerr << "differ for the compressed vector!" << std::endl;
      ++ citer;
    } else {
      if (iter->second != citer->second)
	std::cerr << "differ for value" << std::endl;
      
      ++ iter;
      ++ citer;
    }
  }
  
  for (/**/; iter != iter_end; ++ iter)  
    std::cerr << "differ for the original vector!" << std::endl;
  
  for (/**/; citer != citer_end; ++ citer)
    std::cerr << "differ for the compressed vector!" << std::endl;
}

void check_compact(const feature_set_type& features)
{
  feature_compact_type feats(features);
  feature_compact_type feats1(features.begin(), features.end());
  feature_compact_type feats2(features.begin(), features.end(), true);
  
  check_compact(features, feats);
  check_compact(features, feats1);
  check_compact(features, feats2);
}

int main(int argc, char** argv)
{
  srandom(utils::random_seed());

  feature_set_type features1;
  feature_set_type features2;
  feature_set_type features3;
  
  features1["ngram1"] = 1.0;
  features1["ngram2"] = 2.0;
  features1["ngram3"] = 3.0;

  features2["ngram2"] = 2.0;
  features2["ngram3"] = 1.0;
  
  features3 = features1;

  check_compact(features1);
  check_compact(features2);
  check_compact(features3);
  
  for (int i = 0; i != 8; ++ i) {
    std::string feat = "bad:" + utils::lexical_cast<std::string>(i);
    
    features1[feat] = i;
  }

  check_compact(features1);
  check_compact(features2);
  check_compact(features3);
  
  std::cout << "feature1" << std::endl;
  std::cout << features1;

  std::cout << "feature2" << std::endl;
  std::cout << features2;

  std::cout << "feature3" << std::endl;
  std::cout << features3;

  std::cout << "feature1 + feature2" << std::endl;
  std::cout << features1 + features2;
  
  std::cout << "feature1 * feature2" << std::endl;
  std::cout << features1 * features2;

  std::cout << "dot_product(feature1, feature2)" << std::endl;
  std::cout << dot_product(features1, features2) << std::endl;

  check_compact(features1);
  check_compact(features2);
  check_compact(features3);
  
  std::cout << "feature1 - feature3" << std::endl;
  std::cout << features1 - features3;

  std::cout << "feature3 - feature1" << std::endl;
  std::cout << features3 - features1;

  features1.erase_prefix(std::string("bad"));
  std::cout << "erased prefix for bad" << std::endl;
  std::cout << features1;

  check_compact(features1);
  check_compact(features2);
  check_compact(features3);

  std::cout << "inserted many" << std::endl;
  for (int i = 0; i != 32; ++ i) {
    std::string feat = "good:" + utils::lexical_cast<std::string>(i);
    
    features1[feat] = i;
  }
  
  check_compact(features1);
  check_compact(features2);
  check_compact(features3);

  std::cout << "feature1" << std::endl;
  std::cout << features1;
  
  std::cout << "feature1 + feature2" << std::endl;
  std::cout << features1 + features2;
  
  std::cout << "feature1 * feature2" << std::endl;
  std::cout << features1 * features2;
  
  std::cout << "dot_product(feature1, feature2)" << std::endl;
  std::cout << dot_product(features1, features2) << std::endl;
  
  std::cout << "feature1 - feature3" << std::endl;
  std::cout << features1 - features3;
  
  std::cout << "feature3 - feature1" << std::endl;
  std::cout << features3 - features1;
  
  check_compact(features1);
  check_compact(features2);
  check_compact(features3);
  
  for (int iter = 0; iter != 64; ++ iter) {

    for (int i = 0; i != 32; ++ i) {
      std::string feat = "double:" + utils::lexical_cast<std::string>(random());
      
      features1[feat] = (1.0 * random()) / random();
    }

    for (int i = 0; i != 32; ++ i) {
      std::string feat = "double:" + utils::lexical_cast<std::string>(random());
      
      features1[feat] = (- 1.0 * random()) / random();
    }

    for (int i = 0; i != 32; ++ i) {
      std::string feat = "int:" + utils::lexical_cast<std::string>(random());
      
      features1[feat] = random() % 16;
    }

    for (int i = 0; i != 32; ++ i) {
      std::string feat = "int:" + utils::lexical_cast<std::string>(random());
      
      features1[feat] = - (int(random()) % 16);
    }
    
    check_compact(features1);
  }
}
