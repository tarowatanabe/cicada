//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cluster.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " cluster-file" << std::endl;
    return 1;
  }
  
  cicada::Cluster cluster(argv[1]);
  
  std::string word;
  while (std::cin >> word)
    std::cout << "word: " << word << " cluster: " << cluster[word] << std::endl;
  
}
