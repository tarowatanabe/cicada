//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "parameter.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Parameter parameter_type;
  
  parameter_type parameter("ngram:file=tmptmp,order=good,name=ngram,name=bad,good=\"bad morning\"");
  
  std::cout << parameter.name();
  for (parameter_type::const_iterator iter = parameter.begin(); iter != parameter.end(); ++ iter)
    std::cout << ' ' << iter->first << ':' << iter->second;
  std::cout << std::endl;

  std::cout << "generated: " << parameter << std::endl;

  parameter.erase("name");
  
  std::cout << "erased: " << parameter << std::endl;
  
  std::string line;
  while (std::getline(std::cin, line)) {
    parameter_type parameter(line);
    
    std::cout << parameter.name();
    for (parameter_type::const_iterator iter = parameter.begin(); iter != parameter.end(); ++ iter)
      std::cout << ' ' << iter->first << ':' << iter->second;
    std::cout << std::endl;
    
    std::cout << "generated: " << parameter << std::endl;
  }

}
