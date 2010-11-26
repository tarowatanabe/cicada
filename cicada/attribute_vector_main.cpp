//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "attribute_vector.hpp"

int main(int argc, char** argv)
{
  typedef cicada::AttributeVector attribute_set_type;
  
  attribute_set_type attr1("{\"good\":1,\"bad\":4.5,\"bad2\":1e-5, \"neutral\":\"big\"}");

  std::cout << "attr1 size: " << attr1.size() << std::endl;
  
  attribute_set_type::const_iterator giter = attr1.find("good");
  if (giter != attr1.end())
    std::cout << giter->first << " : " << giter->second << std::endl;
  
  attribute_set_type::const_iterator biter = attr1.find("bad");
  if (biter != attr1.end())
    std::cout << biter->first << " : " << biter->second << std::endl;

  attribute_set_type::const_iterator niter = attr1.find("neutral");
  if (niter != attr1.end())
    std::cout << niter->first << " : " << niter->second << std::endl;
  
  std::cout << "attr1: " << attr1 << std::endl;
  
}
