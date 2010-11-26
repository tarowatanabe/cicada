//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "attribute_vector.hpp"

void find_attrs(const cicada::AttributeVector& attrs, const std::string& attr)
{
  typedef cicada::AttributeVector attribute_set_type;
  
  attribute_set_type::const_iterator iter = attrs.find(attr);
  if (iter != attrs.end())
    std::cout << "found: " << iter->first << " : " << iter->second << std::endl;
  
  {
    attribute_set_type::const_iterator iter = attrs.find_int(attr);
    if (iter != attrs.end())
      std::cout << "int: " << iter->first << " : " << iter->second << std::endl;
  }

  {
    attribute_set_type::const_iterator iter = attrs.find_float(attr);
    if (iter != attrs.end())
      std::cout << "float: " << iter->first << " : " << iter->second << std::endl;
  }
  
  {
    attribute_set_type::const_iterator iter = attrs.find_string(attr);
    if (iter != attrs.end())
      std::cout << "string: " << iter->first << " : " << iter->second << std::endl;
  }
}

int main(int argc, char** argv)
{
  typedef cicada::AttributeVector attribute_set_type;
  
  attribute_set_type attr1("{\"good\":1,\"bad\":4.5,\"bad2\":1e-5, \"neutral\":\"big\"}");

  std::cout << "attr1 size: " << attr1.size() << std::endl;
  
  find_attrs(attr1, "good");
  find_attrs(attr1, "bad");
  find_attrs(attr1, "bad2");
  find_attrs(attr1, "neutral");
  
  std::cout << "attr1: " << attr1 << std::endl;
  
}
