//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "attribute_vector.hpp"
#include "attribute_vector_codec.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"

typedef cicada::AttributeVector attribute_set_type;

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

void check_codec(const attribute_set_type& attributes)
{
  std::vector<char> buffer;
  
  cicada::attribute_vector_encode(attributes, std::back_inserter(buffer));
  
  attribute_set_type decoded;
  cicada::attribute_vector_decode(buffer.begin(), buffer.end(), decoded);

  if (attributes != decoded)
    std::cerr << "differ: " << attributes << " codec: "  << decoded << std::endl;
}

int main(int argc, char** argv)
{
  typedef cicada::AttributeVector attribute_set_type;

  srandom(utils::random_seed());
  
  std::cout << "size: " << sizeof(attribute_set_type) << " value: " << sizeof(attribute_set_type::value_type) << std::endl;
  
  attribute_set_type attr1("{\"good\":1,\"bad\":4.5,\"bad2\":1e-5, \"neutral\":\"bi\\u0020g\"}");

  std::cout << "attr1 size: " << attr1.size() << std::endl;
  
  find_attrs(attr1, "good");
  find_attrs(attr1, "bad");
  find_attrs(attr1, "bad2");
  find_attrs(attr1, "neutral");

  std::cout << "sizeof attr-value: " << sizeof(attribute_set_type::data_type) << std::endl;
  
  std::cout << "attr1: " << attr1 << std::endl;

  attr1["bad"] = attribute_set_type::int_type(5);
  std::cout << "attr1: " << attr1 << std::endl;
  attr1.erase("bad");
  std::cout << "attr1: " << attr1 << std::endl;
  attr1["bad"] = 4.5;
  std::cout << "attr1: " << attr1 << std::endl;

  attribute_set_type attributes;
  
  for (int iter = 0; iter != 32; ++ iter) {
    for (int i = 0; i != 16; ++ i) {
      std::string attr = "double:" + utils::lexical_cast<std::string>(random());
      attributes[attr] = (1.0 * random()) / random();
    }

    check_codec(attributes);

    for (int i = 0; i != 16; ++ i) {
      std::string attr = "int:" + utils::lexical_cast<std::string>(random());
      attributes[attr] = attribute_set_type::int_type(random());
    }

    check_codec(attributes);
    
    for (int i = 0; i != 16; ++ i) {
      std::string attr = "string:" + utils::lexical_cast<std::string>(random());
      attributes[attr] = (utils::lexical_cast<std::string>(random())
			  + utils::lexical_cast<std::string>(random())
			  + utils::lexical_cast<std::string>(random())
			  + utils::lexical_cast<std::string>(random()));
    }

    check_codec(attributes);
  }
}
