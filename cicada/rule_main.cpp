
#include <iostream>

#include "rule.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Rule rule_type;
  

  std::cout << "rule: " << rule_type("good boy ||| bad boy ||| feature1=1") << std::endl;
  std::cout << "rule: " << rule_type("good boy ||| bad boy") << std::endl;
  
  std::cout << "rule: " << rule_type("good boy ||| ||| feature1=1") << std::endl;
  std::cout << "rule: " << rule_type("good boy ||| ") << std::endl;
  
  std::cout << "rule: " << rule_type("[s] ||| good boy ||| bad boy ||| feature1=1") << std::endl;
  std::cout << "rule: " << rule_type("[s] ||| good boy ||| bad boy") << std::endl;
  
  std::cout << "rule: " << rule_type("[s] ||| good boy ||| ||| feature1=1") << std::endl;
  std::cout << "rule: " << rule_type("[s] ||| good boy ||| ") << std::endl;
  
  std::cout << "rule: " << rule_type("[s] ||| good boy [x,1] ||| bad boy [x,1] ||| feature1=1") << std::endl;
  std::cout << "rule: " << rule_type("[s] ||| good boy [x,1] ||| bad boy [x,1]") << std::endl;
  
  std::cout << "rule: " << rule_type("[s] ||| good boy [x,1] ||| ||| feature1=1") << std::endl;
  std::cout << "rule: " << rule_type("[s] ||| good boy [x,1] ||| ") << std::endl;
}
