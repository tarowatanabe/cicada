//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "rule.hpp"

#include <cicada/msgpack/rule.hpp>

#include "msgpack_main_impl.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Rule rule_type;
 
  std::cout << "rule: " << rule_type("good boy") << std::endl;
  msgpack_test(rule_type("good boy"));

  std::cout << "rule: " << rule_type("[s] ||| good boy") << std::endl;
  std::cout << "rhs: " << rule_type("[s] ||| good boy").rhs.size() << std::endl;
  msgpack_test(rule_type("[s] ||| good boy"));
  
  std::cout << "rule: " << rule_type("[s] ||| good boy [x,1]") << std::endl;
  std::cout << "rhs: " << rule_type("[s] ||| good boy [x,1]").rhs.size() << std::endl;
  msgpack_test(rule_type("[s] ||| good boy [x,1]"));
}
