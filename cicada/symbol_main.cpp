//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "symbol.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Symbol symbol_type;
  
  std::cout << symbol_type("good") << " is a non-terminal? " << symbol_type("good").is_non_terminal() << std::endl;
  
  std::cout << symbol_type("[x]") << " is a non-terminal? " << symbol_type("[x]").is_non_terminal() << std::endl;

  std::cout << symbol_type("[x]") << " non-terminal-index? " << symbol_type("[x]").non_terminal_index() << std::endl;

  std::cout << symbol_type("[x,3]") << " non-terminal-index? " << symbol_type("[x,3]").non_terminal_index() << std::endl;
  std::cout << symbol_type("[x,3]") << " non-terminal-index? " << symbol_type("[x,3]").non_terminal_index() << std::endl;
  std::cout << symbol_type("[x,5]") << " non-terminal-index? " << symbol_type("[x,5]").non_terminal_index() << std::endl;
  std::cout << symbol_type("[x,3]") << " non-terminal-index? " << symbol_type("[x,3]").non_terminal(4) << std::endl;
  std::cout << symbol_type("[x,3]") << " non-terminal-index? " << symbol_type("[x,3]").non_terminal(5) << std::endl;
  std::cout << symbol_type("[x,3]") << " non-terminal-index? " << symbol_type("[x,3]").non_terminal(0) << std::endl;
  std::cout << symbol_type("[x,3]") << " non-terminal-strip? " << symbol_type("[x,3]").non_terminal(5).non_terminal_strip() << std::endl;
  
  std::cout << symbol_type("[x,3]") << " non-terminal? " << symbol_type("[x,3]").non_terminal() << std::endl;
  std::cout << symbol_type("[x]") << " non-terminal? " << symbol_type("[x]").non_terminal() << std::endl;

  std::cout << symbol_type("[x,3g]") << " non-terminal-index? " << symbol_type("[x,3g]").non_terminal_index() << std::endl;
  std::cout << symbol_type("[x,3g]") << " non-terminal-index? " << symbol_type("[x,3g]").non_terminal() << std::endl;

}
