//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include <iostream>

#include <boost/lexical_cast.hpp>

#include "symbol.hpp"

void process(const cicada::Symbol& x, int index)
{
  std::cout << x << ' ' << index << std::endl
	    << "is a non-terminal? " << x.is_non_terminal() << std::endl
	    << "non-terminal? " << x.non_terminal() << std::endl
	    << "non-terminal index? " << x.non_terminal_index() << std::endl
	    << "non-terminal index? " << x.non_terminal(index) << std::endl
	    << "non-terminal strip? " << x.non_terminal_strip() << std::endl
	    << "annotation 1 true?" << x.annotation(1, true) << std::endl
	    << "annotation 1 false?" << x.annotation(1, false) << std::endl
	    << "annotation 2 true?" << x.annotation(2, true) << std::endl
	    << "annotation 2 false?" << x.annotation(2, false) << std::endl
	    << "coarse 1? " << x.coarse(1) << std::endl
	    << "coarse 2? " << x.coarse(2) << std::endl
	    << "coarse? " << x.coarse() << std::endl;
}

int main(int argc, char** argv)
{
  typedef cicada::Symbol symbol_type;

  process("good", 1);
  process("[,]", 2);
  process("[x]", 1);
  process("[x]", 2);
  process("[x,3]", 2);
  process("[x,3]", 3);
  process("[x,4]", 3);
  process("[x,100]", 5000);
  process("[x,5g]", 4);
  process("[BAD,4]", 3);
  process("[GOOD,5]", 3);

  process("[BAD@1]", 3);
  process("[BAD@1,5]", 3);
  process("[BAD@2]", 3);
  process("[BAD@3]", 3);
  process("[BAD@3,2]", 3);
  process("[BAD@1,4]", 3);
  
  std::cerr << "pos: " << symbol_type("Good/[ADJ]").pos() << " terminal: " << symbol_type("Good/[ADJ]").terminal() << std::endl;
  std::cerr << "pos: " << symbol_type("Good|[ADJ]").pos() << " terminal: " << symbol_type("Good|[ADJ]").terminal() << std::endl;
  std::cerr << "pos: " << symbol_type("Good\\[ADJ]").pos() << " terminal: " << symbol_type("Good\\[ADJ]").terminal() << std::endl;

  srandom(time(0) * getpid());

  for (int i = 0; i != 1024 * 4; ++ i) {
    const std::string rnd = boost::lexical_cast<std::string>(random());
    
    symbol_type symbol(rnd);
  }
}
