
#include <iostream>

#include "parameter.hpp"

int main(int argc, char** argv)
{
  typedef cicada::Parameter parameter_type;
  
  parameter_type parameter("ngram,file=tmptmp,order=good,name=ngram");
  
  std::cout << parameter.name();
  for (parameter_type::const_iterator iter = parameter.begin(); iter != parameter.end(); ++ iter)
    std::cout << ' ' << iter->first << ':' << iter->second;
  std::cout << std::endl;
}
