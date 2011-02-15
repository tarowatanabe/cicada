
#include <iostream>

#include "utils/lexical_cast.hpp"
#include "utils/resource.hpp"

#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv)
{
  std::cout << utils::lexical_cast<bool>("yes") << std::endl;
  std::cout << utils::lexical_cast<bool>("Yes ") << std::endl;
  std::cout << utils::lexical_cast<bool>(" true") << std::endl;
  std::cout << utils::lexical_cast<bool>(" tRue ") << std::endl;
  std::cout << utils::lexical_cast<bool>(" 1 ") << std::endl;
  std::cout << utils::lexical_cast<bool>(" 2 ") << std::endl;

  std::cout << utils::lexical_cast<bool>("no") << std::endl;
  std::cout << utils::lexical_cast<bool>("No") << std::endl;
  std::cout << utils::lexical_cast<bool>("nil ") << std::endl;
  std::cout << utils::lexical_cast<bool>("NIL ") << std::endl;
  std::cout << utils::lexical_cast<bool>("false") << std::endl;
  std::cout << utils::lexical_cast<bool>("False") << std::endl;
  std::cout << utils::lexical_cast<bool>(" -1 ") << std::endl;
  std::cout << utils::lexical_cast<bool>(" 0 ") << std::endl;
  std::cout << utils::lexical_cast<bool>(" ") << std::endl;

  std::cout << "integers" << std::endl;
  std::cout << utils::lexical_cast<int>("567") << std::endl;
  std::cout << utils::lexical_cast<size_t>("567") << std::endl;
  std::cout << utils::lexical_cast<double>(" inf ") << std::endl;
  std::cout << utils::lexical_cast<double>("nan") << std::endl;

  std::cout << "generator" << std::endl;
  std::cout << utils::lexical_cast<std::string>(6.77) << std::endl;
  std::cout << utils::lexical_cast<std::string>(6) << std::endl;
  std::cout << utils::lexical_cast<std::string>(1e-60) << std::endl;
  std::cout << utils::lexical_cast<std::string>(std::numeric_limits<double>::infinity()) << std::endl;
  std::cout << utils::lexical_cast<std::string>(utils::lexical_cast<double>("nan")) << std::endl;

  ::srandom(time(0));
  
  {
    utils::resource boost_start;
    
    for (int i = 0; i != 1024 * 64; ++ i)
      boost::lexical_cast<int>(boost::lexical_cast<std::string>(random()));

    utils::resource boost_end;

    utils::resource utils_start;

    for (int i = 0; i != 1024 * 64; ++ i)
      utils::lexical_cast<int>(utils::lexical_cast<std::string>(random()));

    utils::resource utils_end;

    utils::resource std_start;
    char buffer[256];
    int integer;
    
    for (int i = 0; i != 1024 * 64; ++ i) {
      sprintf(buffer, "%d", static_cast<int>(random()));
      sscanf(buffer, "%d", &integer);
    }
    
    utils::resource std_end;


    std::cout << "boost cpu time: " << (boost_end.cpu_time() - boost_start.cpu_time())
	      << " user time: " << (boost_end.user_time() - boost_start.user_time())
	      << std::endl;
    std::cout << "utils cpu time: " << (utils_end.cpu_time() - utils_start.cpu_time())
	      << " user time: " << (utils_end.user_time() - utils_start.user_time())
	      << std::endl;
    std::cout << "std cpu time: " << (std_end.cpu_time() - std_start.cpu_time())
	      << " user time: " << (std_end.user_time() - std_start.user_time())
	      << std::endl;
  }

  {
    utils::resource boost_start;
    
    for (int i = 0; i != 1024 * 64; ++ i)
      boost::lexical_cast<double>(boost::lexical_cast<std::string>(double(random()) / random()));
    
    utils::resource boost_end;
    
    utils::resource utils_start;
    
    for (int i = 0; i != 1024 * 64; ++ i)
      utils::lexical_cast<double>(utils::lexical_cast<std::string>(double(random()) / random()));
    
    utils::resource utils_end;
    
    utils::resource std_start;
    char buffer[256];
    double integer;
    
    for (int i = 0; i != 1024 * 64; ++ i) {
      sprintf(buffer, "%g", double(random()) / random());
      std::string tmp(buffer);
      sscanf(tmp.c_str(), "%lg", &integer);
    }
    
    utils::resource std_end;

    
    std::cout << "boost cpu time: " << (boost_end.cpu_time() - boost_start.cpu_time())
	      << " user time: " << (boost_end.user_time() - boost_start.user_time())
	      << std::endl;
    std::cout << "utils cpu time: " << (utils_end.cpu_time() - utils_start.cpu_time())
	      << " user time: " << (utils_end.user_time() - utils_start.user_time())
	      << std::endl;
    std::cout << "std cpu time: " << (std_end.cpu_time() - std_start.cpu_time())
	      << " user time: " << (std_end.user_time() - std_start.user_time())
	      << std::endl;
  }

}
