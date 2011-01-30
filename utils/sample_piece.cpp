
#include "piece.hpp"
#include "space_separator.hpp"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

int main(int argc, char** argv)
{
  std::string str1("good thing");
  const char* str2("bad thing");
  
  utils::piece pie1(str1);
  utils::piece pie2(str2);
  
  std::cout << "piece1: " << pie1 << " base: " << (void*) pie1.c_str() << std::endl
	    << "piece2: " << pie2 << " base: " << (void*) pie2.c_str() << std::endl;
  
  std::cout << "piece1 substr: " << pie1.substr(5) << " base: " << (void*) pie1.substr(5).c_str() << std::endl
	    << "piece2 substr: " << pie2.substr(4) << " base: " << (void*) pie2.substr(4).c_str() << std::endl
	    << "equal? " << (pie1.substr(5) == pie2.substr(4)) << std::endl;
  
  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
  
  tokenizer_type tokens1(pie1);
  tokenizer_type tokens2(pie2);
  
  for (tokenizer_type::iterator iter = tokens1.begin(); iter != tokens1.end(); ++ iter)
    std::cout << "token1: " << *iter << " base: " << (void*) (*iter).c_str() << " cast: " << static_cast<std::string>(*iter) << std::endl;
  for (tokenizer_type::iterator iter = tokens2.begin(); iter != tokens2.end(); ++ iter)
    std::cout << "token2: " << *iter << " base: " << (void*) (*iter).c_str() << " cast: " << static_cast<std::string>(*iter) << std::endl;

  const char* str3 = "Good";
  const char* str4 = "GOOD";
  const char* str5 = "Good-bye";
  
  std::cout << "less: " << (utils::piece(str3) < str4) << std::endl
	    << "greater: " << (utils::piece(str3) > str4) << std::endl
	    << "equal: " << (utils::piece(str3) == str4) << std::endl;
  
  std::cout << "less: " << (utils::ipiece(str3) < str4) << std::endl
	    << "greater: " << (utils::ipiece(str3) > str4) << std::endl
	    << "equal: " << (utils::ipiece(str3) == str4) << std::endl;
  
  std::cout << "less: " << (utils::ipiece(str3) < str5) << std::endl
	    << "greater: " << (utils::ipiece(str3) > str5) << std::endl
	    << "equal: " << (utils::ipiece(str3) == str5) << std::endl;
  
  std::cout << "cast: " << boost::lexical_cast<int>(utils::ipiece("500")) << std::endl;
}
