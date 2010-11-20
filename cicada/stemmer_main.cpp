#include "stemmer.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cout << argv[0] << " stemmer-spec" << std::endl;
    return 1;
  }

  std::cout << cicada::Stemmer::lists();
  
  cicada::Stemmer& stemmer(cicada::Stemmer::create(argv[1]));
  
  std::string word;
  while (std::cin >> word)
    std::cout << "word: " << word << " stemmed: " << stemmer[word] << std::endl;

  
}
