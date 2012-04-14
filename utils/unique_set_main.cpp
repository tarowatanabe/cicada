
#include "unique_set.hpp"

#include <iostream>
#include <string>
#include <set>

int main (int argc, char** argv)
{
  typedef utils::unique_set<std::string> unique_set_type;
  typedef std::set<std::string> word_set_type;
  
  word_set_type  words;
  unique_set_type uniques;
  std::string word;
  
  while (std::cin >> word) {
    uniques[word];
    words.insert(word);
  }
  
  std::cerr << "unique size: " << uniques.size() << std::endl
	    << "set sizes: " << words.size() << std::endl;
  
  word_set_type::const_iterator witer_end = words.end();
  for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer)
    if (uniques.find(*witer) == uniques.end())
      std::cerr << "not inserted: " << *witer <<std::endl;
}
