#include <stdexcept>
#include <iostream>

#include "wordnet.hpp"

int main(int argc, char** argv)
{
  try {
    wn::WordNet  wordnet(argc > 1 ? std::string(argv[1]) : std::string());
    
    std::string word;
    wn::WordNet::synset_set_type synsets;
    wn::WordNet::morph_set_type morphs;
    
    while (std::getline(std::cin, word)) {
      wordnet(word, morphs);

      if (! morphs.empty()) {
	wn::WordNet::morph_set_type::const_iterator miter_end = morphs.end();
	for (wn::WordNet::morph_set_type::const_iterator miter = morphs.begin(); miter != miter_end; ++ miter)
	  std::cout << "morph: " << *miter << std::endl;
      }
      
      wordnet(word, synsets);
      if (! synsets.empty()) {
	std::cout << "word: " << word << std::endl;
	
	wn::WordNet::synset_set_type::const_iterator siter_end = synsets.end();
	for (wn::WordNet::synset_set_type::const_iterator siter = synsets.begin(); siter != siter_end; ++ siter)
	  std::cout << "pos: " << siter->pos
		    << " word: " << siter->word
		    << " sense: " << siter->sense
		    << std::endl;
	
      }
    }

  } 
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
