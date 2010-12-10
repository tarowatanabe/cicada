// -*- mode: c++ -*-

#ifndef __WN__WORDNET__HPP__
#define __WN__WORDNET__HPP__ 1

#include <string>
#include <vector>

namespace wn
{
  class WordNet
  {
  public:
    struct SynSet
    {
      std::string pos;
      std::string word;
      int sense;
      
      SynSet() {}
    };
    typedef SynSet synset_type;
    typedef std::vector<synset_type, std::allocator<synset_type> > synset_set_type;
    
  public:
    WordNet()  { initialize(""); }
    WordNet(const std::string& path) { initialize(path); }
    
  public:
    void operator()(const std::string& word, synset_set_type& synsets) const;
    
  private:
    static void initialize(const std::string& path);
  };
};

#endif
