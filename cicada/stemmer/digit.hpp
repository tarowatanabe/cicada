// -*- mode: c++ -*-

#ifndef __CICADA__STEMMER_DIGIT__HPP__
#define __CICADA__STEMMER_DIGIT__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Digit : public Stemmer
    {
    private:
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    
    public:
      Digit() {}
    
    public:
      symbol_type operator[](const symbol_type& x) const;
    
    private:
      symbol_set_type cache;
    };
  };
};

#endif
