// -*- mode: c++ -*-

#ifndef __CICADA__STEMMER_PREFIX__HPP__
#define __CICADA__STEMMER_PREFIX__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    class Prefix : public Stemmer
    {
    private:
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    
    public:
      Prefix(const size_type __size) : size(__size) {}
    
    public:
      symbol_type operator[](const symbol_type& x) const;
    
    private:
      symbol_set_type cache;
      size_type       size;
    };
  };
};

#endif
