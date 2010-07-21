// -*- mode: c++ -*-

#ifndef __CICADA__STEMMER_SIMPLE__HPP__
#define __CICADA__STEMMER_SIMPLE__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  class StemmerPrefix : public Stemmer
  {
  private:
    typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    
  public:
    StemmerPrefix(const size_type __size) : size(__size) {}
    
  public:
    symbol_type operator[](const symbol_type& x) const;
    
  private:
    symbol_set_type cache;
    size_type       size;
  };

  class StemmerSuffix : public Stemmer
  {
  private:
    typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    
  public:
    StemmerSuffix(const size_type __size) : size(__size) {}
    
  public:
    symbol_type operator[](const symbol_type& x) const;
    
  private:
    symbol_set_type cache;
    size_type       size;
  };
  
};

#endif
