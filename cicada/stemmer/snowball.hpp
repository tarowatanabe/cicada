// -*- mode: c++ -*-

#ifndef __CICADA__STEMMER_SNOWBALL__HPP__
#define __CICADA__STEMMER_SNOWBALL__HPP__ 1

#include <vector>

#include <cicada/stemmer.hpp>

namespace cicada
{
  namespace stemmer
  {
    struct SnowballImpl;

    class Snowball : public Stemmer
    {
    private:
      typedef SnowballImpl impl_type;
      typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
      
    public:
      Snowball(const std::string& language);
      ~Snowball();
    
    public:
      symbol_type operator[](const symbol_type& x) const;
    
    private:
      impl_type*      pimpl;
      symbol_set_type cache;
    };
  };
};

#endif
