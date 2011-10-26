// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_CACHE__HPP__
#define __CICADA__NGRAM_CACHE__HPP__ 1

#include <cicada/symbol.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/vector2.hpp>
#include <utils/array_power2.hpp>

namespace cicada
{
  class NGramCache : public utils::hashmurmur<size_t>
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol symbol_type;
    
    typedef utils::hashmurmur<size_t> hasher_type;

    static const size_type cache_size = 1024 * 64;
    
    typedef utils::vector2<symbol_type::id_type, std::allocator<symbol_type::id_type> > context_set_type;
    typedef utils::array_power2<double, cache_size, std::allocator<double> >            score_set_type;
    typedef utils::array_power2<uint8_t, cache_size, std::allocator<uint8_t> >          length_set_type;
    
  public:
    NGramCache(const int order=3)
      : contexts(cache_size, order), scores(), length() {}
    
  public:
    void clear()
    {
      scores.clear();
      length.clear();
    }
    
    template <typename Iterator>
    size_type operator()(Iterator first, Iterator last) const
    {
      return hasher_type::operator()(first, last, 0) & (cache_size - 1);
    }
    
    inline       double& score(size_type pos) { return scores[pos]; }
    inline const double& score(size_type pos) const { return scores[pos]; }
    
    template <typename Iterator>
    bool equal_to(size_type pos, Iterator first, Iterator last) const
    {
      return std::distance(first, last) == length[pos] && std::equal(first, last, contexts.begin(pos));
    }
    
    template <typename Iterator>
    void assign(size_type pos, Iterator first, Iterator last)
    {
      std::copy(first, last, contexts.begin(pos));
      length[pos] = last - first;
    }
    
  private:
    context_set_type contexts;
    score_set_type   scores;
    length_set_type  length;
  };
};

#endif
