// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MATCHER__WORDNET__HPP__
#define __CICADA__MATCHER__WORDNET__HPP__ 1

#include <string>

#include <cicada/matcher.hpp>

#include <wn/wordnet.hpp>

#include <utils/array_power2.hpp>


namespace cicada
{
  namespace matcher
  {
    class WordNet : public cicada::Matcher
    {
    private:
      typedef wn::WordNet wordnet_type;
      
      typedef wordnet_type::synset_type     synset_type;
      typedef wordnet_type::synset_set_type synset_set_type;
      
    public:
      WordNet() : wordnet() {}
      WordNet(const std::string& path) : wordnet(path) {}
      
    public:
      bool operator()(const symbol_type& x, const symbol_type& y) const;

    private:
      struct cache_type
      {
	symbol_type     word;
	synset_set_type synsets;
	
	cache_type() : word(), synsets() {}
      };
      
      typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;
      
    private:
      wordnet_type wordnet;
      cache_set_type caches;
    };
  };
};

#endif
