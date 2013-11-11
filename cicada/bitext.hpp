// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__BITEXT__HPP__
#define __CICADA__BITEXT__HPP__ 1

#include <cicada/sentence.hpp>

namespace cicada
{
  class Bitext
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Sentence sentence_type;
    
    typedef sentence_type::symbol_type symbol_type;
    typedef sentence_type::word_type   word_type;

  public:
    Bitext() : source_(), target_() {}
    Bitext(const std::pair<sentence_type, sentence_type>& x) : source_(x.first), target_(x.second) {}
    Bitext(const sentence_type& source, const sentence_type& target) : source_(source), target_(target) {}
    
    void clear()
    {
      source_.clear();
      target_.clear();
    }
    
    void swap(Bitext& x)
    {
      source_.swap(x.source_);
      target_.swap(x.target_);
    }
    
  public:
    sentence_type source_;
    sentence_type target_;
  };
  
};

namespace std
{
  inline
  void swap(cicada::Bitext& x, cicada::Bitext& y)
  {
    x.swap(y);
  }
};

#endif
