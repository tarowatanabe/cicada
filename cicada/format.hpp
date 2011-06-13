// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__HPP__
#define __CICADA__FORMAT__HPP__ 1

#include <string>
#include <vector>

#include <utils/piece.hpp>

namespace cicada
{
  class Format
  {
  public:
    typedef std::string phrase_type;
    typedef std::vector<phrase_type, std::allocator<phrase_type> > phrase_set_type;

  public:
    virtual ~Format() {}

  public:
    static Format&  create(const utils::piece& parameter);
    static const char* lists();
    
  public:  
    virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const = 0;

  private:
    std::string __algorithm;
  };
};

#endif

