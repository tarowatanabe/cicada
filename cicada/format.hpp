// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__HPP__
#define __CICADA__FORMAT__HPP__ 1

#include <string>
#include <vector>

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
    const phrase_set_type& operator[](const phrase_type& phrase) const { return operator()(phrase); }
    virtual const phrase_set_type& operator()(const phrase_type& phrase) const = 0;
  };
};

#endif

