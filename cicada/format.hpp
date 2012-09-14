// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
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
    typedef std::string tag_type;
    
    struct phrase_tag_type
    {
      phrase_type phrase;
      tag_type    tag;
      
      phrase_tag_type()
	: phrase(), tag() {}
      phrase_tag_type(const phrase_type& __phrase, const tag_type& __tag)
	: phrase(__phrase), tag(__tag) {}
      phrase_tag_type(const std::pair<phrase_type, tag_type>& x)
	: phrase(x.first), tag(x.second) {}
    };
    typedef std::vector<phrase_tag_type, std::allocator<phrase_tag_type> > phrase_set_type;

  public:
    virtual ~Format() {}

  public:
    static Format&  create(const utils::piece& parameter);
    static const char* lists();
    
  public:  
    virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const = 0;
    const std::string& algorithm() const { return __algorithm; }

  private:
    std::string __algorithm;
  };
};

#endif

