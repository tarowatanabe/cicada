// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MATCHER__HPP__
#define __CICADA__MATCHER__HPP__ 1

#include <string>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

namespace cicada
{
  class Matcher
  {
  public:
    typedef Symbol    symbol_type;
    typedef Vocab     vocab_type;
    
    typedef symbol_type          word_type;
    typedef symbol_type::id_type id_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  public:
    Matcher() {}
    virtual ~Matcher() {}
    
  private:
    // we do not allow copy/construct
    Matcher& operator=(const Matcher& x) { return *this; }
    Matcher(const Matcher& x) {}
    
  public:
    static Matcher&    create(const std::string& parameter);
    static const char* lists();
    
  public:
    virtual bool operator()(const symbol_type& x, const symbol_type& y) const = 0;
    
    const std::string& algorithm() const { return __algorithm; }

  private:
    std::string __algorithm;
  };
};

#endif
