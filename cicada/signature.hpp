// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SIGNATURE__HPP__
#define __CICADA__SIGNATURE__HPP__ 1

#include <string>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/piece.hpp>

namespace cicada
{
  class Signature
  {
  public:
    typedef Symbol    symbol_type;
    typedef Vocab     vocab_type;
    
    typedef symbol_type          word_type;
    typedef symbol_type::id_type id_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  public:
    Signature() {}
    virtual ~Signature() {}
    
  private:
    // we do not allow copy/construct
    Signature& operator=(const Signature& x) { return *this; }
    Signature(const Signature& x) {}
    
  public:
    static Signature&    create(const utils::piece& parameter);
    static const char* lists();
    
  public:
    symbol_type operator()(const symbol_type& x) const { return operator[](x); }
    virtual symbol_type operator[](const symbol_type& x) const = 0;
    const std::string& algorithm() const { return __algorithm; }

  private:
    std::string __algorithm;
  };
};

#endif
