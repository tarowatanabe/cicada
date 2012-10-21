// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TRAITS__HPP__
#define __CICADA__TRAITS__HPP__ 1

#include <cicada/symbol.hpp>
#include <cicada/feature.hpp>
#include <cicada/attribute.hpp>

#include <utils/compact_func.hpp>

namespace utils
{
  template <>
  struct unassigned<cicada::Symbol>
  {
    cicada::Symbol operator()() const { return cicada::Symbol(cicada::Symbol::id_type(-1)); }
  };
  
  template <>
  struct deleted<cicada::Symbol>
  {
    cicada::Symbol operator()() const { return cicada::Symbol(cicada::Symbol::id_type(-2)); }
  };

  template <>
  struct unassigned<cicada::Attribute>
  {
    cicada::Attribute operator()() const { return cicada::Attribute(cicada::Attribute::id_type(-1)); }
  };
  
  template <>
  struct deleted<cicada::Attribute>
  {
    cicada::Attribute operator()() const { return cicada::Attribute(cicada::Attribute::id_type(-2)); }
  };

  template <>
  struct unassigned<cicada::Feature>
  {
    cicada::Feature operator()() const { return cicada::Feature(cicada::Feature::id_type(-1)); }
  };
  
  template <>
  struct deleted<cicada::Feature>
  {
    cicada::Feature operator()() const { return cicada::Feature(cicada::Feature::id_type(-2)); }
  };

  
};

#endif
