// -*- mode: c++ -*-

#ifndef __CICADA__RULE__HPP__
#define __CICADA__RULE__HPP__ 1

#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol_vector.hpp>
#include <cicada/feature_vector.hpp>

namespace cicada
{
  class Rule
  {
  public:
    typedef cicada::Symbol       symbol_type;
    typedef cicada::Symbol       word_type;
    typedef cicada::Vocab        vocab_type;
    typedef cicada::SymbolVector symbol_set_type;
    typedef cicada::FeatureVector<double, std::allocator<double> > feature_set_type;
    
  public:
    Rule(const std::string& x) { assign(x); }
    Rule() : arity(0) {}
    Rule(const symbol_type& x_lhs,
	 const symbol_set_type& x_source)
      : lhs(x_lhs), source(x_source), arity(0) {}
    
    Rule(const symbol_type& x_lhs,
	 const symbol_set_type& x_source,
	 const symbol_set_type& x_target)
      : lhs(x_lhs), source(x_source), target(x_target), arity(0) {}

    Rule(const symbol_type& x_lhs,
	 const symbol_set_type& x_source,
	 const symbol_set_type& x_target,
	 const int x_arity)
      : lhs(x_lhs), source(x_source), target(x_target), arity(x_arity) {}

    void assign(const std::string& x);

    void clear()
    {
      lhs = symbol_type();
      source.clear();
      target.clear();
      features.clear();
      arity = 0;
    }
    
    // sort non-terminal index wrt source-side or target-side
    void sort_source_index();
    void sort_target_index();

    friend
    std::ostream& operator<<(std::ostream& os, const Rule& x);
    friend
    std::istream& operator>>(std::istream& is, Rule& x);
    
  public:
    symbol_type      lhs;
    symbol_set_type  source;
    symbol_set_type  target;
    feature_set_type features;
    int              arity;
  };
};

#endif
