// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_SHARED__HPP__
#define __CICADA__GRAMMAR_SHARED__HPP__ 1

// very siple shared grammar class..

#include <string>

#include <cicada/transducer.hpp>
#include <cicada/grammar_mutable.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  
  class GrammarShared : public Transducer
  {
  private:
    typedef GrammarMutable impl_type;
    
  public:
    GrammarShared(const std::string& parameter) : pimpl(new impl_type(parameter)) { }
    GrammarShared(const GrammarShared& x) : pimpl(x.pimpl) {}
    GrammarShared& operator=(const GrammarShared& x)
    {
      pimpl = x.pimpl;
      return *this;
    }
    
  private:
    GrammarShared(const int __max_span=0) {}
    
  public:
    // virtual members
    transducer_ptr_type clone() const { return transducer_ptr_type(new GrammarShared(*this)); }
    
    bool valid_span(int first, int last, int distance) const { return pimpl->valid_span(first, last, distance); }
    id_type root() const { return pimpl->root(); }
    id_type next(const id_type& node, const symbol_type& symbol) const { return pimpl->next(node, symbol); }
    bool has_next(const id_type& node) const { return pimpl->has_next(node); }
    const rule_pair_set_type& rules(const id_type& node) const { return pimpl->rules(node); }
    
  private:
    boost::shared_ptr<impl_type> pimpl;
  };
  
};

#endif
