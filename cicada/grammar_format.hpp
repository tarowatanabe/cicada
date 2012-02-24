// -*- mode: c++ -*-
//
//  Copyright(C) 2011-12 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR_FORMAT__HPP__
#define __CICADA__GRAMMAR_FORMAT__HPP__ 1

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/format.hpp>

#include <utils/hashmurmur.hpp>

namespace cicada
{
  class GrammarFormat : public GrammarMutable
  {
  private:
    typedef GrammarMutable base_type;
    
    typedef int32_t uchar_type;
    
    typedef feature_set_type::feature_type feature_type;
    typedef Format format_type;
    
    typedef std::string prefix_type;
    typedef std::vector<prefix_type, std::allocator<prefix_type> > prefix_set_type;
    
  public:
    GrammarFormat(const symbol_type& __non_terminal,
		  const std::string& param_formatter,
		  const bool __remove_space)
      : base_type(),
	non_terminal(__non_terminal),
	format(&format_type::create(param_formatter)),
	remove_space(__remove_space),
	prefix(),
	feature("format-penalty")
    { }
    
    transducer_ptr_type clone() const
    {
      std::auto_ptr<GrammarFormat> __tmp(new GrammarFormat(*this));
      __tmp->format = &format_type::create(format->algorithm());
      return transducer_ptr_type(__tmp.release());
    }

    void assign(const hypergraph_type& graph) { base_type::clear(); prefix.clear(); }
    void assign(const lattice_type& lattice) { base_type::clear(); prefix.clear(); }
    
    id_type next(const id_type& node, const symbol_type& symbol) const;
    
    // has_next is always true...!
    bool has_next(const id_type& node) const { return true; }
    
  private:
    // we will keep actual rules in base_type + queried types via "node-id + word"
    symbol_type non_terminal;
    const format_type* format;
    bool remove_space;
    
    prefix_set_type prefix;
    
    feature_type feature;
  };
};

#endif
