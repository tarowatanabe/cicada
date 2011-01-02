// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HEAD__CHINESE__HPP__
#define __CICADA__HEAD__CHINESE__HPP__ 1

#include <vector>

#include <cicada/head_finder.hpp>

namespace cicada
{
  namespace head
  {
    
    class Chinese : public cicada::HeadFinder
    {
    public:
      Chinese();
      
      size_type find_marked_head(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const { return size_type(-1); }
      size_type find_head(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const
      {
	typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
	typedef std::vector<int, std::allocator<int> > index_set_type;
	
	const symbol_type& mother = edge.rule->lhs;
	
	symbol_set_type rhs;
	index_set_type  index;
	
	{
	  rule_type::symbol_set_type::const_iterator riter_begin = edge.rule->rhs.begin();
	  rule_type::symbol_set_type::const_iterator riter_end   = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter)
	    if (riter->is_non_terminal()) {
	      rhs.push_back(*riter);
	      index.push_back(riter - riter_begin);
	    }
	}
	
	symbol_set_type::const_iterator riter_begin = rhs.begin();
	symbol_set_type::const_iterator riter_end   = rhs.end();

	// default to right most...
	
	category_info_type::const_iterator citer = categories.find(mother);
	if (citer == categories.end())
	  return index.back();
	
	category_map_type::const_iterator iter_begin = citer->second.begin();
	category_map_type::const_iterator iter_end   = citer->second.end();
	for (category_map_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	  const bool fallback = (iter == iter_end - 1);

	  symbol_set_type::const_iterator riter = riter_end;
	  
	  switch (iter->first) {
	  case left:        riter = traverse_left(iter->second, riter_begin, riter_end); break;
	  case leftdis:     riter = traverse_leftdis(iter->second, riter_begin, riter_end); break;
	  case leftexcept:  riter = traverse_leftexcept(iter->second, riter_begin, riter_end); break;
	  case right:       riter = traverse_right(iter->second, riter_begin, riter_end); break;
	  case rightdis:    riter = traverse_rightdis(iter->second, riter_begin, riter_end); break;
	  case rightexcept: riter = traverse_rightexcept(iter->second, riter_begin, riter_end); break;
	  }

	  if (riter == riter_end && fallback)
	    switch (iter->first) {
	    case left:
	    case leftdis:
	    case leftexcept:
	      riter = riter_begin;
	      break;
	    case right:
	    case rightdis:
	    case rightexcept:
	      riter = riter_end - 1;
	      break;
	    }
	  
	  if (riter != riter_end)
	    return index[riter - riter_begin];
	}
	
	return size_type(-1);
      }
      
    private:
      category_info_type categories;
    };
  };
};

#endif
