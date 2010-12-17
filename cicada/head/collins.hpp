// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HEAD__COLLINS__HPP__
#define __CICADA__HEAD__COLLINS__HPP__ 1

#include <cicada/head_finder.hpp>

namespace cicada
{
  namespace head
  {
    
    class Collins : public cicada::HeadFinder
    {
    public:
      Collins();
      
      size_type find_marked_head(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const { return size_type(-1); }
      size_type find_head(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const
      {
	const symbol_type& mother = edge.rule->lhs;
	
	rule_type::symbol_set_type::const_iterator riter_begin = edge.rule->rhs.begin();
	rule_type::symbol_set_type::const_iterator riter_end   = edge.rule->rhs.end();
	
	category_info_type::const_iterator citer = categories.find(mother);
	if (citer == categories.end())
	  return size_type(-1);
	
	category_map_type::const_iterator iter_begin = citer->second.begin();
	category_map_type::const_iterator iter_end = citer->second.end();
	for (category_map_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
	  const bool fallback = (iter == iter_end - 1);

	  rule_type::symbol_set_type::const_iterator riter = riter_end;
	  
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
	      while (riter != riter_end && ! riter->is_non_terminal()) ++ riter;
	      break;
	    case right:
	    case rightdis:
	    case rightexcept:
	      riter = riter_end;
	      while (riter != riter_begin && ! (riter - 1)->is_non_terminal()) -- riter;
	      -- riter;
	      break;
	    }
	  
	  if (riter != riter_end) {
	    typedef std::vector<rule_type::symbol_set_type::const_iterator, std::allocator<rule_type::symbol_set_type::const_iterator> > iterator_set_type;
	    iterator_set_type iterators;
	    for (rule_type::symbol_set_type::const_iterator iiter = riter_begin; iiter != riter_end; ++ iiter)
	      if (iiter->is_non_terminal())
		iterators.push_back(iiter);
	    
	    iterator_set_type::const_iterator iiter_begin = iterators.begin();
	    iterator_set_type::const_iterator iiter_end   = iterators.end();
	    iterator_set_type::const_iterator iiter = std::find(iiter_begin, iiter_end, riter);
	    if (iiter >= iterators.begin() + 2) {
	      const symbol_type cat_prev = (*(iiter - 1))->non_terminal();
	      if (cat_prev == "[CC]" || cat_prev == "[CONJP]") {
		// skip punctuation that is pre-terminal...
		
		iterator_set_type::const_iterator iiter_new = skip_pre_terminals(graph, edge.tails, iiter_begin, iiter - 1, iiter_end);
		
		if (iiter_new != iiter_end)
		  iiter = iiter_new;
	      }
	    }
	    
	    return *iiter - riter_begin;
	  }
	}
	
	return size_type(-1);
      }
      
    private:
      
      template <typename Iterator, typename Tails>
      Iterator skip_pre_terminals(const hypergraph_type& graph, const Tails& tails, Iterator first, Iterator iter, Iterator last) const
      {
	for (/**/; iter != first; -- iter) {
	  const symbol_type& cat = *(*(iter - 1));
	  
	  if (! is_punctuation(cat.non_terminal())) return iter - 1;
	  
	  // cat is punctuation...
	  int pos = cat.non_terminal_index() - 1;
	  if (pos < 0)
	    pos = (iter - 1) - first;
	  
	  if (pos >= static_cast<int>(tails.size()))
	    throw std::runtime_error("invalid tails");
	  
	  if (graph.nodes[tails[pos]].edges.empty())
	    throw std::runtime_error("no edges");
	  
	  const edge_type& edge = graph.edges[graph.nodes[tails[pos]].edges.front()];
	  
	  // we have tail, meaning that cat is NOT pre-terminal
	  if (! edge.tails.empty()) return iter - 1;
	}
	
	return last;
      }
      
      bool is_punctuation(const symbol_type& cat) const
      {
	category_set_type::const_iterator citer = std::lower_bound(punctuations.begin(), punctuations.end(), cat);
	return citer != punctuations.end() && *citer == cat;
      }
      
      
    private:
      category_info_type categories;
      category_set_type  punctuations;
    };
  };
};

#endif
