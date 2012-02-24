// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HEAD_FINDER__HPP__
#define __CICADA__HEAD_FINDER__HPP__ 1

// HeadFinder mostly taken from StanfordParser

#include <string>
#include <vector>
#include <algorithm>

#include <cicada/hypergraph.hpp>

#include <utils/unordered_map.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  class HeadFinder
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::edge_type        edge_type;
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef enum {
      left,
      leftdis,
      leftexcept,
      right,
      rightdis,
      rightexcept,
    } direction_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > category_set_type;
    typedef std::pair<direction_type, category_set_type> category_list_type;
    typedef std::vector<category_list_type, std::allocator<category_list_type> > category_map_type;

    typedef utils::unordered_map<symbol_type, category_map_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				 std::allocator<std::pair<const symbol_type, category_map_type> > >::type category_info_type;
    
  public:    
    HeadFinder() : __algorithm() {}
    HeadFinder(const std::string& x) : __algorithm(x) {}
    virtual ~HeadFinder() {}
    
  public:
    size_type operator()(const hypergraph_type& graph, const edge_type& edge) const
    {
      return operator()(graph, edge, symbol_type());
    }

    size_type operator()(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const
    {
      size_type num_tails = 0;
      size_type pos_tail = 0;
      rule_type::symbol_set_type::const_iterator riter_begin = edge.rule->rhs.begin();
      rule_type::symbol_set_type::const_iterator riter_end   = edge.rule->rhs.end();
      for (rule_type::symbol_set_type::const_iterator riter = riter_begin; riter != riter_end; ++ riter)
	if (riter->is_non_terminal()) {
	  ++ num_tails;
	  pos_tail = riter - riter_begin;
	}
      
      switch (num_tails) {
      case 0: return size_type(-1); // leaf ... no-head
      case 1: return pos_tail;      // single tail ... this is the head
      default:
	const size_type head_marked = find_marked_head(graph, edge, parent);
	if (head_marked != size_type(-1))
	  return head_marked;
	else
	  return find_head(graph, edge, parent);
      }
    }
    
    void operator()(const hypergraph_type& source, hypergraph_type& target) const
    {
      target = source;
      operator()(target);
    }
    
    void operator()(hypergraph_type& graph) const
    {
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;
      typedef std::vector<edge_set_type, std::allocator<edge_set_type> > node_map_type;

      // first, collect out-edges...
      node_map_type out_edges(graph.nodes.size());
      {
	hypergraph_type::edge_set_type::iterator eiter_end = graph.edges.end();
	for (hypergraph_type::edge_set_type::iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = *eiter;
	  
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_begin = edge.tails.begin();
	  hypergraph_type::edge_type::node_set_type::const_iterator titer_end   = edge.tails.end();
	  
	  for (hypergraph_type::edge_type::node_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    out_edges[*titer].push_back(edge.id);
	}
      }
      
      const attribute_type attribute("head");
      
      hypergraph_type::edge_set_type::iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) {
	hypergraph_type::edge_type& edge = *eiter;
	
	const size_type index = (! out_edges[edge.head].empty()
				 ? operator()(graph, edge, graph.edges[out_edges[edge.head].front()].rule->lhs)
				 : operator()(graph, edge));
	
	if (index == size_type(-1))
	  edge.attributes.erase(attribute);
	else
	  edge.attributes[attribute] = attribute_set_type::int_type(index);
      }
    }
    
  protected:
    // actual members called by operator()(const rule_type& rule)
    
    virtual size_type find_marked_head(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const=0;
    virtual size_type find_head(const hypergraph_type& graph, const edge_type& edge, const symbol_type& parent) const=0;

  protected:
    
    template <typename Iterator>
    Iterator traverse_left(const category_set_type& categories, Iterator first, Iterator last) const
    {
      category_set_type::const_iterator citer_end = categories.end();
      for (category_set_type::const_iterator citer = categories.begin(); citer != citer_end; ++ citer) {
	Iterator iter = std::find(first, last, *citer);
	if (iter != last)
	  return iter;
      }
      return last;
    }
    
    
    template <typename Iterator>
    Iterator traverse_leftdis(const category_set_type& categories, Iterator first, Iterator last) const
    {
      for (/**/; first != last; ++ first) 
	if (first->is_non_terminal()) {
	  category_set_type::const_iterator citer = std::find(categories.begin(), categories.end(), *first);
	  if (citer != categories.end())
	    return first;
	}
      return last;
    }

    template <typename Iterator>
    Iterator traverse_leftexcept(const category_set_type& categories, Iterator first, Iterator last) const
    {
      for (/**/; first != last; ++ first) 
	if (first->is_non_terminal()) {
	  category_set_type::const_iterator citer = std::find(categories.begin(), categories.end(), *first);
	  if (citer == categories.end())
	    return first;
	}
      return last;
    }

    template <typename Iterator>
    Iterator traverse_right(const category_set_type& categories, Iterator first, Iterator last) const
    {
      category_set_type::const_iterator citer_end = categories.end();
      for (category_set_type::const_iterator citer = categories.begin(); citer != citer_end; ++ citer) {
	for (Iterator iter = last; iter != first; -- iter)
	  if (*citer == *(iter - 1))
	    return iter - 1;
      }
      return last;
    }

    template <typename Iterator>
    Iterator traverse_rightdis(const category_set_type& categories, Iterator first, Iterator last) const
    {
      
      for (Iterator iter = last; iter != first; -- iter)
	if (first->is_non_terminal()) {
	  category_set_type::const_iterator citer = std::find(categories.begin(), categories.end(), *(iter - 1));
	  if (citer != categories.end())
	    return iter - 1;
	}
      return last;
    }
    
    template <typename Iterator>
    Iterator traverse_rightexcept(const category_set_type& categories, Iterator first, Iterator last) const
    {
      for (Iterator iter = last; iter != first; -- iter) 
	if (first->is_non_terminal()) {
	  category_set_type::const_iterator citer = std::find(categories.begin(), categories.end(), *(iter - 1));
	  if (citer == categories.end())
	    return iter - 1;
	}
      return last;
    }

    template <size_t N>
    category_set_type assign_category(const char* (&nt)[N]) const
    {
      category_set_type categories;
      for (size_t i = 0; i != N; ++ i)
	categories.push_back('[' + std::string(nt[i]) + ']');
      return categories;
    }

    category_set_type assign_category(const char* (&nt)[0])
    {
      return category_set_type();
    }
    
    
  public:
    static HeadFinder& create(const utils::piece& parameter);
    static const char* lists();

    const std::string& algorithm() const { return __algorithm; }
    
  private:
    std::string __algorithm;
  };
};

#endif
