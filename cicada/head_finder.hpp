// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HEAD_FINDER__HPP__
#define __CICADA__HEAD_FINDER__HPP__ 1

// HeadFinder mostly taken from StanfordParser

#include <string>
#include <vector>

#include <cicada/hypergraph.hpp>

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>


namespace cicada
{
  class HeadFinder
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::symbol_type      symbol_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef enum {
      left,
      leftdis,
      right,
      rightdis,
      leftexcept,
      rightexcept,
    } direction_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > category_set_type;
    typedef std::pair<direction_type, category_set_type> category_list_type;
    typedef std::vector<category_list_type, std::allocator<category_list_type> > category_map_type;

#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<symbol_type, category_map_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				    std::allocator<std::pair<const symbol_type, category_map_type> > > category_info_type;
#else
    typedef sgi::hash_map<symbol_type, category_map_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			  std::allocator<std::pair<const symbol_type, category_map_type> > > category_info_type;
#endif    
    
  public:    
    HeadFinder() : __algorithm() {}
    HeadFinder(const std::string& x) : __algorithm(x) {}
    virtual ~HeadFinder() {}
    
  public:
    size_type operator()(const rule_type& rule) const
    {
      return operator()(rule, symbol_type());
    }

    size_type operator()(const rule_type& rule, const symbol_type& parent) const
    {
      size_type num_tails = 0;
      size_type pos_tail = 0;
      rule_type::symbol_set_type::const_iterator riter_begin = rule.rhs.begin();
      rule_type::symbol_set_type::const_iterator riter_end   = rule.rhs.end();
      for (rule_type::symbol_set_type::const_iterator riter = rule.rhs.begin(); riter != riter_end; ++ riter)
	if (riter->is_non_terminal()) {
	  ++ num_tails;
	  pos_tail = riter - riter_begin;
	}
      
      switch (num_tails) {
      case 0: return size_type(-1); // leaf ... no-head
      case 1: return pos_tail;      // single tail ... this is the head
      default:
	const size_type head_marked = find_marked_head(rule, parent);
	if (head_marked != size_type(-1))
	  return head_marked;
	else
	  return find_head(rule, parent);
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
				 ? operator()(*edge.rule, graph.edges[out_edges[edge.head].front()].rule->lhs)
				 : operator()(*edge.rule));
	
	if (index == size_type(-1))
	  edge.attributes.erase(attribute);
	else
	  edge.attributes[attribute] = attribute_set_type::int_type(index);
      }
    }
    
  protected:
    // actual members called by operator()(const rule_type& rule)
    
    virtual size_type find_marked_head(const rule_type& rule, const symbol_type& parent) const=0;
    virtual size_type find_head(const rule_type& rule, const symbol_type& parent) const=0;
    
  public:
    static HeadFinder& create(const std::string& parameter);
    static const char* lists();

    const std::string& algorithm() const { return __algorithm; }
    
  private:
    std::string __algorithm;
  };
};

#endif
