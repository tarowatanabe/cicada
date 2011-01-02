// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// hungarian algorithm to solve assignment problem...
//
// do we use boost graph? currently, no... since there exists no bipartite graph interface...
//

#ifndef __KUHN__MUNKRES__HPP__
#define __KUHN__MUNKRES__HPP__ 1

#include <stdexcept>
#include <vector>
#include <utility>

#include <boost/numeric/conversion/bounds.hpp>

namespace detail
{
  template <typename CostMatrix>
  struct KuhnMunkresAssignment
  {
    typedef typename CostMatrix::value_type value_type;

    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef difference_type index_type;

    typedef std::vector<value_type, std::allocator<value_type> > label_type;
    
    typedef std::vector<char, std::allocator<char> > covered_type;
    typedef std::vector<index_type, std::allocator<index_type> > tree_type;
    
    typedef std::vector<index_type, std::allocator<index_type> > matched_type;
    
    typedef std::pair<value_type, index_type> index_value_type;
    typedef std::vector<index_value_type, std::allocator<index_value_type> > slack_type;
    
    KuhnMunkresAssignment(const CostMatrix& __costs)
      : costs(__costs),
	label_u(__costs.size1(), boost::numeric::bounds<value_type>::lowest()),
	label_v(__costs.size1(), value_type(0)),
	matched_size(0),
	matched_u(__costs.size1(), index_type(-1)),
	matched_v(__costs.size1(), index_type(-1)),
	covered(__costs.size1()),
	tree(__costs.size1()),
	min_slack(__costs.size1()) {}
    
    void improve_labels(const value_type value)
    {
      const size_type matrix_size = costs.size1();
      
      for (index_type v = 0; v != matrix_size; ++ v) {
	const bool has_parent = (tree[v] >= 0);
	
	label_u[v]         -= value * covered[v];
	label_v[v]         += value * has_parent;
	min_slack[v].first -= value * (! has_parent);
      }
    }
    
    void improve_matching(const index_type v)
    {
      if (tree[v] < 0)
	throw std::runtime_error("no parent?");
      
      const index_type u = tree[v];
      
      if (matched_u[u] >= 0)
	improve_matching(matched_u[u]);
      
      matched_size += bool(matched_u[u] < 0);
      
      matched_u[u] = v;
      matched_v[v] = u;
    }
    
    value_type slack(index_type u, const index_type v) 
    {
      return label_u[u] + label_v[v] - costs(u, v);
    }
    
    void augment()
    {
      const size_type matrix_size = costs.size1();

      for (;;) {
	value_type value = boost::numeric::bounds<value_type>::highest();
	index_type u = 0;
	index_type v = 0;
	
	for (index_type i = 0; i != matrix_size; ++ i)
	  if (tree[i] < 0 && min_slack[i].first < value) {
	    value = min_slack[i].first;
	    u     = min_slack[i].second;
	    v     = i;
	  }
	
	if (! covered[u])
	  throw std::runtime_error("not covered???");
	
	if (value > value_type(0))
	  improve_labels(value);
	
	// (u, v) is saturated! meaning that slack(u, v) == value_type(0)!
	
	tree[v] = u;

	if (matched_v[v] >= 0) {
	  
	  const index_type u1 = matched_v[v];
	  
	  if (covered[u1])
	    throw std::runtime_error("already covered???");
	  
	  covered[u1] = true;
	  
	  for (index_type v = 0; v != matrix_size; ++ v)
	    if (tree[v] < 0 && min_slack[v].first > slack(u1, v))
	      min_slack[v] = std::make_pair(slack(u1, v), u1);
	  
	} else {
	  // v is a free vertex...
	  improve_matching(v);
	  return;
	}
      }
    }
    
    template <typename Iterator>
    void operator()(Iterator result)
    {
      if (costs.size1() != costs.size2())
	throw std::runtime_error("the size must be equal!");
      
      const size_type matrix_size = costs.size1();
      
      for (index_type u = 0; u != matrix_size; ++ u)
	for (index_type v = 0; v != matrix_size; ++ v)
	  label_u[u] = std::max(label_u[u], costs(u, v));
      
      while (matched_size < matrix_size) {
	// choose free vertex: u0
	index_type u0 = 0;
	for (/**/; u0 != matrix_size && matched_u[u0] >= 0; ++ u0);
	
	// initialize tree structure...
	std::fill(tree.begin(), tree.end(), index_type(-1));
	std::fill(covered.begin(), covered.end(), false);
	covered[u0] = true;
	
	// compute min-slack...
	for (index_type v = 0; v < matrix_size; ++ v)
	  min_slack[v] = std::make_pair(slack(u0, v), u0);
	
	// augment!
	augment();
      }
      
      // output results!
      for (index_type u = 0; u != matrix_size; ++ u, ++ result)
	*result = std::make_pair(u, matched_u[u]);
    }
    
    const CostMatrix& costs;
    
    label_type label_u;
    label_type label_v;
    
    size_type     matched_size;
    matched_type  matched_u;
    matched_type  matched_v;
    
    tree_type    tree;
    covered_type covered;
    
    slack_type min_slack;
  };
};

template <typename CostMatrix, typename Iterator>
inline
void kuhn_munkres_assignment(const CostMatrix& costs,
			     Iterator result)
{
  detail::KuhnMunkresAssignment<CostMatrix> assigner(costs);
  
  assigner(result);
}

#endif

