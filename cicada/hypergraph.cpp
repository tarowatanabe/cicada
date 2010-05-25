#include "hypergraph.hpp"

#define BOOST_SPIRIT_THREADSAFE



namespace cicada
{

  // topologically sort...
  
  struct TopologicallySortImpl
  {
    typedef HyperGraph hypergraph_type;

    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
 
    enum color_type {
      white,
      gray,
      black
    };

    struct dfs_type
    {
      id_type node;
      int edge;
      int tail;
      
      dfs_type(const id_type& _node, const int& _edge, const int& _tail) 
	: node(_node), edge(_edge), tail(_tail) {}
    };
    
    typedef std::vector<int, std::allocator<int> > reloc_set_type;
    typedef std::vector<color_type, std::allocator<color_type> > color_set_type;
    typedef std::vector<dfs_type, std::allocator<dfs_type> > stack_type;

    struct no_filter_edge
    {
      bool operator()(const edge_type& edge) const
      {
	return false;
      }
    };
    
    template <typename Filter>
    void operator()(hypergraph_type& x, Filter filter)
    {
      reloc_set_type reloc_node(x.nodes.size(), -1);
      reloc_set_type reloc_edge(x.edges.size(), -1);
      color_set_type color(x.nodes.size(), white);
      stack_type stack;
      
      stack.reserve(x.nodes.size());
      stack.push_back(dfs_type(x.goal, 0, 0));
      int node_count = 0;
      int edge_count = 0;
      
      while (! stack.empty()) {
	const dfs_type& dfs = stack.back();
	id_type node_id = dfs.node;
	int pos_edge = dfs.edge;
	int pos_tail = dfs.tail;
	
	stack.pop_back();
	
	const node_type* curr_node = &(x.nodes[node_id]);
	while (pos_edge != curr_node->in_edges.size()) {
	  const edge_type& curr_edge = x.edges[curr_node->in_edges[pos_edge]];
	  
	  if (pos_tail == curr_edge.tail_nodes.size() || filter(curr_edge)) {
	    // reach end: proceed to the next edge with pos_tail initialized
	    ++ pos_edge;
	    pos_tail = 0;
	    continue;
	  }
	  
	  const id_type tail_node = curr_edge.tail_nodes[pos_tail];
	  const color_type tail_color = color[tail_node];
	  
	  switch (tail_color) {
	  case white:
	    ++ pos_tail;
	    stack.push_back(dfs_type(node_id, pos_edge, pos_tail));
	    
	    node_id = tail_node;
	    curr_node = &(x.nodes[node_id]);
	    
	    color[node_id] = gray;
	    pos_edge = 0;
	    pos_tail = 0;
	    
	    break;
	  case black:
	    ++ pos_tail;
	    break;
	  case gray:
	    throw std::runtime_error("detected cycle!");
	    break;
	  }
	}
	
	color[node_id] = black;
	reloc_node[node_id] = node_count ++;
	for (int i = 0; i < curr_node->in_edges.size(); ++ i)
	  if (! filter(x.edges[curr_node->in_edges[i]]))
	    reloc_edge[curr_node->in_edges[i]] = edge_count ++;
      }
      
      // sorted graph!
      HyperGraph sorted;
      
      // construct edges...
      for (int i = 0; i < reloc_edge.size(); ++ i)
	if (reloc_edge[i] >= 0) {
	  const edge_type& edge_old = x.edges[i];

	  const id_type edge_id = sorted.edges.size();
	  
	  sorted.edges.push_back(edge_old);
	  
	  edge_type& edge_new = sorted.edges.back();
	  
	  edge_new.id = edge_id;
	  
	  edge_new.head_node = reloc_node[edge_new.head_node];
	  edge_type::node_set_type::iterator niter_end = edge_new.tail_nodes.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tail_nodes.begin(); niter != niter_end; ++ niter)
	    *niter = reloc_node[*niter];
	  
	  reloc_edge[i] = edge_id;
	}
      
      
      // construct reverse node-map ...
      reloc_set_type reloc_map_node(node_count, -1);
      for (int i = 0; i < x.nodes.size(); ++ i)
	if (reloc_node[i] >= 0)
	  reloc_map_node[reloc_node[i]] = i;
      
      // construct nodes
      for (int i = 0; i < reloc_map_node.size(); ++ i) {
	const node_type& node_old = x.nodes[reloc_map_node[i]];
	node_type& node_new = sorted.add_node();
	
	node_type::edge_set_type::const_iterator eiter_end = node_old.in_edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node_old.in_edges.begin(); eiter != eiter_end; ++ eiter)
	  if (reloc_edge[*eiter] >= 0)
	    node_new.in_edges.push_back(reloc_edge[*eiter]);
	
	eiter_end = node_old.out_edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node_old.out_edges.begin(); eiter != eiter_end; ++ eiter)
	  if (reloc_edge[*eiter] >= 0)
	    node_new.in_edges.push_back(reloc_edge[*eiter]);
      }
      
      sorted.goal = sorted.nodes.size() - 1;
      sorted.is_sorted = true;
      
      sorted.swap(x);
    }
  };
  
  void HyperGraph::topologically_sort()
  {
    TopologicallySortImpl()(*this, TopologicallySortImpl::no_filter_edge());
  }
  
  // comput union
  struct UnionImpl
  {
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::id_type id_type;
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    typedef Rule   rule_type;
    
    void operator()(hypergraph_type& x, const hypergraph_type& y)
    {
      if (&x == &y || y.nodes.empty() || y.goal == hypergraph_type::invalid) return;
      if (x.nodes.empty() || x.goal == hypergraph_type::invalid) {
	x = y;
	return;
      }
      
      // check if we share the same goal... otherwise, we will create new edge and fill...
      //
      // 
      const symbol_type& goal_x = x.edges[x.nodes[x.goal].in_edges.front()].rule->lhs;
      const symbol_type& goal_y = y.edges[y.nodes[y.goal].in_edges.front()].rule->lhs;

      if (goal_x == goal_y) {
	const int y_node_offset = x.nodes.size();
	const int y_edge_offset = x.edges.size();
	
	// -1 to adjust goal
	x.nodes.resize(x.nodes.size() + y.nodes.size() - 1);
	x.edges.resize(x.edges.size() + y.edges.size());
	
	x.is_sorted = false;
	
	// merge nodes
	for (id_type id = 0; id < y.nodes.size(); ++ id)
	  if (id != y.goal) {
	    const id_type id_new = id + y_node_offset - (id > y.goal);
	    
	    const node_type& node_old = y.nodes[id];
	    node_type& node_new = x.nodes[id_new];
	    
	    node_new = node_old;
	    
	    node_new.id =  id_new;
	    node_type::edge_set_type::iterator eiter_end = node_new.in_edges.end();
	    for (node_type::edge_set_type::iterator eiter = node_new.in_edges.begin(); eiter != eiter_end; ++ eiter)
	      *eiter += y_edge_offset;
	    
	    eiter_end = node_new.out_edges.end();
	    for (node_type::edge_set_type::iterator eiter = node_new.out_edges.begin(); eiter != eiter_end; ++ eiter)
	      *eiter += y_edge_offset;
	  }
	
	// merge edges
	for (id_type id = 0; id < y.edges.size(); ++ id) {
	  const edge_type& edge_old = y.edges[id];
	  edge_type& edge_new = x.edges[id + y_edge_offset];
	  
	  edge_new = edge_old;
	  edge_new.id = id + y_edge_offset;
	  
	  if (edge_new.head_node == y.goal) {
	    edge_new.head_node = x.goal;
	    x.nodes[x.goal].in_edges.push_back(edge_new.id);
	  } else
	    edge_new.head_node += y_node_offset - (edge_new.head_node > y.goal);
	  
	  edge_type::node_set_type::iterator niter_end = edge_new.tail_nodes.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tail_nodes.begin(); niter != niter_end; ++ niter)
	    *niter += y_node_offset - (*niter > y.goal);
	}
      } else {
	const int y_node_offset = x.nodes.size();
	const int y_edge_offset = x.edges.size();
	
	// +1 to adjust goal
	x.nodes.resize(x.nodes.size() + y.nodes.size() + 1);
	// +2 to adjust edges toward goal
	x.edges.resize(x.edges.size() + y.edges.size() + 2);
	
	x.is_sorted = false;
	
	// merge nodes
	for (id_type id = 0; id < y.nodes.size(); ++ id) {
	  const id_type id_new = id + y_node_offset;
	  
	  const node_type& node_old = y.nodes[id];
	  node_type& node_new = x.nodes[id_new];
	  
	  node_new = node_old;
	    
	  node_new.id =  id_new;
	  node_type::edge_set_type::iterator eiter_end = node_new.in_edges.end();
	  for (node_type::edge_set_type::iterator eiter = node_new.in_edges.begin(); eiter != eiter_end; ++ eiter)
	    *eiter += y_edge_offset;
	  
	  eiter_end = node_new.out_edges.end();
	  for (node_type::edge_set_type::iterator eiter = node_new.out_edges.begin(); eiter != eiter_end; ++ eiter)
	    *eiter += y_edge_offset;
	}
	
	// merge edges
	for (id_type id = 0; id < y.edges.size(); ++ id) {
	  const edge_type& edge_old = y.edges[id];
	  edge_type& edge_new = x.edges[id + y_edge_offset];
	  
	  edge_new = edge_old;
	  edge_new.id = id + y_edge_offset;
	  
	  edge_new.head_node += y_node_offset;
	  edge_type::node_set_type::iterator niter_end = edge_new.tail_nodes.end();
	  for (edge_type::node_set_type::iterator niter = edge_new.tail_nodes.begin(); niter != niter_end; ++ niter)
	    *niter += y_node_offset;
	}
	
	// create new node...
	node_type& goal_node = x.nodes.back();
	edge_type& goal_edge_x = x.edges[x.edges.size() - 2];
	edge_type& goal_edge_y = x.edges[x.edges.size() - 1];
	
	goal_node.id = x.nodes.size() - 1;
	goal_edge_x.id = x.edges.size() - 2;
	goal_edge_y.id = x.edges.size() - 1;

	goal_edge_x.rule.reset(new rule_type(vocab_type::GOAL,
					     rule_type::symbol_set_type(1, goal_x.non_terminal(1)),
					     rule_type::symbol_set_type(1, goal_x.non_terminal(1)),
					     1));
	goal_edge_y.rule.reset(new rule_type(vocab_type::GOAL,
					     rule_type::symbol_set_type(1, goal_y.non_terminal(1)),
					     rule_type::symbol_set_type(1, goal_y.non_terminal(1)),
					     1));
	
	goal_edge_x.head_node = goal_node.id;
	goal_edge_y.head_node = goal_node.id;
	
	goal_edge_x.tail_nodes = edge_type::node_set_type(1, x.goal);
	goal_edge_y.tail_nodes = edge_type::node_set_type(1, y.goal + y_node_offset);
	
	x.nodes[x.goal].out_edges.push_back(goal_edge_x.id);
	x.nodes[y.goal + y_node_offset].out_edges.push_back(goal_edge_y.id);
	
	goal_node.in_edges.push_back(goal_edge_x.id);
	goal_node.in_edges.push_back(goal_edge_y.id);
	
      }
    }
  };
  
  void HyperGraph::unite(const HyperGraph& x)
  {
    UnionImpl()(*this, x);
  }
    
  

  
};
