// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HYPERGRAPH__HPP__
#define __CICADA__HYPERGRAPH__HPP__ 1

#include <iostream>
#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/rule.hpp>

#include <utils/simple_vector.hpp>
#include <utils/chunk_vector.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  
  class HyperGraph
  {
  public:
    typedef uint32_t       id_type;
    typedef cicada::Symbol symbol_type;
    typedef cicada::Vocab  vocab_type;
    
    typedef cicada::Rule                 rule_type;
    typedef boost::shared_ptr<rule_type> rule_ptr_type;

    typedef cicada::FeatureVector<double, std::allocator<double> > feature_set_type;
    
  public:
    static const id_type invalid;
    
  public:
    HyperGraph() : goal(invalid) {}
    
  public:
    struct Node
    {
      typedef std::vector<id_type, std::allocator<id_type> > edge_set_type;
      
      Node() : id(invalid) {}
      
      edge_set_type edges;
      
      id_type id;
    };
    
    struct Edge
    {
      typedef utils::simple_vector<id_type, std::allocator<id_type> > node_set_type;
      typedef cicada::Rule rule_type;
      
      Edge()
	: head(invalid), tails(), rule(), first(0), last(0) {}
      
      template <typename Iterator>
      Edge(Iterator first, Iterator last)
	: head(invalid), tails(first, last), rule(), first(0), last(0) {}
      
      id_type       head;
      node_set_type tails;
      
      feature_set_type features;
      
      rule_ptr_type rule;
      
      // meta data...
      int first;
      int last;
      
      id_type id;
    };
    
    typedef Node node_type;
    typedef Edge edge_type;
    
  public:
    typedef utils::chunk_vector<node_type, 4096 / sizeof(node_type), std::allocator<node_type> > node_set_type;
    typedef utils::chunk_vector<edge_type, 4096 / sizeof(edge_type), std::allocator<edge_type> > edge_set_type;

  public:

    edge_type& add_edge(const edge_type& edge)
    {
      const id_type edge_id = edges.size();
      
      edges.push_back(edge);
      edges.back().id = edge_id;
      
      return edges.back();
    }
    

    edge_type& add_edge()
    {
      const id_type edge_id = edges.size();
      
      edges.push_back(edge_type());
      edges.back().id = edge_id;
      
      return edges.back();
    }

    template <typename Iterator>
    edge_type& add_edge(Iterator first, Iterator last)
    {
      const id_type edge_id = edges.size();
      
      edges.push_back(edge_type(first, last));
      edges.back().id = edge_id;
            
      return edges.back();
    }
    
    node_type& add_node()
    {
      const id_type node_id = nodes.size();
      
      nodes.push_back(node_type());
      nodes.back().id = node_id;
      
      return nodes.back();
    }
    
    void connect_edge(const id_type edge, const id_type head)
    {
      edges[edge].head = head;
      nodes[head].edges.push_back(edge);
    };
    
    void clear()
    {
      edges.clear();
      nodes.clear();
      
      goal = invalid;
    }
    
    void swap(HyperGraph& x)
    {
      nodes.swap(x.nodes);
      edges.swap(x.edges);
      std::swap(goal, x.goal);
    }

    bool is_valid() const
    {
      return goal != invalid;
    }

  public:
    // algorithms...
    
    void topologically_sort();
    
    void unite(const HyperGraph& x);

  public:
    void assign(const std::string& x);
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);

    friend
    std::ostream& operator<<(std::ostream& os, const HyperGraph& x);
    friend
    std::istream& operator>>(std::istream& is, HyperGraph& x);
    
  public:
    node_set_type nodes;
    edge_set_type edges;
    
    id_type goal;
  };
  
};

namespace std
{
  inline
  void swap(cicada::HyperGraph& x, cicada::HyperGraph& y)
  {
    x.swap(y);
  }
  
};

#endif
