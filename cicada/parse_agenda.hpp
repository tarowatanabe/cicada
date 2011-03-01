// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__PARSE_AGENDA__HPP__
#define __CICADA__PARSE_AGENDA__HPP__ 1

#include <vector>
#include <algorithm>
#include <set>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/lattice.hpp>
#include <cicada/grammar.hpp>
#include <cicada/transducer.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/semiring.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/b_heap.hpp>
#include <utils/std_heap.hpp>
#include <utils/bithack.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

namespace cicada
{
  
  // an agenda-based chart parsing algorithm processing in bottom-up fashion
  //
  // @INPROCEEDINGS{Klein01parsingand,
  //  author = {Dan Klein and Christopher D. Manning},
  //  title = {Parsing and Hypergraphs},
  //  booktitle = {IN IWPT},
  //  year = {2001},
  //  pages = {123--134},
  //  publisher = {}
  // }
  // 
  // One of the main reason for "bottom-up" is the grammar encoding represented in transducer class
  // which do not support Earlye-style lhs-first encoding and differentiating terminals/non-terminals
  //

  template <typename Semiring, typename Function>
  struct ParseAgenda
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;
    
    typedef Lattice    lattice_type;
    typedef Grammar    grammar_type;
    typedef Transducer transducer_type;
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::feature_set_type   feature_set_type;
    typedef hypergraph_type::attribute_set_type attribute_set_type;

    typedef attribute_set_type::attribute_type attribute_type;
    
    typedef hypergraph_type::rule_type     rule_type;
    typedef hypergraph_type::rule_ptr_type rule_ptr_type;

    typedef Semiring semiring_type;
    typedef Semiring score_type;
    
    typedef Function function_type;
    
    // since we need to differentiate by lhs, we need to check "pos"
    struct Dot
    {
      size_type                table;
      transducer_type::id_type node;
      size_type                pos;
      
      Dot() : table(size_type(-1)), node(), pos(size_type(-1)) {}
      Dot(const size_type& __table, const transducer_type::id_type& __node)
	: table(__table), node(__node), pos(size_type(-1)) {}
      Dot(const size_type& __table, const transducer_type::id_type& __node, const size_type& __pos)
	: table(__table), node(__node), pos(__pos) {}
    };
    
    typedef Dot dot_type;

    struct Span
    {
      int first;
      int last;
      
      Span() : first(0), last(0) {}
      Span(const int& __first, const int& __last) : first(__first), last(__last) {}
    };

    typedef Span span_type;

    struct Edge
    {
      typedef Edge edge_type;
      
      const edge_type* active;
      const edge_type* passive;
      
      // active edge
      dot_type                                  dot;
      rule_ptr_type                             rule;
      hypergraph_type::edge_type::node_set_type tails;
      feature_set_type                          features;
      attribute_set_type                        attributes;
      
      // passive edge
      symbol_type lhs;
      span_type   span;
      int         level;
      
    public:
      bool is_passive() const { return dot.pos != size_type(-1); }
      bool is_active() const { return dot.pos == size_type(-1); }
      
      bool is_scanned() const { return active && ! passive; }
      bool is_predicted() const { return span.first == span.last; }
      bool is_completed() const { return active && passive; }
    };
    
    typedef Edge edge_type;
    
    struct Traversal
    {
      const edge_type* active;
      const edge_type* passive;
      bool is_active;
      
      Traversal(const edge_type* __active, const edge_type* __passive, const bool __is_active)
	: active(__active), passive(__passive), is_active(__is_active) {}
      Traversal()
	: active(0), passive(0), is_active(0) {}
    };
    typedef Traversal traversal_type;
    
    typedef utils::hashmurmur<size_type> traversal_hash_type;
    struct traversal_equal_type
    {
      bool operator()(const traversal_type& x, const traversal_type& y) const
      {
	return x.passive == y.passive && x.active == y.active && x.is_active == y.is_active;
      }
    };
    
    typedef google::dense_hash_set<traversal_type, traversal_hash_type, traversal_equal_type > traversal_set_type;

    
    void complete_active(const edge_type& active)
    {
      const transducer_type& transducer = grammar[active.dot.table];

      const edge_type query(active.span.last, active.span.last);
      
      std::pair<edge_set_passive_type::const_iterator, edge_set_passive_type::const_iterator> result = edges_passive.equal_range(&query);
      for (edge_set_passive_type::const_iterator piter = result.first; piter != result.second; ++ piter) {
	const edge_type& passive = *(*piter);
	
	const transducer_type::id_type node = transducer.next(active.dot.node, non_terminal);
	if (node == transducer.root()) continue;
	
	const transducer_type::rule_pair_set_type& rules = transducer.rules(node);

	size_type pos = 0;
	transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	  
	  
	}
	
	if (transducer.has_next(node)) {
	  
	}
      }
    }

  };
};

#endif
