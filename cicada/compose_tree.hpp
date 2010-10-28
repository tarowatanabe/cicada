// -*- mode: c++ -*-

#ifndef __CICADA__COMPOSE_TREE__HPP__
#define __CICADA__COMPOSE_TREE__HPP__ 1

#include <vector>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/tree_grammar.hpp>
#include <cicada/tree_transducer.hpp>
#include <cicada/hypergraph.hpp>

#include <utils/chunk_vector.hpp>
#include <utils/chart.hpp>
#include <utils/hashmurmur.hpp>

#include <google/dense_hash_map>

//
// hypergraph-to-hypergraph transduction by
//
//@InProceedings{zhang-EtAl:2009:EMNLP1,
//    author    = {Zhang, Hui  and  Zhang, Min  and  Li, Haizhou  and  Tan, Chew Lim},
//    title     = {Fast Translation Rule Matching for Syntax-based Statistical Machine Translation},
//    booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//    month     = {August},
//    year      = {2009},
//    address   = {Singapore},
//    publisher = {Association for Computational Linguistics},
//    pages     = {1037--1045},
//    url       = {http://www.aclweb.org/anthology/D/D09/D09-1108}
//}
//
// The terminologies used in the above paper is somewhat confusing:
//   input-forest, input-hypergraph, encoded hyper-path etc.
//
// The algorithm computes:
//  for each node, try matching with tree fragments
//  when matched, the match is represented by a set of hypergraph-node-id of input-hypergraph
//  and a set of matched rules.
//  We can uncover output hypergraph by enumerating matched rules.
//  Book-keep the matched input-hypergraph in a chart so that we can construct translational packed forest.
// 

namespace cicada
{
  
  struct ComposeTree
  {
    typedef Symbol symbol_type;
    typedef Vocab  vocab_type;

    typedef TreeGrammar    grammar_type;
    typedef TreeTransducer transducer_type;
    typedef HyperGraph     hypergraph_type;
    
    typedef hypergraph_type::feature_set_type feature_set_type;
    typedef hypergraph_type::rule_type        rule_type;
    typedef hypergraph_type::rule_ptr_type    rule_ptr_type;

    typedef transducer_type::rule_pair_set_type tree_rule_pair_set_type;
    typedef transducer_type::rule_pair_type     tree_rule_pair_type;
    typedef transducer_type::rule_type          tree_rule_type;
    typedef transducer_type::rule_ptr_type      tree_rule_ptr_type;
    
    ComposeTree(const grammar_type& __grammar)
      : grammar(__grammar) {}
    
    void operator()(const hypergraph_type& graph_in, hypergraph_type& graph_out)
    {
      
      
      
    }
    
    const grammar_type& grammar;
  };
  
  
  inline
  void compose_tree(const Grammar& grammar, const HyperGraph& graph_in, HyperGraph& graph_out)
  {
    ComposeCKY __composer(grammar);
    __composer(graph_in, graph_out);
  }
};

#endif
