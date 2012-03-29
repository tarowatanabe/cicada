//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-SynAlign based on GHKM work:

// @InProceedings{cohn-blunsom:2009:EMNLP,
//   author    = {Cohn, Trevor  and  Blunsom, Phil},
//   title     = {A {Bayesian} Model of Syntax-Directed Tree to String Grammar Induction},
//   booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//   month     = {August},
//   year      = {2009},
//   address   = {Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {352--361},
//   url       = {http://www.aclweb.org/anthology/D/D09/D09-1037}
// }
//
//
// We do not perform composition, and each node is mapped to corresponding span in the target-side.
//

#include <map>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/tree_rule.hpp>
#include <cicada/tree_rule_compact.hpp>
#include <cicada/rule.hpp>

#include "utils/chunk_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/restaurant.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/dense_hash_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

typedef cicada::Vocab      vocab_type;
typedef cicada::Sentence   sentence_type;
typedef cicada::Symbol     word_type;
typedef cicada::HyperGraph hypergraph_type;

typedef cicada::Rule               rule_type;
typedef rule_type::rule_ptr_type   rule_ptr_type;
typedef rule_type::symbol_set_type symbol_set_type;

//
// syntactic alignment
//
// We assume one side is syntactic, and the other is hiero-like string
//
// Thus, it is possible to induce a synchronous-CFG by mapping non-temrinal symbols to the string-side
// Or, keep them separate like a synchronous-TSG.
//

struct PYPSynAlign
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  struct rule_pair_type
  {
    const rule_type* lhs;
    symbol_set_type  rhs;
    
    rule_pair_type() : lhs(0), rhs() {}
    rule_pair_type(const rule_ptr_type& __lhs, const symbol_set_type& __rhs) : lhs(&(*__lhs)), rhs(__rhs) {}
  };
  
  
  
  
};
