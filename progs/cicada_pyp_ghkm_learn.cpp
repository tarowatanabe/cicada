//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-GHKM!

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
#include "utils/chinese_restaurant_process.hpp"
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

