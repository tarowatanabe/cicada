// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/dependency.hpp>

#include <utils/compact_trie.hpp>
#include <utils/dense_hash_map.hpp>

struct DependencyCounts
{
  typedef size_t             size_type;
  typedef ptrdiff_t          difference_type;
  
  typedef uint64_t           count_type;
  typedef cicada::Symbol     word_type;
  typedef cicada::Vocab      vocab_type;
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Dependency dependency_type;
  
  typedef utils::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>,
				std::allocator<std::pair<const word_type, count_type> > > root_count_set_type;
  
  typedef utils::compact_trie<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>,
  			      std::allocator<std::pair<const word_type, count_type> > > trie_count_set_type;
  
  struct count_set_type
  {
    root_count_set_type root;
    trie_count_set_type left;
    trie_count_set_type right;
    
    count_set_type() { root.set_empty_key(word_type()); }
    
    void clear()
    {
      root.clear();
      left.clear();
      right.clear();
    }
  };
};
