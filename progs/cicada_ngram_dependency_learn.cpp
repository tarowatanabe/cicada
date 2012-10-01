//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// dependency ngram language model learning...
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/functional/hash.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <sstream>

#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/dependency.hpp>

#include <utils/compact_trie.hpp>
#include <utils/dense_hash_map.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/compress_stream.hpp>
#include <utils/spinlock.hpp>

// dependency counts format:
// 
// counts/
//   unigram/
//     vocab.gz
//     vocab_cs.gz
//
//   root/
//     counts-xxx.gz
//     ...
//     root.idx
//   left/
//     counts-xxx.gz
//     ...
//     left.idx
//
//   left-head/
//     counts-xxx.gz
//     ...
//     left-head.idx
//
//   right/
//     counts-xxx.gz
//     ...
//     right.idx
//
//   right-head/
//     counts-xxx.gz
//     ...
//     right-head.idx
//

// after model construction, we will re-sort so that
// ngrams are stored in suffix order

struct DependencyModel
{
  // we will learn model from counts...

  typedef size_t             size_type;
  typedef ptrdiff_t          difference_type;
  typedef uint32_t           id_type;
  
  typedef uint64_t           count_type;
  typedef float              logprob_type;
  typedef double             prob_type;
  typedef cicada::Symbol     word_type;
  typedef cicada::Vocab      vocab_type;
  typedef cicada::Sentence   ngram_type;

  struct IndexModel
  {
    struct parameter_type
    {
      count_type   count;
      count_type   types;
      logprob_type logprob;
      logprob_type backoff;
      
      parameter_type() : count(0), types(0), logrpob(0), backoff(0) {}
    };

    typedef utils::dense_hash_map<word_type, parameter_type, boost::hash<word_type>, std::equal_to<word_type>,
				  std::allocator<std::pair<word_type, parameter_type> > > count_set_type;

    struct node_type
    {
      count_set_type counts;
      id_type parent;
      
      node_type() : counts(), parent(id_type(-1)) {}
    };
    
    typedef utils::compact_trie_dense<word_type, node_type, boost::hash<word_type>, std::equal_to<word_type>,
				      std::allocator<std::pair<const word_type, node_type> > > count_trie_type;



    struct shard_type
    {
      typedef utils::spinlock spinlock_type;
      typedef spinlock_type::scoped_lock     lock_type;
      typedef spinlock_type::scoped_try_lock trylock_type;
      
      count_trie_type counts;
      spinlock_type   mutex;
    };
    
    typedef std::vector<shard_type, std::allocator<shard_type> > shard_set_type;
    
    count_set_type root;
    shard_set_type left;
    shard_set_type right;
    shard_set_type left_head;
    shard_set_type right_head;
    
    void clear()
    {
      root.clear();
      left.clear();
      right.clear();
      left_head.clear();
      right_head.clear();
    }
    
    static
    void increment(shard_type& sahrd, const ngram_type& ngram, const count_type& count)
    {
      if (ngram.empty()) return;
      
      count_trie_type::id_type id = shard.counts.root();
      
      ngram_type::const_reverse_iterator niter_end(ngram.rend());
      for (ngram_type::const_reverse_iterator niter(ngram.end() - 1); niter != niter_end; ++ niter) {
	if (id == shard.counts.root()) {
	  shard_type::lock_type lock(shard.mutex);
	  
	  const count_trie_type::id_type id_next = shard.counts.insert(id, *niter);
	  
	  shard.counts[id_next].parent = id;
	  id = id_next;
	} else {
	  const count_trie_type::id_type id_next = shard.counts.insert(id, *niter);
	  
	  shard.counts[id_next].parent = id;
	  id = id_next;
	}
      }
      
      shard.counts[id].counts[ngram.back()].count += count;
    }
  };
};

void dependency_ngram_read(const path_type& path, const char* name, model_type::shard_set_type& counts)
{
  
  
}

void dependency_ngram_read(const path_type& path, model_type& model)
{
  {
    // first, read vocabulary...
    
    const path_type path_vocab = path / "unigram" / "vocab_cs.gz";
    
    if (! boost::filesystem::exists(path_vocab))
      throw std::runtime_error("no sorted vocabulary file? " + path_vocab.string());

    utils::compress_istream is(path_vocab, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
    
    qi::uint_parser<count_type> count_parser;
    
    std::string word;
    count_type  count;
    while (iter != iter_end) {
      word.clear();
      if (! qi::phrase_parse(iter, iter_end,
			     qi::lexeme[+(standard::char - standard::blank)] >> count_parser >> (qi::eol || qi::eoi),
			     standard::blank,
			     word,
			     count))
	throw std::runtime_error("parsing failed");
      
      const word_type tmptmp(word);
    }
  }
  
  {
    // second, read root
    const path_type path_root = path / "root";
    const path_type path_root_idx = path_root / "root.idx";
    
    if (! boost::filesystem::exists(path_root))
      throw std::runtime_error("no root path? " + path_root.string());
    
    if (! boost::filesystem::exists(path_root_idx))
      throw std::runtime_error("no root index path? " + path_root_idx.string());

    qi::uint_parser<count_type> count_parser;
    
    utils::compress_istream is(path_root_idx);
    std::string name;
    
    while (std::getline(is, name)) {
      const path_type path_counts = path_root / name;
      
      if (! boost::filesystem::exists(path))
	throw std::runtime_error("no root counts? " + path_counts.string());
      
      utils::compress_istream is(path_counts, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iter_type iter(is);
      iter_type iter_end;
      
      std::string word;
      count_type  count;
      while (iter != iter_end) {
	word.clear();
	if (! qi::phrase_parse(iter, iter_end,
			       qi::lexeme[+(standard::char - standard::blank)] >> count_parser >> (qi::eol || qi::eoi),
			       standard::blank,
			       word,
			       count))
	  throw std::runtime_error("parsing failed");
	
	model.root[word].count += count;
      }
    }
  }
  
  // thrid, collect ngram counts...
  
  dependency_ngram_read(path, "left",  model.left);
  dependency_ngram_read(path, "right", model.right);
  dependency_ngram_read(path, "left-head",  model.left_head);
  dependency_ngram_read(path, "right-head", model.right_head);
}
