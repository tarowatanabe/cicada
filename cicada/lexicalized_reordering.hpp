// -*- mode: c++ -*-

#ifndef __CICADA__LEXICALIZED_REORDERING__HPP__
#define __CICADA__LEXICALIZED_REORDERING__HPP__ 1

#include <stdint.h>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <succinct_db/succinct_trie_db.hpp>

#include <boost/filesystem.hpp>

#include <utils/array_power2.hpp>
#include <utils/mathop.hpp>


namespace cicada
{

  class LexicalizedReordering : public utils::hashmurmur<uint64_t>
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Symbol             word_type;
    typedef Vocab              vocab_type;
    typedef word_type::id_type word_id_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef float weight_type;
    
    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;
    
  private:
    typedef word_id_type key_type;
    typedef weight_type  mapped_type;
    
    typedef std::allocator<std::pair<key_type, mapped_type> >  lexicon_alloc_type;
    typedef succinctdb::succinct_trie_db<key_type, mapped_type, lexicon_alloc_type > lexicon_type;
    
    struct NodeCache
    {
      word_type::id_type      word;
      lexicon_type::size_type node;
      
      NodeCache()
	: word(word_type::id_type(-1)), node(0) {}
    };
    typedef NodeCache node_cache_type;
    
    typedef utils::array_power2<node_cache_type, 1024 * 8, std::allocator<node_cache_type> > node_cache_set_type;
    
    struct LexiconCache
    {
      word_type::id_type word;
      lexicon_type::size_type prev;
      lexicon_type::size_type next;
      
      LexiconCache()
	: word(word_type::id_type(-1)), prev(0), next(0) {}
    };
    typedef LexiconCache lexicon_cache_type;
    typedef utils::array_power2<lexicon_cache_type, 1024 * 64, std::allocator<lexicon_cache_type> > lexicon_cache_set_type;

  public:
    
    LexicalizedReordering() { clear(); }
    LexicalizedReordering(const path_type& path) { open(path); }
    LexicalizedReordering(const LexicalizedReordering& x)
      : lexicon(x.lexicon),
	vocab(x.vocab) {}
    
    LexicalizedReordering& operator=(const LexicalizedReordering& x)
    {
      clear();
      
      lexicon = x.lexicon;
      vocab = x.vocab;
      return *this;
    }
    
  public:
    
    template <typename Iterator>
    weight_type operator()(const word_type& target, Iterator first, Iterator last) const
    {
      if (empty()) return 0.0;
      
      const lexicon_type::size_type node = find(target);
      
      double weight = 0.0;
      
      if (lexicon.is_valid(node)) {
	// bias term...
	const lexicon_type::size_type node_next = find(node, vocab_type::NONE);
	if (lexicon.is_valid(node_next) && lexicon.exists(node_next))
	  weight += lexicon[node_next];
	
	// others...
	for (/**/; first != last; ++ first) {
	  const lexicon_type::size_type node_next = find(node, *first);
	  
	  if (lexicon.is_valid(node_next) && lexicon.exists(node_next))
	    weight += lexicon[node_next];
	}
      }
      
      return - boost::math::log1p(std::exp(- weight));
    }

  private:
    lexicon_type::size_type find(const word_type& word) const
    {
      const size_type cache_pos = hash_value(word) & (cache_node.size() - 1);
      node_cache_type& cache = const_cast<node_cache_type&>(cache_node[cache_pos]);
      
      if (cache.word != word.id()) {
	const word_id_type word_id = vocab[word];
	
	cache.word = word.id();
	cache.node = lexicon.find(&word_id, 1);
      }
      
      return cache.node;
    }
    
    lexicon_type::size_type find(const lexicon_type::size_type& node, const word_type& word) const
    {
      const size_type cache_pos = hasher_type::operator()(word.id(), node) & (cache_lexicon.size() - 1);
      lexicon_cache_type& cache = const_cast<lexicon_cache_type&>(cache_lexicon[cache_pos]);
      if (cache.word != word.id() || cache.prev != node) {
	const word_id_type word_id = vocab[word];
	
	cache.word = word.id();
	cache.prev = node;
	cache.next = lexicon.find(&word_id, 1, node);
      }
      return cache.next;
    }

  public:
    void open(const path_type& path);
    void write(const path_type& file) const;
    
    void close() { clear(); }
    void clear()
    {
      cache_lexicon.clear();
      cache_node.clear();
      lexicon.clear();
      vocab.clear();
    }
    
    path_type path() const { return lexicon.path().parent_path(); }
    bool empty() const { return lexicon.empty(); }

  private:
    lexicon_cache_set_type cache_lexicon;
    node_cache_set_type    cache_node;

    lexicon_type lexicon;
    vocab_type   vocab;
  };
  
};

#endif
