// -*- mode: c++ -*-

#ifndef __CICADA__LEXICALIZED_REORDERING__HPP__
#define __CICADA__LEXICALIZED_REORDERING__HPP__ 1

#include <stdint.h>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/symbol_vector.hpp>

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
    
    typedef Symbol  word_type;
    typedef Vocab   vocab_type;
    
    typedef SymbolVector                        phrase_type;
    typedef std::pair<phrase_type, phrase_type> phrase_pair_type;
    typedef std::pair<phrase_type, size_type>   phrase_node_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;

    typedef float          score_type;
    typedef uint8_t        quantized_type;
    typedef uint32_t       id_type;

    typedef utils::simple_vector<score_type, std::allocator<score_type> > feature_set_type;
    
  private:
    typedef word_type::id_type key_type;
    typedef id_type            mapped_type;
    
    typedef std::allocator<std::pair<key_type, mapped_type> >  phrase_alloc_type;
    typedef succinctdb::succinct_trie_db<key_type, mapped_type, phrase_alloc_type > phrase_db_type;
    
    class ScoreSet
    {
    public:
      typedef utils::map_file<score_type, std::allocator<score_type> >                     score_set_type;
      typedef utils::packed_vector_mapped<quantized_type, std::allocator<quantized_type> > quantized_set_type;
      typedef boost::array<score_type, 256>                                                score_map_type;
      
      ScoreSet() {}
      ScoreSet(const path_type& path) { read(path); }
      
      void read(const path_type& path);
      void write(const path_type& file) const;
      
      void clear()
      {
	score.clear();
	quantized.clear();
      }
      void close() { clear(); }
      
      score_type operator[](size_type pos) const
      {
	return (quantized.is_open()
		? maps[quantized[pos]]
		: score[pos]);
      }
      
      path_type path() const
      {
	return (quantized.is_open()
		? quantized.path().parent_path()
		: score.path().parent_path());
      }
      bool empty() const { return quantized.empty() && score.empty(); }
      size_type size() const
      {
	return (quantized.is_open()
		? quantized.size()
		: score.size());
      }
      
      score_set_type     score;
      quantized_set_type quantized;
      score_map_type     maps;
    };
    
    typedef ScoreSet score_set_type;
    typedef std::vector<score_set_type, std::allocator<score_set_type> > score_db_type;

    template <typename Tp>
    struct __cache
    {
      Tp        value;
      size_type index;
      
      __cache() : value(), index(size_type(-1)) {}
    };

    typedef __cache<feature_set_type> cache_feature_type;
    typedef __cache<phrase_type>      cache_phrase_type;
    typedef __cache<phrase_node_type> cache_phrase_node_type;

    typedef utils::array_power2<cache_feature_type,     1024 * 16, std::allocator<cache_feature_type> >     cache_feature_set_type;
    typedef utils::array_power2<cache_phrase_type,      1024 * 8,  std::allocator<cache_phrase_type> >      cache_phrase_set_type;
    typedef utils::array_power2<cache_phrase_node_type, 1024 * 16, std::allocator<cache_phrase_node_type> > cache_phrase_node_set_type;
    
  public:
    
    LexicalizedReordering() { clear(); }
    LexicalizedReordering(const std::string& parameter) { open(parameter); }
    LexicalizedReordering(const LexicalizedReordering& x)
      : phrases(x.phrases),
	scores(x.scores),
	vocab(x.vocab),
	fe(x.fe),
	bidirectional(x.bidirectional),
	monotonicity(x.monotonicity) {}
    
    LexicalizedReordering& operator=(const LexicalizedReordering& x)
    {
      clear();
      
      phrases = x.phrases;
      scores  = x.scores;
      vocab   = x.vocab;
      
      fe            = x.fe;
      bidirectional = x.bidirectional;
      monotonicity  = x.monotonicity;

      return *this;
    }
    
  public:
    const feature_set_type& operator[](size_type pos) const
    {
      if (! phrases.is_valid(pos) || ! phrases.exists(pos)) {
	if (none_features.empty())
	  const_cast<feature_set_type&>(none_features).resize(scores.size(), 0.0);
	
	return none_features;
      } else {
	const size_type cache_pos = hasher_type::operator()(pos) & (cache_features.size() - 1);
	cache_feature_type& cache = const_cast<cache_feature_type&>(cache_features[cache_pos]);
	
	if (cache.index != pos) {
	  cache.value.resize(scores.size());
	  for (size_t feature = 0; feature < scores.size(); ++ feature)
	    cache.value[feature] = scores[feature][pos];
	  
	  cache.index = pos;
	}
	
	return cache.value;
      }
    }

    template <typename Iterator>
    size_type find(Iterator first, Iterator last) const
    {
      const size_type cache_pos = hasher_type::operator()(first, last, 0) & (cache_phrases.size() - 1);
      cache_phrase_type& cache = const_cast<cache_phrase_type&>(cache_phrases[cache_pos]);
      
      if (! equal_phrase(first, last, cache.value)) {
	cache.value.assign(first, last);
	
	size_type node = 0;
	for (Iterator iter = first; iter != last && phrases.is_valid(node); ++ iter) {
	  const word_type::id_type id = vocab[*iter];
	  node = phrases.find(&id, 1, node);
	}
	
	cache.index = node;
      }
      
      return cache.index;
    }
    
    template <typename IteratorSource, typename IteratorTarget>
    size_type find(IteratorSource first_source, IteratorSource last_source,
		   IteratorTarget first_target, IteratorTarget last_target) const
    {
      size_type node = find(first_source, last_source);
      
      if (! phrases.is_valid(node))
	return node;
      
      const size_type cache_pos = hasher_type::operator()(first_target, last_target, node) & (cache_phrase_nodes.size() - 1);
      cache_phrase_node_type& cache = const_cast<cache_phrase_node_type&>(cache_phrase_nodes[cache_pos]);
      if (cache.value.second != node || ! equal_phrase(first_target, first_target, cache.value.first)) {
	
	cache.value.first.assign(first_target, last_target);
	cache.value.second = node;
	
	const word_type::id_type id = vocab[vocab_type::EMPTY];
	node = phrases.find(&id, 1, node);
	
	if (phrases.is_valid(node)) {
	  for (IteratorTarget iter = first_target; iter != last_target && phrases.is_valid(node); ++ iter) {
	    const word_type::id_type id = vocab[*iter];
	    node = phrases.find(&id, 1, node);
	  }
	}
	
	cache.index = node;
      }
      
      return cache.index;
    }
    
  private:
    
    template <typename Iterator, typename __Phrase>
    inline
    bool equal_phrase(Iterator first, Iterator last, const __Phrase& x) const
    {
      return static_cast<int>(x.size()) == std::distance(first, last) && std::equal(first, last, x.begin());
    }

  public:
    void open(const std::string& parameter);
    void write(const path_type& file) const;
    void quantize();
    
    void close() { clear(); }
    void clear()
    {
      none_features.clear();
      cache_features.clear();
      cache_phrases.clear();
      cache_phrase_nodes.clear();
      
      phrases.clear();
      scores.clear();
		    
      vocab.clear();

      fe = false;
      bidirectional = false;
      monotonicity = false;
    }
    
    path_type path() const { return phrases.path().parent_path(); }
    bool empty() const { return phrases.empty(); }

  private:
    feature_set_type           none_features;
    cache_feature_set_type     cache_features;
    cache_phrase_set_type      cache_phrases;
    cache_phrase_node_set_type cache_phrase_nodes;

    phrase_db_type phrases;
    score_db_type  scores;

    vocab_type vocab;
    
    bool fe;
    bool bidirectional;
    bool monotonicity;
  };
  
};

#endif
