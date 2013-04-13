// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_INDEX__HPP__
#define __CICADA__NGRAM_INDEX__HPP__ 1

// NGramIndex structure shared by NGram and NGramCounts
// Actually, the difference is the data associated with index:
// NGramCounts has conts and counts-modified
// NGram has logprob, backoff, lobbound


#include <stdint.h>

#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/packed_vector.hpp>
#include <utils/succinct_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/hashmurmur3.hpp>
#include <utils/spinlock.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  
  class NGramIndex
  {
  public:
    typedef Symbol                  word_type;
    typedef Vocab                   vocab_type;
    
    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef word_type::id_type      id_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;

  public:
    // root:    is_root_shard() && is_root_node()
    // unigram: is_root_shard() && ! is_root_node()
    // bigram, trigram etc: ! is_root_shard() && ! is_root_node()

    struct State
    {
      typedef uint64_t state_type;
      
      State() : state(state_type(-1)) {}
      State(const state_type shard)
	: state(((shard & 0xffff) << 48) | (state_type(-1) & 0xffffffffffffll)) {}
      State(const state_type shard, const state_type node)
	: state(((shard & 0xffff) << 48) | (node & 0xffffffffffffll)) {}
      
      size_type shard() const
      {
	return utils::bithack::branch(((state >> 48) & 0xffff) == 0xffff, size_type(-1), size_type((state >> 48) & 0xffff));
      }
      
      size_type node() const
      {
	return utils::bithack::branch((state & 0xffffffffffffll) == 0xffffffffffffll, size_type(-1), size_type(state & 0xffffffffffffll));
      }
      
      bool is_root() const { return state == state_type(-1); }
      bool is_root_shard() const { return ((state >> 48) & 0xffff) == 0xffff; }
      bool is_root_node() const { return (state & 0xffffffffffffll) == 0xffffffffffffll; } 

      const state_type& value() const { return state; }
      
      friend
      bool operator==(const State& x, const State& y) { return x.state == y.state; }
      friend
      bool operator!=(const State& x, const State& y) { return x.state != y.state; }
      friend
      bool operator<(const State& x, const State& y) { return x.state < y.state; }
      friend
      bool operator>(const State& x, const State& y) { return x.state > y.state; }
      friend
      bool operator<=(const State& x, const State& y) { return x.state <= y.state; }
      friend
      bool operator>=(const State& x, const State& y) { return x.state >= y.state; }
      
      friend
      size_t  hash_value(State const& x) { return utils::hashmurmur3<size_t>()(x.state); }
      
    private:
      state_type state;
    };

    typedef State state_type;

    struct Shard : public utils::hashmurmur3<size_t>
    {
    private:
      typedef utils::hashmurmur3<size_t> hasher_type;
      
    public:
      typedef utils::packed_vector_mapped<id_type, std::allocator<id_type> >   id_set_type;
      typedef utils::succinct_vector_mapped<std::allocator<int32_t> >          position_set_type;
      typedef std::vector<size_type, std::allocator<size_type> >               off_set_type;

    public:
      typedef utils::spinlock spinlock_type;
      typedef spinlock_type::scoped_lock     lock_type;
      typedef spinlock_type::scoped_try_lock trylock_type;
      
      struct cache_pos_type
      {
	size_type pos;
	size_type pos_next;
	id_type id;
        
	cache_pos_type() : pos(size_type(-1)), pos_next(size_type(-1)), id(id_type(-1)) {}
      };
      
      struct cache_suffix_type
      {
	state_type state;
	state_type suffix;
	
	cache_suffix_type() : state(), suffix() {}
      };
      
      typedef std::vector<cache_pos_type,    std::allocator<cache_pos_type> >    cache_pos_set_type;
      typedef std::vector<cache_suffix_type, std::allocator<cache_suffix_type> > cache_suffix_set_type;
      
    public:
      Shard() {}
      Shard(const path_type& path) { open(path); }
      
      Shard(const Shard& x) : ids(x.ids), positions(x.positions), offsets(x.offsets) { clear_cache(); }
      Shard& operator=(const Shard& x)
      {
	ids = x.ids;
	positions = x.positions;
	offsets = x.offsets;
	
	clear_cache();
	
	return *this;
      }
      
    public:
      void close() { clear(); }
      void clear()
      {
	ids.clear();
	positions.clear();
	offsets.clear();
	
	clear_cache();
      };
      
      void clear_cache();

      void open(const path_type& path);
      
      void populate()
      {
	ids.populate();
	positions.populate();
      }

    public:
      id_type operator[](size_type pos) const { return index(pos); }
      id_type index(size_type pos) const { return (pos < offsets[1] ? id_type(pos) : ids[pos - offsets[1]]); }
      size_type position_size() const { return offsets[offsets.size() - 2]; }
      size_type size() const { return offsets.back(); }
      bool empty() const { return offsets.empty(); }
      path_type path() const { return ids.path().parent_path(); }
      
      size_type parent(size_type pos) const
      {
	const size_type offset = offsets[1];
	
	return (pos < offset ? size_type(-1) : positions.select(pos + 1 - offset, true) + (offset + 1) - pos - 1);
      }

      bool has_child(size_type pos) const
      {
	return children_first(pos) != children_last(pos);
      }
      
      size_type children_first(size_type pos) const
      {
	if (pos == size_type(-1) || pos == 0)
	  return (~pos) & offsets[1];
	else
	  return children_last(pos - 1);
      }
      
      size_type children_last(size_type pos) const
      {
	const size_type offset = offsets[1];
	
	if (pos == size_type(-1))
	  return offset;
	else if (pos >= position_size())
	  return size();
	else {
	  position_set_type::size_type last = positions.select(pos + 2 - 1, false);
	  
	  return utils::bithack::branch(last == position_set_type::size_type(-1), size(), (last + 1 + offset + 1) - (pos + 2));
	  
	  //return (last == position_set_type::size_type(-1) ? size() : (last + 1 + offset + 1) - (pos + 2));
	}
      }
      
      template <typename Iterator>
      std::pair<Iterator, size_type> traverse(Iterator first, Iterator last, const vocab_type& vocab) const
      {
	typedef typename std::iterator_traits<Iterator>::value_type value_type;
	return __traverse_dispatch(first, last, vocab, value_type());
      }
      
      size_type __find(size_type pos, const id_type& id) const
      {
	const size_type pos_first = children_first(pos);
	const size_type pos_last  = children_last(pos);
	
	if (pos_first == pos_last) return size_type(-1);
	
	const size_type child = lower_bound(pos_first, pos_last, id);
	return utils::bithack::branch(child != pos_last && !(id < operator[](child)), child, size_type(-1));
      }
            
      size_type find(size_type pos, const id_type& id) const
      {
	if (id == id_type(-1))
	  return size_type(-1);
	else if (pos == size_type(-1))
	  return utils::bithack::branch(id < offsets[1], size_type(id), size_type(-1));
	else {
	  // if id <= 21bit range and last <= 21bit range
	  //   use simple caching
	  // else
	  //   lock-base caching
	  
	  // lock here!
	  trylock_type lock(const_cast<spinlock_type&>(spinlock_pos));
	  
	  if (lock) {
	    const size_type cache_pos = hasher_type::operator()(id, pos) & (caches_pos.size() - 1);
	    cache_pos_type& cache = const_cast<cache_pos_type&>(caches_pos[cache_pos]);
	    if (cache.pos != pos || cache.id != id) {
	      cache.pos      = pos;
	      cache.id       = id;
	      cache.pos_next = __find(pos, id);
	    }
	    return cache.pos_next;
	  } else
	    return __find(pos, id);
        }
      }
      
      size_type lower_bound(size_type first, size_type last, const id_type& id) const
      {
	const size_type offset = offsets[1];

	if (last <= offset)
	  return utils::bithack::min(size_type(id), last); // unigram!
	else {
	  // otherwise...
	  size_type length = last - first;
	  first -= offset;
	  last  -= offset;
	  
	  if (length <= 32) {
	    for (/**/; first != last && ids[first] < id; ++ first);
	    return first + offset;
	  } else {
	    while (length) {
	      const size_t half  = length >> 1;
	      const size_t middle = first + half;
	      
	      const bool is_less = ids[middle] < id;
	      
	      first  = utils::bithack::branch(is_less, middle + 1, first);
	      length = utils::bithack::branch(is_less, length - half - 1, half);
	    }
	    return first + offset;
	  }
	}
      }
      
      template <typename Iterator, typename _Word>
      std::pair<Iterator, size_type> __traverse_dispatch(Iterator first, Iterator last, const vocab_type& vocab, _Word) const
      {
	size_type pos = size_type(-1);
	for (/**/; first != last; ++ first) {
	  const size_type node = find(pos, vocab[word_type(*first)]);
	  
	  if (node == size_type(-1))
	    return std::make_pair(first, pos);
	  pos = node;
	}
	return std::make_pair(first, pos);
      }
      
      template <typename Iterator>
      std::pair<Iterator, size_type> __traverse_dispatch(Iterator first, Iterator last, const vocab_type& vocab, id_type) const
      {
	size_type pos = size_type(-1);
	for (/**/; first != last; ++ first) {
	  const size_type node = find(pos, *first);
	  
	  if (node == size_type(-1))
	    return std::make_pair(first, pos);
	  pos = node;
	}
	return std::make_pair(first, pos);
      }
      
    public:
      id_set_type        ids;
      position_set_type  positions;
      off_set_type       offsets;
      
      spinlock_type          spinlock_pos;
      spinlock_type          spinlock_suffix;
      
      cache_pos_set_type     caches_pos;
      cache_suffix_set_type  caches_suffix;
    };
    
    typedef Shard shard_type;
    typedef std::vector<shard_type, std::allocator<shard_type> > shard_set_type;
    
    typedef shard_set_type::const_iterator  const_iterator;
    typedef shard_set_type::iterator              iterator;
    
    typedef shard_set_type::const_reference const_reference;
    typedef shard_set_type::reference             reference;
    
  public:
    NGramIndex() {}
    NGramIndex(const path_type& path) { open(path); }
    
  public:
    state_type root() const { return state_type(); }
    
    template <typename Iterator>
    std::pair<state_type, Iterator> next(state_type state, Iterator first, Iterator last) const
    {
      for (/**/; first != last; ++ first) {
	const state_type state_next = next(state, *first);
	if (state_next.is_root_node())
	  return std::make_pair(state, first);
	
	state = state_next;
      }
      
      return std::make_pair(state, first);
    }

    template <typename _Word>
    state_type next(const state_type& state, const _Word& word) const
    {
      return next(state, __vocab[word]);
    }
    
    state_type next(const state_type& state, const id_type& word) const
    {
      if (state.is_root())
	return state_type(state.shard(), __shards[0].find(state.node(), word));
      else {
	if (state.is_root_node())
	  throw std::runtime_error("invalid state");
	
	if (! state.is_root_shard())
	  return state_type(state.shard(), __shards[state.shard()].find(state.node(), word));
	else {
	  // state.node() is equal to unigram's id
	  const size_type shard = shard_index(state.node(), word);
	  
	  return state_type(shard, __shards[shard].find(state.node(), word));
	}
      }
    }

    int order(const state_type& state) const
    {
      if (state.is_root())
	return 0;
      else if (state.is_root_shard())
	return 1;
      else {
	const shard_type& shard = __shards[state.shard()];
	const size_type node = state.node();
	
	size_type order = 2;
	for (/**/; order < shard.offsets.size(); ++ order)
	  if (node < shard.offsets[order])
	    return order;
	return order;
      }
    }

    template <typename Iterator>
    std::pair<Iterator, Iterator> prefix(Iterator first, Iterator last) const
    {
      last = std::min(last, first + order() - 1);
      
      if (first == last || first + 1 == last) return std::make_pair(first, last);
      
      return std::make_pair(first, next(state_type(), first, last - 1).second + 1);
    }

    template <typename Iterator>
    state_type suffix(Iterator first, Iterator last) const 
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      
      return __suffix_dispatch(first, last, value_type());
    }
    
    template <typename Iterator, typename __Word>
    state_type __suffix_dispatch(Iterator first, Iterator last, __Word) const
    {
      typedef std::vector<id_type, std::allocator<id_type> > buffer_type;
      
      const size_type length = std::distance(first, last);
      
      if (length == 0)
	return state_type();
      else if (length == 1) {
	const state_type state = next(state_type(), *first);
	return (state.is_root_node() ? state_type() : state);
      }
      
      buffer_type buffer(length);
      for (typename buffer_type::iterator biter = buffer.begin(); first != last; ++ first, ++ biter)
	*biter = __vocab[*first];
      
      return __suffix_dispatch(buffer.begin(), buffer.end(), id_type());
    }

    template <typename Iterator>
    state_type __suffix_dispatch(Iterator first, Iterator last, id_type) const
    {
      const size_type length = std::distance(first, last);
      
      if (length == 0)
	return state_type();
      else if (length == 1) {
	const state_type state = next(state_type(), *first);
	return (state.is_root_node() ? state_type() : state);
      }
      
      first = std::max(first, last - order());
      
      state_type state;
      
      while (first != last) {
	std::pair<state_type, Iterator> result = next(state, first, last);
	
	if (result.second == last)
	  return result.first;
	else {
	  state = suffix(result.first);
	  first = result.second;
	}
      }
      
      return state_type();
    }
    
    state_type suffix(const state_type& state) const
    {
      typedef std::vector<id_type, std::allocator<id_type> > context_type;
      
      // root or unigram's suffix is root
      if (state.is_root() || state.is_root_shard())
	return state_type();
      
      shard_type& shard = const_cast<shard_type&>(__shards[state.shard()]);
      
      // if we are bigram, we will simply take root_shard + current id
      if (state.node() < shard.offsets[2])
	return state_type(size_type(-1), shard[state.node()]);
      
      const size_type cache_pos = hash_value(state) & (shard.caches_suffix.size() - 1);
      
      // trylock...
      {
	shard_type::trylock_type lock(shard.spinlock_suffix);
	
	if (lock && shard.caches_suffix[cache_pos].state == state)
	  return shard.caches_suffix[cache_pos].suffix;
      }

      // we will push in reverse order...
      context_type context(order());
      context_type::reverse_iterator riter = context.rbegin();
	
      {
	size_type node = state.node();
	for (/**/; node != size_type(-1); ++ riter) {
	  *riter = shard[node];
	  node = shard.parent(node);
	}
      }
      
      context_type::const_iterator first = riter.base() + 1;
      context_type::const_iterator last  = context.end();

      state_type state_suffix;
      for (/**/; first != last - 1; ++ first) {
	std::pair<state_type, context_type::const_iterator> result = next(state_type(), first, last);
	
	if (result.second == last) {
	  state_suffix = result.first;
	  break;
	}
      }

      // trylock...
      {
	shard_type::trylock_type lock(shard.spinlock_suffix);
	
	if (lock) {
	  shard.caches_suffix[cache_pos].state  = state;
	  shard.caches_suffix[cache_pos].suffix = state_suffix;
	}
      }
      
      return state_suffix;
    }

    size_type shard_index(const id_type& first, const id_type& second) const
    {
      return __hasher(first, __hasher(second, 0)) % __shards.size();
    }

    template <typename Iterator>
    size_type shard_index(Iterator first, Iterator last) const
    {
      if (std::distance(first, last) <= 1) return 0;
      
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      return __shard_index_dispatch(first, last, value_type());
    }
    
    template <typename Iterator>
    std::pair<Iterator, size_type> traverse(size_type shard, Iterator first, Iterator last) const
    {
      return __shards[shard].traverse(first, last, __vocab);
    }
        
    template <typename Iterator>
    std::pair<Iterator, size_type> traverse(Iterator first, Iterator last) const
    {
      return __shards[shard_index(first, last)].traverse(first, last, __vocab);
    }
    
    bool is_bos(const id_type& id) const
    {
      return __vocab[vocab_type::BOS] == id;
    }
    
    bool is_bos(const word_type& word) const
    {
      return vocab_type::BOS == word;
    }
    
    inline const_reference operator[](size_type pos) const { return __shards[pos]; }
    inline       reference operator[](size_type pos)       { return __shards[pos]; }
    
    inline const_iterator begin() const { return __shards.begin(); }
    inline       iterator begin()       { return __shards.begin(); }
    
    inline const_iterator end() const { return __shards.end(); }
    inline       iterator end()       { return __shards.end(); }

    inline const vocab_type& vocab() const { return __vocab; }
    inline       vocab_type& vocab()       { return __vocab; }
    
    size_type size() const { return __shards.size(); }
    bool empty() const { return __shards.empty(); }
    
    void reserve(size_type n) { __shards.reserve(n); }
    void resize(size_type n) { __shards.resize(n); }
    void clear()
    {
      __shards.clear();
      __vocab.clear();
      __order = 0;
      __path = path_type();
    }
    void close() { clear(); }
    
    void open(const path_type& path);

    void populate()
    {
      shard_set_type::iterator siter_end = __shards.end();
      for (shard_set_type::iterator siter = __shards.begin(); siter != siter_end; ++ siter)
	siter->populate();
      __vocab.populate();
    }

    bool is_open() const { return ! __shards.empty() && ! __path.empty(); }
    path_type path() const { return __path; }
    
    inline const int& order() const { return __order; }
        
    size_type ngram_size(int order) const
    {
      switch (order) {
      case 0: return 0;
      case 1: return __shards.front().offsets[1];
      default:
	size_type sum = 0;
	for (size_t shard = 0; shard < __shards.size(); ++ shard)
	  sum += __shards[shard].offsets[order] - __shards[shard].offsets[order - 1];
	return sum;
      }
    }
      
  private:
    
    template <typename Iterator, typename _Word>
    size_type __shard_index_dispatch(Iterator first, Iterator last, _Word) const
    {
      return shard_index(__vocab[word_type(*first)], __vocab[word_type(*(first + 1))]);
    }
    
    template <typename Iterator>
    size_type __shard_index_dispatch(Iterator first, Iterator last, id_type) const
    {
      return shard_index(*first, *(first + 1));
    }
        
  private:
    shard_set_type __shards;
    vocab_type     __vocab;
    hasher_type    __hasher;
    
    int            __order;
    path_type      __path;
  };
  
};

#endif
