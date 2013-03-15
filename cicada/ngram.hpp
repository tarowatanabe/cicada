// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM__HPP__
#define __CICADA__NGRAM__HPP__ 1

#include <stdint.h>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/ngram_index.hpp>

#include <boost/array.hpp>

#include <utils/packed_vector.hpp>
#include <utils/succinct_vector.hpp>
#include <utils/map_file.hpp>
#include <utils/mathop.hpp>

namespace cicada
{
  class NGram
  {
  public:
    typedef Symbol             word_type;
    typedef Vocab              vocab_type;
    
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    
    typedef word_type::id_type id_type;
    typedef uint64_t           count_type;
    typedef float              logprob_type;
    typedef double             prob_type;
    typedef uint8_t            quantized_type;
    
    typedef boost::filesystem::path path_type;

  public:
    struct ShardData
    {
    public:
      typedef utils::map_file<logprob_type, std::allocator<logprob_type> >                 logprob_set_type;
      typedef utils::packed_vector_mapped<quantized_type, std::allocator<quantized_type> > quantized_set_type;
      typedef boost::array<logprob_type, 256>                                              logprob_map_type;
      typedef std::vector<logprob_map_type, std::allocator<logprob_map_type> >             logprob_map_set_type;
      
      ShardData()
	: logprobs(), quantized(), maps(), offset(0) {}
      ShardData(const path_type& path)
      	: logprobs(), quantized(), maps(), offset(0) { open(path); }
      
      void open(const path_type& path);
      
      void populate()
      {
	logprobs.populate();
	quantized.populate();
      }

      path_type path() const { return (quantized.is_open() ? quantized.path().parent_path() : logprobs.path().parent_path()); }
      
      void close() { clear(); }
      void clear()
      {
	logprobs.clear();
	quantized.clear();
	maps.clear();
	offset = 0;
      }
      
      logprob_type operator()(size_type pos, int order) const
      {
	return (quantized.is_open() ? maps[order][quantized[pos - offset]] : logprobs[pos - offset]);
      }
      
      size_type size() const
      {
	return (quantized.is_open() ? quantized.size() + offset : logprobs.size() + offset);
      }

      bool is_quantized() const { return quantized.is_open(); }
            
      logprob_set_type     logprobs;
      quantized_set_type   quantized;
      logprob_map_set_type maps;
      size_type            offset;
    };
    
    typedef ShardData  shard_data_type;
    typedef std::vector<shard_data_type, std::allocator<shard_data_type> > shard_data_set_type;
    typedef NGramIndex shard_index_type;

    typedef shard_index_type::state_type state_type;
    
  public:
    NGram(const int _debug=0) : debug(_debug) { clear(); }
    NGram(const path_type& path, const int _debug=0) : debug(_debug) { open(path); }
    
  public:
    static const logprob_type logprob_min() { return boost::numeric::bounds<logprob_type>::lowest(); }
    static const logprob_type logprob_bos() { return -99.0 * M_LN10; }
    
  public:
    
    state_type root() const { return index.root(); }

    template <typename Iterator>
    std::pair<Iterator, Iterator> prefix(Iterator first, Iterator last) const
    {
      return index.prefix(first, last);
    }
    
    template <typename Iterator>
    state_type suffix(Iterator first, Iterator last) const
    {
      return index.suffix(first, last);
    }
    
    template <typename _Word>
    std::pair<state_type, logprob_type> logbound(state_type state, const _Word& word, bool backoffed=false, int max_order=0) const
    {
      return logbound(state, index.vocab()[word], backoffed, max_order);
    }
    
    std::pair<state_type, logprob_type> logbound(state_type state, const id_type& word, bool backoffed=false, int max_order=0) const
    {
      // returned state... maximum suffix of state + word, since we may forced backoff :-)
      max_order = utils::bithack::branch(max_order <= 0, index.order(), utils::bithack::min(index.order(), max_order));
      
      int order = index.order(state) + 1;
      for (/**/; order > max_order; -- order)
	state = index.prev(state);
      
      logprob_type backoff = 0.0;
      for (/**/; order; -- order) {
	const state_type result = index.next(state, word);
	
	if (! result.is_root_node()) {
	  const size_type shard_index = index.shard_index(result);
	  const logprob_type score = (! backoffed && result.node() < logbounds[shard_index].size()
				      ? logbounds[shard_index](result.node(), order)
				      : logprobs[shard_index](result.node(), order));
	  
	  if (score != logprob_min()) {
	    state = index.suffix(result);
	    for (int order_state = index.order(state) + 1; order_state > max_order; -- order_state)
	      state = index.prev(state);
	    
	    return std::make_pair(state, score + backoff);
	  }
	}
	
	if (state.is_root())
	  return std::make_pair(state_type(), index.is_bos(word) ? logprob_bos() : smooth + backoff);
	
	backoffed = true;
	backoff += backoffs[index.shard_index(state)](state.node(), order - 1);
	state = index.prev(state);
      }
      
      // we will not reach here...
      return std::make_pair(state_type(), index.is_bos(word) ? logprob_bos() : smooth + backoff);
    }
    
    template <typename _Word>
    std::pair<state_type, logprob_type> logprob(state_type state, const _Word& word, bool backoffed = false, int max_order = 0) const
    {
      return logprob(state, index.vocab()[word], backoffed, max_order);
    }
    
    std::pair<state_type, logprob_type> logprob(state_type state, const id_type& word, bool backoffed = false, int max_order = 0) const
    {
      // returned state... maximum suffix of state + word, since we may forced backoff :-)
      max_order = utils::bithack::branch(max_order <= 0, index.order(), utils::bithack::min(index.order(), max_order));
      
      int order = index.order(state) + 1;
      for (/**/; order > max_order; -- order)
	state = index.prev(state);
      
      logprob_type backoff = 0.0;
      for (/**/; order; -- order) {
	const state_type result = index.next(state, word);
	
	if (! result.is_root_node()) {
	  const logprob_type score = logprobs[index.shard_index(result)](result.node(), order);
	  
	  if (score != logprob_min()) {
	    state = index.suffix(result);
	    for (int order_state = index.order(state) + 1; order_state > max_order; -- order_state)
	      state = index.prev(state);
	    
	    return std::make_pair(state, score + backoff);
	  }
	}
	
	if (state.is_root())
	  return std::make_pair(state_type(), index.is_bos(word) ? logprob_bos() : smooth + backoff);
	
	backoffed = true;
	backoff += backoffs[index.shard_index(state)](state.node(), order - 1);
	state = index.prev(state);
      }
      
      // we will not reach here...
      return std::make_pair(state_type(), index.is_bos(word) ? logprob_bos() : smooth + backoff);
    }

    
    template <typename Iterator>
    logprob_type logbound(Iterator first, Iterator last, bool smooth_smallest=false) const
    {
      if (first == last) return 0.0;
      
      const int order = std::distance(first, last);
      
      if (order >= index.order())
	return logprob(first, last, smooth_smallest);
      
      if (order == 1) {
	const state_type result = index.next(state_type(), *(last - 1));
	
	if (! result.is_root_node()) {
	  const logprob_type logprob = (result.node() < logbounds[0].size()
					? logbounds[0](result.node(), 1)
					: logprobs[0](result.node(), 1));
	  
	  if (logprob != logprob_min())
	    return logprob;
	}
	
	return (index.is_bos(*(last - 1)) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth));
      } else {
	// search for the "suffix" history...
	int order_trie = 1;
	state_type state(index.shard_index_backward(first, last));
	
	for (Iterator iter = last - 1; iter != first; -- iter, ++ order_trie) {
	  const state_type result = index.next(state, *(iter - 1));
	  
	  if (result.is_root_node()) break;
	  
	  state = result;
	}
	
	if (order_trie <= 2)
	  state = state_type(size_type(-1), state.node());

	bool backoffed = (order_trie != order);
	logprob_type backoff = 0.0;
	
	for (/**/; order_trie; -- order_trie) {
	  const state_type result = index.next(state, *(last - 1));
	  
	  if (! result.is_root_node()) {
	    const size_type shard_index = index.shard_index(result);
	    const logprob_type score = (! backoffed && result.node() < logbounds[shard_index].size()
					? logbounds[shard_index](result.node(), order_trie)
					: logprobs[shard_index](result.node(), order_trie));
	    
	    if (score != logprob_min())
	      return score + backoff;
	  }
	  
	  if (state.is_root())
	    return (index.is_bos(*(last - 1)) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth + backoff));
	  
	  backoffed = true;
	  backoff += backoffs[index.shard_index(state)](state.node(), order_trie - 1);
	  state = index.prev(state);
	}
	
	return (index.is_bos(*(last - 1)) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth + backoff));
      }
    }
    
    template <typename Iterator>
    logprob_type operator()(Iterator first, Iterator last, bool smooth_smallest=false) const
    {
      return logprob(first, last, smooth_smallest);
    }
    
    template <typename Iterator>
    logprob_type logprob(Iterator first, Iterator last, bool smooth_smallest=false) const
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      
      return __logprob_dispatch(first, last, value_type(), smooth_smallest);
    }

  private:
    
    template <typename Iterator, typename _Word>
    logprob_type __logprob_dispatch(Iterator first, Iterator last, _Word __word, bool smooth_smallest=false) const
    {
      if (first == last) return 0.0;

      first = std::max(first, last - index.order());
      
      const int order = std::distance(first, last);
      const id_type word_id = index.vocab()[*(last - 1)];
      
      if (order == 1) {
	const state_type result = index.next(state_type(), word_id);

	if (! result.is_root_node()) {
	  const logprob_type logprob = logprobs[0](result.node(), 1);

	  if (logprob != logprob_min())
	    return logprob;
	}
	
	return (index.is_bos(word_id) ? logprob_bos() : smooth);
      } else {
	// search for the "suffix" history...
	int order_trie = 1;
	state_type state(index.shard_index_backward(first, last));
	
	for (Iterator iter = last - 1; iter != first; -- iter, ++ order_trie) {
	  const state_type result = index.next(state, *(iter - 1));

	  if (result.is_root_node()) break;
	  
	  state = result;
	}
	
	if (order_trie <= 2)
	  state = state_type(size_type(-1), state.node());

	logprob_type backoff = 0.0;
	for (/**/; order_trie; -- order_trie) {
	  // modify this for bigram case!
	  const state_type result = index.next(state, word_id);

	  if (! result.is_root_node()) {
	    const logprob_type score = logprobs[index.shard_index(result)](result.node(), order_trie);
	    
	    if (score != logprob_min())
	      return score + backoff;
	  }
	  
	  if (state.is_root())
	    return (index.is_bos(word_id) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth + backoff));
	  
	  backoff += backoffs[index.shard_index(state)](state.node(), order_trie - 1);
	  state = index.prev(state);
	}
	
	return (index.is_bos(word_id) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth + backoff));
      }
    }

    template <typename Iterator>
    logprob_type __logprob_dispatch(Iterator first, Iterator last, const id_type __word, bool smooth_smallest=false) const
    {
      if (first == last) return 0.0;

      first = std::max(first, last - index.order());
      
      const int order = std::distance(first, last);
      const id_type word_id = *(last - 1);
      
      if (order == 1) {
	const state_type result = index.next(state_type(), word_id);
	
	if (! result.is_root_node()) {
	  const logprob_type logprob = logprobs[0](result.node(), 1);

	  if (logprob != logprob_min())
	    return logprob;
	}
	
	return (index.is_bos(word_id) ? logprob_bos() : smooth);
      } else {
	// search for the "suffix" history...
	int order_trie = 1;
	state_type state(index.shard_index_backward(first, last));
	
	for (Iterator iter = last - 1; iter != first; -- iter, ++ order_trie) {
	  const state_type result = index.next(state, *(iter - 1));
	  
	  if (result.is_root_node()) break;
	  
	  state = result;
	}

	if (order_trie <= 2)
	  state = state_type(size_type(-1), state.node());

	logprob_type backoff = 0.0;
	for (/**/; order_trie; -- order_trie) {
	  const state_type result = index.next(state, word_id);
	  
	  if (! result.is_root_node()) {
	    const logprob_type score = logprobs[index.shard_index(result)](result.node(), order_trie);
	    
	    if (score != logprob_min())
	      return score + backoff;
	  }
	  
	  if (state.is_root())
	    return (index.is_bos(word_id) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth + backoff));
	  
	  backoff += backoffs[index.shard_index(state)](state.node(), order_trie - 1);
	  state = index.prev(state);
	}
	
	return (index.is_bos(word_id) ? logprob_bos() : (smooth_smallest ? logprob_min() : smooth + backoff));
      }
    }
    
    
  public:
    path_type path() const { return index.path().parent_path(); }
    size_type size() const { return index.size(); }
    bool empty() const { return index.empty(); }
    
    void open(const path_type& path);

    void populate()
    {
      index.populate();
      
      populate(logprobs.begin(), logprobs.end());
      populate(backoffs.begin(), backoffs.end());
      populate(logbounds.begin(), logbounds.end());
    }

    template <typename Iterator>
    void populate(Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first)
	first->populate();
    }
    
    void close() { clear(); }
    void clear()
    {
      index.clear();
      logprobs.clear();
      backoffs.clear();
      logbounds.clear();
      smooth = utils::mathop::log(1e-7);
    }
    
    bool is_open() const { return index.is_open(); }
    bool has_bounds() const { return ! logbounds.empty(); }

  public:
    static NGram& create(const path_type& path);
    
  public:
    shard_index_type    index;
    shard_data_set_type logprobs;
    shard_data_set_type backoffs;
    shard_data_set_type logbounds;
    
    logprob_type   smooth;
    int debug;
  };
  
};

#endif
