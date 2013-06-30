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
#include <cicada/ngram_state.hpp>
#include <cicada/ngram_state_chart.hpp>

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
    
    // this is a returned result... POD
    struct Result
    {
      typedef NGram::logprob_type    logprob_type;
      typedef NGram::state_type      state_type;
      typedef NGram::size_type       size_type;
      typedef NGram::difference_type difference_type;
      
      state_type state;  // current state
      
      logprob_type prob;  // probability w/o backoff
      logprob_type bound; // upper bound w/o backoff
      
      uint32_t length;   // mathed ngram length
      uint32_t complete; // complete implies a complete ngram w/o further context

      Result() : state(), prob(0), bound(0), length(0), complete(false) {}
    };
    typedef Result result_type;

  public:
    NGram(const int _debug=0) : debug(_debug) { clear(); }
    NGram(const path_type& path, const int _debug=0) : debug(_debug) { open(path); }
    
  public:
    static const logprob_type logprob_min() { return boost::numeric::bounds<logprob_type>::lowest(); }
    static const logprob_type logprob_bos() { return -99.0 * M_LN10; }
    
  public:
    
    template <typename Word_>
    result_type ngram_score(const void* buffer_in, const Word_& word, void* buffer_out) const
    {
      NGramState ngram_state(index.order());
      
      result_type result = lookup(buffer_in, word, buffer_out);
      
      const size_type context_length = ngram_state.size(buffer_in);
      
      const logprob_type* biter     = ngram_state.backoff(buffer_in) + result.length + (result.length == 0) - 1;
      const logprob_type* biter_end = ngram_state.backoff(buffer_in) + context_length;
      
      for (/**/; biter < biter_end; ++ biter) {
	result.prob  += *biter;
      }
      
      return result;
    }
    
    result_type ngram_partial_score(const void* buffer_in, const state_type& state, const int order, void* buffer_out) const
    {
      NGramState ngram_state(index.order());
      
      result_type result = lookup_partial(buffer_in, state, order, buffer_out);
      
      const size_type context_length = ngram_state.size(buffer_in);
      
      const logprob_type* biter     = ngram_state.backoff(buffer_in) + utils::bithack::max(static_cast<int>(result.length + (result.length == 0) - order), 0);
      const logprob_type* biter_end = ngram_state.backoff(buffer_in) + context_length;
      
      for (/**/; biter < biter_end; ++ biter) {
	result.prob  += *biter;
      }
      
      return result;
    }

    template <typename Iterator>
    double ngram_score_update(Iterator first, Iterator last, int order) const
    {
      double adjust = 0.0;

      for (/**/; first != last; ++ first, ++ order) {
	if (first->is_root_node()) continue;
	
	const size_type shard_index = utils::bithack::branch(first->is_root_shard(), size_type(0), first->shard());
	const size_type shard_node = first->node();

	const logprob_type logprob = logprobs[shard_index](shard_node, order);
	
	adjust += logprob;
	adjust -= (! logbounds.empty() && shard_node < logbounds[shard_index].size()
		   ? logbounds[shard_index](shard_node, order)
		   : logprob);
      }
      
      return adjust;
    }
    
    template <typename Word_>
    result_type lookup(const void* buffer_in, const Word_& word, void* buffer_out) const
    {
      NGramState ngram_state(index.order());

      word_type::id_type* output         = ngram_state.context(buffer_out);
      logprob_type*       output_backoff = ngram_state.backoff(buffer_out);
      
      const result_type result = lookup(ngram_state.context(buffer_in),
					ngram_state.context(buffer_in) + ngram_state.size(buffer_in),
					word,
					output,
					output_backoff);
      
      ngram_state.size(buffer_out) = output - ngram_state.context(buffer_out);
      
      return result;
    }
    
    template <typename Iterator, typename Temp, typename Output, typename OutputBackoff>
    result_type lookup(Iterator rfirst, Iterator rend,
		       const Temp& word,
		       Output& output,
		       OutputBackoff& output_backoff) const
    {
      return lookup(rfirst, rend, index.vocab()[word], output, output_backoff);
    }
    
    template <typename Iterator, typename Output, typename OutputBackoff>
    result_type lookup(Iterator rfirst, Iterator rend,
		       const word_type::id_type& word,
		       Output& output,
		       OutputBackoff& output_backoff) const
    {
      result_type result;
      
      state_type state = index.next(state_type(), word);
      
      // non-unigram match...
      if (state.is_root_node()) {
	result.state    = state;
	result.prob     = smooth;
	result.bound    = smooth;
	result.length   = 0;
	result.complete = true;
	
	return result;
      }
      
      // we need a special handling for bigram/unigram here...!
      
      int order = 1;
      
      *output = word;
      ++ output;
      *output_backoff = backoffs[0](state.node(), order);
      ++ output_backoff;
      
      // at least we have unigram...
      
      // we will try to find out the longest matches...
      // Here, we do not check .shard(), since we already know we are working with bigram and higher...
      for (/**/; rfirst != rend; ++ rfirst) {
	const state_type state_next = index.next(state, *rfirst);
	
	if (state_next.is_root_node()) {
	  result.complete = true;
	  break;
	}
	
	++ order;
	
	// do we need to check whether it is possible to extend further...?
	if (order < index.order()) {
	  *output = *rfirst;
	  ++ output;
	  *output_backoff = backoffs[state_next.shard()](state_next.node(), order);
	  ++ output_backoff;
	}
	
	state = state_next;
      }
      
      const size_type shard_index = utils::bithack::branch(state.is_root_shard(), size_type(0), state.shard());
      const size_type shard_node = state.node();
      
      if (! result.complete)
	result.complete = (order == index.order() || ! index[shard_index].has_child(shard_node));
      
      lookup_result(state, order, result);
      
      return result;
    }

    // lookup ngram from partial state
    
    result_type lookup_partial(const void* buffer_in,
			       state_type state,
			       int order,
			       void* buffer_out) const
    {
      NGramState ngram_state(index.order());
      
      word_type::id_type* output         = ngram_state.context(buffer_out);
      logprob_type*       output_backoff = ngram_state.backoff(buffer_out);
      
      const result_type result = lookup_partial(ngram_state.context(buffer_in),
						ngram_state.context(buffer_in) + ngram_state.size(buffer_in),
						state,
						order,
						output,
						output_backoff);
      
      ngram_state.size(buffer_out) = output - ngram_state.context(buffer_out);
      
      return result;
    }
    
    template <typename Iterator, typename Output, typename OutputBackoff>
    result_type lookup_partial(Iterator rfirst, Iterator rend,
			       state_type state,
			       int order,
			       Output& output,
			       OutputBackoff& output_backoff) const
    {
      result_type result;
      
      if (order <= 0)
	throw std::runtime_error("invalid ngram state/length for partial lookup!");
      
      lookup_result(state, order, result);

      if (result.complete)
	throw std::runtime_error("we assume non-completed context");
            
      // we use bound score...
      //
      // TODO: we need to different lobounds and logprobs....
      //
      const logprob_type adjust = result.bound;
      
      // at least we have unigram...
      for (/**/; rfirst != rend; ++ rfirst) {
	const state_type state_next = index.next(state, *rfirst);
	
	if (state_next.is_root_node()) {
	  result.complete = true;
	  break;
	}
	
	++ order;
	
	// do we need to check whether it is possible to extend further...?
	if (order < index.order()) {
	  *output = *rfirst;
	  ++ output;
	  *output_backoff = backoffs[state_next.shard()](state_next.node(), order);
	  ++ output_backoff;
	}
	
	state = state_next;
      }
      
      const size_type shard_index = utils::bithack::branch(state.is_root_shard(), size_type(0), state.shard());
      const size_type shard_node = state.node();
      
      if (! result.complete)
	result.complete = (order == index.order() || ! index[shard_index].has_child(shard_node));
      
      lookup_result(state, order, result);
      
      // make an adjustment to the score...
      result.prob  -= adjust;
      result.bound -= adjust;
      
      return result;
    }
    
    template <typename Iterator>
    state_type lookup_context(Iterator first, Iterator last, void* buffer_out) const
    {
      NGramState ngram_state(index.order());
      
      word_type::id_type* output         = ngram_state.context(buffer_out);
      logprob_type*       output_backoff = ngram_state.backoff(buffer_out);
      
      const state_type state = lookup_context(first, last, output, output_backoff);
      
      ngram_state.size(buffer_out) = output - ngram_state.context(buffer_out);
      
      return state;
    }

    template <typename Iterator, typename Output, typename OutputBackoff>
    state_type lookup_context(Iterator first, Iterator last,
			      Output& output,
			      OutputBackoff& output_backoff) const
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      
      return lookup_context(first, last, output, output_backoff, value_type());
    }
    
    template <typename Iterator>
    logprob_type logbound(Iterator first, Iterator last) const
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      
      if (first == last) return 0.0;
      
      logbound(first, last, value_type());
    }    
    
    template <typename Iterator>
    logprob_type logprob(Iterator first, Iterator last) const
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      
      if (first == last) return 0.0;
      
      return logprob(first, last, value_type());
    }

  private:
    
    void lookup_result(const state_type& state, int order, result_type& result) const
    {
      if (state.is_root()) {
	result.state = state;
	result.prob     = smooth;
	result.bound    = smooth;
	result.length   = 0;
	
	return;
      }
      
      size_type shard_index = utils::bithack::branch(state.is_root_shard(), size_type(0), state.shard());
      size_type shard_node = state.node();
      
      result.state = state;
            
      while (1) {
	result.prob     = logprobs[shard_index](shard_node, order);
	result.bound    = result.prob;
	result.bound    = (! logbounds.empty() && shard_node < logbounds[shard_index].size()
			   ? logbounds[shard_index](shard_node, order)
			   : result.prob);
	result.length   = order;
	
	if (result.prob != logprob_min()) break;

	if (order == 1) {
	  // very strange, though...
	  result.prob   = smooth;
	  result.bound  = smooth;
	  result.length = 0;
	  break;
	}
	
	shard_node = index[shard_index].parent(shard_node);
	-- order;
	if (order == 1)
	  shard_index = 0;
      }
    }

    struct ExtractVocab
    {
      const vocab_type& vocab_;
      ExtractVocab(const vocab_type& vocab)  : vocab_(vocab) {}

      template <typename Word_>
      word_type::id_type operator()(const Word_& word) const { return vocab_[word]; }
    };

    struct ExtractId
    {
      word_type::id_type operator()(const word_type::id_type& id) const { return id; }
    };
    
    template <typename Iterator, typename Word_>
    logprob_type logbound(Iterator first, Iterator last, Word_) const
    {
      return logbound_dispatch(first, last, ExtractVocab(index.vocab()));
    }
    
    template <typename Iterator>
    logprob_type logbound(Iterator first, Iterator last, word_type::id_type) const
    {
      return logbound_dispatch(first, last, ExtractId());
    }

    template <typename Iterator, typename Word_>
    logprob_type logprob(Iterator first, Iterator last, Word_) const
    {
      return logprob_dispatch(first, last, ExtractVocab(index.vocab()));
    }
    
    template <typename Iterator>
    logprob_type logprob(Iterator first, Iterator last, word_type::id_type) const
    {
      return logprob_dispatch(first, last, ExtractId());
    }

    template <typename Iterator, typename Output, typename OutputBackoff, typename Word_>
    state_type lookup_context(Iterator first, Iterator last,
			      Output& output,
			      OutputBackoff& output_backoff,
			      Word_) const
    {
      return lookup_context_dispatch(first, last, output, output_backoff, ExtractVocab(index.vocab()));
    }
    
    template <typename Iterator, typename Output, typename OutputBackoff>
    state_type lookup_context(Iterator first, Iterator last,
			      Output& output,
			      OutputBackoff& output_backoff,
			      word_type::id_type) const
    {
      return lookup_context_dispatch(first, last, output, output_backoff, ExtractId());
    }

    template <typename Iterator, typename Output, typename OutputBackoff, typename Extract>
    state_type lookup_context_dispatch(Iterator first, Iterator last,
				       Output& output,
				       OutputBackoff& output_backoff,
				       Extract extractor) const
    {
      // clip context size...
      first = std::max(first, last - (index.order() - 1));
      
      state_type state;
      int order = 1;
      
      for (/**/; last != first; -- last, ++ order) {
	const word_type::id_type word = extractor(*(last - 1));
	const state_type state_next = index.next(state, word);
	
	if (state_next.is_root_node()) break;

	const size_type shard_index = utils::bithack::branch(state_next.is_root_shard(), size_type(0), state_next.shard());
	
	*output = word;
	++ output;
	*output_backoff = backoffs[shard_index](state_next.node(), order);
	++ output_backoff;
	
	state = state_next;
      }
      
      return state;
    }

    struct IgnoreIterator
    {
      IgnoreIterator& operator=(const logprob_type& x) { return *this; }
      IgnoreIterator& operator=(const word_type::id_type& x) { return *this; }
      IgnoreIterator& operator*()  { return *this; }
      IgnoreIterator& operator++() { return *this; }
    };

    template <typename Iterator, typename Extract>
    logprob_type logbound_dispatch(Iterator first, Iterator last, Extract extractor) const
    {
      typedef std::vector<char, std::allocator<char> > buffer_type;
      
      NGramState ngram_state(index.order());
      
      buffer_type __buffer(ngram_state.buffer_size());
      void* buffer = &(*__buffer.begin());
      
      // first, lookup from last - 1 to first and fill-in the buffer...
      word_type::id_type* citer = ngram_state.context(buffer);
      logprob_type*       biter = ngram_state.backoff(buffer);

      lookup_context(first, last - 1, citer, biter);

      const size_type context_length = citer - ngram_state.context(buffer);
      
      // second, use lookup(rfirst, rend, word, output, output_backoff) with output/output-backoff ignored!
      IgnoreIterator output;
      IgnoreIterator output_backoff;
      
      const result_type result = lookup(ngram_state.context(buffer), ngram_state.context(buffer) + context_length,
					extractor(*(last - 1)),
					output,
					output_backoff);
      {
	const logprob_type* backoff = ngram_state.backoff(buffer);
	const logprob_type* biter = backoff + result.length - (result.length != 0);
	const logprob_type* blast = backoff + context_length;

	// if this is a complete context or backoff, use probability, not upper bound!
	if (biter < blast || result.complete) {
	  for (/**/; biter < blast; ++ biter)
	    result.prob += *biter;
	  return result.prob;
	} else
	  return result.bound;
      }
    }
    
    template <typename Iterator, typename Extract>
    logprob_type logprob_dispatch(Iterator first, Iterator last, Extract extractor) const
    {
      typedef std::vector<char, std::allocator<char> > buffer_type;
      
      NGramState ngram_state(index.order());
      
      buffer_type __buffer(ngram_state.buffer_size());
      void* buffer = &(*__buffer.begin());
      
      // first, lookup from last - 1 to first and fill-in the buffer...
      
      word_type::id_type* citer = ngram_state.context(buffer);
      logprob_type*       biter = ngram_state.backoff(buffer);

      lookup_context(first, last - 1, citer, biter);

      const size_type context_length = citer - ngram_state.context(buffer);
      
      // second, use lookup(rfirst, rend, word, output, output_backoff) with output/output-backoff ignored!
      IgnoreIterator output;
      IgnoreIterator output_backoff;
      
      result_type result = lookup(ngram_state.context(buffer), ngram_state.context(buffer) + context_length,
				  extractor(*(last - 1)),
				  output,
				  output_backoff);
      
      const logprob_type* backoff = ngram_state.backoff(buffer);
      for (const logprob_type* biter = backoff + result.length - (result.length != 0); biter < backoff + context_length; ++ biter)
	result.prob += *biter;
      
      return result.prob;
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
