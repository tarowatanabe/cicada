// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_STATE_CHART__HPP__
#define __CICADA__NGRAM_STATE_CHART__HPP__ 1

//
// ngram state manager for chart
//

#include <algorithm>

#include <cicada/ngram_index.hpp>
#include <cicada/ngram_state.hpp>

namespace cicada
{
  struct NGramStateChart
  {
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    
    typedef Symbol             symbol_type;
    typedef Symbol             word_type;
    typedef float              logprob_type;

    typedef NGramIndex::state_type state_type;
    typedef NGramState suffix_state_type;
    
    NGramStateChart(const size_type order=1) : suffix_(order) {}
    
    size_type buffer_size() const
    {
      return offset_suffix() + suffix_.buffer_size();
    }

    size_type offset_prefix() const
    {
      return 0;
    }

    size_type offset_suffix() const
    {
      return sizeof(size_type) * 2 + sizeof(state_type) * (suffix_.order_ - 1);
    }

    void* suffix(void* buffer) const
    {
      return (char*) buffer + offset_suffix();
    }

    const void* suffix(const void* buffer) const
    {
      return (const char*) buffer + offset_suffix();
    }
    
    size_type& size_suffix(void* buffer) const
    {
      return suffix_.size(suffix(buffer));
    }
    
    const size_type& size_suffix(const void* buffer) const
    {
      return suffix_.size(suffix(buffer));
    }

    size_type& size_prefix(void* buffer) const
    {
      return *reinterpret_cast<size_type*>(buffer);
    }
    
    const size_type& size_prefix(const void* buffer) const
    {
      return *reinterpret_cast<const size_type*>(buffer);
    }
    
    size_type& complete(void* buffer) const
    {
      return *reinterpret_cast<size_type*>((char*) buffer + sizeof(size_type));
    }

    const size_type& complete(const void* buffer) const
    {
      return *reinterpret_cast<const size_type*>((const char*) buffer + sizeof(size_type));
    }

    state_type* state(void* buffer) const
    {
      return reinterpret_cast<state_type*>((char*) buffer + sizeof(size_type) * 2);
    }

    const state_type* state(const void* buffer) const
    {
      return reinterpret_cast<const state_type*>((const char*) buffer + sizeof(size_type) * 2);
    }
    
    word_type::id_type* context(void* buffer) const
    {
      return suffix_.context(suffix(buffer));
    }
    
    const word_type::id_type* context(const void* buffer) const
    {
      return suffix_.context(suffix(buffer));
    }

    logprob_type* backoff(void* buffer) const
    {
      return suffix_.backoff(suffix(buffer));
    }

    const logprob_type* backoff(const void* buffer) const
    {
      return suffix_.backoff(suffix(buffer));
    }
    
    void fill(void* buffer) const
    {
      const size_type len = size_prefix(buffer);
      
      std::fill(reinterpret_cast<char*>(state(buffer) + len), reinterpret_cast<char*>(state(buffer) + suffix_.order_ - 1), 0);
      
      suffix_.fill(suffix(buffer));
    }

    void copy(const void* buffer, void* copied)
    {
      std::copy((char*) buffer, ((char*) buffer) + buffer_size(), (char*) copied);
    }

    void copy_prefix(const void* buffer, void* copied)
    {
      std::copy((char*) buffer, ((char*) buffer) + offset_suffix(), (char*) copied);
    }

    void copy_suffix(const void* buffer, void* copied)
    {
      suffix_.copy(suffix(buffer), suffix(copied));
    }
    
    suffix_state_type suffix_;
  };

};

#endif
