// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_STATE__HPP__
#define __CICADA__NGRAM_STATE__HPP__ 1

//
// ngram state manager
//

#include <algorithm>

#include <cicada/symbol.hpp>

namespace cicada
{
  struct NGramState
  {
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    
    typedef Symbol             symbol_type;
    typedef Symbol             word_type;
    typedef float              logprob_type;
    
    NGramState(const size_type order) : order_(order) {}
    
    size_type order_;
    
    size_type buffer_size() const
    {
      return sizeof(size_type) + (sizeof(word_type::id_type) + sizeof(logprob_type)) * (order_ - 1);
    }
    
    size_type& size(void* buffer) const
    {
      return *reinterpret_cast<size_type*>(buffer);
    }
    
    const size_type& size(const void* buffer) const
    {
      return *reinterpret_cast<const size_type*>(buffer);
    }
    
    word_type::id_type* context(void* buffer) const
    {
      return reinterpret_cast<word_type::id_type*>((char*) buffer + sizeof(size_type));
    }
    
    const word_type::id_type* context(const void* buffer) const
    {
      return reinterpret_cast<const word_type::id_type*>((const char*) buffer + sizeof(size_type));
    }
    
    logprob_type* backoff(void* buffer) const
    {
      return reinterpret_cast<logprob_type*>((char*) buffer + sizeof(size_type) + sizeof(word_type::id_type) * (order_ - 1));
    }
    
    const logprob_type* backoff(const void* buffer) const
    {
      return reinterpret_cast<const logprob_type*>((const char*) buffer + sizeof(size_type) + sizeof(word_type::id_type) * (order_ - 1));
    }
    
    void fill(void* buffer) const
    {
      const size_type len = size(buffer);
      
      std::fill(reinterpret_cast<char*>(context(buffer) + len), reinterpret_cast<char*>(context(buffer) + order_ - 1), 0);
      std::fill(reinterpret_cast<char*>(backoff(buffer) + len), reinterpret_cast<char*>(backoff(buffer) + order_ - 1), 0);
    }

    void copy(const void* buffer, void* copied)
    {
      std::copy((char*) buffer, ((char*) buffer) + buffer_size(), (char*) copied);
    }
    

    void append(const void* buffer1, const void* buffer2, void* appended)
    {
      const size_type len1 = size(buffer1);
      const size_type len2 = size(buffer2);

      // merge context
      std::copy(context(buffer1), context(buffer1) + len1, context(appended));
      std::copy(context(buffer2), context(buffer2) + len2, context(appended) + len1);
      
      // merge backoff
      std::copy(backoff(buffer1), backoff(buffer1) + len1, backoff(appended));
      std::copy(backoff(buffer2), backoff(buffer2) + len2, backoff(appended) + len1);
      
      // merge size
      
      size(appended) = len1 + len2;
    }
  };
};

#endif
