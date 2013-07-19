// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_STATE__HPP__
#define __CICADA__NGRAM_STATE__HPP__ 1

//
// ngram state representation based on
//
// @InProceedings{heafield:2011:WMT,
//  author    = {Heafield, Kenneth},
//  title     = {KenLM: Faster and Smaller Language Model Queries},
//  booktitle = {Proceedings of the Sixth Workshop on Statistical Machine Translation},
//  month     = {July},
//  year      = {2011},
//  address   = {Edinburgh, Scotland},
//  publisher = {Association for Computational Linguistics},
//  pages     = {187--197},
//  url       = {http://www.aclweb.org/anthology/W11-2123}
// }
//
// and
//
// @InProceedings{heafield-koehn-lavie:2012:EMNLP-CoNLL,
//   author    = {Heafield, Kenneth  and  Koehn, Philipp  and  Lavie, Alon},
//   title     = {Language Model Rest Costs and Space-Efficient Storage},
//   booktitle = {Proceedings of the 2012 Joint Conference on Empirical Methods in Natural Language Processing and Computational Natural Language Learning},
//   month     = {July},
//   year      = {2012},
//   address   = {Jeju Island, Korea},
//   publisher = {Association for Computational Linguistics},
//   pages     = {1169--1178},
//   url       = {http://www.aclweb.org/anthology/D12-1107}
// }

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
    
    NGramState(const size_type order=1) : order_(order) {}
    
    size_type order_;
    
    size_type buffer_size() const
    {
      return sizeof(size_type) + (sizeof(word_type::id_type) + sizeof(logprob_type)) * (order_ - 1);
    }

    bool empty(const void* buffer) const
    {
      return size(buffer) == 0;
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
    
    void clear(void* buffer) const
    {
      size(buffer) = 0;
      
      std::fill(reinterpret_cast<char*>(context(buffer)), reinterpret_cast<char*>(context(buffer) + order_ - 1), 0);
      std::fill(reinterpret_cast<char*>(backoff(buffer)), reinterpret_cast<char*>(backoff(buffer) + order_ - 1), 0);
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
