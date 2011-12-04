// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_RULE_CODEC__HPP__
#define __CICADA__TREE_RULE_CODEC__HPP__ 1

#include <cicada/tree_rule.hpp>

#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>

namespace cicada
{
  struct TreeRuleCODEC
  {
    typedef cicada::Symbol       symbol_type;
    typedef cicada::Symbol       label_type;
    
    typedef cicada::TreeRule rule_type;
    
    typedef uint8_t  byte_type;
    typedef uint64_t off_type;

    static size_t encode(byte_type* buffer, const off_type& value)
    {
      return utils::byte_aligned_encode(value, reinterpret_cast<char*>(buffer));
    }
      
    static size_t encode(byte_type* buffer, const symbol_type::id_type& value)
    {
      return utils::byte_aligned_encode(value, reinterpret_cast<char*>(buffer));
    }
      
    static size_t encode(byte_type* buffer, const symbol_type& value)
    {
      return encode(buffer, value.id());
    }
      
    static size_t decode(const byte_type* buffer, off_type& value)
    {
      return utils::byte_aligned_decode(value, reinterpret_cast<const char*>(buffer));
    }

    static size_t decode(const byte_type* buffer, symbol_type::id_type& value)
    {
      return utils::byte_aligned_decode(value, reinterpret_cast<const char*>(buffer));
    }
      
    static size_t decode(const byte_type* buffer, symbol_type& value)
    {
      symbol_type::id_type value_id = 0;
      const size_t ret = utils::byte_aligned_decode(value_id, reinterpret_cast<const char*>(buffer));
      value = symbol_type(value_id);
      return ret;
    }
    
    template <typename Iterator>
    static inline
    void encode(const rule_type& rule, Iterator result)
    {
      byte_type buf[32];
      byte_type* begin = buf;
      byte_type* iter  = buf;
      
      std::advance(iter, encode(&(*iter), rule.label));
      std::advance(iter, encode(&(*iter), off_type(rule.antecedents.size())));
      
      std::copy(reinterpret_cast<char*>(begin), reinterpret_cast<char*>(iter), result);
      
      for (size_t i = 0; i != rule.antecedents.size(); ++ i)
	encode(rule.antecedents[i], result);
    }
    
    template <typename Iterator>
    static inline
    void decode(Iterator& iter, TreeRule& rule)
    {
      off_type size = 0;
      std::advance(iter, decode(reinterpret_cast<const byte_type*>(&(*iter)), rule.label));
      std::advance(iter, decode(reinterpret_cast<const byte_type*>(&(*iter)), size));
      
      rule.antecedents.resize(size);
      for (off_type i = 0; i != size; ++ i)
	decode(iter, rule.antecedents[i]);
    }

    template <typename Iterator, typename _Vocab>
    static inline
    void decode(Iterator& iter, const _Vocab& vocab, TreeRule& rule)
    {
      symbol_type::id_type id = 0;
      off_type size = 0;
      std::advance(iter, decode(reinterpret_cast<const byte_type*>(&(*iter)), id));
      std::advance(iter, decode(reinterpret_cast<const byte_type*>(&(*iter)), size));
      
      rule.label = vocab[id];
      rule.antecedents.resize(size);
      for (off_type i = 0; i != size; ++ i)
	decode(iter, vocab, rule.antecedents[i]);
    }
  };
  
  template <typename Iterator>
  inline
  void tree_rule_encode(const TreeRule& x, Iterator result)
  {
    TreeRuleCODEC codec;
    codec.encode(x, result);
  }
  
  template <typename Iterator>
  inline
  void tree_rule_decode(Iterator first, Iterator last, TreeRule& x)
  {
    TreeRuleCODEC codec;
    codec.decode(first, x);
  }
  
  template <typename Iterator, typename _Vocab>
  inline
  void tree_rule_decode(Iterator first, Iterator last, const _Vocab& vocab, TreeRule& x)
  {
    TreeRuleCODEC codec;
    codec.decode(first, vocab, x);
  }
  
};


#endif
