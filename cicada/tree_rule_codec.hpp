// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_RULE_CODEC__HPP__
#define __CICADA__TREE_RULE_CODEC__HPP__ 1

//
// a compact tree-rule encoder/decoder
// based on the representation from:
//
// @InProceedings{ghodke-bird-zhang:2011:IJCNLP-2011,
//   author    = {Ghodke, Sumukh  and  Bird, Steven  and  Zhang, Rui},
//   title     = {A Breadth-First Representation for Tree Matching in Large Scale Forest-Based Translation},
//   booktitle = {Proceedings of 5th International Joint Conference on Natural Language Processing},
//   month     = {November},
//   year      = {2011},
//   address   = {Chiang Mai, Thailand},
//   publisher = {Asian Federation of Natural Language Processing},
//   pages     = {785--793},
//   url       = {http://www.aclweb.org/anthology/I11-1088}
// }
//



#include <cicada/tree_rule.hpp>

#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>

namespace cicada
{
  struct TreeRuleCODEC
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef cicada::Symbol       symbol_type;
    typedef cicada::Symbol       label_type;
    
    typedef cicada::TreeRule rule_type;
    
    typedef uint8_t  byte_type;
    typedef uint64_t off_type;
    
    struct bit_ostream_type
    {
      typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;
      
      bit_ostream_type() : size(0), buffer() {}

      void clear()
      {
	size = 0;
	buffer.clear();
      }

      void push(size_type value)
      {
	const size_type pos = size + value;
	const size_type buf_pos  = pos >> 3;
	const size_type mask_pos = pos & size_type(7);
	
	if (buf_pos >= buffer.size())
	  buffer.resize(buf_pos + 1, 0);
	
	buffer[buf_pos] |= (1 << mask_pos);
	size = pos + 1;
      }
      
      size_type   size;
      buffer_type buffer;
    };

    struct bit_istream_type
    {
      bit_istream_type() : size(0), buffer(0) {}

      void clear()
      {
	size = 0;
	buffer = 0;
      }

      template <typename Iterator>
      size_type pop(Iterator& iter)
      {
	size_type value = 0;
	
	for (;;) {
	  const size_type bit_pos = size & size_type(7);
	  if (bit_pos == 0) {
	    buffer = *iter;
	    ++ iter;
	    
	    // optimize...!
	    while (! buffer) {
	      buffer = *iter;
	      ++ iter;
	      size += 8;
	      value += 8;
	    }
	  }
	  
	  const bool test = buffer & (1 << bit_pos);
	  ++ size;
	  
	  if (test)
	    return value;
	  else
	    ++ value;
	}
      }

      size_type   size;
      byte_type   buffer;
    };
    

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

    inline
    void encode_tree(const rule_type& rule, bit_ostream_type& os)
    {
      os.push(rule.antecedents.size());
      
      for (size_t i = 0; i != rule.antecedents.size(); ++ i)
	encode_tree(rule.antecedents[i], os);
    }
    
    template <typename Iterator>
    inline
    void decode_tree(Iterator& iter, bit_istream_type& is, rule_type& rule)
    {
      const size_type size = is.pop(iter);
      
      rule.antecedents.resize(size);
      for (off_type i = 0; i != size; ++ i)
	decode_tree(iter, is, rule.antecedents[i]);
    }

    template <typename Iterator>
    inline
    void decode_label(Iterator& iter, TreeRule& rule)
    {
      std::advance(iter, decode(reinterpret_cast<const byte_type*>(&(*iter)), rule.label));
      
      for (off_type i = 0; i != rule.antecedents.size(); ++ i)
	decode_label(iter, rule.antecedents[i]);
    }

    template <typename Iterator, typename _Vocab>
    inline
    void decode_label(Iterator& iter, const _Vocab& vocab, TreeRule& rule)
    {
      symbol_type::id_type id = 0;
      std::advance(iter, decode(reinterpret_cast<const byte_type*>(&(*iter)), id));
      rule.label = vocab[id];
      
      for (off_type i = 0; i != rule.antecedents.size(); ++ i)
	decode_label(iter, vocab, rule.antecedents[i]);
    }
    
    
    template <typename Iterator>
    inline
    void encode_label(const rule_type& rule, Iterator result)
    {
      byte_type buf[8];
      byte_type* begin = buf;
      byte_type* iter  = buf;
      
      std::advance(iter, encode(&(*iter), rule.label));
      std::copy(reinterpret_cast<char*>(begin), reinterpret_cast<char*>(iter), result);
      
      for (size_t i = 0; i != rule.antecedents.size(); ++ i)
	encode_label(rule.antecedents[i], result);
    }
    
    template <typename Iterator>
    inline
    void encode(const rule_type& rule, Iterator result)
    {
      os.clear();
      
      encode_tree(rule, os);
      std::copy(os.buffer.begin(), os.buffer.end(), result);
      
      encode_label(rule, result);
    }
    
    template <typename Iterator>
    inline
    void decode(Iterator iter, TreeRule& rule)
    {
      is.clear();
      decode_tree(iter, is, rule);
      decode_label(iter, rule);
    }

    template <typename Iterator, typename _Vocab>
    inline
    void decode(Iterator iter, const _Vocab& vocab, TreeRule& rule)
    {
      is.clear();
      decode_tree(iter, is, rule);
      decode_label(iter, vocab, rule);
    }

    bit_ostream_type os;
    bit_istream_type is;
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
