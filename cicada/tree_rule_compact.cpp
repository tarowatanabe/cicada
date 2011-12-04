
#include "tree_rule_compact.hpp"

#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>

namespace cicada
{
  
  // for each label, we encode # of antecedents.
  // encode data in pre-order

  struct TreeRuleCODEC
  {
    typedef cicada::Symbol symbol_type;
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
  };
  
  template <typename Buffer>
  static inline
  void __encode(const TreeRule& rule, Buffer& buffer)
  {
    typedef TreeRuleCODEC codec_type;
    typedef codec_type::off_type off_type;
    typedef codec_type::byte_type byte_type;
    
    byte_type buf[32];
    byte_type* begin = buf;
    byte_type* iter = buf;
    
    std::advance(iter, codec_type::encode(&(*iter), rule.label));
    std::advance(iter, codec_type::encode(&(*iter), off_type(rule.antecedents.size())));
    
    buffer.insert(buffer.end(), begin, iter);
    
    for (size_t i = 0; i != rule.antecedents.size(); ++ i)
      __encode(rule.antecedents[i], buffer);
  }
  

  template <typename Iterator>
  static inline
  void __decode(TreeRule& rule, Iterator& iter)
  {
    typedef TreeRuleCODEC codec_type;
    typedef codec_type::off_type off_type;
    
    off_type size = 0;
    std::advance(iter, codec_type::decode(&(*iter), rule.label));
    std::advance(iter, codec_type::decode(&(*iter), size));
    
    rule.antecedents.resize(size);
    for (off_type i = 0; i != size; ++ i)
      __decode(rule.antecedents[i], iter);
  }
  
  void TreeRuleCompact::encode(const TreeRule& x)
  {
    typedef TreeRuleCODEC codec_type;
    typedef codec_type::off_type off_type;
    typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;
    
    buffer_type buffer;
    __encode(x, buffer);
    impl.assign(buffer.begin(), buffer.end());
  }
  
  
  TreeRule TreeRuleCompact::decode() const
  {
    typedef TreeRuleCODEC codec_type;
    typedef codec_type::off_type off_type;
    
    TreeRule tree_rule;
    impl_type::const_iterator iter = impl.begin();
    __decode(tree_rule, iter);
    return tree_rule;
  }
};
