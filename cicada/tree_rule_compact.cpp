
#include "tree_rule_compact.hpp"
#include "tree_rule_codec.hpp"

#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>

namespace cicada
{
  void TreeRuleCompact::encode(const TreeRule& x)
  {
    typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;
    
    buffer_type buffer;
    tree_rule_encode(x, std::back_inserter(buffer));
    impl.assign(buffer.begin(), buffer.end());
  }
  
  
  TreeRule TreeRuleCompact::decode() const
  {
    TreeRule tree_rule;
    tree_rule_decode(impl.begin(), impl.end(), tree_rule);
    return tree_rule;
  }
};
