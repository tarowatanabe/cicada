// -*- mode: c++ -*-

#ifndef __CICADA__TREE_TRANSDUCER__HPP__
#define __CICADA__TREE_TRANSDUCER__HPP__ 1

// base class for tree-transducer implementation...
// we will derive all the transducer based storage from here:

#include <stdint.h>

#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/rule.hpp>
#include <cicada/tree_rule.hpp>
#include <cicada/hypergraph.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  
  class TreeTransducer
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol   symbol_type;
    typedef cicada::Vocab    vocab_type;
    typedef cicada::TreeRule rule_type;
    
    typedef cicada::HyperGraph::feature_set_type feature_set_type;
    typedef boost::shared_ptr<rule_type> rule_ptr_type;
    
    struct RulePair
    {
      rule_ptr_type    source;
      rule_ptr_type    target;
      feature_set_type features;
      
      RulePair() : source(), target(), features() {}
    };
    
    typedef RulePair rule_pair_type;
    typedef std::vector<rule_pair_type, std::allocator<rule_pair_type> > rule_pair_set_type;
    
    // 64-bit id type!
    typedef uint64_t id_type;

    typedef TreeTransducer transducer_type;
    typedef boost::shared_ptr<transducer_type> transducer_ptr_type;
    
  public:
    virtual ~TreeTransducer() {}
    virtual transducer_ptr_type clone() const = 0;

    virtual id_type root() const = 0;
    virtual id_type next(const id_type& node, const symbol_type& symbol) const = 0;
    virtual id_type next_epsilon(const id_type& node) const = 0;
    virtual id_type next_comma(const id_type& node) const = 0;
    virtual id_type next_delimitter(const id_type& node) const = 0;
    virtual bool has_next(const id_type& node) const = 0;
    virtual const rule_pair_set_type& rules(const id_type& node) const = 0;
  };
  
  
};

#endif


