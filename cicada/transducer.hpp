// -*- mode: c++ -*-

#ifndef __CICADA__TRANSDUCER__HPP__
#define __CICADA__TRANSDUCER__HPP__ 1

// base class for transducer implementation...
// we will derive all the transducer based storage from here:

#include <stdint.h>

#include <vector>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/rule.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  
  class Transducer
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol symbol_type;
    typedef cicada::Vocab  vocab_type;
    typedef cicada::Rule   rule_type;
    
    typedef boost::shared_ptr<rule_type> rule_ptr_type;
    typedef std::vector<rule_ptr_type, std::allocator<rule_ptr_type> > rule_set_type;
    
    // 64-bit id type!
    typedef uint64_t id_type;
    
  public:
    virtual ~Transducer() {}
    virtual bool valid_span(int first, int last, int distance) const = 0;
    virtual id_type root() const = 0;
    virtual id_type next(const id_type& node, const symbol_type& symbol) const = 0;
    virtual bool has_next(const id_type& node) const = 0;
    virtual const rule_set_type& rules(const id_type& node) const = 0;
  };
  
  
};

#endif


