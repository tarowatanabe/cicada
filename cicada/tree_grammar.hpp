// -*- mode: c++ -*-

#ifndef __CICADA__TREE_GRAMMAR__HPP__
#define __CICADA__TREE_GRAMMAR__HPP__ 1

#include <cicada/tree_transducer.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  
  class TreeGrammar
  {
  public:
    typedef TreeTransducer transducer_type;
    typedef boost::shared_ptr<transducer_type> transducer_ptr_type;

    typedef transducer_type::rule_type          rule_type;
    typedef transducer_type::rule_ptr_type      rule_ptr_type;
    typedef transducer_type::rule_pair_type     rule_pair_type;
    typedef transducer_type::rule_pair_set_type rule_pair_set_type;

  private:
    typedef std::vector<transducer_ptr_type, std::allocator<transducer_ptr_type> > transducer_ptr_set_type;

  public:
    typedef transducer_ptr_set_type::size_type       size_type;
    typedef transducer_ptr_set_type::difference_type difference_type;
    typedef transducer_ptr_set_type::value_type      value_type;
    
    typedef transducer_ptr_set_type::reference       reference;
    typedef transducer_ptr_set_type::const_reference const_reference;
    typedef transducer_ptr_set_type::iterator        iterator;
    typedef transducer_ptr_set_type::const_iterator  const_iterator;
    
  public:
    TreeGrammar() : transducers() {}
    
    size_type size() const { return transducers.size(); }
    bool empty() const { return transducers.empty(); }
    
    const transducer_type& operator[](size_type x) const { return *transducers[x]; }
    transducer_type& operator[](size_type x) { return *transducers[x]; }
        
    const_iterator begin() const { return transducers.begin(); }
    iterator begin() { return transducers.begin(); }
    
    const_iterator end() const { return transducers.end(); }
    iterator end() { return transducers.end(); }
    
    void push_back(const transducer_ptr_type& x) { transducers.push_back(x); }
    void pop_back() { transducers.pop_back(); }
    void clear() { transducers.clear(); }
    
    TreeGrammar clone() const
    {
      TreeGrammar __grammar;
      
      transducer_ptr_set_type::const_iterator iter_end = transducers.end();
      for (transducer_ptr_set_type::const_iterator iter = transducers.begin(); iter != iter_end; ++ iter)
	__grammar.push_back((*iter)->clone());
      
      return __grammar;
    }
    
  private:
    transducer_ptr_set_type transducers;
  };
  
};

#endif
