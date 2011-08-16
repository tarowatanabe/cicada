// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__GRAMMAR__HPP__
#define __CICADA__GRAMMAR__HPP__ 1

#include <cicada/transducer.hpp>

#include <boost/shared_ptr.hpp>

#include <utils/piece.hpp>

namespace cicada
{
  
  class Grammar
  {
  public:
    typedef Transducer transducer_type;
    typedef boost::shared_ptr<transducer_type> transducer_ptr_type;

    typedef transducer_type::rule_type          rule_type;
    typedef transducer_type::rule_ptr_type      rule_ptr_type;
    typedef transducer_type::rule_pair_type     rule_pair_type;
    typedef transducer_type::rule_pair_set_type rule_pair_set_type;

    typedef transducer_type::hypergraph_type hypergraph_type;
    typedef transducer_type::lattice_type    lattice_type;

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
    Grammar() : transducers() {}
    template <typename Iterator>
    Grammar(Iterator first, Iterator last)
    {
      for (/**/; first != last; ++ first)
	push_back(*first);
    }

    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      clear();
      for (/**/; first != last; ++ first)
	push_back(*first);
    }
    
    size_type size() const { return transducers.size(); }
    bool empty() const { return transducers.empty(); }
    
    const transducer_type& operator[](size_type x) const { return *transducers[x]; }
    transducer_type& operator[](size_type x) { return *transducers[x]; }
        
    const_iterator begin() const { return transducers.begin(); }
    iterator begin() { return transducers.begin(); }
    
    const_iterator end() const { return transducers.end(); }
    iterator end() { return transducers.end(); }
    
    void push_back(const transducer_ptr_type& x) { transducers.push_back(x); }
    void push_back(const utils::piece& parameter) { transducers.push_back(transducer_type::create(parameter)); }
    void pop_back() { transducers.pop_back(); }
    void clear() { transducers.clear(); }

    void swap(Grammar& x)
    {
      transducers.swap(x.transducers);
    }
    
    Grammar clone() const
    {
      Grammar __grammar;
      
      transducer_ptr_set_type::const_iterator iter_end = transducers.end();
      for (transducer_ptr_set_type::const_iterator iter = transducers.begin(); iter != iter_end; ++ iter)
	__grammar.push_back((*iter)->clone());
      
      return __grammar;
    }

    void assign(const lattice_type& lattice) const
    {
      transducer_ptr_set_type::const_iterator iter_end = transducers.end();
      for (transducer_ptr_set_type::const_iterator iter = transducers.begin(); iter != iter_end; ++ iter)
	const_cast<transducer_ptr_type&>(*iter)->assign(lattice);
    }
    
    void assign(const lattice_type& lattice, const lattice_type& lattice2) const
    {
      transducer_ptr_set_type::const_iterator iter_end = transducers.end();
      for (transducer_ptr_set_type::const_iterator iter = transducers.begin(); iter != iter_end; ++ iter)
	const_cast<transducer_ptr_type&>(*iter)->assign(lattice, lattice2);
    }
    
    void assign(const hypergraph_type& hypergraph) const
    {
      transducer_ptr_set_type::const_iterator iter_end = transducers.end();
      for (transducer_ptr_set_type::const_iterator iter = transducers.begin(); iter != iter_end; ++ iter)
	const_cast<transducer_ptr_type&>(*iter)->assign(hypergraph);
    }
    
    void assign(const hypergraph_type& hypergraph, const lattice_type& lattice) const
    {
      transducer_ptr_set_type::const_iterator iter_end = transducers.end();
      for (transducer_ptr_set_type::const_iterator iter = transducers.begin(); iter != iter_end; ++ iter)
	const_cast<transducer_ptr_type&>(*iter)->assign(hypergraph, lattice);
    }

  public:    
    static const char* lists() { return transducer_type::lists(); }
    static transducer_ptr_type create(const utils::piece& parameter) { return transducer_type::create(parameter); }
    
  private:
    transducer_ptr_set_type transducers;
  };
  
};

namespace std
{
  inline
  void swap(cicada::Grammar& x, cicada::Grammar& y)
  {
    x.swap(y);
  }
};

#endif
