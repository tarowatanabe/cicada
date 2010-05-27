// -*- mode: c++ -*-

#ifndef __CICADA__LATTICE__HPP__
#define __CICADA__LATTICE__HPP__ 1

#include <vector>
#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/sentence.hpp>
#include <cicada/feature_vector.hpp>

#include <utils/vector2.hpp>

namespace cicada
{
  // lattice class.. actually, this will also serve as base class for sentences......
  // currently, we will use this term for future extention to actual lattice..
  
  // lattice format:
  // JSON lattice format (recommended format!)
  // arc = ["double-quoted label", {"double-quoted-feature-name": 0.4, "another": 0.5}, 1] %%% label, feature(s), distance
  // arc-set = [arc, arc, ..., arc]
  // lattice = [arc-set, arc-set, ..., arc-set]
  
  // Python lattice format:
  // arc = ('single-quoted label', 0.4, 1) %%% label, cost, distance. cost is mapped to "lattice-cost" feature
  // arc-set = (arc, arc, arc, ..., arc,)
  // lattice = (arc-set, arc-set, ..., arc-set,)
  
  // The use of "shortest path" is motivated by:
  //
  // @InProceedings{dyer-muresan-resnik:2008:ACLMain,
  //   author    = {Dyer, Christopher  and  Muresan, Smaranda  and  Resnik, Philip},
  //   title     = {Generalizing Word Lattice Translation},
  //   booktitle = {Proceedings of ACL-08: HLT},
  //   month     = {June},
  //   year      = {2008},
  //   address   = {Columbus, Ohio},
  //   publisher = {Association for Computational Linguistics},
  //   pages     = {1012--1020},
  //   url       = {http://www.aclweb.org/anthology/P/P08/P08-1115}
  //  }


  
  class Lattice
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol symbol_type;
    typedef cicada::FeatureVector<double, std::allocator<double> > feature_set_type;
    
    struct arc_type
    {
      arc_type() : features(), label(), distance(0) {}
      arc_type(const symbol_type& __label,
	       const feature_set_type& __features,
	       const int& __distance)
	: features(__features), label(__label), distance(__distance) {}
      
      feature_set_type features;
      symbol_type      label;
      int              distance;
    };
    
    typedef std::vector<arc_type, std::allocator<arc_type> > arc_set_type;
    typedef std::vector<arc_set_type, std::allocator<arc_set_type> >  lattice_type;

  private:
    typedef utils::vector2<difference_type, std::allocator<difference_type> >     distance_type;

  public:
    typedef lattice_type::value_type      value_type;
    
    typedef lattice_type::reference       reference;
    typedef lattice_type::const_reference const_reference;
    
    typedef lattice_type::iterator       iterator;
    typedef lattice_type::const_iterator const_iterator;
    
  public:
    Lattice() : lattice(), dist() {}
    Lattice(const std::string& x) { assign(x); }
    Lattice(const Sentence& x) : lattice(x.size()), dist()
    {
      iterator liter = lattice.begin();
      Sentence::const_iterator iter_end = x.end();
      for (Sentence::const_iterator iter = x.begin(); iter != iter_end; ++ iter, ++ liter)
	*liter = arc_set_type(1, arc_type(*iter, feature_set_type(), 1));
    }
    
    void assign(const std::string& x);
    
    size_type size() const { return lattice.size(); }
    bool empty() const { return lattice.empty(); }
    
    const_reference operator[](size_type x) const { return lattice[x]; }
    reference operator[](size_type x) { return lattice[x]; }
    
    difference_type shortest_distance(difference_type begin, difference_type end) const 
    {
      return (dist.empty() ? end - begin : dist(begin, end));
    }
    
    void swap(Lattice& x)
    {
      lattice.swap(x.lattice);
      dist.swap(x.dist);
    }

    void clear() { lattice.clear(); dist.clear(); }

    friend
    std::ostream& operator<<(std::ostream& os, const Lattice& x);
    friend
    std::istream& operator>>(std::istream& is, Lattice& x);
    
  private:
    void initialize_distance();
    
  private:
    lattice_type lattice;
    distance_type dist;
  };
  
};

namespace std
{
  inline
  void swap(cicada::Lattice& x, cicada::Lattice& y)
  {
    x.swap(y);
  }
  
};

#endif
