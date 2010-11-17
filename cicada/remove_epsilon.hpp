// -*- mode: c++ -*-

#ifndef __CICADA__REMOVE_EPSILON__HPP__
#define __CICADA__REMOVE_EPSILON__HPP__ 1

#include <vector>
#include <queue>
#include <deque>

#include <cicada/lattice.hpp>
#include <cicada/vocab.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/mathop.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  
  struct RemoveEpsilon
  {
    // epsilon removal, but mainly concerns confusion-network lattice...
    typedef Lattice lattice_type;
    typedef Vocab   vocab_type;
    
    typedef lattice_type::feature_set_type feature_set_type;
    
    typedef std::pair<int, feature_set_type> epsilon_type;
    typedef std::vector<epsilon_type, std::allocator<epsilon_type> > closure_type;
    typedef std::vector<closure_type, std::allocator<closure_type> > closure_set_type;
    
    
    void operator()(const lattice_type& source, lattice_type& target)
    {
      target.clear();
      
      closure_set_type closure(source.size() + 1);
      
      // initial closure...
      for (size_t first = 0; first != closure.size(); ++ first)
	closure[first].push_back(epsilon_type(first, feature_set_type()));
      
      // then...
      for (int first = source.size() - 1; first >= 0; -- first) {
	lattice_type::arc_set_type::const_iterator aiter_end = source[first].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[first].begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label == vocab_type::EPSILON) {
	    const int last = first + aiter->distance;
	    
	    closure[first].push_back(epsilon_type(last, aiter->features));
	    
	    closure_type::const_iterator citer_end = closure[last].end();
	    for (closure_type::const_iterator citer = closure[last].begin() + 1; citer != citer_end; ++ citer)
	      closure[first].push_back(epsilon_type(citer->first, citer->second + aiter->features));
	  }
      }
      
      for (size_t i = 0; i != source.size(); ++ i)
	target.push_back(lattice_type::arc_set_type());
      
      for (size_t state = 0; state != source.size(); ++ state) {
	
	lattice_type::arc_set_type::const_iterator aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter) {
	  if (aiter->label == vocab_type::EPSILON) {
	    const size_t state_next = state + aiter->distance;
	    
	    closure_type::const_iterator citer_end = closure[state_next].end();
	    for (closure_type::const_iterator citer = closure[state_next].begin(); citer != citer_end; ++ citer) {
	      
	      lattice_type::arc_set_type::const_iterator niter_end = source[citer->first].end();
	      for (lattice_type::arc_set_type::const_iterator niter = source[citer->first].begin(); niter != niter_end; ++ niter) 
		if (niter->label != vocab_type::EPSILON) {
		  // label: niter->label, features: aiter->features + citer->second + niter->features, distance = citer->first + niter->distance
		  
		  target[state].push_back(lattice_type::arc_type(niter->label, aiter->features + citer->second + niter->features, (citer->first + niter->distance) - state));
		}
	    }
	    
	  } else 
	    target[state].push_back(*aiter);
	}
      }
      
    }
  };
  
  inline
  void remove_epsilon(Lattice& lattice)
  {
    RemoveEpsilon __remover;
    Lattice __lattice;
    __remover(lattice, __lattice);
    lattice.swap(__lattice);
  }
  

  inline
  void remove_epsilon(const Lattice& lattice, Lattice& removed)
  {
    RemoveEpsilon __remover;
    __remover(lattice, removed);
  }
  
};


#endif
