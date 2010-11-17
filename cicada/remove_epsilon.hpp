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
    
    // @inproceedings{Mohri:2000:GER:647267.721706,
    //   author = {Mohri, Mehryar},
    //   title = {Generic epsilon -Removal Algorithm for Weighted Automata},
    //   booktitle = {Revised Papers from the 5th International Conference on Implementation and Application of Automata},
    //   series = {CIAA '00},
    //   year = {2001},
    //   isbn = {3-540-42491-1},
    //   pages = {230--242},
    //   numpages = {13},
    //   url = {http://portal.acm.org/citation.cfm?id=647267.721706},
    //   acmid = {721706},
    //   publisher = {Springer-Verlag},
    //   address = {London, UK},
    // }

    typedef Lattice lattice_type;
    typedef Vocab   vocab_type;
    
    typedef lattice_type::feature_set_type feature_set_type;

    typedef std::vector<bool, std::allocator<bool> > covered_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > distance_set_type;

    
    void operator()(const lattice_type& source, lattice_type& target)
    {
      typedef std::pair<size_t, feature_set_type > epsilon_type;
      typedef std::set<epsilon_type, std::less<epsilon_type>, std::allocator<epsilon_type> > epsilon_path_type;
      typedef std::vector<epsilon_path_type, std::allocator<epsilon_path_type> > closure_type;
      
      
      target.clear();
      
      closure_type closure(source.size() + 1);

      covered_type      done(source.size() + 1);
      distance_set_type distances;
      
      // calculate epsilon-closure via single-source shortest distance algorithm
      
      for (size_t state = 0; state != source.size(); ++ state) {
	done.clear();
	done.resize(source.size() + 1, false);
	distances.clear();
	
	bool distances_calculated = false;
	
	std::deque<int, std::allocator<int> > stack;
	stack.push(state);
	done[state] = true;
	
	while (! stack.empty()) {
	  const size_t state_id = stack.back();
	  stack.pop_back();
	  
	  if (state_id == source.size()) continue;
	  
	  lattice_type::arc_set_type::const_iterator aiter_end = source[state_id].end();
	  for (lattice_type::arc_set_type::const_iterator aiter = source[state_id].begin(); aiter != aiter_end; ++ aiter) 
	    if (aiter->label == vocab_type::EPSILON && ! done[state_id + aiter->distance]) {
	      const size_t state_next = state_id + aiter->distance;

	      if (! distance_calculated) {
		single_source_shortest_paths(source, state, distances);
		distance_calculated = true;
	      } 
	      done[state_next] = true;
	      
	      closure[state].insert(epsilon_type(state_next, distances[state_next]));
	      
	      stack.push(state_next);
	    }
	}
      }
      
      // replacing epsilon arcs...
      target.clear();
      target.resize(source.size());
      
      for (size_t state = 0; state != source.size(); ++ state) {
	
	lattice_type::arc_set_type::const_iterator aiter_end = aiter_end = source[state].end();
	for (lattice_type::arc_set_type::const_iterator aiter = source[state].begin(); aiter != aiter_end; ++ aiter)
	  if (aiter->label != vocab_type::EPSILON)
	    target[state].push_back(*aiter);
	
	epsilon_path_type::const_iterator citer_end = closure[state].end();
	for (epsilon_path_type::const_iterator citer = closure[state].begin(); citer != citer_end; ++ citer) {
	  
	  lattice_type::arc_set_type::const_iterator niter_end = source[citer->first].end();
	  for (lattice_type::arc_set_type::const_iterator niter = source[citer->first].begin(); niter != niter_end; ++ niter) 
	    if (niter->label != vocab_type::EPSILON) {
	      
	      
	    }
	}
	
      }
      
      
    }
    
    // here is some trick related feature-set...
    //    we use logsum, since feature-set are expected log-linear combination...
    //
    void single_source_shortest_paths(const lattice_type& lattice, const int start, distance_set_type& distances)
    {
      typedef std::set<int, std::less<int>, std::allocator<int> > queue_type;
      
      distances.resize(lattice.size() + 1);

      distance_set_type& d = distances;
      covered_type       d_zero(lattice.size() + 1, true);
      
      distance_set_type r(lattice.size() + 1);
      covered_type      r_zero(lattice.size() + 1, true);
      
      d_zero[start] = false;
      r_zero[start] = false;
      
      queue_type S;
      S.push(start);
      
      while (! S.empty()) {
	const int q = *S.begin();
	S.erase(S.begin());
	
	feature_set_type R = r[q];
	bool R_zero = r_zero[q];
	
	r[q].clear();
	r_zero[q] = true;
	
	if (q == lattice.size()) continue;
	
	lattice_type::arc_set_type::const_iterator aiter_end = lattice[q].end();
	for (lattice_type::arc_set_type::const_iterator aiter = lattice[q].begin(); aiter != aiter_end; ++ aiter) 
	  if (aiter->label == vocab_type::EPSILON) {
	    const size_t next = q + aiter->distance;
	    
	    // r \times w[e]
	    feature_set_type W = R;
	    bool W_zero = r_zero;
	    multiply(W, W_zero, aiter->featuers, false);
	    
	    // d[next] \plus (r \times w[e])
	    feature_set_type D = W;
	    feature_set_type D_zero = W_zero;
	    add(D, D_zero, d[next], d_zero[next]);
	    
	    if (d[next] != D && d_zero[next] != D_zero) {
	      d[next] = D;
	      d_zero[next] = D_zero;
	      
	      // r[next] \ plus (r \times w[e])
	      add(r[next], r_zero[next], W, W_zero);
	      
	      S.insert(next);
	    }
	  }
      }
      
      d[start].clear();
    }
  };

  void multiply(feature_set_type& x, bool& x_zero, const feature_set_type& y, const bool y_zero)
  {
    if (x_zero || y_zero) {
      x.clear();
      x_zero = true;
    }
    x += y;
    x_zero = false;
  }
  
  void add(feature_set_type& x, bool& x_zero, const feature_set_type& y, const bool y_zero)
  {
    if (y_zero) return;
    
    if (x_zero) {
      x = y;
      x_zero = y_zero;
      return;
    }
    
    feature_set_type features;
    
    feature_set_type::iterator iter1 = x.begin();
    feature_set_type::iterator iter1_end = x.end();
    
    feature_set_type::const_iterator iter2 = y.begin();
    feature_set_type::const_iterator iter2_end = y.end();
    
    while (iter1 != iter1_end && iter2 != iter2_end) {
      if (iter1->first < iter2->first) {
	features[iter1->first] = utils::mathop::logsum(iter1->second, 0.0);
	++ iter1;
      } else if (iter2->first < iter1->first) {
	features[iter2->first] = utils::mathop::logsum(iter2->second, 0.0);
	++ iter2;
      } else {
	features[iter1->first] = utils::mathop::logsum(iter1->second, iter2->second);
	++ iter1;
	++ iter2;
      }
    }
    
    for (/**/; iter1 != iter1_end; ++ iter1)
      features[iter1->first] = utils::mathop::logsum(iter1->second, 0.0);
    
    for (/**/; iter2 != iter2_end; ++ iter2)
      features[iter2->first] = utils::mathop::logsum(iter2->second, 0.0);
    
    features.swap(x);
    x_zero = false;
  }
  
  inline
  void remove_epsilin(Lattice& lattice)
  {
    RemoveEpsilon __remover;
    Lattice __lattice;
    __remover(lattice, __lattice);
    lattice.swap(__lattice);
  }
  

  inline
  void remove_epsilin(const Lattice& lattice, Lattice& removed)
  {
    RemoveEpsilon __remover;
    __remover(lattice, removed);
  }
  
};


#endif
