// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __DEPENDENCY__DEGREE2__HPP__
#define __DEPENDENCY__DEGREE2__HPP__ 1

#include <vector>
#include <utility>

#include <utils/vector3.hpp>
#include <utils/chart.hpp>

// degree-2 non-projective parser based on the deduction system presented in
//
// @InProceedings{cohen-gomezrodriguez-satta:2011:EMNLP,
//   author    = {Cohen, Shay B.  and  G\'{o}mez-Rodr\'{i}guez, Carlos  and  Satta, Giorgio},
//   title     = {Exact Inference for Generative Probabilistic Non-Projective Dependency Parsing},
//   booktitle = {Proceedings of the 2011 Conference on Empirical Methods in Natural Language Processing},
//   month     = {July},
//   year      = {2011},
//   address   = {Edinburgh, Scotland, UK.},
//   publisher = {Association for Computational Linguistics},
//   pages     = {1234--1245},
//   url       = {http://www.aclweb.org/anthology/D11-1114}
// }
//

struct DependencyDegree2
{
  typedef std::pair<int, int> dep_type;
  typedef std::vector<dep_type, std::allocator<dep_type> > dep_set_type;
  
  struct Hypothesis
  {
    double       score;
    dep_set_type deps;
    
    Hypothesis() : score(lowest()), deps() {}
    
    bool valid() const { return score > lowest(); }
    
    void clear() { score = lowest(); deps.clear(); }

    Hypothesis& operator+=(const dep_type& x)
    {
      deps.push_back(x);
      return *this;
    }
    
    Hypothesis& operator+=(const Hypothesis& x)
    {
      score += x.score;
      deps.insert(deps.end(), x.deps.begin(), x.deps.end());
      return *this;
    }
    
    static inline double lowest()  { return - std::numeric_limits<double>::infinity(); }
  };
  
  typedef Hypothesis hypothesis_type;
  
  // hypothesis pair: first for root-assigned hypothesis, second for root-non-assigned hypothesis
  typedef std::pair<hypothesis_type, hypothesis_type> hypothesis_pair_type;

  typedef utils::vector3<hypothesis_pair_type, std::allocator<hypothesis_pair_type> > hypothesis_set_type;
  
  typedef utils::chart<hypothesis_set_type, std::allocator<hypothesis_set_type> >  hypothesis_chart_type;
  
  hypothesis_chart_type actives;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    
    actives.clear();
    actives.resize(sentence_size + 2, hypothesis_set_type(sentence_size + 2, sentence_size + 2, sentence_size + 2, hypothesis_pair_type()));
    
    // initialize by axioms...
    for (size_t pos = 0; pos != sentence_size; ++ pos) {
      // we need to shift + 1 for correct indexing...
      // [h3, j, h3 j, j + 1] where j == pos + 1 and h3 should starts from -1
      
      const int j = pos + 1;
      for (int h3 = -1; h3 != j; ++ h3)
	actives(j, j + 1)(h3 + 1, h3 + 1, j + 1).second.score = 0.0;
    }
    
    const int last_max = sentence_size + 1;
    for (int last = 2; last <= last_max; ++ last)
      for (int length = 2; last - length >= 0; ++ length)  {
	const int first = last - length;
	
	//
	// is it correct???
	// [h1, first, h2 h3, middle] [h3, middle, h4 h5, last]
	//
	// first  <= h3 < middle
	// middle <= h5 < last
	// h3 <= h4 < h5
	// -1 <= h1 < h3 (or first??? given h3 < middle...?)
	// h1 <= h2 < h3
	
	hypothesis_set_type& cells = actives(first, last);
	
	for (int middle = first + 1; middle < last; ++ middle) {
	  const hypothesis_set_type& lefts  = actives(first, middle);
	  const hypothesis_set_type& rights = actives(middle, last);
	  
	  for (int h3 = first; h3 < middle; ++ h3)
	    for (int h5 = middle; h5 < last; ++ h5)
	      for (int h4 = h3; h4 < h5; ++ h4)
		for (int h1 = -1; h1 < first; ++ h1)
		  for (int h2 = h1; h2 < h3; ++ h2) {
		    const hypothesis_pair_type& left  = lefts(h1 + 1, h2 + 1, h3 + 1);
		    const hypothesis_pair_type& right = rights(h3 + 1, h4 + 1, h5 + 1);
		    
		    // [h1, i, h2 h5, j] (la1; h5 -> h4)
		    if (h4 > 0) {
		      // left attachment
		      hypothesis_pair_type& cell = cells(h1 + 1, h2 + 1, h5 + 1);
		      const double& score_edge = scores(h5, h4);
		      
		      enumerate(cell.first,  left.first,  right.second, h5, h4, score_edge);
		      enumerate(cell.first,  left.second, right.first,  h5, h4, score_edge);
		      enumerate(cell.second, left.second, right.second, h5, h4, score_edge);
		    }
		    
		    // [h1, i, h2 h4, j] (ra1; h4 -> h5)
		    if (h4 >= 0) {
		      // right attachment
		      hypothesis_pair_type& cell = cells(h1 + 1, h2 + 1, h4 + 1);
		      const double& score_edge = scores(h4, h5);
		      
		      if (h4 == 0)
			enumerate(cell.first, left.second, right.second, h4, h5, score_edge);
		      else {
			enumerate(cell.first,  left.first,  right.second, h4, h5, score_edge);
			enumerate(cell.first,  left.second, right.first,  h4, h5, score_edge);
			enumerate(cell.second, left.second, right.second, h4, h5, score_edge);
		      }
		    }
		    
		    // [h1, i, h4 h5, j] (la2; h5 -> h2)
		    if (h2 > 0) {
		      // left attachment
		      hypothesis_pair_type& cell = cells(h1 + 1, h4 + 1, h5 + 1);
		      const double score_edge = scores(h5, h2);
		      
		      enumerate(cell.first,  left.first,  right.second, h5, h2, score_edge);
		      enumerate(cell.first,  left.second, right.first,  h5, h2, score_edge);
		      enumerate(cell.second, left.second, right.second, h5, h2, score_edge);
		    }
		    
		    // [h1, i, h2 h4, j] (ra2; h2 -> h5)
		    if (h2 >= 0) {
		      // right attachment
		      hypothesis_pair_type& cell = cells(h1 + 1, h2 + 1, h4 + 1);
		      const double score_edge = scores(h2, h5);
		      
		      if (h2 == 0)
			enumerate(cell.first, left.second, right.second, h2, h5, score_edge);
		      else {
			enumerate(cell.first,  left.first,  right.second, h2, h5, score_edge);
			enumerate(cell.first,  left.second, right.first,  h2, h5, score_edge);
			enumerate(cell.second, left.second, right.second, h2, h5, score_edge);
		      }
		    }
		  }
	}
      }
    
    const dep_set_type& deps = actives(0, last_max)(-1 + 1, -1 + 1, 0 + 1).first.deps;
    
    dep_set_type::const_iterator diter_end = deps.end();
    for (dep_set_type::const_iterator diter = deps.begin(); diter != diter_end; ++ diter)
      dependency[diter->second - 1] = diter->first;
  }
  
  void enumerate(hypothesis_type& cell, const hypothesis_type& left, const hypothesis_type& right, const int& head, const int& dep, const double& score_edge)
  {
    if (! left.valid() || ! right.valid()) return;
    
    if (score_edge + left.score + right.score <= cell.score) return;
    
    cell.clear();
    cell.score = score_edge;
    
    cell += dep_type(head, dep);
    cell += left;
    cell += right;
  }
  
};

#endif
