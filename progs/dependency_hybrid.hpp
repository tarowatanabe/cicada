// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __DEPENDENCY__HYBRID__HPP__
#define __DEPENDENCY__HYBRID__HPP__ 1

#include <vector>
#include <utility>

#include <utils/chart.hpp>

// hybrid parser based on the deduction system presented in
//
// @InProceedings{kuhlmann-gomezrodriguez-satta:2011:ACL-HLT2011,
//   author    = {Kuhlmann, Marco  and  G\'{o}mez-Rodr\'{i}guez, Carlos  and  Satta, Giorgio},
//   title     = {Dynamic Programming Algorithms for Transition-Based Dependency Parsers},
//   booktitle = {Proceedings of the 49th Annual Meeting of the Association for Computational Linguistics: Human Language Technologies},
//   month     = {June},
//   year      = {2011},
//   address   = {Portland, Oregon, USA},
//   publisher = {Association for Computational Linguistics},
//   pages     = {673--682},
//   url       = {http://www.aclweb.org/anthology/P11-1068}
// }
//
// which is originally presented in
//
// @INPROCEEDINGS{Yamada03statisticaldependency,
//   author = {Hiroyasu Yamada and Yuji Matsumoto},
//   title = {Statistical Dependency Analysis with Support Vector Machines},
//   booktitle = {In Proceedings of IWPT},
//   year = {2003},
//   pages = {195--206}
// }

struct DependencyHybrid
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
  
  typedef utils::chart<hypothesis_pair_type, std::allocator<hypothesis_pair_type> >  hypothesis_chart_type;
  
  hypothesis_chart_type actives;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();

    actives.clear();
    actives.resize(dependency.size() + 2, hypothesis_pair_type());
    
    // initialize by axioms...
    for (size_t pos = 0; pos != sentence_size; ++ pos)
      actives(pos + 1, pos + 2).second.score = 0.0;
    
    const int last_max = sentence_size + 1;
    for (int last = 2; last <= last_max; ++ last)
      for (int length = 2; last - length >= 0; ++ length)  {
	const int first = last - length;
	
	hypothesis_pair_type& cell = actives(first, last);
	
	for (int middle = first + 1; middle < last; ++ middle) {
	  const hypothesis_pair_type& left = actives(first, middle);
	  const hypothesis_pair_type& right = actives(middle, last);
	  
	  if (last < last_max) {
	    // left attachment
	    
	    // head is last and middle is dependent
	    // we will compute three cases, first and second, second and first and second and second

	    const double score_edge = scores(last, middle);
	    
	    if (left.first.valid() && right.second.valid() && score_edge + left.first.score + right.second.score > cell.first.score) {
	      cell.first.clear();

	      cell.first.score = score_edge;
	      
	      cell.first += dep_type(last, middle);
	      cell.first += left.first;
	      cell.first += right.second;
	    }
	    
	    if (left.second.valid() && right.first.valid() && score_edge + left.second.score + right.first.score > cell.first.score) {
	      cell.first.clear();
	      
	      cell.first.score = score_edge;
	      
	      cell.first += dep_type(last, middle);
	      cell.first += left.second;
	      cell.first += right.first;
	    }
	    
	    if (left.second.valid() && right.second.valid() && score_edge + left.second.score + right.second.score > cell.second.score) {
	      cell.second.clear();
	      
	      cell.second.score = score_edge;
	      
	      cell.second += dep_type(last, middle);
	      cell.second += left.second;
	      cell.second += right.second;
	    }
	  }
	  
	  // right attachment
	  if (first == 0) {
	    const double score_edge = scores(first, middle);
	    
	    if (left.second.valid() && right.second.valid() && score_edge + left.second.score + right.second.score > cell.first.score) {
	      cell.first.clear();
	      
	      cell.first.score = score_edge;
	      
	      cell.first += dep_type(first, middle);
	      cell.first += left.second;
	      cell.first += right.second;
	    }
	  } else {
	    const double score_edge = scores(first, middle);
	    
	    if (left.first.valid() && right.second.valid() && score_edge + left.first.score + right.second.score > cell.first.score) {
	      cell.first.clear();
	      
	      cell.first.score = score_edge;
	      
	      cell.first += dep_type(first, middle);
	      cell.first += left.first;
	      cell.first += right.second;
	    }
	    
	    if (left.second.valid() && right.first.valid() && score_edge + left.second.score + right.first.score > cell.first.score) {
	      cell.first.clear();
	      
	      cell.first.score = score_edge;
	      
	      cell.first += dep_type(first, middle);
	      cell.first += left.second;
	      cell.first += right.first;
	    }
	    
	    if (left.second.valid() && right.second.valid() && score_edge + left.second.score + right.second.score > cell.second.score) {
	      cell.second.clear();
	      
	      cell.second.score = score_edge;
	      
	      cell.second += dep_type(first, middle);
	      cell.second += left.second;
	      cell.second += right.second;
	    }
	  }
	}
      }
    
    const dep_set_type& deps = actives(0, last_max).first.deps;
    
    dep_set_type::const_iterator diter_end = deps.end();
    for (dep_set_type::const_iterator diter = deps.begin(); diter != diter_end; ++ diter)
      dependency[diter->second - 1] = diter->first;
  }
};

#endif
