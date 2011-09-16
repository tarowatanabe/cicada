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
  
  typedef utils::chart<hypothesis_type, std::allocator<hypothesis_type> >  hypothesis_chart_type;
  
  hypothesis_chart_type actives;
  
  template <typename Scores, typename Dependency>
  void operator()(const Scores& scores,
		  Dependency& dependency)
  {
    const size_t sentence_size = dependency.size();
    
    actives.clear();
    actives.resize(dependency.size() + 2, hypothesis_type());
    
    // initialize by axioms...
    for (size_t pos = 0; pos != sentence_size; ++ pos)
      actives(pos + 1, pos + 2).score = 0.0;
    
    const int last_max = sentence_size + 1;
    for (int last = 2; last <= last_max; ++ last)
      for (int length = 2; last - length >= 0; ++ length)  {
	const int first = last - length;
	
	hypothesis_type& cell = actives(first, last);
	
	for (int middle = first + 1; middle < last; ++ middle) {
	  const hypothesis_type& left = actives(first, middle);
	  const hypothesis_type& right = actives(middle, last);
	  
	  // left attachment
	  if (last < last_max)
	    enumerate(cell, left, right, last, middle, scores(last, middle));
	  
	  // right attachment
	  enumerate(cell, left, right, first, middle, scores(first, middle));
	}
      }
    
    const dep_set_type& deps = actives(0, last_max).deps;
    
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
