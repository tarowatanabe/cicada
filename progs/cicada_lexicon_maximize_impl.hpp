//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_MAXIMIZE_IMPL__HPP__
#define __CICADA_LEXICON_MAXIMIZE_IMPL__HPP__ 1

#include "utils/mathop.hpp"

struct Maximize
{
  template <typename Counts, typename Probs>
  void operator()(const Counts& counts, const Probs& probs, Probs& probs_new, const double& prior, const double& smooth)
  {
    probs_new.clear();

    double sum = 0.0;
    typename Counts::const_iterator iter_end = counts.end();
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior;
    
    const double factor = 1.0 / sum;
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      probs_new[iter->first] = (iter->second + prior) * factor;
  }
};

struct MaximizeBayes
{
  template <typename Counts, typename Probs>
  void operator()(const Counts& counts, const Probs& probs, Probs& probs_new, const double& prior, const double& smooth)
  {
    probs_new.clear();
    
    double sum = 0.0;
    typename Counts::const_iterator iter_end = counts.end();
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior;
    
    const double sum_digamma = utils::mathop::digamma(sum);
    for (typename Counts::const_iterator iter = counts.begin(); iter != iter_end; ++ iter)
      probs_new[iter->first] = utils::mathop::exp(utils::mathop::digamma(iter->second + prior) - sum_digamma);
  }
};

struct MaximizeL0
{
  typedef std::vector<double, std::allocator<double> > prob_set_type;

  MaximizeL0(const double& __alpha, const double& __beta)
    : alpha(__alpha), beta(__beta) {}

  template <typename Counts, typename Probs>
  void operator()(const Counts& counts, const Probs& probs, Probs& probs_new, const double& prior, const double& smooth)
  {
    probs_new.clear();
    
    
  }  
  
  void project()
  {
    
    
  }
  
  double alpha;
  double beta;
};

#endif
