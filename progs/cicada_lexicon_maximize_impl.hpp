//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_MAXIMIZE_IMPL__HPP__
#define __CICADA_LEXICON_MAXIMIZE_IMPL__HPP__ 1

#include "utils/mathop.hpp"

struct Maximize
{
  template <typename Counts>
  void operator()(Counts& counts)
  {
    double sum = 0.0;
    typename Counts::iterator iter_end = counts.end();
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior_lexicon;
    
    const double factor = 1.0 / sum;
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      iter->second = (iter->second + prior_lexicon) * factor;
  }
};

struct MaximizeBayes
{
  template <typename Counts>
  void operator()(Counts& counts)
  {
    double sum = 0.0;
    typename Counts::iterator iter_end = counts.end();
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      sum += iter->second + prior_lexicon;
    
    const double sum_digamma = utils::mathop::digamma(sum);
    for (typename Counts::iterator iter = counts.begin(); iter != iter_end; ++ iter)
      iter->second = utils::mathop::exp(utils::mathop::digamma(iter->second + prior_lexicon) - sum_digamma);
  }
};

#endif
