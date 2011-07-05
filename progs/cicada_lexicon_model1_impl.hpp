//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_MODEL1_IMPL__HPP__
#define __CICADA_LEXICON_MODEL1_IMPL__HPP__ 1

#include <set>

#include "cicada_lexicon_impl.hpp"

#include "utils/vector2.hpp"
#include "utils/mathop.hpp"

#include "kuhn_munkres.hpp"
#include "itg_alignment.hpp"

struct LearnModel1 : public LearnBase
{
  LearnModel1(const LearnBase& __base)
    : LearnBase(__base) {}
  
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  
  typedef std::set<int, std::less<int>, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > point_map_type;
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const alignment_type& alignment,
	     const bool inverse,
	     const ttable_type& ttable,
	     ttable_type& counts,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum = 0.0;
    
    probs.resize(source_size + 1);
    
    points.clear();
    points.resize(target_size);
    
    if (inverse) {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter)
	points[aiter->source].insert(aiter->target);
    } else {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter)
	points[aiter->target].insert(aiter->source);
    }
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      if (! points[trg].empty()) {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
	
	double prob_max    = - std::numeric_limits<double>::infinity();
	word_type word_max = vocab_type::NONE;
	
	point_set_type::const_iterator iter_end = points[trg].end();
	for (point_set_type::const_iterator iter = points[trg].begin(); iter != iter_end; ++ iter) {
	  const int src = *iter;
	  
	  probs[src] = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += probs[src];
	  
	  if (probs[src] > prob_max) {
	    prob_max = probs[src];
	    word_max = source[src];
	  }
	}
	
	logsum += utils::mathop::log(prob_sum);
	
	const double factor = 1.0 / prob_sum;
	for (point_set_type::const_iterator iter = points[trg].begin(); iter != iter_end; ++ iter)
	  counts[source[*iter]][target[trg]] += probs[*iter] * factor;
	
	aligned[word_max].insert(target[trg]);
	
      } else {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
      
	prob_set_type::iterator piter = probs.begin();
	*piter = ttable(vocab_type::NONE, target[trg]) * prob_null;
	prob_sum += *piter;
	
	double prob_max    = *piter;
	word_type word_max = vocab_type::NONE;
	
	++ piter;
	
	for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	  *piter = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += *piter;
	  
	  if (*piter > prob_max) {
	    prob_max = *piter;
	    word_max = source[src];
	  }
	}
	
	logsum += utils::mathop::log(prob_sum);
	
	const double factor = 1.0 / prob_sum;
	piter = probs.begin();
	counts[vocab_type::NONE][target[trg]] += (*piter) * factor;
	++ piter;
	
	for (size_type src = 0; src != source_size; ++ src, ++ piter)
	  counts[source[src]][target[trg]] += (*piter) * factor;
	
	aligned[word_max].insert(target[trg]);
      }
    }
    
    objective += logsum / target_size;    
  }
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const ttable_type& ttable,
	     ttable_type& counts,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum = 0.0;
    
    probs.resize(source_size + 1);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      
      prob_set_type::iterator piter = probs.begin();
      *piter = ttable(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;
      
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      piter = probs.begin();
      counts[vocab_type::NONE][target[trg]] += (*piter) * factor;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter)
	counts[source[src]][target[trg]] += (*piter) * factor;
      
      aligned[word_max].insert(target[trg]);
    }
    
    objective += logsum / target_size;
  }
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    learn(source, target, ttable_source_target, ttable_counts_source_target, aligned_source_target, objective_source_target);
    learn(target, source, ttable_target_source, ttable_counts_target_source, aligned_target_source, objective_target_source);
  }

  void operator()(const sentence_type& source, const sentence_type& target, const alignment_type& alignment)
  {
    learn(source, target, alignment, false, ttable_source_target, ttable_counts_source_target, aligned_source_target, objective_source_target);
    learn(target, source, alignment, true,  ttable_target_source, ttable_counts_target_source, aligned_target_source, objective_target_source);
  }

  prob_set_type  probs;
  point_map_type points;
};

struct LearnModel1Posterior : public LearnBase
{
  LearnModel1Posterior(const LearnBase& __base)
    : LearnBase(__base) {}  

  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  typedef std::vector<prob_type, std::allocator<prob_type> > prob_set_type;

  typedef std::set<int, std::less<int>, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > point_map_type;
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const alignment_type& alignment,
	     const bool inverse,
	     const ttable_type& ttable,
	     ttable_type& counts,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum = 0.0;
    
    posterior.clear();
    posterior.reserve(target_size + 1, source_size + 1);
    posterior.resize(target_size + 1, source_size + 1, 0.0);
    
    probs.clear();
    probs.reserve(target_size + 1, source_size + 1);
    probs.resize(target_size + 1, source_size + 1, 0.0);
    
    phi.clear();
    phi.resize(source_size + 1, 0.0);
    
    exp_phi.clear();
    exp_phi.resize(source_size + 1, 1.0);
    
    points.clear();
    points.resize(target_size);

    if (inverse) {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter)
	points[aiter->source].insert(aiter->target);
    } else {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter)
	points[aiter->target].insert(aiter->source);
    }
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      if (! points[trg].empty()) {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
	
	double prob_max    = - std::numeric_limits<double>::infinity();
	word_type word_max = vocab_type::NONE;
	
	point_set_type::const_iterator iter_end = points[trg].end();
	for (point_set_type::const_iterator iter = points[trg].begin(); iter != iter_end; ++ iter) {
	  const int src = *iter;
	  
	  double& prob = probs(trg + 1, src + 1);
	  
	  prob = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += prob;
	  
	  if (prob > prob_max) {
	    prob_max = prob;
	    word_max = source[src];
	  }
	}
	
	logsum += utils::mathop::log(prob_sum);
	
	const double factor = 1.0 / prob_sum;
	posterior_set_type::const_iterator piter     = probs.begin(trg + 1);
	posterior_set_type::const_iterator piter_end = probs.end(trg + 1);
	posterior_set_type::iterator siter = posterior.begin(trg + 1);
	for (/**/; piter != piter_end; ++ piter, ++ siter)
	  (*siter) = (*piter) * factor;
	
	aligned[word_max].insert(target[trg]);
      } else {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
      
	posterior_set_type::iterator piter     = probs.begin(trg + 1);
	posterior_set_type::iterator piter_end = probs.end(trg + 1);
	*piter = ttable(vocab_type::NONE, target[trg]) * prob_null;
	prob_sum += *piter;
      
	double prob_max    = *piter;
	word_type word_max = vocab_type::NONE;
      
	++ piter;
      
	for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	  *piter = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += *piter;
	
	  if (*piter > prob_max) {
	    prob_max = *piter;
	    word_max = source[src];
	  }
	}
      
	logsum += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
	piter = probs.begin(trg + 1);
	posterior_set_type::iterator siter = posterior.begin(trg + 1);
	for (/**/; piter != piter_end; ++ piter, ++ siter)
	  (*siter) = (*piter) * factor;
      
	aligned[word_max].insert(target[trg]);
      }
    }
    
    objective += logsum / target_size;
    
    for (int iter = 0; iter < 5; ++ iter) {
      // update phi.. but ignore NULL...
      
      bool updated = false;
      for (size_type src = 1; src <= source_size; ++ src) {
	double sum = 0.0;
	for (size_type trg = 1; trg <= target_size; ++ trg)
	  sum += posterior(trg, src);
	
	phi[src] += 1.0 - sum;
	if (phi[src] > 0.0)
	  phi[src] = 0.0;
	
	updated |= (phi[src] != 0.0);
	exp_phi[src] = utils::mathop::exp(phi[src]);
      }
      
      if (! updated) break;
      
      for (size_type trg = 1; trg <= target_size; ++ trg) {
	double sum = 0.0;
	for (size_type src = 0; src <= source_size; ++ src)
	  sum += probs(trg, src) * exp_phi[src];
	
	const double factor = 1.0 / sum;
	for (size_type src = 0; src <= source_size; ++ src)
	  posterior(trg, src) = probs(trg, src) * factor * exp_phi[src];
      }
    }
    
    // update...
    for (size_type trg = 1; trg <= target_size; ++ trg)
      for (size_type src = 0; src <= source_size; ++ src)
	counts[src == 0 ? vocab_type::NONE : source[src - 1]][target[trg - 1]] += posterior(trg, src);
  }
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const ttable_type& ttable,
	     ttable_type& counts,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    double logsum = 0.0;
    
    posterior.reserve(target_size + 1, source_size + 1);
    posterior.resize(target_size + 1, source_size + 1, 0.0);
    
    probs.reserve(target_size + 1, source_size + 1);
    probs.resize(target_size + 1, source_size + 1, 0.0);
    
    phi.clear();
    phi.resize(source_size + 1, 0.0);
    
    exp_phi.clear();
    exp_phi.resize(source_size + 1, 1.0);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      posterior_set_type::iterator piter     = probs.begin(trg + 1);
      posterior_set_type::iterator piter_end = probs.end(trg + 1);
      *piter = ttable(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;
      
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      piter = probs.begin(trg + 1);
      posterior_set_type::iterator siter = posterior.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
      
      aligned[word_max].insert(target[trg]);
    }
    
    objective += logsum / target_size;
    
    for (int iter = 0; iter < 5; ++ iter) {
      // update phi.. but ignore NULL...
      
      bool updated = false;
      for (size_type src = 1; src <= source_size; ++ src) {
	double sum = 0.0;
	for (size_type trg = 1; trg <= target_size; ++ trg)
	  sum += posterior(trg, src);
	
	phi[src] += 1.0 - sum;
	if (phi[src] > 0.0)
	  phi[src] = 0.0;
	
	updated |= (phi[src] != 0.0);
	exp_phi[src] = utils::mathop::exp(phi[src]);
      }
      
      if (! updated) break;
      
      for (size_type trg = 1; trg <= target_size; ++ trg) {
	double sum = 0.0;
	for (size_type src = 0; src <= source_size; ++ src)
	  sum += probs(trg, src) * exp_phi[src];
	
	const double factor = 1.0 / sum;
	for (size_type src = 0; src <= source_size; ++ src)
	  posterior(trg, src) = probs(trg, src) * factor * exp_phi[src];
      }
    }
    
    // update...
    for (size_type trg = 1; trg <= target_size; ++ trg)
      for (size_type src = 0; src <= source_size; ++ src)
	counts[src == 0 ? vocab_type::NONE : source[src - 1]][target[trg - 1]] += posterior(trg, src);
  }

  void operator()(const sentence_type& source, const sentence_type& target)
  {
    learn(source, target, ttable_source_target, ttable_counts_source_target, aligned_source_target, objective_source_target);
    learn(target, source, ttable_target_source, ttable_counts_target_source, aligned_target_source, objective_target_source);
  }

  void operator()(const sentence_type& source, const sentence_type& target, const alignment_type& alignment)
  {
    learn(source, target, alignment, false, ttable_source_target, ttable_counts_source_target, aligned_source_target, objective_source_target);
    learn(target, source, alignment, true,  ttable_target_source, ttable_counts_target_source, aligned_target_source, objective_target_source);
  }

  posterior_set_type posterior;
  posterior_set_type probs;
  point_map_type     points;
  
  prob_set_type      phi;
  prob_set_type      exp_phi;
};

struct LearnModel1Symmetric : public LearnBase
{
  LearnModel1Symmetric(const LearnBase& __base)
    : LearnBase(__base) {}
  
  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;

  typedef std::set<int, std::less<int>, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > point_map_type;

  void operator()(const sentence_type& source, const sentence_type& target, const alignment_type& alignment)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!

    points_source_target.clear();
    points_target_source.clear();
    
    points_source_target.resize(target_size);
    points_target_source.resize(source_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      points_source_target[aiter->target].insert(aiter->source);
      points_target_source[aiter->source].insert(aiter->target);
    }
    
    posterior_source_target.clear();
    posterior_target_source.clear();

    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1, 0.0);
    posterior_target_source.resize(source_size + 1, target_size + 1, 0.0);
    
    prob_source_target.reserve(source_size + 1);
    prob_target_source.reserve(target_size + 1);
    
    prob_source_target.resize(source_size + 1);
    prob_target_source.resize(target_size + 1);

    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      if (! points_source_target[trg].empty()) {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
	
	double prob_max    = - std::numeric_limits<double>::infinity();
	word_type word_max = vocab_type::NONE;
	
	std::fill(prob_source_target.begin(), prob_source_target.end(), 0.0);
	
	point_set_type::const_iterator iter_end = points_source_target[trg].end();
	for (point_set_type::const_iterator iter = points_source_target[trg].begin(); iter != iter_end; ++ iter) {
	  const int src = *iter;
	  
	  double& prob = prob_source_target[src + 1];
	  prob = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += prob;
	  
	  if (prob > prob_max) {
	    prob_max = prob;
	    word_max = source[src];
	  }
	}
	
	logsum_source_target += utils::mathop::log(prob_sum);
	
	const double factor = 1.0 / prob_sum;
	
	prob_set_type::const_iterator piter      = prob_source_target.begin();
	prob_set_type::const_iterator piter_end  = prob_source_target.end();
	posterior_set_type::iterator  siter = posterior_source_target.begin(trg + 1);
	for (/**/; piter != piter_end; ++ piter, ++ siter)
	  (*siter) = (*piter) * factor;
	
	aligned_source_target[word_max].insert(target[trg]);
      } else {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
      
	prob_set_type::iterator piter     = prob_source_target.begin();
	prob_set_type::iterator piter_end = prob_source_target.end();
	*piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
	prob_sum += *piter;
      
	double prob_max    = *piter;
	word_type word_max = vocab_type::NONE;

	++ piter;
      
	for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	  *piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += *piter;
	
	  if (*piter > prob_max) {
	    prob_max = *piter;
	    word_max = source[src];
	  }
	}
      
	logsum_source_target += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
      
	piter = prob_source_target.begin();
	posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
	for (/**/; piter != piter_end; ++ piter, ++ siter)
	  (*siter) = (*piter) * factor;
	
	aligned_source_target[word_max].insert(target[trg]);
      }
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      if (! points_target_source[src].empty()) {
	const double prob_align_norm = 1.0 / target_size;
	double prob_sum = 0.0;
	
	double prob_max    = - std::numeric_limits<double>::infinity();
	word_type word_max = vocab_type::NONE;

	std::fill(prob_target_source.begin(), prob_target_source.end(), 0.0);
	
	point_set_type::const_iterator iter_end = points_target_source[src].end();
	for (point_set_type::const_iterator iter = points_target_source[src].begin(); iter != iter_end; ++ iter) {
	  const int trg = *iter;
	  
	  double& prob = prob_target_source[trg + 1];
	  prob = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	  prob_sum += prob;
	  
	  if (prob > prob_max) {
	    prob_max = prob;
	    word_max = target[trg];
	  }
	}
	
	logsum_target_source += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
      
	prob_set_type::const_iterator piter     = prob_target_source.begin();
	prob_set_type::const_iterator piter_end = prob_target_source.end();
	posterior_set_type::iterator  titer     = posterior_target_source.begin(src + 1);
	for (/**/; piter != piter_end; ++ piter, ++ titer)
	  (*titer) = (*piter) * factor;
	
	aligned_target_source[word_max].insert(source[src]);
      } else {
	const double prob_align_norm = 1.0 / target_size;
	double prob_sum = 0.0;
      
	prob_set_type::iterator piter     = prob_target_source.begin();
	prob_set_type::iterator piter_end = prob_target_source.end();
	*piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
	prob_sum += *piter;
      
	double prob_max    = *piter;
	word_type word_max = vocab_type::NONE;

	++ piter;
      
	for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	  *piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	  prob_sum += *piter;

	  if (*piter > prob_max) {
	    prob_max = *piter;
	    word_max = target[trg];
	  }
	}
      
	logsum_target_source += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
      
	piter = prob_target_source.begin();
	posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
	for (/**/; piter != piter_end; ++ piter, ++ titer)
	  (*titer) = (*piter) * factor;
	
	aligned_target_source[word_max].insert(source[src]);
      }
    }
    
    // accumulate!
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	double count = ((trg == 0 ? 1.0 : posterior_source_target(trg, src))
			* (src == 0 ? 1.0 : posterior_target_source(src, trg)));
	
	if (src != 0 && trg != 0)
	  count = utils::mathop::sqrt(count);
	
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  ttable_counts_source_target[source_word][target_word] += count;
	
	if (src != 0)
	  ttable_counts_target_source[target_word][source_word] += count;
      }
    
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
  } 
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!
    
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1);
    posterior_target_source.resize(source_size + 1, target_size + 1);
    
    prob_source_target.reserve(source_size + 1);
    prob_target_source.reserve(target_size + 1);
    
    prob_source_target.resize(source_size + 1);
    prob_target_source.resize(target_size + 1);

    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_source_target.begin();
      prob_set_type::iterator piter_end = prob_source_target.end();
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin();
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
      
      aligned_source_target[word_max].insert(target[trg]);
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_target_source.begin();
      prob_set_type::iterator piter_end = prob_target_source.end();
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;

	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = target[trg];
	}
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin();
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
      
      aligned_target_source[word_max].insert(source[src]);
    }
    
    // accumulate!
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	double count = ((trg == 0 ? 1.0 : posterior_source_target(trg, src))
			* (src == 0 ? 1.0 : posterior_target_source(src, trg)));
	
	if (src != 0 && trg != 0)
	  count = utils::mathop::sqrt(count);
	
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  ttable_counts_source_target[source_word][target_word] += count;
	
	if (src != 0)
	  ttable_counts_target_source[target_word][source_word] += count;
      }
    
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
  }

  prob_set_type      prob_source_target;
  prob_set_type      prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;

  point_map_type points_source_target;
  point_map_type points_target_source;
  
};

struct LearnModel1SymmetricPosterior : public LearnBase
{
  LearnModel1SymmetricPosterior(const LearnBase& __base)
    : LearnBase(__base) {}
  
  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;

  typedef std::set<int, std::less<int>, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > point_map_type;

  void operator()(const sentence_type& source, const sentence_type& target, const alignment_type& alignment)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    points_source_target.clear();
    points_target_source.clear();
    
    points_source_target.resize(target_size);
    points_target_source.resize(source_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      points_source_target[aiter->target].insert(aiter->source);
      points_target_source[aiter->source].insert(aiter->target);
    }
    
    // we do not have to clearn!
    posterior_source_target.clear();
    posterior_target_source.clear();
    
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1, 0.0);
    posterior_target_source.resize(source_size + 1, target_size + 1, 0.0);
    
    prob_source_target.clear();
    prob_target_source.clear();
    
    prob_source_target.reserve(target_size + 1, source_size + 1);
    prob_target_source.reserve(source_size + 1, target_size + 1);
    
    prob_source_target.resize(target_size + 1, source_size + 1, 0.0);
    prob_target_source.resize(source_size + 1, target_size + 1, 0.0);
    
    phi.clear();
    exp_phi.clear();
    
    phi.resize(target_size + 1, source_size + 1, 0.0);
    exp_phi.resize(target_size + 1, source_size + 1, 1.0);
    
    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      if (! points_source_target[trg].empty()) {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
	
	double prob_max    = - std::numeric_limits<double>::infinity();
	word_type word_max = vocab_type::NONE;
	
	point_set_type::const_iterator iter_end = points_source_target[trg].end();
	for (point_set_type::const_iterator iter = points_source_target[trg].begin(); iter != iter_end; ++ iter) {
	  const int src = *iter;
	  
	  double& prob = prob_source_target(trg + 1, src + 1);
	  prob = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += prob;
	  
	  if (prob > prob_max) {
	    prob_max = prob;
	    word_max = source[src];
	  }
	}
	
	logsum_source_target += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
	posterior_set_type::const_iterator piter     = prob_source_target.begin(trg + 1);
	posterior_set_type::const_iterator piter_end = prob_source_target.end(trg + 1);
	posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
	for (/**/; piter != piter_end; ++ piter, ++ siter)
	  (*siter) = (*piter) * factor;

	aligned_source_target[word_max].insert(target[trg]);
      } else {
	const double prob_align_norm = 1.0 / source_size;
	double prob_sum = 0.0;
      
	posterior_set_type::iterator piter     = prob_source_target.begin(trg + 1);
	posterior_set_type::iterator piter_end = prob_source_target.end(trg + 1);
	*piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
	prob_sum += *piter;
      
	double prob_max    = *piter;
	word_type word_max = vocab_type::NONE;

	++ piter;
      
	for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	  *piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	  prob_sum += *piter;
	
	  if (*piter > prob_max) {
	    prob_max = *piter;
	    word_max = source[src];
	  }
	}
      
	logsum_source_target += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
      
	piter = prob_source_target.begin(trg + 1);
	posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
	for (/**/; piter != piter_end; ++ piter, ++ siter)
	  (*siter) = (*piter) * factor;

	aligned_source_target[word_max].insert(target[trg]);
      }
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      if (! points_target_source[src].empty()) {
	const double prob_align_norm = 1.0 / target_size;
	double prob_sum = 0.0;
	
	double prob_max    = - std::numeric_limits<double>::infinity();
	word_type word_max = vocab_type::NONE;
	
	point_set_type::const_iterator iter_end = points_target_source[src].end();
	for (point_set_type::const_iterator iter = points_target_source[src].begin(); iter != iter_end; ++ iter) {
	  const int trg = *iter;
	  
	  double& prob = prob_target_source(src + 1, trg + 1);
	  
	  prob = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	  prob_sum += prob;
	  
	  if (prob > prob_max) {
	    prob_max = prob;
	    word_max = target[trg];
	  }
	}
	
	logsum_target_source += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
      
	posterior_set_type::const_iterator piter     = prob_target_source.begin(src + 1);
	posterior_set_type::const_iterator piter_end = prob_target_source.end(src + 1);
	posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
	for (/**/; piter != piter_end; ++ piter, ++ titer)
	  (*titer) = (*piter) * factor;
	
	aligned_target_source[word_max].insert(source[src]);
      } else {
	const double prob_align_norm = 1.0 / target_size;
	double prob_sum = 0.0;
      
	posterior_set_type::iterator piter     = prob_target_source.begin(src + 1);
	posterior_set_type::iterator piter_end = prob_target_source.end(src + 1);
	*piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
	prob_sum += *piter;
      
	double prob_max    = *piter;
	word_type word_max = vocab_type::NONE;

	++ piter;
      
	for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	  *piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	  prob_sum += *piter;
	
	  if (*piter > prob_max) {
	    prob_max = *piter;
	    word_max = target[trg];
	  }
	}
      
	logsum_target_source += utils::mathop::log(prob_sum);
      
	const double factor = 1.0 / prob_sum;
      
	piter = prob_target_source.begin(src + 1);
	posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
	for (/**/; piter != piter_end; ++ piter, ++ titer)
	  (*titer) = (*piter) * factor;
      
	aligned_target_source[word_max].insert(source[src]);
      }
    }
    
    // perplexity..
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
    
    // now we will adjust posterior..
    
    for (int iter = 0; iter != 5; ++ iter) {
      
      bool updated = false;
      
      // update phi... we do not consider null alignment!
      for (size_type src = 1; src <= source_size; ++ src)
	for (size_type trg = 1; trg <= target_size; ++ trg) {
	  const double epsi = posterior_source_target(trg, src) - posterior_target_source(src, trg);
	  const double update = - epsi;
	  
	  phi(trg, src) += update;
	  
	  updated |= (phi(trg, src) != 0.0);
	  exp_phi(trg, src) = utils::mathop::exp(phi(trg, src));
	}
      
      if (! updated) break;
      
      // recompute...
      for (size_type trg = 1; trg <= target_size; ++ trg) {
	double prob_sum = 0.0;
	for (size_type src = 0; src <= source_size; ++ src)
	  prob_sum += prob_source_target(trg, src) * exp_phi(trg, src);
	
	const double factor = 1.0 / prob_sum;
	for (size_type src = 0; src <= source_size; ++ src)
	  posterior_source_target(trg, src) = prob_source_target(trg, src) * factor *  exp_phi(trg, src);
      }
      
      for (size_type src = 1; src <= source_size; ++ src) {
	double prob_sum = 0.0;
	for (size_type trg = 0; trg <= target_size; ++ trg)
	  prob_sum += prob_target_source(src, trg) / exp_phi(trg, src);
	
	const double factor = 1.0 / prob_sum;
	for (size_type trg = 0; trg <= target_size; ++ trg)
	  posterior_target_source(src, trg) = prob_target_source(src, trg) * factor / exp_phi(trg, src);
      }
    }
    
    // since we have already adjusted posterior, we simply accumulate individual counts...
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  ttable_counts_source_target[source_word][target_word] += posterior_source_target(trg, src);
	
	if (src != 0)
	  ttable_counts_target_source[target_word][source_word] += posterior_target_source(src, trg);
      }
  }

  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!
    
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1);
    posterior_target_source.resize(source_size + 1, target_size + 1);

    prob_source_target.reserve(target_size + 1, source_size + 1);
    prob_target_source.reserve(source_size + 1, target_size + 1);
    
    prob_source_target.resize(target_size + 1, source_size + 1);
    prob_target_source.resize(source_size + 1, target_size + 1);
    
    phi.clear();
    exp_phi.clear();
    
    phi.resize(target_size + 1, source_size + 1, 0.0);
    exp_phi.resize(target_size + 1, source_size + 1, 1.0);
    
    double logsum_source_target = 0.0;
    double logsum_target_source = 0.0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      posterior_set_type::iterator piter     = prob_source_target.begin(trg + 1);
      posterior_set_type::iterator piter_end = prob_source_target.end(trg + 1);
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = source[src];
	}
      }
      
      logsum_source_target += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin(trg + 1);
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;

      aligned_source_target[word_max].insert(target[trg]);
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      posterior_set_type::iterator piter     = prob_target_source.begin(src + 1);
      posterior_set_type::iterator piter_end = prob_target_source.end(src + 1);
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      
      double prob_max    = *piter;
      word_type word_max = vocab_type::NONE;

      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
	
	if (*piter > prob_max) {
	  prob_max = *piter;
	  word_max = target[trg];
	}
      }
      
      logsum_target_source += utils::mathop::log(prob_sum);
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin(src + 1);
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
      
      aligned_target_source[word_max].insert(source[src]);
    }
    
    // perplexity..
    objective_source_target += logsum_source_target / target_size;
    objective_target_source += logsum_target_source / source_size;
    
    // now we will adjust posterior..
    
    for (int iter = 0; iter != 5; ++ iter) {
      
      bool updated = false;
      
      // update phi... we do not consider null alignment!
      for (size_type src = 1; src <= source_size; ++ src)
	for (size_type trg = 1; trg <= target_size; ++ trg) {
	  const double epsi = posterior_source_target(trg, src) - posterior_target_source(src, trg);
	  const double update = - epsi;
	  
	  phi(trg, src) += update;
	  
	  updated |= (phi(trg, src) != 0.0);
	  exp_phi(trg, src) = utils::mathop::exp(phi(trg, src));
	}
      
      if (! updated) break;
      
      // recompute...
      for (size_type trg = 1; trg <= target_size; ++ trg) {
	double prob_sum = 0.0;
	for (size_type src = 0; src <= source_size; ++ src)
	  prob_sum += prob_source_target(trg, src) * exp_phi(trg, src);
	
	const double factor = 1.0 / prob_sum;
	for (size_type src = 0; src <= source_size; ++ src)
	  posterior_source_target(trg, src) = prob_source_target(trg, src) * factor *  exp_phi(trg, src);
      }
      
      for (size_type src = 1; src <= source_size; ++ src) {
	double prob_sum = 0.0;
	for (size_type trg = 0; trg <= target_size; ++ trg)
	  prob_sum += prob_target_source(src, trg) / exp_phi(trg, src);
	
	const double factor = 1.0 / prob_sum;
	for (size_type trg = 0; trg <= target_size; ++ trg)
	  posterior_target_source(src, trg) = prob_target_source(src, trg) * factor / exp_phi(trg, src);
      }
    }
    
    // since we have already adjusted posterior, we simply accumulate individual counts...
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	const word_type& source_word = (src == 0 ? vocab_type::NONE : source[src - 1]);
	const word_type& target_word = (trg == 0 ? vocab_type::NONE : target[trg - 1]);
	
	if (trg != 0)
	  ttable_counts_source_target[source_word][target_word] += posterior_source_target(trg, src);
	
	if (src != 0)
	  ttable_counts_target_source[target_word][source_word] += posterior_target_source(src, trg);
      }
  }

  posterior_set_type prob_source_target;
  posterior_set_type prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;
  
  posterior_set_type phi;
  posterior_set_type exp_phi;
  
  point_map_type points_source_target;
  point_map_type points_target_source;
};

struct ViterbiModel1 : public ViterbiBase
{
  ViterbiModel1(const ttable_type& __ttable_source_target,
		const ttable_type& __ttable_target_source)
    : ViterbiBase(__ttable_source_target, __ttable_target_source) {}

  void viterbi(const sentence_type& source,
	       const sentence_type& target,
	       const ttable_type& ttable,
	       alignment_type& alignment)
  {
    alignment.clear();
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      
      double prob_max = ttable(vocab_type::NONE, target[trg]) * prob_null;
      int    align_max = -1;
      
      for (size_type src = 0; src != source_size; ++ src) {
	const double prob = ttable(source[src], target[trg]) * prob_align * prob_align_norm;
	
	if (prob > prob_max) {
	  prob_max = prob;
	  align_max = src;
	}
      }
      
      if (align_max >= 0)
	alignment.push_back(std::make_pair(align_max, static_cast<int>(trg)));
    }
  }

  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  const span_set_type& span_source,
		  const span_set_type& span_target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    operator()(source, target, alignment_source_target, alignment_target_source);
  }
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    viterbi(source, target, ttable_source_target, alignment_source_target);
    viterbi(target, source, ttable_target_source, alignment_target_source);
  }
};

struct ITGModel1 : public ViterbiBase
{
  typedef utils::vector2<double, std::allocator<double> > matrix_type;
  typedef utils::vector2<double, std::allocator<double> > posterior_set_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  
  ITGModel1(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source)
    : ViterbiBase(__ttable_source_target, __ttable_target_source) {}

  class insert_align
  {
    alignment_type& alignment_source_target;
    alignment_type& alignment_target_source;
    
  public:
    insert_align(alignment_type& __alignment_source_target,
		 alignment_type& __alignment_target_source)
      : alignment_source_target(__alignment_source_target),
	alignment_target_source(__alignment_target_source) {}
    
    template <typename Edge>
    insert_align& operator=(const Edge& edge)
    {	
      alignment_source_target.push_back(edge);
      alignment_target_source.push_back(std::make_pair(edge.second, edge.first));
      
      return *this;
    }
    
    insert_align& operator*() { return *this; }
    insert_align& operator++() { return *this; }
    insert_align operator++(int) { return *this; }
  };

  const span_set_type __span_source;
  const span_set_type __span_target;
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    operator()(source, target, __span_source, __span_target, alignment_source_target, alignment_target_source);
  }
  
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  const span_set_type& span_source,
		  const span_set_type& span_target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1);
    posterior_target_source.resize(source_size + 1, target_size + 1);
    
    prob_source_target.reserve(source_size + 1);
    prob_target_source.reserve(target_size + 1);
    
    prob_source_target.resize(source_size + 1);
    prob_target_source.resize(target_size + 1);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_source_target.begin();
      prob_set_type::iterator piter_end = prob_source_target.end();
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin();
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_target_source.begin();
      prob_set_type::iterator piter_end = prob_target_source.end();
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin();
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
    }

    costs.clear();
    costs.resize(source_size + 1, target_size + 1, boost::numeric::bounds<double>::lowest());
    
    for (size_type src = 1; src <= source_size; ++ src)
      for (size_type trg = 1; trg <= target_size; ++ trg)
	costs(src, trg) = 0.5 * (utils::mathop::log(posterior_source_target(trg, src)) 
				 + utils::mathop::log(posterior_target_source(src, trg)));
    
    for (size_type trg = 1; trg <= target_size; ++ trg)
      costs(0, trg) = utils::mathop::log(posterior_source_target(trg, 0));
    
    for (size_type src = 1; src <= source_size; ++ src)
      costs(src, 0) = utils::mathop::log(posterior_target_source(src, 0));
    
    alignment_source_target.clear();
    alignment_target_source.clear();
    
    if (span_source.empty() && span_target.empty())
      aligner(costs, insert_align(alignment_source_target, alignment_target_source));
    else
      aligner(costs, span_source, span_target, insert_align(alignment_source_target, alignment_target_source));
    
    std::sort(alignment_source_target.begin(), alignment_source_target.end());
    std::sort(alignment_target_source.begin(), alignment_target_source.end());
  }
  
  matrix_type costs;
  
  prob_set_type      prob_source_target;
  prob_set_type      prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;
  
  detail::ITGAlignment aligner;
};

struct MaxMatchModel1 : public ViterbiBase
{
  typedef utils::vector2<double, std::allocator<double> > matrix_type;
  typedef utils::vector2<double, std::allocator<double> > posterior_set_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;

  MaxMatchModel1(const ttable_type& __ttable_source_target,
		 const ttable_type& __ttable_target_source)
    : ViterbiBase(__ttable_source_target, __ttable_target_source) {}
  
  class insert_align
  {
    int source_size;
    int target_size;
    
    alignment_type& alignment_source_target;
    alignment_type& alignment_target_source;
    
  public:
    insert_align(const int& _source_size,
		 const int& _target_size,
		 alignment_type& __alignment_source_target,
		 alignment_type& __alignment_target_source)
      : source_size(_source_size), target_size(_target_size),
	alignment_source_target(__alignment_source_target),
	alignment_target_source(__alignment_target_source) {}
    
    template <typename Edge>
    insert_align& operator=(const Edge& edge)
    {	
      if (edge.first < source_size && edge.second < target_size) {
	alignment_source_target.push_back(edge);
	alignment_target_source.push_back(std::make_pair(edge.second, edge.first));
      }
      
      return *this;
    }
    
    insert_align& operator*() { return *this; }
    insert_align& operator++() { return *this; }
    insert_align operator++(int) { return *this; }
  };
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  const span_set_type& span_source,
		  const span_set_type& span_target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    operator()(source, target, alignment_source_target, alignment_target_source);
  }

  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const double prob_null  = p0;
    const double prob_align = 1.0 - p0;
    
    // we do not have to clearn!
    posterior_source_target.reserve(target_size + 1, source_size + 1);
    posterior_target_source.reserve(source_size + 1, target_size + 1);
    
    posterior_source_target.resize(target_size + 1, source_size + 1);
    posterior_target_source.resize(source_size + 1, target_size + 1);
    
    prob_source_target.reserve(source_size + 1);
    prob_target_source.reserve(target_size + 1);
    
    prob_source_target.resize(source_size + 1);
    prob_target_source.resize(target_size + 1);
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      const double prob_align_norm = 1.0 / source_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_source_target.begin();
      prob_set_type::iterator piter_end = prob_source_target.end();
      *piter = ttable_source_target(vocab_type::NONE, target[trg]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type src = 0; src != source_size; ++ src, ++ piter) {
	*piter = ttable_source_target(source[src], target[trg]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_source_target.begin();
      posterior_set_type::iterator siter = posterior_source_target.begin(trg + 1);
      for (/**/; piter != piter_end; ++ piter, ++ siter)
	(*siter) = (*piter) * factor;
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      const double prob_align_norm = 1.0 / target_size;
      double prob_sum = 0.0;
      
      prob_set_type::iterator piter     = prob_target_source.begin();
      prob_set_type::iterator piter_end = prob_target_source.end();
      *piter = ttable_target_source(vocab_type::NONE, target[src]) * prob_null;
      prob_sum += *piter;
      ++ piter;
      
      for (size_type trg = 0; trg != target_size; ++ trg, ++ piter) {
	*piter = ttable_target_source(target[trg], source[src]) * prob_align * prob_align_norm;
	prob_sum += *piter;
      }
      
      const double factor = 1.0 / prob_sum;
      
      piter = prob_target_source.begin();
      posterior_set_type::iterator titer = posterior_target_source.begin(src + 1);
      for (/**/; piter != piter_end; ++ piter, ++ titer)
	(*titer) = (*piter) * factor;
    }
    
    costs.clear();
    costs.resize(source_size + target_size, target_size + source_size, 0.0);
    
    for (size_type src = 0; src != source_size; ++ src)
      for (size_type trg = 0; trg != target_size; ++ trg) {
	costs(src, trg) = 0.5 * (utils::mathop::log(posterior_source_target(trg + 1, src + 1))
				 + utils::mathop::log(posterior_target_source(src + 1, trg + 1)));
	
	costs(src, trg + source_size) = utils::mathop::log(posterior_target_source(src + 1, 0));
	costs(src + target_size, trg) = utils::mathop::log(posterior_source_target(trg + 1, 0));
      }
    
    alignment_source_target.clear();
    alignment_target_source.clear();
    
    kuhn_munkres_assignment(costs, insert_align(source_size, target_size, alignment_source_target, alignment_target_source));
    
    std::sort(alignment_source_target.begin(), alignment_source_target.end());
    std::sort(alignment_target_source.begin(), alignment_target_source.end());
  }
  
  matrix_type costs;

  prob_set_type      prob_source_target;
  prob_set_type      prob_target_source;
  posterior_set_type posterior_source_target;
  posterior_set_type posterior_target_source;
};

#endif
