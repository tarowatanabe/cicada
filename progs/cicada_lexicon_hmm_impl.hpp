//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_HMM_IMPL__HPP__
#define __CICADA_LEXICON_HMM_IMPL__HPP__ 1

#include "cicada_lexicon_impl.hpp"

#include "utils/vector2_aligned.hpp"
#include "utils/vector3_aligned.hpp"
#include "utils/mathop.hpp"
#include "utils/aligned_allocator.hpp"
#include "utils/config.hpp"

#include "kuhn_munkres.hpp"
#include "itg_alignment.hpp"

struct LearnHMM : public LearnBase
{
  LearnHMM(const ttable_type& __ttable_source_target,
	   const ttable_type& __ttable_target_source,
	   const atable_type& __atable_source_target,
	   const atable_type& __atable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source, __atable_source_target, __atable_target_source) {}
  

  struct HMMData
  {
    typedef utils::vector2_aligned<prob_type, utils::aligned_allocator<prob_type> > forward_type;
    typedef utils::vector2_aligned<prob_type, utils::aligned_allocator<prob_type> > backward_type;
    
    typedef utils::vector2_aligned<prob_type, utils::aligned_allocator<prob_type> > emission_type;
    typedef utils::vector3_aligned<prob_type, utils::aligned_allocator<prob_type> > transition_type;
    
    typedef std::vector<prob_type, utils::aligned_allocator<prob_type> > scale_type;
    
    typedef utils::vector2_aligned<prob_type, utils::aligned_allocator<prob_type> > posterior_type;
    
    typedef std::vector<int, std::allocator<int> > point_set_type;
    typedef std::vector<point_set_type, std::allocator<point_set_type> > point_map_type;
    
    forward_type    forward;
    backward_type   backward;
    
    emission_type   emission;
    transition_type transition;
    
    scale_type scale;
    
    posterior_type posterior;
    
    sentence_type source;
    sentence_type target;
    sentence_type source_class;
    sentence_type target_class;
    
    point_map_type points;
    
    void prepare(const sentence_type& __source,
		 const sentence_type& __target,
		 const ttable_type& ttable,
		 const atable_type& atable,
		 const classes_type& classes_source,
		 const classes_type& classes_target)
    {
      const size_type source_size = __source.size();
      const size_type target_size = __target.size();
      
      const double prob_null  = p0;
      const double prob_align = 1.0 - p0;

      source.clear();
      target.clear();
      
      source.reserve((source_size + 2) * 2);
      target.reserve((target_size + 2) * 2);
      source.resize((source_size + 2) * 2, vocab_type::NONE);
      target.resize((target_size + 2) * 2, vocab_type::NONE);
      
      source[0] = vocab_type::BOS;
      target[0] = vocab_type::BOS;
      source[(source_size + 2) - 1] = vocab_type::EOS;
      target[(target_size + 2) - 1] = vocab_type::EOS;
      
      std::copy(__source.begin(), __source.end(), source.begin() + 1);
      std::copy(__target.begin(), __target.end(), target.begin() + 1);
      
      source_class.reserve(source_size + 2);
      target_class.reserve(target_size + 2);
      source_class.resize(source_size + 2);
      target_class.resize(target_size + 2);
      
      source_class[0] = vocab_type::BOS;
      target_class[0] = vocab_type::BOS;
      source_class[(source_size + 2) - 1] = vocab_type::EOS;
      target_class[(target_size + 2) - 1] = vocab_type::EOS;
      
      sentence_type::iterator csiter = source_class.begin() + 1;
      for (sentence_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter, ++ csiter)
	*csiter = classes_source[*siter];
      
      sentence_type::iterator ctiter = target_class.begin() + 1;
      for (sentence_type::const_iterator titer = target.begin(); titer != target.end(); ++ titer, ++ ctiter)
	*ctiter = classes_target[*titer];
      
      
      emission.clear();
      transition.clear();
      
      emission.reserve(target_size + 2, (source_size + 2) * 2);
      transition.reserve(target_size + 2, (source_size + 2) * 2, (source_size + 2) * 2);
      
      emission.resize(target_size + 2, (source_size + 2) * 2, 0.0);
      transition.resize(target_size + 2, (source_size + 2) * 2, (source_size + 2) * 2, 0.0);
      
      // compute emission table...
      emission(0, 0) = 1.0;
      emission(target_size + 2 - 1, source_size + 2 - 1) = 1.0;
      for (int trg = 1; trg <= target_size; ++ trg) {
	
	// translation into non-NULL word
	prob_type* eiter = &(*emission.begin(trg)) + 1;
	for (int src = 1; src <= source_size; ++ src, ++ eiter)
	  (*eiter) = ttable(source[src], target[trg]);
	
	prob_type* eiter_first = &(*emission.begin(trg)) + source_size + 2;
	prob_type* eiter_last  = eiter_first + source_size + 2 - 1; // -1 to exclude EOS
	
	std::fill(eiter_first, eiter_last, ttable(vocab_type::NONE, target[trg]));
      }

      
      // compute transition table...
      // we start from 1, since there exists no previously aligned word before BOS...
      for (int trg = 1; trg < target_size + 2; ++ trg) {
	// alignment into non-null
	// start from 1 to exlude transition into <s>
	for (int next = 1; next < source_size + 2; ++ next) {
	  prob_type* titer1 = &(*transition.begin(trg, next)); // from word
	  prob_type* titer2 = titer1 + (source_size + 2);      // from NONE
	  // - 1 to exclude previous </s>...
	  for (int prev = 0; prev < source_size + 2 - 1; ++ prev, ++ titer1, ++ titer2) {
	    const prob_type prob = prob_align * atable(source_class[prev], target_class[trg],
						       source_size, target_size,
						       prev - 1, next - 1);
	    
	    *titer1 = prob;
	    *titer2 = prob;
	  }
	}
	
	// null transition
	// we will exclude EOS
	for (int next = 0; next < (source_size + 2) - 1; ++ next) {
	  const int next_none = next + (source_size + 2);
	  const int prev1_none = next;
	  const int prev2_none = next + (source_size + 2);
	  
	  transition(trg, next_none, prev1_none) = prob_null;
	  transition(trg, next_none, prev2_none) = prob_null;
	}
      }
    }
    
    void forward_backward(const sentence_type& __source,
			  const sentence_type& __target)
    {
      const size_type source_size = __source.size();
      const size_type target_size = __target.size();
      
      forward.clear();
      backward.clear();
      scale.clear();
      
      forward.reserve(target_size + 2, (source_size + 2) * 2);
      backward.reserve(target_size + 2, (source_size + 2) * 2);
      scale.reserve(target_size + 2);

      forward.resize(target_size + 2, (source_size + 2) * 2, 0.0);
      backward.resize(target_size + 2, (source_size + 2) * 2, 0.0);
      scale.resize(target_size + 2, 1.0);
      
      
      forward(0, 0) = 1.0;
      for (int trg = 1; trg < target_size + 2; ++ trg) {
	
	// +1 to exclude BOS
	prob_type*       niter = &(*forward.begin(trg)) + 1;
	const prob_type* eiter = &(*emission.begin(trg)) + 1;
	for (int next = 1; next < source_size + 2; ++ next, ++ niter, ++ eiter) {
	  
	  const prob_type* piter = &(*forward.begin(trg - 1));
	  const prob_type* titer = &(*transition.begin(trg, next));
	  
#if defined(__SSE2__) && defined(HAVE_EMMINTRIN_H)
	  const prob_type* piter_end = piter + (source_size + 2) * 2;
	  
	  const double factor = *eiter;
	  const __m128d factors = _mm_set_pd(factor, factor);
	  __m128d accum = _mm_setzero_pd();
	  
	  while (piter != piter_end) {
	    __m128d prob_prevs = _mm_load_pd(piter);
	    __m128d prob_trans = _mm_load_pd(titer);
	    
	    __m128d tmp = _mm_mul_pd(prob_prevs, factors);
	    tmp = _mm_mul_pd(tmp, prob_trans);
	    
	    accum = _mm_add_pd(accum, tmp);
	    
	    piter += 2;
	    titer += 2;
	  }
	  
	  double result[2] __attribute__((aligned(16)));
	  _mm_store_pd(result, accum);
	  *niter += result[0] + result[1];
#else
	  const double factor = *eiter;
	  if (factor > 0.0) {
	    
	    double accum[4] = {0.0, 0.0, 0.0, 0.0};
	    const int loop_size = (source_size + 2) * 2;
	    for (int prev = 0; prev < loop_size - 3; prev += 4) {
	      accum[0] += piter[prev + 0] * titer[prev + 0] * factor;
	      accum[1] += piter[prev + 1] * titer[prev + 1] * factor;
	      accum[2] += piter[prev + 2] * titer[prev + 2] * factor;
	      accum[3] += piter[prev + 3] * titer[prev + 3] * factor;
	    }
	    
	    switch (loop_size & 0x03) {
	    case 3: accum[4 - 3] += piter[loop_size - 3] * titer[loop_size - 3] * factor;
	    case 2: accum[4 - 2] += piter[loop_size - 2] * titer[loop_size - 2] * factor;
	    case 1: accum[4 - 1] += piter[loop_size - 1] * titer[loop_size - 1] * factor;
	    }
	    *niter += accum[0] + accum[1] + accum[2] + accum[3];
	  }
#endif
	  
#if 0
	  // -1 to exclude EOS...
	  double next_word = 0.0;
	  double next_none = 0.0;
	  for (int prev = 0; prev < (source_size + 2) - 1; ++ prev) {
	    next_word += piter[prev] * titer[prev] * factor;
	    next_none += piter[prev + source_size + 2] * titer[prev + source_size + 2] * factor;
	  }
	  *niter += next_word + next_none;
#endif
	  
#if 0
	  for (int prev = 0; prev < (source_size + 2) * 2; ++ prev, ++ piter, ++ titer)
	    *niter += (*piter) * (*eiter) * (*titer);
#endif
	}
	
	prob_type*       niter_none = &(*forward.begin(trg)) + (source_size + 2);
	const prob_type* eiter_none = &(*emission.begin(trg)) + (source_size + 2);
	const prob_type* piter_none1 = &(*forward.begin(trg - 1));
	const prob_type* piter_none2 = piter_none1 + (source_size + 2);
	for (int next = 0; next < source_size + 2 - 1; ++ next, ++ niter_none, ++ eiter_none, ++ piter_none1, ++ piter_none2) {
	  const int next_none = next + source_size + 2;
	  const int prev_none1 = next;
	  const int prev_none2 = next + source_size + 2;
	  
	  *niter_none += (*piter_none1) * (*eiter_none) * transition(trg, next_none, prev_none1);
	  *niter_none += (*piter_none2) * (*eiter_none) * transition(trg, next_none, prev_none2);
	  
#if 0
	  forward(trg, next_none) += forward(trg - 1, prev_none1) * emission(trg, next_none) * transition(trg, next_none, prev_none1);
	  forward(trg, next_none) += forward(trg - 1, prev_none2) * emission(trg, next_none) * transition(trg, next_none, prev_none2);
#endif
	}
	
#if defined(__SSE2__) && defined(HAVE_EMMINTRIN_H)
	
	__m128d accum = _mm_setzero_pd();
	prob_type* piter_start = &(*forward.begin(trg));
	prob_type* piter_end   = &(*forward.end(trg));
	
	for (prob_type* piter = piter_start; piter != piter_end; piter += 2) {
	  __m128d probs = _mm_load_pd(piter);
	  accum = _mm_add_pd(accum, probs);
	}
	
	double results[2] __attribute__((aligned(16)));
	_mm_store_pd(results, accum);
	double result = results[0] + results[1];
	result = (result == 0.0 ? 1.0 : 1.0 / result);
	
	if (result != 1.0) {
	  const __m128d factors = _mm_set_pd(result, result);
	  for (prob_type* piter = piter_start; piter != piter_end; piter += 2) {
	    __m128d tmp = _mm_load_pd(piter);
	    tmp = _mm_mul_pd(tmp, factors);
	    _mm_store_pd(piter, tmp);
	  }
	}
	
	scale[trg] = result;
#else
	scale[trg] = std::accumulate(forward.begin(trg), forward.end(trg), 0.0);
	scale[trg] = (scale[trg] == 0.0 ? 1.0 : 1.0 / scale[trg]);
	if (scale[trg] != 1.0)
	  std::transform(forward.begin(trg), forward.end(trg), forward.begin(trg), std::bind2nd(std::multiplies<double>(), scale[trg]));
#endif
      }
      
      backward(target_size + 2 - 1, source_size + 2 - 1) = 1.0;
      for (int trg = target_size + 2 - 2; trg >= 0; -- trg) {
	
	const prob_type scale = scale[trg];
	
	// +1 to exclude BOS
	const prob_type* niter = &(*backward.begin(trg + 1)) + 1;
	const prob_type* eiter = &(*emission.begin(trg + 1)) + 1;
	for (int next = 1; next < source_size + 2; ++ next, ++ niter, ++ eiter) {
	  
	  prob_type*       piter = &(*backward.begin(trg));
	  const prob_type* titer = &(*transition.begin(trg + 1, next));
	  
#if defined(__SSE2__) && defined(HAVE_EMMINTRIN_H)
	  prob_type* piter_end = piter + (source_size + 2) * 2;
	  
	  const double factor = (*eiter) * (*niter) * scale;
	  
	  const __m128d factors = _mm_set_pd(factor, factor);
	  
	  while (piter != piter_end) {
	    __m128d prob_trans = _mm_load_pd(titer);
	    __m128d prob_prevs = _mm_load_pd(piter);
	    
	    __m128d tmp = _mm_mul_pd(prob_trans, factors);
	    tmp = _mm_add_pd(tmp, prob_prevs);
	    _mm_store_pd(piter, tmp);
	    
	    piter += 2;
	    titer += 2;
	  }
#else
	  const double factor = (*eiter) * (*niter) * scale;
	  if (factor > 0.0) {
	  
	    const int loop_size = (source_size + 2) * 2;
	    for (int prev = 0; prev < loop_size - 3; prev += 4) {
	      piter[prev + 0] += titer[prev + 0] * factor;
	      piter[prev + 1] += titer[prev + 1] * factor;
	      piter[prev + 2] += titer[prev + 2] * factor;
	      piter[prev + 3] += titer[prev + 3] * factor;
	    }
	    switch (loop_size & 0x03) {
	    case 3: piter[loop_size - 3] += titer[loop_size - 3] * factor;
	    case 2: piter[loop_size - 2] += titer[loop_size - 2] * factor;
	    case 1: piter[loop_size - 1] += titer[loop_size - 1] * factor;
	    }
	  }
#endif
#if 0
	  // -1 to exclude EOS...
	  for (int prev = 0; prev < (source_size + 2) - 1; ++ prev) {
	    piter[prev] += titer[prev] * factor;
	    piter[prev + source_size + 2] += titer[prev + source_size + 2] * factor;
	  }
#endif
	  
#if 0
	  for (int prev = 0; prev < (source_size + 2) * 2; ++ prev)
	    backward(trg, prev) += backward(trg + 1, next) * emission(trg + 1, next) * transition(trg + 1, next, prev) * scale;
#endif
	}
	
	const prob_type* niter_none = &(*backward.begin(trg + 1)) + (source_size + 2);
	const prob_type* eiter_none = &(*emission.begin(trg + 1)) + (source_size + 2);
	prob_type*       piter_none1 = &(*backward.begin(trg));
	prob_type*       piter_none2 = piter_none1 + (source_size + 2);
	for (int next = 0; next < source_size + 2 - 1; ++ next, ++ niter_none, ++ eiter_none, ++ piter_none1, ++ piter_none2) {
	  const int next_none = next + source_size + 2;
	  const int prev_none1 = next;
	  const int prev_none2 = next + source_size + 2;
	  
	  *piter_none1 += (*niter_none) * (*eiter_none) * transition(trg + 1, next_none, prev_none1) * scale;
	  *piter_none2 += (*niter_none) * (*eiter_none) * transition(trg + 1, next_none, prev_none2) * scale;
	  
#if 0
	  backward(trg, prev_none1) += backward(trg + 1, next_none) * emission(trg + 1, next_none) * transition(trg + 1, next_none, prev_none1) * scale;
	  backward(trg, prev_none2) += backward(trg + 1, next_none) * emission(trg + 1, next_none) * transition(trg + 1, next_none, prev_none2) * scale;
#endif
	}
      }
    }

    void estimate_posterior(const sentence_type& __source,
			    const sentence_type& __target)
    {
      const size_type source_size = __source.size();
      const size_type target_size = __target.size();
      
      posterior.clear();
      posterior.reserve(target_size + 1, source_size + 1);
      posterior.resize(target_size + 1, source_size + 1, 0.0);
      
      const prob_type sum = forward(target_size + 2 - 1, source_size + 2 - 1);
      
      for (int trg = 1; trg <= target_size; ++ trg) {
	const double factor = 1.0 / (scale[trg] * sum);
	
	// + 1 to skip BOS
	const prob_type* fiter = &(*forward.begin(trg)) + 1;
	const prob_type* biter = &(*backward.begin(trg)) + 1;
	prob_type* piter = &(*posterior.begin(trg)) + 1;
	for (int src = 1; src <= source_size; ++ src, ++ fiter, ++ biter, ++ piter) {
	  const prob_type count = (*fiter) * (*biter) * factor;
	  
	  if (std::isfinite(count) && count > 0.0)
	    (*piter) += count;
	}
	
	fiter = &(*forward.begin(trg)) + source_size + 2;
	biter = &(*backward.begin(trg)) + source_size + 2;
	prob_type count_none = 0.0;
	for (int src = 0; src < source_size + 2; ++ src, ++ fiter, ++ biter) {
	  const prob_type count = (*fiter) * (*biter) * factor;
	  
	  if (std::isfinite(count) && count > 0.0)
	    count_none += count;
	}
	posterior(trg, 0) = count_none;
      }
    }

    struct logminus
    {
      double operator()(const double& logsum, const double& prob) const
      {
	return logsum - utils::mathop::log(prob);
      }
    };
    
    double objective() const
    {
      return std::accumulate(scale.begin(), scale.end(), 0.0, logminus());
    }

    void accumulate(const sentence_type& __source,
		    const sentence_type& __target,
		    ttable_type& counts)
    {
      const size_type source_size = __source.size();
      const size_type target_size = __target.size();
      
      const prob_type sum = forward(target_size + 2 - 1, source_size + 2 - 1);
      
      // accumulate lexcion
      for (int trg = 1; trg <= target_size; ++ trg) {
	const double factor = 1.0 / (scale[trg] * sum);
      
	// + 1 to skip BOS
	const prob_type* fiter = &(*forward.begin(trg)) + 1;
	const prob_type* biter = &(*backward.begin(trg)) + 1;
	for (int src = 1; src < source_size + 2; ++ src, ++ fiter, ++ biter) {
	  const prob_type count = (*fiter) * (*biter) * factor;
	  
	  if (std::isfinite(count) && count > 0.0)
	    counts[source[src]][target[trg]] += count;
	}
      
	// null alignment...
	double count_none = 0.0;
	for (int src = 0; src < source_size + 2; ++ src, ++ fiter, ++ biter) {
	  const prob_type count = (*fiter) * (*biter) * factor;
	  if (std::isfinite(count) && count > 0.0)
	    count_none += count;
	}
      
	counts[vocab_type::NONE][target[trg]] += count_none;
      }
    }
    
    void accumulate(const sentence_type& __source,
		    const sentence_type& __target,
		    atable_type& counts)
    {
      typedef std::vector<atable_type::mapped_type*, std::allocator<atable_type::mapped_type*> > mapped_type;

      const size_type source_size = __source.size();
      const size_type target_size = __target.size();
      
      const double sum = forward(target_size + 2 - 1, source_size + 2 - 1);
      const double factor = 1.0 / sum;

      mapped_type mapped((source_size + 2) - 1);
      
      for (int trg = 1; trg < target_size + 2; ++ trg) {
	const prob_type* biter = &(*data.backward.begin(trg)) + 1;
	const prob_type* eiter = &(*data.emission.begin(trg)) + 1;
	
	for (int prev = 0; prev < (source_size + 2) - 1; ++ prev)
	  mapped[prev] = &(counts[std::make_pair(source_class[prev], target_class[trg])]);
	
	// + 1, we will exclude <s>, since we will never aligned to <s>
	for (int next = 1; next < source_size + 2; ++ next, ++ biter, ++ eiter) {
	  const prob_type factor_backward = (*biter) * (*eiter);
	  
	  if (factor_backward > 0.0) {
	    const prob_type* fiter_word = &(*data.forward.begin(trg - 1));
	    const prob_type* titer_word = &(*data.transition.begin(trg, next));
	    const prob_type* fiter_none = &(*data.forward.begin(trg - 1)) + (source_size + 2);
	    const prob_type* titer_none = &(*data.transition.begin(trg, next)) + (source_size + 2);
	    
	    // - 1 to exlude EOS
	    for (int prev = 0; prev < (source_size + 2) - 1; ++ prev, ++ fiter_word, ++ titer_word, ++ fiter_none, ++ titer_none) {
	      const double count_word = (*fiter_word) * factor_backward * (*titer_word) * factor;
	      const double count_none = (*fiter_none) * factor_backward * (*titer_none) * factor;
	      
	      mapped[prev]->operator[](next - prev) += ((count_word > 0.0 && std::isfinite(count_word) ? count_word : 0.0)
							+ (count_none > 0.0 && std::isfinite(count_none) ? count_none : 0.0));
	    }
	  }
	}
      }
    }
  };

  typedef HMMData hmm_data_type;
  
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const ttable_type& ttable,
	     const atable_type& atable,
	     const classes_type& classes_source,
	     const classes_type& classes_target,
	     ttable_type& counts_ttable,
	     atable_type& counts_atable,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    hmm.prepare(source, target, ttable, atable, classes_source, classes_target);
    
    hmm.forward_backward(source, target);
    
    objective += hmm.objective() / target_size;
    
    hmm.accumulate(source, target, counts_ttable);
    
    hmm.accumulate(source, target, counts_atable);
  }
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    learn(source,
	  target,
	  ttable_source_target,
	  atable_source_target,
	  classes_source,
	  classes_target,
	  ttable_counts_source_target,
	  atable_counts_source_target,
	  aligned_source_target,
	  objective_source_target);
    learn(target,
	  source,
	  ttable_target_source,
	  atable_target_source,
	  classes_target,
	  classes_source,
	  ttable_counts_target_source,
	  atable_counts_target_source,
	  aligned_target_source,
	  objective_target_source);
  }

  hmm_data_type hmm;
};

struct LearnHMMPosterior : public LearnBase
{
  LearnHMMPosterior(const ttable_type& __ttable_source_target,
		    const ttable_type& __ttable_target_source,
		    const atable_type& __atable_source_target,
		    const atable_type& __atable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source, __atable_source_target, __atable_target_source) {}

  typedef LearnHMM::hmm_data_type hmm_data_type;
  
  typedef std::vector<prob_type, std::allocator<prob_type> > prob_set_type;
  
  void learn(const sentence_type& source,
	     const sentence_type& target,
	     const ttable_type& ttable,
	     const atable_type& atable,
	     const classes_type& classes_source,
	     const classes_type& classes_target,
	     ttable_type& counts_ttable,
	     atable_type& counts_atable,
	     aligned_type& aligned,
	     double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    hmm.prepare(source, target, ttable, atable, classes_source, classes_target);
    
    hmm.forward_backward(source, target);
    
    objective += hmm.objective() / target_size;
    
    phi.clear();
    phi.resize(source_size + 1, 0.0);
    
    exp_phi_old.clear();
    exp_phi_old.resize(source_size + 1, 1.0);
    
    for (int iter = 0; iter < 5; ++ iter) {
      hmm.compute_posterior(source, target);
      
      exp_phi.clear();
      exp_phi.resize(source_size + 1, 1.0);
      
      size_type count_zero = 0;
      for (int src = 1; src <= source_size; ++ src) {
	double sum_posterior = 0.0;
	for (int trg = 1; trg <= target_size; ++ trg)
	  sum_posterior += hmm.posterior(trg, src);
	
	phi[src] += 1.0 - sum_posterior;
	if (phi[src] > 0.0)
	  phi[src] = 0.0;
	
	count_zero += (phi[src] == 0.0);
	exp_phi[src] = utils::mathop::exp(phi[src]);
      }
      
      if (count_zero == source_size) break;
      
      // rescale emission table...
      for (int trg = 1; trg <= target_size; ++ trg) {
	// translation into non-NULL word
	prob_type* eiter = &(*data.emission.begin(trg)) + 1;
	for (int src = 1; src <= source_size; ++ src, ++ eiter)
	  (*eiter) *= exp_phi[src] / exp_phi_old[src];
      }
      // swap...
      exp_phi_old.swap(exp_phi);
      
      hmm.forward_backward(source, target);
    }
    
    hmm.accumulate(source, target, counts_ttable);
    
    hmm.accumulate(source, target, counts_atable);
  }

  void operator()(const sentence_type& source, const sentence_type& target)
  {
    learn(source,
	  target,
	  ttable_source_target,
	  atable_source_target,
	  classes_source,
	  classes_target,
	  ttable_counts_source_target,
	  atable_counts_source_target,
	  aligned_source_target,
	  objective_source_target);
    learn(target,
	  source,
	  ttable_target_source,
	  atable_target_source,
	  classes_target,
	  classes_source,
	  ttable_counts_target_source,
	  atable_counts_target_source,
	  aligned_target_source,
	  objective_target_source);
  }

  hmm_data_type hmm;
  
  prob_set_type phi;
  prob_set_type exp_phi;
  prob_set_type exp_phi_old;
};

struct LearnModel1Symmetric : public LearnBase
{
  LearnModel1Symmetric(const ttable_type& __ttable_source_target,
		       const ttable_type& __ttable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source) {}
      
  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  
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
};

struct LearnModel1SymmetricPosterior : public LearnBase
{
  LearnModel1SymmetricPosterior(const ttable_type& __ttable_source_target,
				const ttable_type& __ttable_target_source)
    : LearnBase(__ttable_source_target, __ttable_target_source) {}

  typedef utils::vector2<prob_type, std::allocator<prob_type> > posterior_set_type;
  
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
