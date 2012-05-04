// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __ITG_ALIGNMENT__HPP__
#define __ITG_ALIGNMENT__HPP__ 1

#include <vector>

#include <cicada/span_vector.hpp>

#include <utils/chart.hpp>
#include <utils/bichart.hpp>
#include <utils/vector2.hpp>
#include <utils/mathop.hpp>
#include <utils/bithack.hpp>

namespace detail
{
#if 1
  struct ITGAlignment
  {
    static const double threshold = 1e-2;
    
    
    typedef cicada::SpanVector       span_set_type;
    typedef span_set_type::span_type span_type;
    
    struct span_pair_type
    {
      span_type source;
      span_type target;
      
      span_pair_type() : source(), target() {}
      span_pair_type(const span_type& __source, const span_type& __target) : source(__source), target(__target) {}
      span_pair_type(const int source_first, const int source_last,
		     const int target_first, const int target_last)
	: source(source_first, source_last), target(target_first, target_last) {}

      size_t size() const { return source.size() + target.size(); }
      bool empty() const { return source.empty() && target.empty(); }
      
      friend
      bool operator==(const span_pair_type& x, const span_pair_type& y)
      {
	return x.source == y.source && x.target == y.target;
      }
    };
    
    typedef double logprob_type;
    
    typedef utils::bichart<logprob_type, std::allocator<logprob_type> > chart_type;

    typedef utils::chart<logprob_type, std::allocator<logprob_type> > chart_mono_type;
    typedef std::vector<logprob_type, std::allocator<logprob_type> > alpha_type;
    typedef std::vector<logprob_type, std::allocator<logprob_type> > beta_type;
    
    typedef std::pair<span_pair_type, span_pair_type> backptr_type;
    typedef utils::bichart<backptr_type, std::allocator<backptr_type> > backptr_chart_type;
    
    typedef std::vector<span_pair_type, std::allocator<span_pair_type> > stack_type;
    
    typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
    typedef std::vector<span_pair_set_type, std::allocator<span_pair_set_type> > agenda_type;

    typedef std::pair<logprob_type, span_pair_type> score_span_pair_type;
    typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;

    // sort by less so that we can pop from a greater score item.
    struct heap_compare
    {
      bool operator()(const score_span_pair_type& x, const score_span_pair_type& y) const
      {
	return x.first < y.first;
      }
    };
    
    
    struct PruneNone
    {
      bool operator()(const int source_first, const int source_last,
		      const int target_first, const int target_last) const
      {
	return false;
      }
    };
    
    struct PruneSpan
    {
      typedef utils::chart<bool, std::allocator<bool> > span_prune_type;

      template <typename Costs>
      PruneSpan(const Costs& costs, const span_set_type& spans_source, const span_set_type& spans_target)
      {
	const int source_size = costs.size1() - 1;
	const int target_size = costs.size2() - 1;
	
	prune_source.clear();
	prune_target.clear();
	
	prune_source.reserve(source_size + 1);
	prune_target.reserve(target_size + 1);
	
	prune_source.resize(source_size + 1, ! spans_source.empty());
	prune_target.resize(target_size + 1, ! spans_target.empty());
	
	if (! spans_source.empty()) {
	  span_set_type::const_iterator siter_end = spans_source.end();
	  for (span_set_type::const_iterator siter = spans_source.begin(); siter != siter_end; ++ siter)
	    prune_source(siter->first, siter->last) = false;
	}
	
	if (! spans_target.empty()) {
	  span_set_type::const_iterator siter_end = spans_target.end();
	  for (span_set_type::const_iterator siter = spans_target.begin(); siter != siter_end; ++ siter)
	    prune_target(siter->first, siter->last) = false;
	}
      }
      
      bool operator()(const int source_first, const int source_last,
		      const int target_first, const int target_last) const
      {
	return prune_source(source_first, source_last) || prune_target(target_first, target_last);
      }
      
      span_prune_type prune_source;
      span_prune_type prune_target;
    };

    template <typename Costs>
    void initialize(const Costs& costs)
    {
      const int source_size = costs.size1() - 1;
      const int target_size = costs.size2() - 1;
      
      const double infinity = - std::numeric_limits<double>::infinity();
      
      const backptr_type backptr_invalid(span_pair_type(span_type(-1, -1), span_type(-1, -1)),
					 span_pair_type(span_type(-1, -1), span_type(-1, -1)));

      agenda.clear();
      agenda.resize(source_size + target_size + 1);
      
      chart.clear();
      chart.reserve(source_size + 1, target_size + 1);
      chart.resize(source_size + 1, target_size + 1, infinity);
      
      backptr.clear();
      backptr.reserve(source_size + 1, target_size + 1);
      backptr.resize(source_size + 1, target_size + 1, backptr_invalid);
      
      chart_source.clear();
      chart_source.reserve(source_size + 1);
      chart_source.resize(source_size + 1, infinity);
      
      chart_target.clear();
      chart_target.reserve(target_size + 1);
      chart_target.resize(target_size + 1, infinity);
      
      alpha_source.clear();
      alpha_target.clear();
      beta_source.clear();
      beta_target.clear();
      
      alpha_source.resize(source_size + 1, infinity);
      alpha_target.resize(target_size + 1, infinity);
      beta_source.resize(source_size + 1, infinity);
      beta_target.resize(target_size + 1, infinity);
      
      for (int src = 0; src <= source_size; ++ src)
	for (int trg = 0; trg <= target_size; ++ trg) {
	  
	  if (src < source_size && trg < target_size) {
	    // one-to-one alignment
	    const logprob_type score = costs(src + 1, trg + 1);
	    
	    chart(src, src + 1, trg, trg + 1) = score;
	    chart_source(src, src + 1) = std::max(chart_source(src, src + 1), score);
	    chart_target(trg, trg + 1) = std::max(chart_target(trg, trg + 1), score);
	    
	    agenda[2].push_back(span_pair_type(src, src + 1, trg, trg + 1));
	    
	    // one-to-many alignment...
	    for (int src_last = src + 2; src_last <= source_size; ++ src_last) {
	      const int trg_last = trg + 1;
	      
	      const logprob_type score = chart(src, src_last - 1, trg, trg_last) + costs(src_last, trg_last);
	      
	      chart(src, src_last, trg, trg_last) = score;
	      
	      chart_source(src, src_last) = std::max(chart_source(src, src_last), score);
	      chart_target(trg, trg_last) = std::max(chart_target(trg, trg_last), score);
	      
	      agenda[src_last - src + 1].push_back(span_pair_type(src, src_last, trg, trg_last));
	    }
	    
	    for (int trg_last = trg + 2; trg_last <= target_size; ++ trg_last) {
	      const int src_last = src + 1;
	      
	      const logprob_type score = chart(src, src_last, trg, trg_last - 1) + costs(src_last, trg_last);
	      
	      chart(src, src_last, trg, trg_last) = score;
	      
	      chart_source(src, src_last) = std::max(chart_source(src, src_last), score);
	      chart_target(trg, trg_last) = std::max(chart_target(trg, trg_last), score);
	      
	      agenda[1 + trg_last - trg].push_back(span_pair_type(src, src_last, trg, trg_last));
	    }
	  }
	  
	  // null-to-many
	  if (src < source_size) {
	    logprob_type score = 0.0;
	    
	    const int src_first = src;
	    for (int src_last = src_first + 1; src_last <= source_size; ++ src_last) {
	      score += costs(src_last, 0);
	      
	      chart(src_first, src_last, trg, trg) = score;
	      
	      chart_source(src_first, src_last) = std::max(chart_source(src_first, src_last), score);
	      
	      agenda[src_last - src_first].push_back(span_pair_type(src_first, src_last, trg, trg));
	    }
	  }
	  
	  if (trg < target_size) {
	    logprob_type score = 0.0;
	    
	    const int trg_first = trg;
	    for (int trg_last = trg_first + 1; trg_last <= target_size; ++ trg_last) {
	      score += costs(0, trg_last);
	      
	      chart(src, src, trg_first, trg_last) = score;
	      
	      chart_target(trg_first, trg_last) = std::max(chart_target(trg_first, trg_last), score);
	      
	      agenda[trg_last - trg_first].push_back(span_pair_type(src, src, trg_first, trg_last));
	    }
	  }
	}

      // forward-backward to compute estiamtes...
      forward_backward(chart_source, alpha_source, beta_source);
      forward_backward(chart_target, alpha_target, beta_target);
    }

    void forward_backward(const chart_mono_type& chart, alpha_type& alpha, beta_type& beta)
    {
      const size_t sentence_size = chart.size() - 1;
      
      // forward...
      alpha[0] = 0;
      for (size_t last = 1; last <= sentence_size; ++ last)
	for (size_t first = 0; first != last; ++ first)
	  alpha[last] = std::max(alpha[last], alpha[first] + chart(first, last));
      
      // backward...
      beta[sentence_size] = 0;
      for (ptrdiff_t first = sentence_size - 1; first >= 0; -- first)
	for (size_t last = first + 1; last <= sentence_size; ++ last)
	  beta[first] = std::max(beta[first], chart(first, last) + beta[last]);
    }
    
    template <typename Costs, typename Prune>
    void construct(const Costs& costs, const Prune& prune)
    {
      const int source_size = costs.size1() - 1;
      const int target_size = costs.size2() - 1;
      
      const double infinity = - std::numeric_limits<double>::infinity();
      
      // use of the new parsing strategy without span-pruning!
      
      const int T = source_size;
      const int V = target_size;
      const int L = T + V;
      
      double beam = std::log(threshold);
      do {
	
	for (int l = 1; l != L; ++ l) {
	  span_pair_set_type& spans = agenda[l];
	  
	  heap.clear();
	  heap.reserve(spans.size());
	  
	  span_pair_set_type::const_iterator siter_end = spans.end();
	  for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)  {
	    const double score = (chart(siter->source.first, siter->source.last, siter->target.first, siter->target.last)
				  + std::min(alpha_source[siter->source.first] + beta_source[siter->source.last],
					     alpha_target[siter->target.first] + beta_target[siter->target.last]));
	    
	    heap.push_back(score_span_pair_type(score, *siter));
	    std::push_heap(heap.begin(), heap.end(), heap_compare());
	  }

	  heap_type::iterator hiter_begin = heap.begin();
	  heap_type::iterator hiter       = heap.end();
	  heap_type::iterator hiter_end   = heap.end();
	  
	  const double score_threshold = hiter_begin->first + beam;
	  for (/**/; hiter_begin != hiter && hiter_begin->first > score_threshold; -- hiter)
	    std::pop_heap(hiter_begin, hiter, heap_compare());
	  
	  for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	    const span_pair_type& span_pair = iter->second;
	    
	    const int s = span_pair.source.first;
	    const int t = span_pair.source.last;
	    const int u = span_pair.target.first;
	    const int v = span_pair.target.last;
	    
	    for (int S = utils::bithack::max(s - l, 0); S <= s; ++ S) {
	      const int L = l - (s - S);
	      
	      // straight
	      for (int U = utils::bithack::max(u - L, 0); U <= u - (S == s); ++ U) {
		// parent span: StUv
		// span1: SsUu
		// span2: stuv
		
		if (chart(S, s, U, u) == infinity || prune(S, t, U, v)) continue;
		
		const span_pair_type  span1(S, s, U, u);
		const span_pair_type& span2(span_pair);
		const span_pair_type  span_head(S, t, U, v);
		
		logprob_type& value = chart(S, t, U, v);
		backptr_type& back  = backptr(S, t, U, v);

		if (value == infinity)
		  agenda[span_head.size()].push_back(span_head);
		
		const logprob_type score = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		
		if (score > value) {
		  value = score;
		  
		  back.first =  span1;
		  back.second = span2;
		}
	      }
	      
	      // inversion
	      for (int U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
		// parent span: StuU
		// span1: SsvU
		// span2: stuv
		
		if (chart(S, s, v, U) == infinity || prune(S, t, u, U)) continue;
		
		const span_pair_type  span1(S, s, v, U);
		const span_pair_type& span2(span_pair);
		const span_pair_type  span_head(S, t, u, U);
		
		logprob_type& value = chart(S, t, u, U);
		backptr_type& back  = backptr(S, t, u, U);
		
		if (value == infinity)
		  agenda[span_head.size()].push_back(span_head);
		
		const logprob_type score = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		
		if (score > value) {
		  value = score;
		  
		  back.first =  span1;
		  back.second = span2;
		}
	      }
	    }
	    
	    for (int S = t; S <= utils::bithack::min(t + l, T); ++ S) {
	      const int L = l - (S - t);
		
	      // inversion
	      for (int U = utils::bithack::max(u - L, 0); U <= u - (S == t); ++ U) {
		// parent span: sSUv
		// span1: stuv
		// span2: tSUu
		
		if (chart(t, S, U, u) == infinity || prune(s, S, U, v)) continue;
	      
		const span_pair_type& span1(span_pair);
		const span_pair_type  span2(t, S, U, u);
		const span_pair_type  span_head(s, S, U, v);
		
		logprob_type& value = chart(s, S, U, v);
		backptr_type& back  = backptr(s, S, U, v);
		
		if (value == infinity)
		  agenda[span_head.size()].push_back(span_head);
		
		const logprob_type score = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		
		if (score > value) {
		  value = score;
		  
		  back.first =  span1;
		  back.second = span2;
		}
	      }
		
	      // straight
	      for (int U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
		// parent span: sSuU
		// span1: stuv
		// span2: tSvU
		
		if (chart(t, S, v, U) == infinity || prune(s, S, u, U)) continue;
		  
		const span_pair_type& span1(span_pair);
		const span_pair_type  span2(t, S, v, U);
		const span_pair_type  span_head(s, S, u, U);
		
		logprob_type& value = chart(s, S, u, U);
		backptr_type& back  = backptr(s, S, u, U);
		
		if (value == infinity)
		  agenda[span_head.size()].push_back(span_head);
		
		const logprob_type score = (chart(span1.source.first, span1.source.last, span1.target.first, span1.target.last)
					    + chart(span2.source.first, span2.source.last, span2.target.first, span2.target.last));
		
		if (score > value) {
		  value = score;
		  
		  back.first =  span1;
		  back.second = span2;
		}
	      }
	    }
	  }
	}
	
	beam += std::log(threshold);
      } while (chart(0, T, 0, V) == infinity);
    }

    template <typename Costs, typename Iterator>
    void viterbi(const Costs& costs, Iterator result)
    {
      const int source_size = costs.size1() - 1;
      const int target_size = costs.size2() - 1;
      
      const backptr_type backptr_invalid(span_pair_type(span_type(-1, -1), span_type(-1, -1)),
					 span_pair_type(span_type(-1, -1), span_type(-1, -1)));
      
      
      stack_type stack;
      stack.push_back(span_pair_type(span_type(0, source_size), span_type(0, target_size)));
      
      while (! stack.empty()) {
	const span_pair_type span = stack.back();
	
	stack.pop_back();
	
	const backptr_type& back = backptr(span.source.first, span.source.last, span.target.first, span.target.last);
	
	if (back == backptr_invalid) {
	  // terminal!
	  if (! span.source.empty() && ! span.target.empty()) {
	    for (int src = span.source.first; src != span.source.last; ++ src)
	      for (int trg = span.target.first; trg != span.target.last; ++ trg) {
		*result = std::make_pair(src, trg);
		++ result;
	      }
	  }
	} else {
	  // non-terminal!
	  stack.push_back(back.second);
	  stack.push_back(back.first);
	}
      }
    }
    
    template <typename Costs, typename Iterator>
    void operator()(const Costs& costs, Iterator result)
    {
      initialize(costs);
      
      construct(costs, PruneNone());
      
      viterbi(costs, result);
    }
    
    template <typename Costs, typename Iterator>
    void operator()(const Costs& costs, const span_set_type& spans_source, const span_set_type& spans_target, Iterator result)
    {
      initialize(costs);
      
      construct(costs, PruneSpan(costs, spans_source, spans_target));
      
      viterbi(costs, result);
    }

    void shrink()
    {
      chart.clear();
      backptr.clear();
      
      chart_type(chart).swap(chart);
      backptr_chart_type(backptr).swap(backptr);
    }
    
    heap_type       heap;
    agenda_type     agenda;
    chart_mono_type chart_source;
    chart_mono_type chart_target;
    
    alpha_type alpha_source;
    alpha_type alpha_target;
    beta_type  beta_source;
    beta_type  beta_target;

    chart_type         chart;
    backptr_chart_type backptr;
  };
  
#endif
  
#if 0
  struct ITGAlignment
  {
    static const int max_fertility = 8;
    static const double threshold_global = 1e-5;
    static const double threshold_local  = 1e-2;
    
    typedef double logprob_type;
    
    typedef utils::chart<logprob_type, std::allocator<logprob_type> > chart_type;
    typedef utils::chart<chart_type, std::allocator<chart_type> > chart_set_type;

    typedef cicada::SpanVector       span_set_type;
    typedef span_set_type::span_type span_type;
    
    struct span_pair_type
    {
      span_type source;
      span_type target;
      
      span_pair_type() : source(), target() {}
      span_pair_type(const span_type& __source, const span_type& __target) : source(__source), target(__target) {}

      friend
      bool operator==(const span_pair_type& x, const span_pair_type& y)
      {
	return x.source == y.source && x.target == y.target;
      }
    };
    
    typedef std::pair<span_pair_type, span_pair_type> backptr_type;
    typedef utils::chart<backptr_type, std::allocator<backptr_type> > backptr_set_type;
    typedef utils::chart<backptr_set_type, std::allocator<backptr_set_type> > backptr_chart_type;

    typedef std::vector<span_pair_type, std::allocator<span_pair_type> > stack_type;
    
    ITGAlignment() {}
    
    struct SpanBeam
    {
      typedef utils::chart<bool, std::allocator<bool> > span_prune_type;
      typedef utils::chart<span_prune_type, std::allocator<span_prune_type> > span_prune_chart_type;

      typedef utils::vector2<double, std::allocator<double> > prob_set_type;
      typedef utils::vector2<double, std::allocator<double> > prefix_set_type;
      typedef utils::vector2<double, std::allocator<double> > suffix_set_type;
      typedef utils::vector2<double, std::allocator<double> > delta_set_type;
      typedef std::vector<double, std::allocator<double> >    inside_prune_type;
      
      template <typename Costs>
      SpanBeam(const Costs& costs)
      {
	prune_span(costs, prune_straight, false);
	prune_span(costs, prune_inverted, true);
      }
      
      bool operator()(const int source_first, const int source_last) const
      {
	return false;
      }
      
      bool operator()(const int source_first, const int source_last,
		      const int target_first, const int target_last) const
      {
	return (prune_straight(source_first, source_last)(target_first, target_last)
		&& prune_inverted(target_first, target_last)(source_first, source_last));
      }

      static const int L = 0;
      static const int M = 1;
      static const int R = 2;
      
      template <typename Costs, typename Spans>
      void prune_span(const Costs& costs,
		      Spans& span_prune,
		      const bool inverse)
      {
	const int source_size = (! inverse ? costs.size1() : costs.size2()) - 1;
	const int target_size = (! inverse ? costs.size2() : costs.size1()) - 1;
	
	span_prune.clear();
	span_prune.reserve(source_size + 1);
	span_prune.resize(source_size + 1, typename Spans::value_type(target_size + 1, true));
	
	probs.clear();
	prefix.clear();
	suffix.clear();
	delta.clear();
	inside.clear();

	probs.reserve(source_size + 2, target_size + 2);
	prefix.reserve(source_size + 2, target_size + 2);
	suffix.reserve(source_size + 2, target_size + 2);
      
	probs.resize(source_size + 2, target_size + 2, 0.0);
	prefix.resize(source_size + 2, target_size + 2, 0.0);
	suffix.resize(source_size + 2, target_size + 2, 0.0);
	
	delta.reserve(target_size + 2, 3);
	inside.reserve(target_size + 2);
	
	delta.resize(target_size + 2, 3, 0.0);
	inside.resize(target_size + 2, 1.0);

	if (inverse) {
	  for (int src = 0; src <= source_size; ++ src)
	    for (int trg = 1; trg <= target_size; ++ trg)
	      probs(src, trg) = utils::mathop::exp(costs(trg, src));
	} else {
	  for (int src = 0; src <= source_size; ++ src)
	    for (int trg = 1; trg <= target_size; ++ trg)
	      probs(src, trg) = utils::mathop::exp(costs(src, trg));
	}
	
	double merit_total = 1.0;
	
	for (int trg = 1; trg <= target_size; ++ trg) {
	  prefix(0, trg) = probs(0, trg);
	  for (int src = 1; src <= source_size; ++ src)
	    prefix(src, trg) = prefix(src - 1, trg) + probs(src, trg);
	  
	  merit_total *= prefix(source_size, trg);
	  
	  for (int src = source_size; src >= 0; -- src)
	    suffix(src, trg) = suffix(src + 1, trg) + probs(src, trg);
	}
            
	for (int i = 1; i <= source_size; ++ i) {
	
	  // fill-in NULL probabilities...
	  for (int a = 1; a <= target_size; ++ a)
	    inside[a] = probs(0, a);
	  
	  for (int j = i; j <= source_size; ++ j) {
	    // [i, j]
	  
	    std::fill(delta.begin(), delta.end(), 0.0);
	    delta(0, 0) = 1.0;
	  
	    prefix(i - 1, target_size + 1) = 0.0;
	    suffix(j + 1, target_size + 1) = 1.0;
	  
	    for (int a = 1; a <= target_size + 1; ++ a) {
	      inside[a] += probs(j, a);
	    
	      const double outside = prefix(i - 1, a) + suffix(j + 1, a);
	    
	      delta(a, L) = delta(a - 1, L) * outside;
	      delta(a, M) = std::max(delta(a - 1, L), delta(a - 1, M)) * inside[a];
	      delta(a, R) = std::max(delta(a - 1, M), delta(a - 1, R)) * outside;
	    }
	  
	    const double merit_span = delta(target_size + 1, R);
	  
	    // the first threshold...
	    if (merit_span / merit_total < threshold_global)
	      continue;
	    
	    // the source span [i, j] survived...
	    // next, step: search for actual span...

	    double outside_accumulated = 1.0;
	    for (int right = target_size; right >= 1; -- right) {
	      const double outside = (prefix(i - 1, right + 1) + suffix(j + 1, right + 1));
	      outside_accumulated *= outside;
	    
	      //const double merit_r = delta(right, R) * outside_accumulated;
	      //const double merit_m = delta(right, M) * outside_accumulated;
	    
	      //if (merit_r / merit_span < threshold_local) break;
	      //if (merit_m / merit_span < threshold_local) continue;
	    
	      // enumerate left boundary...
	      double inside_accumulated = 1.0;
	      for (int left = right; left >= 1; -- left) {
		const double merit_m = delta(left, M) * inside_accumulated * outside_accumulated;
	      
		if (merit_m / merit_span < threshold_local) break;
	      
		inside_accumulated *= inside[left];
		
		const double merit = delta(left - 1, L) * inside_accumulated * outside_accumulated;
		
		// survived! source:[i, j] and target:[left, right]
		if (merit / merit_span >= threshold_local)
		  span_prune(i - 1, j)(left - 1, right) = false;
	      }
	    }
	  }
	}
	
	span_prune(0, source_size)(0, target_size) = false;
      }
      
      span_prune_chart_type prune_straight;
      span_prune_chart_type prune_inverted;

      prob_set_type     probs;
      prefix_set_type   prefix;
      suffix_set_type   suffix;
      delta_set_type    delta;
      inside_prune_type inside;
    };

    struct SpanPrune
    {
      typedef utils::chart<bool, std::allocator<bool> > span_prune_type;

      template <typename Costs>
      SpanPrune(const Costs& costs, const span_set_type& spans_source, const span_set_type& spans_target)
      {
	const int source_size = costs.size1() - 1;
	const int target_size = costs.size2() - 1;
	
	prune_source.clear();
	prune_target.clear();
	
	prune_source.reserve(source_size + 1);
	prune_target.reserve(target_size + 1);
	
	prune_source.resize(source_size + 1, ! spans_source.empty());
	prune_target.resize(target_size + 1, ! spans_target.empty());
	
	if (! spans_source.empty()) {
	  span_set_type::const_iterator siter_end = spans_source.end();
	  for (span_set_type::const_iterator siter = spans_source.begin(); siter != siter_end; ++ siter)
	    prune_source(siter->first, siter->last) = false;
	}
	
	if (! spans_target.empty()) {
	  span_set_type::const_iterator siter_end = spans_target.end();
	  for (span_set_type::const_iterator siter = spans_target.begin(); siter != siter_end; ++ siter)
	    prune_target(siter->first, siter->last) = false;
	}
      }
      
      bool operator()(const int source_first, const int source_last) const
      {
	return prune_source(source_first, source_last);
      }
      
      bool operator()(const int source_first, const int source_last,
		      const int target_first, const int target_last) const
      {
	return prune_source(source_first, source_last) || prune_target(target_first, target_last);
      }
      
      span_prune_type prune_source;
      span_prune_type prune_target;
    };
    
    template <typename Costs>
    void initialize(const Costs& costs)
    {
      const int source_size = costs.size1() - 1;
      const int target_size = costs.size2() - 1;
      
      const int sent_fert = utils::bithack::max(source_size, target_size) / utils::bithack::min(source_size, target_size);
      const int max_fert  = utils::bithack::max(max_fertility, sent_fert);
      
      const double lowest = boost::numeric::bounds<double>::lowest();
      
      const backptr_type backptr_invalid(span_pair_type(span_type(-1, -1), span_type(-1, -1)),
					 span_pair_type(span_type(-1, -1), span_type(-1, -1)));
      
      inside.clear();
      inside.reserve(source_size + 1);
      inside.resize(source_size + 1, chart_type(target_size + 1, lowest));
      
      backptr.clear();
      backptr.reserve(source_size + 1);
      backptr.resize(source_size + 1, backptr_set_type(target_size + 1, backptr_invalid));
      
      for (int src = 0; src <= source_size; ++ src)
	for (int trg = 0; trg <= target_size; ++ trg) {
	  
	  if (src < source_size && trg < target_size) {
	    inside(src, src + 1)(trg, trg + 1) = costs(src + 1, trg + 1);
	    
	    // one-to-many
	    for (int pos = src + 1; pos < utils::bithack::min(source_size, src + max_fert); ++ pos)
	      inside(src, pos + 1)(trg, trg + 1) = inside(src, pos)(trg, trg + 1) + costs(pos + 1, trg + 1);
	    
	    for (int pos = trg + 1; pos < utils::bithack::min(target_size, trg + max_fert); ++ pos)
	      inside(src, src + 1)(trg, pos + 1) = inside(src, src + 1)(trg, pos) + costs(src + 1, pos + 1);
	  }
	  
	  // one-to-null
	  if (src < source_size)
	    inside(src, src + 1)(trg, trg) = costs(src + 1, 0);
	  
	  if (trg < target_size)
	    inside(src, src)(trg, trg + 1) = costs(0, trg + 1);
	}
    }
    
    template <typename Costs, typename Prune>
    void construct(const Costs& costs, const Prune& prune)
    {
      const int source_size = costs.size1() - 1;
      const int target_size = costs.size2() - 1;

      const double lowest = boost::numeric::bounds<double>::lowest();
      
      // compute inside probabilities...
      for (int source_length = 1; source_length <= source_size; ++ source_length)
	for (int source_first = 0; source_first + source_length <= source_size; ++ source_first) {
	  const int source_last = source_first + source_length;
	  
	  if (source_length > 1 && prune(source_first, source_last)) continue;
	  
	  chart_type& inside_source = inside(source_first, source_last);
	  
	  for (int target_length = (source_length == 1 ? 2 : 1); target_length <= target_size; ++ target_length)
	    for (int target_first = 0; target_first + target_length <= target_size; ++ target_first) {
	      const int target_last = target_first + target_length;
	      
	      if (target_length > 1 && prune(source_first, source_last, target_first, target_last)) continue;
	      
	      logprob_type& inside_source_target = inside_source(target_first, target_last);
	      backptr_type& back = backptr(source_first, source_last)(target_first, target_last);
	      
	      for (int src = source_first; src <= source_last; ++ src) {
		if (src - source_first > 1 && prune(source_first, src)) continue;
		if (source_last - src > 1 && prune(src, source_last)) continue;
		
		const chart_type& source1 = inside(source_first, src);
		const chart_type& source2 = inside(src, source_last);
		
		for (int trg = target_first; trg <= target_last; ++ trg)
		  if ((src - source_first) * (source_last - src) + (trg - target_first) * (target_last - trg) != 0) {
		    
		    {// straight...
		      const logprob_type& target1 = source1(target_first, trg);
		      const logprob_type& target2 = source2(trg, target_last);
		      
		      if (target1 > lowest && target2 > lowest && target1 + target2 > inside_source_target) {
			inside_source_target = target1 + target2;
			
			back.first  = span_pair_type(span_type(source_first, src), span_type(target_first, trg));
			back.second = span_pair_type(span_type(src, source_last), span_type(trg, target_last));
		      }
		    }
		    
		    {// inverted...
		      const logprob_type& target1 = source1(trg, target_last);
		      const logprob_type& target2 = source2(target_first, trg);
		      
		      if (target1 > lowest && target2 > lowest && target1 + target2 > inside_source_target) {
			inside_source_target = target1 + target2;
			
			back.first  = span_pair_type(span_type(src, source_last), span_type(target_first, trg));
			back.second = span_pair_type(span_type(source_first, src), span_type(trg, target_last));
		      }
		    }
		  }
	      }
	    }
	}
    }

    template <typename Costs, typename Iterator>
    void viterbi(const Costs& costs, Iterator result)
    {
      const int source_size = costs.size1() - 1;
      const int target_size = costs.size2() - 1;
      
      const backptr_type backptr_invalid(span_pair_type(span_type(-1, -1), span_type(-1, -1)),
					 span_pair_type(span_type(-1, -1), span_type(-1, -1)));
      
      
      stack_type stack;
      stack.push_back(span_pair_type(span_type(0, source_size), span_type(0, target_size)));
      
      while (! stack.empty()) {
	const span_pair_type span = stack.back();
	
	stack.pop_back();
	
	const backptr_type& back = backptr(span.source.first, span.source.last)(span.target.first, span.target.last);
	
	if (back == backptr_invalid) {
	  // terminal!
	  if (! span.source.empty() && ! span.target.empty()) {
	    for (int src = span.source.first; src != span.source.last; ++ src)
	      for (int trg = span.target.first; trg != span.target.last; ++ trg) {
		*result = std::make_pair(src, trg);
		++ result;
	      }
	  }
	} else {
	  // non-terminal!
	  stack.push_back(back.second);
	  stack.push_back(back.first);
	}
      }
    }
    
    
    template <typename Costs, typename Iterator>
    void operator()(const Costs& costs, Iterator result)
    {
      initialize(costs);
      
      construct(costs, SpanBeam(costs));
      
      viterbi(costs, result);
    }
    
    template <typename Costs, typename Iterator>
    void operator()(const Costs& costs, const span_set_type& spans_source, const span_set_type& spans_target, Iterator result)
    {
      initialize(costs);
      
      construct(costs, SpanPrune(costs, spans_source, spans_target));
      
      viterbi(costs, result);
    }

    void shrink()
    {
      inside.clear();
      backptr.clear();
      
      chart_set_type(inside).swap(inside);
      backptr_chart_type(backptr).swap(backptr);
    }
    
    chart_set_type     inside;
    backptr_chart_type backptr;
  };

#endif  
};

#endif
