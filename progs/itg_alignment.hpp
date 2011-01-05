// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __ITG_ALIGNMENT__HPP__
#define __ITG_ALIGNMENT__HPP__ 1

#include <vector>

#include <cicada/span_vector.hpp>

#include <utils/chart.hpp>
#include <utils/vector2.hpp>
#include <utils/mathop.hpp>
#include <utils/bithack.hpp>

namespace detail
{
  
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
	const int source_size = (! inverse ? costs.size1() : costs.size2());
	const int target_size = (! inverse ? costs.size2() : costs.size1());
	
	span_prune.clear();
	span_prune.resize(source_size + 1, typename Spans::value_type(target_size + 1, true));
	
	probs.clear();
	prefix.clear();
	suffix.clear();
	delta.clear();
	inside.clear();
      
	probs.resize(source_size + 2, target_size + 2, 0.0);
	prefix.resize(source_size + 2, target_size + 2, 0.0);
	suffix.resize(source_size + 2, target_size + 2, 0.0);
      
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
	    
	      const double merit_r = delta(right, R) * outside_accumulated;
	      const double merit_m = delta(right, M) * outside_accumulated;
	    
	    
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
	const int source_size = costs.size1();
	const int target_size = costs.size2();
	
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
      const int source_size = costs.size1();
      const int target_size = costs.size2();
      
      const int sent_fert = utils::bithack::max(source_size, target_size) / utils::bithack::min(source_size, target_size);
      const int max_fert  = utils::bithack::max(max_fertility, sent_fert);
      
      const double lowest = boost::numeric::bounds<double>::lowest();
      
      const backptr_type backptr_invalid(span_pair_type(span_type(-1, -1), span_type(-1, -1)),
					 span_pair_type(span_type(-1, -1), span_type(-1, -1)));
      
      inside.clear();
      inside.resize(source_size + 1, chart_type(target_size + 1, lowest));
      
      backptr.clear();
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
      const int source_size = costs.size1();
      const int target_size = costs.size2();

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
      const int source_size = costs.size1();
      const int target_size = costs.size2();
      
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
      const double lowest = boost::numeric::bounds<double>::lowest();
      
      initialize(costs);
      
      construct(costs, SpanBeam(costs));
      
      viterbi(costs, result);
    }
    
    template <typename Costs, typename Iterator>
    void operator()(const Costs& costs, const span_set_type& spans_source, const span_set_type& spans_target, Iterator result)
    {
      const double lowest = boost::numeric::bounds<double>::lowest();
      
      initialize(costs);
      
      construct(costs, SpanPrune(costs, spans_source, spans_target));
      
      viterbi(costs, result);
    }
    
    chart_set_type     inside;
    backptr_chart_type backptr;
  };
  
};

#endif
