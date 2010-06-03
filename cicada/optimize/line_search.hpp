// -*- mode: c++ -*-

#ifndef __CICADA__OPTIMIZE__LINE_SEARCH__HPP__
#define __CICADA__OPTIMIZE__LINE_SEARCH__HPP__ 1

#include <vector>
#include <deque>
#include <algorithm>

#include <cicada/hypergraph.hpp>
#include <cicada/weight_vector.hpp>

#include <cicada/semiring/envelope.hpp>
#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace optimize
  {
    struct LineSearch
    {
      typedef HyperGraph hypergraph_type;
      
      typedef hypergraph_type::feature_set_type feature_set_type;
      typedef cicada::WeightVector<double>      weight_set_type;
      typedef feature_set_type::feature_type    feature_type;

      typedef semiring::Envelope       envelope_type;
      
      typedef envelope_type::line_type line_type;
      typedef eval::Score              score_type;

      typedef boost::shared_ptr<line_type>  line_ptr_type;
      typedef boost::shared_ptr<score_type> score_ptr_type;

      typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> >   score_set_type;
      
      typedef std::pair<line_ptr_type, score_ptr_type>                        segment_type;
      typedef std::vector<segment_type, std::allocator<segment_type> >        segment_set_type;
      typedef std::deque<segment_set_type, std::allocator<segment_set_type> > segment_document_type;
      
      static const double interval_min = 1e-4;
      static const double interval_offset_lower = 200.0;
      static const double interval_offset_upper = 0.2;
      static const double value_min = -100.0;
      static const double value_max =  100.0;
      
    public:
      struct Result
      {
	double error;
	double score;
	double lower;
	double upper;
	
	weight_set_type operator()(const weight_set_type& origin,
				   const weight_set_type& direction)
	{
	  weight_set_type weights;
	  
	  for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	    if (! feature_type(id).empty()) {
	      const feature_type feature(id);

	      const double weight = origin[feature];
	      const double dir    = direction[feature];
	      
	      const double feature_lower = weight + lower * dir;
	      const double feature_upper = weight + upper * dir;
	      
	      const double average = (lower + upper) * 0.5;
	      
	      // if we cross minus and plus, force zero...
	      
	      weights[feature] = (weight + average * dir) * double(! (feature_lower * feature_upper < 0.0));
	    }
	  return weights;
	}
	
	
	Result() : error(std::numeric_limits<double>::infinity()), score(), lower(), upper() {}
	
	Result(const double& _error,
	       const double& _score,
	       const double& _lower,
	       const double& _upper)
	  : error(_error), score(_score), lower(_lower), upper(_upper) {}
      };
      typedef Result value_type;

    private:
      struct Item
      {
	int seg;
	
	segment_set_type::const_iterator first;
	segment_set_type::const_iterator last;

	Item() : seg(), first(), last() {}
	Item(const int& __seg,
	     segment_set_type::const_iterator __first,
	     segment_set_type::const_iterator __last)
	  : seg(__seg), first(__first), last(__last) {}
      };
      typedef Item item_type;
      
      struct item_heap_compare_type
      {
	bool operator()(const item_type& x, const item_type& y) const
	{
	  return x.first->first->x > y.first->first->x;
	}
      };
      
      typedef std::vector<item_type, std::allocator<item_type> > heap_type;

      
    public:
      LineSearch(const weight_set_type& __origin,
		 const weight_set_type& __direction,
		 const int __debug = 0)
	: origin(__origin),
	  direction(__direction),
	  bound_lower(),
 	  bound_upper(),
	  debug(__debug) { initialize_bound(); }
      
      LineSearch(const weight_set_type& __origin,
		 const weight_set_type& __direction,
		 const weight_set_type& __bound_lower,
		 const weight_set_type& __bound_upper,
		 const int __debug = 0)
	: origin(__origin),
	  direction(__direction),
	  bound_lower(__bound_lower),
	  bound_upper(__bound_upper),
	  debug(__debug) { initialize_bound(); }

      template <typename Regularizer>
      value_type operator()(segment_document_type& segments,
			    Regularizer regularizer,
			    const bool minimize)
      {
	// we assume a set of line_ptr and score_ptr pair...
	
	const double error_factor = (minimize ? 1.0 : - 1.0);
	const std::pair<double, double> range = valid_range();

	heap_type heap;
	
	if (debug >= 4)
	  std::cerr << "minimum: " << range.first << " maximum: " << range.second << std::endl;
	
	score_ptr_type score;
	score_set_type scores(segments.size());
	
	for (int seg = 0; seg < segments.size(); ++ seg)
	  if (! segments[seg].empty()) {
	    
	    if (! score)
	      score.reset(segments[seg].front().second->zero());
	    
	    scores[seg] = segments[seg].front().second;
	    *score += *segments[seg].front().second;
	    
	    if (segments[seg].size() > 1)
	      heap.push_back(item_type(seg, segments[seg].begin() + 1, segments[seg].end()));
	  }

	if (heap.empty())
	  return value_type();
	
	// priority queue...
	std::make_heap(heap.begin(), heap.end(), item_heap_compare_type());
	
	const double error = error_factor * score->score().first;
	
	double best_lower = lower_bound(heap.front().first->first->x, range.first);
	double best_upper = heap.front().first->first->x;
	
	double best_error = error + regularizer(origin, direction, best_lower, best_upper);
	double best_score = error;
	
	double segment_prev = best_lower;
	double score_prev   = error;

	if (debug >= 4)
	  std::cerr << "lower: " << best_lower
		    << " upper: " << best_upper
		    << " error: " << error
		    << " regularized: " << best_error
		    << std::endl;

	
	while (! heap.empty()) {
	  // next at heap...
	  
	  const double segment_curr = heap.front().first->first->x;
	  
	  while (! heap.empty() && heap.front().first->first->x == segment_curr) { 
	    
	    std::pop_heap(heap.begin(), heap.end(), item_heap_compare_type());
	    
	    *score -= *scores[heap.back().seg];
	    *score += *heap.back().first->second;
	    scores[heap.back().seg] = heap.back().first->second;
	    
	    // pop and push heap...
	    ++ heap.back().first;
	    if (heap.back().first != heap.back().last)
	      std::push_heap(heap.begin(), heap.end(), item_heap_compare_type());
	    else
	      heap.pop_back();
	  }
	
	  const double segment_next = (heap.empty() ? upper_bound(segment_curr, range.second) : heap.front().first->first->x);
	  
	  // we perform merging of ranges if error counts are equal...
	  const double error = error_factor * score->score().first;
	  if (error != score_prev) {
	    segment_prev = segment_curr;
	    score_prev = error;
	  }
	  
	  const double lower = segment_prev;
	  const double upper = segment_next;
	  const double point = (lower + upper) * 0.5;
	  
	  if (point > range.second) break;   // out of range for upper-bound, quit!
	  if (point < range.first) continue; // out of range for lower-bound...
	  if (std::fabs(point) < interval_min) continue; // interval is very small
	  
	  const double error_regularized = error + regularizer(origin, direction, lower, upper);
	  
	  if (debug >= 4)
	    std::cerr << "lower: " << lower
		      << " upper: " << upper
		      << " error: " << error
		      << " regularized: " << error_regularized
		      << std::endl;
	  
	  if (error_regularized < best_error) {
	    best_error = error_regularized;
	    best_score = error;
	    best_lower = lower;
	    best_upper = upper;
	  }
	}
	
	const double point = (best_lower + best_upper) * 0.5;
	if (point < range.first || range.second < point)
	  return value_type();
	else {
	  if (debug >= 2)
	    std::cerr << "minimum error: " << best_error
		      << " score: " << best_score
		      << " lower: " << best_lower
		      << " upper: " << best_upper << std::endl;
	  return value_type(best_error, best_score, best_lower, best_upper);
	}
      }

    private:
      
      std::pair<double, double> valid_range() const
      {
	double minimum = - std::numeric_limits<double>::infinity();
	double maximum =   std::numeric_limits<double>::infinity();
	
	for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	  if (! feature_type(id).empty()) {
	    const feature_type feature(id);
	    
	    const double ori = origin[feature];
	    const double dir = direction[feature];
	    const double low = bound_lower[feature];
	    const double upp = bound_upper[feature];

	    if (dir > 0.0) {
	      maximum = std::min(maximum, (upp - ori) / dir);
	      minimum = std::max(minimum, (low - ori) / dir);
	    } else {
	      maximum = std::min(maximum, (low  - ori) / dir);
	      minimum = std::max(minimum, (upp - ori) / dir);
	    }
	  }
	
	return std::make_pair(minimum, maximum);
      }
      
      double lower_bound(const double& point, const double& bound)
      {
	return (point <= bound ? point - interval_offset_lower : std::max(bound, point - interval_offset_lower));
      }
      
      double upper_bound(const double& point, const double& bound)
      {
	return (point >= bound ? point + interval_offset_upper : std::min(bound, point + interval_offset_upper));
      }
      

      void initialize_bound()
      {
	if (bound_lower.empty())
	  for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	    bound_lower[feature_type(id)] = value_min;
	
	if (bound_upper.empty())
	  for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	    bound_upper[feature_type(id)] = value_max;
	
	// make sure we have enough space...
	bound_lower[feature_type(feature_type::allocated() - 1)];
	bound_upper[feature_type(feature_type::allocated() - 1)];
	
	// checking...
	for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	  if (bound_upper[feature_type(id)] < bound_lower[feature_type(id)])
	    throw std::runtime_error("invalid lower-upper bound for feature: " + static_cast<const std::string&>(feature_type(id)));
      }
      
    private:
      const weight_set_type&      origin;
      const weight_set_type&      direction;

      weight_set_type bound_lower;
      weight_set_type bound_upper;

      const int debug;
    };
  };
};

#endif
