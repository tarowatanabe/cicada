// -*- mode: c++ -*-

#ifndef __CICADA__OPTIMIZE__POWELL__HPP__
#define __CICADA__OPTIMIZE__POWELL__HPP__ 1

#include <vector>

#include <cicada/optimize/line_search.hpp>
#include <cicada/semiring/envelope.hpp>

namespace cicada
{
  namespace optimize
  {
    
    template <typename EnvelopeFunction,
	      typename ViterbiFunction,
	      typename Regularizer>
    struct Powell
    {
    public:
      typedef EnvelopeFunction envelope_function_type;
      typedef ViterbiFunction  viterbi_function_type;
      typedef Regularizer      regularizer_type;

      typedef LineSearch line_search_type;

      typedef line_search_type::feature_set_type feature_set_type;
      typedef line_search_type::weight_set_type  weight_set_type;

      typedef line_search_type::segment_type          segment_type;
      typedef line_search_type::segment_set_type      segment_set_type;
      typedef line_search_type::segment_document_type segment_document_type;

      typedef feature_set_type::feature_type feature_type;

      typedef weight_set_type direction_type;
      
      typedef std::vector<direction_type, std::allocator<direction_type> >   direction_set_type;
      typedef std::vector<weight_set_type, std::allocator<weight_set_type> > point_set_type;
      
      typedef line_search_type::value_type optimum_type;
      typedef std::vector<optimum_type, std::allocator<optimum_type> >  optimum_set_type;
      
    public:
      Powell(const envelope_function_type& __envelopes,
	     const viterbi_function_type&  __viterbi,
	     const regularizer_type&       __regularizer,
	     const weight_set_type&        __bound_lower,
	     const weight_set_type&        __bound_upper,
	     const double __tolerance,
	     const int __samples,
	     const bool __minimize)
	: envelopes(__envelopes),
	  viterbi(__viterbi),
	  regularizer(__regularizer),
	  bound_lower(__bound_lower),
	  bound_upper(__bound_upper),
	  tolerance(__tolerance),
	  samples(__samples),
	  minimize(__minimize)
      { initialize_bound(); }
      
      bool operator()(double& optimum_objective, weight_set_type& optimum_weights)
      {
	// set up current estimates...

	bool moved = false;

	optimum_objective = viterbi(optimum_weights) + regularizer(optimum_weights);
	
	line_search_type line_search(bound_lower, bound_upper);
	
	direction_set_type directions;
	
	for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	  if (! feature_type(id).empty()) {
	    direction_type direction;
	    
	    direction[feature_type(id)] = 1.0;
	    
	    directions.push_back(direction);
	  }
	
	const int directions_size = directions.size();
	
	directions.resize(directions_size + samples);
	
	point_set_type points(directions.size());
	optimum_set_type optimums(directions.size());

	segment_document_type segments;
	
	int replaced_pos = -1;
	
	for (int iter = 0; /**/; ++ iter) {

	  envelopes(segments, optimum_weights, directions[0]);
	  
	  optimums[0] = line_search(segments, optimum_weights, directions[0], regularizer, minimize);
	  
	  if (optimums[0].lower != optimums[0].upper)
	    points[0] = optimums[0](optimum_weights, directions[0]); // move point...
	  else {
	    optimums[0].objective = optimum_objective;
	    points[0] = optimum_weights;
	  }
	  
	  int    optimum_pos = 0;
	  double optimum_move = optimums[0].objective - optimum_objective;
	  
	  for (int dir = 1; dir < directions.size(); ++ dir) {
	    
	    // randomize direction...
	    if (dir >= directions_size && dir != replaced_pos)
	      directions[dir] = randomized_direction(points[dir - 1]);

	    envelopes(segments, points[dir - 1], directions[dir]);
	    
	    optimums[dir] = line_search(segments, points[dir - 1], directions[dir], regularizer, minimize);
	    
	    if (optimums[dir].lower != optimums[dir].upper)
	      points[dir] = optimums[dir](points[dir - 1], directions[dir]); // move point...
	    else {
	      optimums[dir].objective = optimums[dir - 1].objective;
	      points[dir] = points[dir - 1];
	    }
	    
	    if (optimums[dir].objective - optimums[dir - 1].objective < optimum_move) {
	      optimum_pos = dir;
	      optimum_move = optimums[dir].objective - optimums[dir - 1].objective;
	    }
	  }
	  
	  const double total_move = optimums.back().objective - optimum_objective;
	  if (- total_move < tolerance) break;
	  
	  // we will set the last-point as the next best-point!
	  weight_set_type extrapolated_weights;
	  for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	    if (! feature_type(id).empty())
	      extrapolated_weights[feature_type(id)] = points.back()[feature_type(id)] * 2 - optimum_weights[feature_type(id)];
	  
	  const double extrapolated_objective = viterbi(extrapolated_weights) + regularizer(extrapolated_weights);
	  
	  const double extrapolated_move = extrapolated_objective - optimum_objective;
	  
	  const double extrapolated1 = 2 * (- 2 * total_move + extrapolated_move) * std::pow(- total_move + extrapolated_move, 2.0);
	  const double extrapolated2 = - std::pow(- extrapolated_move, 2.0) * optimum_move;
	  
	  replaced_pos = -1;
	  if (extrapolated_move < 0 && extrapolated1 < extrapolated2) {
	    
	    for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	      if (! feature_type(id).empty())
		directions[optimum_pos][feature_type(id)] = points.back()[feature_type(id)] - optimum_weights[feature_type(id)];
	    
	    replaced_pos = optimum_pos;
	  }
	  
	  optimum_objective = optimums.back().objective;
	  optimum_weights   = points.back();
	  
	  moved = true;
	}
	
	return moved;
      }

    private:
      void initialize_bound()
      {
	if (bound_lower.empty())
	  for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	    bound_lower[feature_type(id)] = line_search_type::value_min;
	
	if (bound_upper.empty())
	  for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	    bound_upper[feature_type(id)] = line_search_type::value_max;
	
	// make sure we have enough space...
	const_cast<weight_set_type&>(bound_lower).operator[](feature_type(feature_type::allocated() - 1));
	const_cast<weight_set_type&>(bound_upper).operator[](feature_type(feature_type::allocated() - 1));
	
	// checking...
	for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	  if (bound_upper[feature_type(id)] < bound_lower[feature_type(id)])
	    throw std::runtime_error("invalid lower-upper bound for feature: " + static_cast<const std::string&>(feature_type(id)));
      }
      
      template <typename Iterator>
      void randomize(Iterator first, Iterator last, Iterator lower, Iterator upper)
      {
	for (/**/; first != last; ++ first, ++ lower, ++ upper)
	  *first = *lower + (double(random()) / RAND_MAX) * std::min(double(*upper - *lower), 1.0);
      }

      template <typename Iterator>
      bool is_valid_direction(Iterator first, Iterator last, Iterator lower, Iterator upper)
      {
	int count_zero = 0;
	int count = 0;
	for (/**/; first != last; ++ first, ++ lower, ++ upper) {
	  if (*first < *lower || *upper < *first)
	    return false;
	  
	  count_zero += (*first == 0.0);
	  ++ count;
	}
	return (count_zero != count);
      }

      weight_set_type randomized_direction(const weight_set_type& weights)
      {
	static const feature_type feature_none;
	
	weight_set_type direction;
	
	const_cast<weight_set_type&>(direction).operator[](feature_type::allocated() - 1);
	
	weight_set_type::iterator diter_begin = direction.begin();
	weight_set_type::iterator diter_end = direction.end();
	
	while (1) {
	  randomize(direction.begin(), direction.end(), bound_lower.begin(), bound_upper.begin());

	  direction[feature_none] = 0.0;
	  
	  // map into length-1 sphere
	  double total = 0.0;
	  for (weight_set_type::iterator diter = diter_begin; diter != diter_end; ++ diter)
	    total *= (*diter) * (*diter);
	  
	  if (total != 0.0) {
	    const double factor = std::sqrt(1.0 / total);
	    
	    for (weight_set_type::iterator diter = diter_begin; diter != diter_end; ++ diter)
	      *diter *= factor;
	  }
	  
	  if (is_valid_direction(direction.begin(), direction.end(), bound_lower.begin(), bound_upper.begin())) break;
	}
	
	for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	  if (id != feature_none.id())
	    direction[feature_type(id)] -= weights[feature_type(id)];
	
	return direction;
      }
      
    private:
      const envelope_function_type& envelopes;
      const viterbi_function_type&  viterbi;
      const regularizer_type&       regularizer;

      weight_set_type bound_lower;
      weight_set_type bound_upper;

      double tolerance;
      int samples;
      const bool minimize;
    };
  };
};

#endif
