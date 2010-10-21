
#include <cicada/optimize/line_search.hpp>

struct LineSearch
{
  typedef cicada::eval::Scorer         scorer_type;
  typedef cicada::eval::ScorerDocument scorer_document_type;
  
  typedef scorer_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_type::score_ptr_type  score_ptr_type;

  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;

  typedef line_search_type::value_type optimum_type;
  
  typedef cicada::semiring::Envelope envelope_type;
  typedef std::vector<envelope_type, std::allocator<envelope_type> >  envelope_set_type;

  typedef cicada::Sentence sentence_type;

  LineSearch(const int __debug=0) : line_search(__debug), debug(__debug) {}

  segment_document_type segments;
  envelope_set_type     envelopes;
  sentence_type         yield;
  
  line_search_type line_search;
  int debug;
  
  template <typename HypergraphSet, typename ScorerSet>
  double operator()(const HypergraphSet& graphs,
		    const ScorerSet& scorers,
		    const weight_set_type& origin,
		    const weight_set_type& direction,
		    const double lower,
		    const double upper,
		    const bool yield_source=false)
  {
    segments.clear();
    segments.resize(graphs.size());
    
    for (int seg = 0; seg < graphs.size(); ++ seg) 
      if (graphs[seg].is_valid()) {
	
	if (debug >= 4)
	  std::cerr << "line-search segment: " << seg << std::endl;
	
	envelopes.clear();
	envelopes.resize(graphs[seg].nodes.size());
	
	cicada::inside(graphs[seg], envelopes, cicada::semiring::EnvelopeFunction<weight_set_type>(origin, direction));
	
	const envelope_type& envelope = envelopes[graphs[seg].goal];
	const_cast<envelope_type&>(envelope).sort();
	
	envelope_type::const_iterator eiter_end = envelope.end();
	for (envelope_type::const_iterator eiter = envelope.begin(); eiter != eiter_end; ++ eiter) {
	  const envelope_type::line_ptr_type& line = *eiter;
	  
	  line->yield(yield, yield_source);
	  
	  scorer_type::score_ptr_type score = scorers[seg]->score(yield);
	  
	  segments[seg].push_back(std::make_pair(line->x, score));
	}
      }
    
    // traverse segments...
    const optimum_type optimum = line_search(segments, lower, upper, false);
    return (optimum.lower + optimum.upper) * 0.5;
  }
  
};

struct OptimizerBase
{
  OptimizerBase(const hypergraph_set_type&           __graphs,
		const feature_function_ptr_set_type& __features)
    : graphs(__graphs),
      features(__features),
      weight_scale(1.0),
      weight_norm(0.0)
  {
    // initialize weights and weights bleu...
    for (int i = 0; i < features.size(); ++ i)
      if (features[i]) {
	feature_bleu = features[i]->feature_name();
	break;
      }
    
    weights_bleu[feature_bleu] = loss_scale;
  }
  
  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
  typedef cicada::WeightVector<weight_type, std::allocator<weight_type> >  expectation_type;
  
  typedef std::vector<weight_type, std::allocator<weight_type> > inside_set_type;
  
  typedef cicada::semiring::Tropical<double> bleu_weight_type;
  typedef std::vector<bleu_weight_type, std::allocator<bleu_weight_type> > bleu_weight_set_type;
  
  struct gradient_set_type
  {
    typedef gradient_type value_type;
      
    template <typename Index>
    gradient_type& operator[](Index)
    {
      return gradient;
    }

    void clear() { gradient.clear(); }
      
    gradient_type gradient;
  };

  struct bleu_weight_function
  {
    typedef cicada::semiring::Tropical<double> value_type;
      
    bleu_weight_function(const weight_set_type& __weights) : weights(__weights) {}

    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
    }
      
    const weight_set_type& weights;
  };

  struct weight_set_function
  {
    typedef cicada::semiring::Logprob<double> value_type;

    weight_set_function(const weight_set_type& __weights)
      : weights(__weights) {}

    const weight_set_type& weights;
      
    value_type operator()(const feature_set_type& x) const
    {
      return cicada::semiring::traits<value_type>::log(x.dot(weights));
    }
  };

  struct weight_function
  {
    typedef weight_type value_type;

    weight_function(const weight_set_type& __weights)
      : weights(__weights) {}

    const weight_set_type& weights;

    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e
      return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
    }
  };

  struct weight_max_function
  {
    typedef weight_type value_type;

    weight_max_function(const weight_set_type& __weights, const bleu_weight_set_type& __bleus, const bleu_weight_type& __max_bleu)
      : weights(__weights), bleus(__bleus), max_bleu(__max_bleu) {}

    const weight_set_type&      weights;
    const bleu_weight_set_type& bleus;
    const bleu_weight_type      max_bleu;

    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e
      if (log(max_bleu) - log(bleus[edge.id]) >= 1e-4)
	return value_type();
      else
	return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
    }
  };

  struct feature_function
  {
    typedef gradient_type value_type;

    feature_function(const weight_set_type& __weights) : weights(__weights) {}

    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e r_e
      gradient_type grad;
	
      const weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights));
	
      feature_set_type::const_iterator fiter_end = edge.features.end();
      for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	if (fiter->second != 0.0)
	  grad[fiter->first] = weight_type(fiter->second) * weight;
	
      return grad;
    }
      
    const weight_set_type& weights;
  };
  

  struct feature_max_function
  {
    typedef gradient_type value_type;

    feature_max_function(const weight_set_type& __weights, const bleu_weight_set_type& __bleus, const bleu_weight_type& __max_bleu)
      : weights(__weights), bleus(__bleus), max_bleu(__max_bleu) {}

    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e r_e
      if (log(max_bleu) - log(bleus[edge.id]) >= 1e-4)
	return gradient_type();

      gradient_type grad;
	
      const weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights));
	
      feature_set_type::const_iterator fiter_end = edge.features.end();
      for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	if (fiter->second != 0.0)
	  grad[fiter->first] = weight_type(fiter->second) * weight;
	
      return grad;
    }
    
    
    const weight_set_type&      weights;
    const bleu_weight_set_type& bleus;
    const bleu_weight_type      max_bleu;
  };
  
  bool operator()(const int id)
  {
    if (! graphs[id].is_valid()) return false;
    
    model_type model;
    model.push_back(features[id]);
    
    if (! apply_exact) {
      weights[feature_bleu] = loss_scale / weight_scale;
      cicada::apply_cube_prune(model, graphs[id], graph_reward, weight_set_function(weights), cube_size);
      
      weights[feature_bleu] = - loss_scale / weight_scale;
      cicada::apply_cube_prune(model, graphs[id], graph_penalty, weight_set_function(weights), cube_size);
      
      weights[feature_bleu] = 0.0;
      
      graph_reward.unite(graph_penalty);
      graph_penalty.clear();
    }
    
    const hypergraph_type& graph = (apply_exact ? graphs[id] : graph_reward);
    
    // compute inside/outside by bleu using tropical semiring...
    bleus_inside.clear();
    bleus_inside_outside.clear();
    
    bleus_inside.resize(graph.nodes.size());
    bleus_inside_outside.resize(graph.edges.size());
    
    cicada::inside_outside(graph, bleus_inside, bleus_inside_outside, bleu_weight_function(weights_bleu), bleu_weight_function(weights_bleu));
    
    const bleu_weight_type bleu_max = *std::max_element(bleus_inside_outside.begin(), bleus_inside_outside.end());
    
    // then, inside/outside to collect potentials...
    
    inside.clear();
    inside_correct.clear();
    
    gradients.clear();
    gradients_correct.clear();
    
    inside.resize(graph.nodes.size());
    inside_correct.resize(graph.nodes.size());
    
    if (softmax_margin) {
      weights[feature_bleu] = - loss_scale / weight_scale;
      
      cicada::inside_outside(graph, inside, gradients, weight_function(weights), feature_function(weights));
      
      cicada::inside_outside(graph, inside_correct, gradients_correct,
			     weight_max_function(weights, bleus_inside_outside, bleu_max),
			     feature_max_function(weights, bleus_inside_outside, bleu_max));
      
      weights[feature_bleu] = 0.0;
    } else {
      cicada::inside_outside(graph, inside, gradients, weight_function(weights), feature_function(weights));
      
      cicada::inside_outside(graph, inside_correct, gradients_correct,
			     weight_max_function(weights, bleus_inside_outside, bleu_max),
			     feature_max_function(weights, bleus_inside_outside, bleu_max));
    }
    
    Z = inside.back();
    Z_correct = inside_correct.back();
    
    gradients.gradient /= Z;
    gradients_correct.gradient /= Z_correct;

    const double margin = log(Z_correct) - log(Z);

    if (debug >= 3)
      std::cerr << "id: " << id
		<< " bleu: " << log(bleu_max)
		<< " correct: " << log(Z_correct)
		<< " partition: " << log(Z)
		<< " margin: " << margin
		<< std::endl;

    return true;
  }
  
  const hypergraph_set_type&           graphs;
  const feature_function_ptr_set_type& features;
  
  weight_set_type weights;
  weight_set_type weights_bleu;
  weight_set_type::feature_type feature_bleu;

  double objective;
  
  double weight_scale;
  double weight_norm;
  
  hypergraph_type graph_reward;
  hypergraph_type graph_penalty;

  bleu_weight_set_type bleus_inside;
  bleu_weight_set_type bleus_inside_outside;

  inside_set_type   inside;
  gradient_set_type gradients;

  inside_set_type   inside_correct;
  gradient_set_type gradients_correct;
  
  weight_type Z;
  weight_type Z_correct;
};

struct OptimizerSGDL2 : public OptimizerBase
{
  OptimizerSGDL2(const hypergraph_set_type&           __graphs,
		 const feature_function_ptr_set_type& __features)
    : OptimizerBase(__graphs, __features),
      samples(0),
      epoch(0),
      lambda(C / __graphs.size()) {}
  
  void initialize()
  {
    samples = 0;
    
    weights[feature_bleu] = 0.0;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
    
    objective = 0.0;
  }
  
  void finalize()
  {
    weights[feature_bleu] = 0.0;
    
    weights *= weight_scale;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  
  void operator()(const int segment)
  {
    const double eta = 1.0 / (lambda * (epoch + 2));
    ++ epoch;
    
    // we will minimize...
    if (OptimizerBase::operator()(segment)) {
      // we have gradients_correct/gradients/Z_correct/Z
      
      rescale(1.0 - eta * lambda);
      
      // update wrt correct gradients
      gradient_type::const_iterator citer_end = gradients_correct.gradient.end();
      for (gradient_type::const_iterator citer = gradients_correct.gradient.begin(); citer != citer_end; ++ citer) 
	if (citer->first != feature_bleu) {
	  const double feature = citer->second;
	  update(weights[citer->first], feature * eta);
	}
      
      // update wrt marginal gradients...
      gradient_type::const_iterator miter_end = gradients.gradient.end();
      for (gradient_type::const_iterator miter = gradients.gradient.begin(); miter != miter_end; ++ miter) 
	if (miter->first != feature_bleu) {
	  const double feature = miter->second;
	  update(weights[miter->first], - feature * eta);
	}
      
      // projection...
      if (weight_norm > 1.0 / lambda)
	rescale(std::sqrt(1.0 / (lambda * weight_norm)));
      
      objective += double(log(Z_correct) - log(Z)) * weight_scale;
      
      if (weight_scale < 0.01 || 100 < weight_scale) {
	weights[feature_bleu] = 0.0;
	weights *= weight_scale;
	weight_scale = 1.0;
	weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
      }

      ++ samples;
    }
  }  
  
  
  void update(double& x, const double& alpha)
  {
    weight_norm += 2.0 * x * alpha * weight_scale + alpha * alpha;
    x += alpha / weight_scale;
  }
  
  void rescale(const double scaling)
  {
    weight_norm *= scaling * scaling;
    if (scaling != 0.0)
      weight_scale *= scaling;
    else {
      weight_scale = 1.0;
      std::fill(weights.begin(), weights.end(), 0.0);
    }
  }
  
  size_t samples;
  size_t epoch;
  double lambda;
};

struct OptimizerSGDL1 : public OptimizerBase
{
  typedef cicada::WeightVector<double> penalty_set_type;
  
  OptimizerSGDL1(const hypergraph_set_type&           __graphs,
		 const feature_function_ptr_set_type& __features)
    : OptimizerBase(__graphs, __features),
      samples(0),
      epoch(0),
      lambda(C / __graphs.size()) {}
  
  void initialize()
  {
    samples = 0;

    weights[feature_bleu] = 0.0;
    
    weight_scale = 1.0;
    
    objective = 0.0;
    penalties.clear();
    penalty = 0.0;
  }
  
  void finalize()
  {
    
  }
  
  void operator()(const int segment)
  {
    //const double eta = 1.0 / (1.0 + double(epoch) / graphs.size());
    const double eta = 1.0 / (lambda * (epoch + 2));
    ++ epoch;
    
    // cummulative penalty
    penalty += eta * lambda;
    
    // we will maximize, not minimize...
    if (OptimizerBase::operator()(segment)) {
      // we have gradients_correct/gradients/Z_correct/Z
      
      gradient_type::const_iterator citer_end = gradients_correct.gradient.end();
      for (gradient_type::const_iterator citer = gradients_correct.gradient.begin(); citer != citer_end; ++ citer) 
	if (citer->first != feature_bleu) {
	  const double feature = citer->second;
	  
	  weights[citer->first] += eta * feature;
	  apply(weights[citer->first], penalties[citer->first], penalty);
	}
      
      gradient_type::const_iterator miter_end = gradients.gradient.end();
      for (gradient_type::const_iterator miter = gradients.gradient.begin(); miter != miter_end; ++ miter) 
	if (miter->first != feature_bleu) {
	  const double feature = miter->second;
	  
	  weights[miter->first] -= eta * feature;
	  apply(weights[miter->first], penalties[miter->first], penalty);
	}
      
      objective += double(log(Z_correct) - log(Z)) * weight_scale;
      ++ samples;
    }
  }

  void apply(double& x, double& penalty, const double& cummulative)
  {
    const double x_half = x;
    if (x > 0.0)
      x = std::max(0.0, x - penalty - cummulative);
    else if (x < 0.0)
      x = std::min(0.0, x - penalty + cummulative);
    penalty += x - x_half;
  }

  size_t samples;
  size_t epoch;
  double lambda;
  
  penalty_set_type penalties;
  double penalty;
};
