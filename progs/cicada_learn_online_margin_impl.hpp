//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_ONLINE_MARGIN_IMPL__HPP__
#define __CICADA_LEARN_ONLINE_MARGIN_IMPL__HPP__ 1

#include <vector>

#include <cicada/hypergraph.hpp>
#include <cicada/feature_vector.hpp>
#include <cicada/weight_vector.hpp>
#include <cicada/semiring.hpp>
#include <cicada/viterbi.hpp>
#include <cicada/operation/traversal.hpp>
#include <cicada/operation/functional.hpp>

#include <utils/mulvector2.hpp>
#include <utils/bithack.hpp>
#include <utils/mathop.hpp>
#include <utils/chunk_vector.hpp>

struct Margin
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::HyperGraph hypergraph_type;
  
  typedef hypergraph_type::feature_set_type   feature_set_type;
  typedef hypergraph_type::attribute_set_type attribute_set_type;
  
  typedef feature_set_type::feature_type     feature_type;
  typedef attribute_set_type::attribute_type attribute_type;
  
  typedef std::pair<feature_type, double> feature_value_type;
  
  typedef utils::mulvector2<feature_value_type, std::allocator<feature_value_type> > delta_set_type;
  typedef std::vector<feature_value_type, std::allocator<feature_value_type> > feature_value_set_type;

  struct less_feature_value
  {
    bool operator()(const feature_value_type& x, const feature_value_type& y) const
    {
      return x.first < y.first;
    }
  };

  typedef cicada::WeightVector<double> weight_set_type;
  
  virtual ~Margin() {}

  virtual void encode(const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle) = 0;
  
  void clear() { deltas.clear(); }

  void push_back(const feature_set_type& features)
  {
    delta.clear();
    delta.reserve(features.size());
    
    feature_set_type::const_iterator fiter_end = features.end();
    for (feature_set_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
      if (fiter->second != 0.0)
	delta.push_back(*fiter);
    
    if (delta.empty()) return;
    
    std::sort(delta.begin(), delta.end(), less_feature_value());
    
    deltas.push_back(delta.begin(), delta.end());
  }

  delta_set_type         deltas;
  feature_value_set_type delta;
};

struct MarginDerivation : public Margin
{
  // full-derivation margin

  typedef cicada::HyperGraph hypergraph_type;
  
  typedef hypergraph_type::feature_set_type feature_set_type;
  
  struct traversal
  {
    typedef feature_set_type value_type;

    template <typename Edge, typename Iterator>
    void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
    {
      yield = edge.features;
      for (/**/; first != last; ++ first)
	yield += *first;
    }
  };

  
  void encode(const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle)
  {
    typedef cicada::semiring::Tropical<double>  weight_type;

    if (! forest.is_valid() || ! oracle.is_valid()) return;
    
    weight_type weight_forest;
    weight_type weight_oracle;
    
    feature_set_type features_forest;
    feature_set_type features_oracle;
    
    cicada::viterbi(oracle,
		    features_oracle,
		    weight_oracle,
		    traversal(),
		    cicada::operation::weight_function<weight_type >(weights));
    cicada::viterbi(forest,
		    features_forest,
		    weight_forest,
		    traversal(),
		    cicada::operation::weight_function<weight_type >(weights));

    bool weights_empty = weights.empty();
    if (! weights_empty) {
      int count_non_zero = 0;
      for (feature_type::id_type id = 0; id != weights.size() && ! count_non_zero; ++ id)
	count_non_zero += (weights[id] != 0.0);
      
      weights_empty = (! count_non_zero);
    }
    
    // already achieved optimum...
    if (! weights_empty && weight_oracle >= weight_forest) return;
    
    features_oracle -= features_forest;
    
    push_back(features_oracle);
  }
};

struct MarginViolation : public Margin
{
  // single max-violation node margin
  typedef cicada::semiring::Tropical<double>  weight_type;

  struct traversal
  {
    typedef feature_set_type value_type;

    template <typename Edge, typename Iterator>
    void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
    {
      yield = edge.features;
      for (/**/; first != last; ++ first)
	yield += *first;
    }
  };

  // Initially, I tried to implement on top of the inside/outside framework, but it seems to be easier
  // to implement a special functions...

  struct Inside
  {    
    typedef std::vector<weight_type, std::allocator<weight_type> >           weight_map_type;
    typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_map_type;
    
    typedef std::vector<weight_type, std::allocator<weight_type> >           weight_bin_type;
    typedef utils::chunk_vector<feature_set_type, 4096/sizeof(feature_set_type), std::allocator<feature_set_type> > feature_bin_type;

    struct attribute_int : public boost::static_visitor<attribute_set_type::int_type>
    {
      // we will not throw, but simply return zero. (TODO: return negative?)
      attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
      attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
      attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
    };
    
    Inside() : attr_prune_bin("prune-bin") {}
    
    void operator()(const weight_set_type& weights, const hypergraph_type& forest)
    {
      weights_inside.clear();
      features_inside.clear();
      weights_bin.clear();
      features_bin.clear();

      cicada::operation::weight_function<weight_type > function(weights);
      
      weights_inside.resize(forest.nodes.size());
      features_inside.resize(forest.nodes.size());
      
      hypergraph_type::node_set_type::const_iterator niter_end = forest.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = forest.nodes.begin(); niter != niter_end; ++ niter) {
	typedef hypergraph_type::node_type node_type;
	
	const node_type& node = *niter;
	
	weight_type& weight = weights_inside[node.id];
	feature_set_type& features = features_inside[node.id];
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  typedef hypergraph_type::edge_type edge_type;
	  
	  const edge_type& edge = forest.edges[*eiter];
	  
	  weight_type score = function(edge);
	  edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	    score *= weights_inside[*niter];

	  bool features_computed = false;
	  
	  if (score > weight) {
	    weight = score;
	    
	    features = edge.features;
	    edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	    for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	      features += features_inside[*niter];

	    features_computed = true;
	  }
	  
	  attribute_set_type::const_iterator piter = edge.attributes.find(attr_prune_bin);
	  if (piter == edge.attributes.end()) continue;
	  
	  const int bin_pos = boost::apply_visitor(attribute_int(), piter->second);
	  
	  if (bin_pos < 0) continue;
	  
	  if (bin_pos >= weights_bin.size())
	    weights_bin.resize(bin_pos + 1);
	  if (bin_pos >= features_bin.size())
	    features_bin.resize(bin_pos + 1);
	  
	  if (score > weights_bin[bin_pos])  {
	    weights_bin[bin_pos] = score;

	    if (features_computed)
	      features_bin[bin_pos] = features;
	    else {
	      features_bin[bin_pos] = edge.features;
	      edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	      for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
		features_bin[bin_pos] += features_inside[*niter];
	    }
	  }
	}
      }
    }
    
    weight_map_type  weights_inside;
    feature_map_type features_inside;
    
    weight_bin_type  weights_bin;
    feature_bin_type features_bin;

    const attribute_type attr_prune_bin;
  };
  

  Inside inside_forest;
  Inside inside_oracle;
};

struct MarginViolationSingle : public MarginViolation
{
  // single max-violation node margin
  
  MarginViolationSingle()  {}
  
  void encode(const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle)
  {
    if (! forest.is_valid() || ! oracle.is_valid()) return;

    bool weights_empty = weights.empty();
    if (! weights_empty) {
      int count_non_zero = 0;
      for (feature_type::id_type id = 0; id != weights.size() && ! count_non_zero; ++ id)
	count_non_zero += (weights[id] != 0.0);
      
      weights_empty = (! count_non_zero);
    }
    
    if (weights_empty) {
      weight_type weight_forest;
      weight_type weight_oracle;
      
      feature_set_type features_forest;
      feature_set_type features_oracle;
      
      cicada::viterbi(oracle,
		      features_oracle,
		      weight_oracle,
		      traversal(),
		      cicada::operation::weight_function<weight_type >(weights));
      cicada::viterbi(forest,
		      features_forest,
		      weight_forest,
		      traversal(),
		      cicada::operation::weight_function<weight_type >(weights));
      
      features_oracle -= features_forest;
      
      
      push_back(features_oracle);
    } else {
      // First, take maximum from oracle wrt weights
      inside_oracle(weights, oracle);
      
      // Second, take maximum from forest wrt weights
      inside_forest(weights, forest);
      
      // Third, compute the largest margin
      
      const size_type bin_max = utils::bithack::min(inside_forest.weights_bin.size(), inside_oracle.weights_bin.size());
      
      if (! bin_max) return;
      
      size_type margin_pos(-1);
      double    margin_max(- std::numeric_limits<double>::infinity());
      
      for (size_type bin = 0; bin != bin_max; ++ bin) 
	if (! inside_forest.features_bin[bin].empty() && ! inside_oracle.features_bin[bin].empty()) {
	  const double margin = (cicada::semiring::log(inside_forest.weights_bin[bin])
				 - cicada::semiring::log(inside_oracle.weights_bin[bin]));
	  
	  if (margin > 0.0 && margin > margin_max) {
	    margin_pos = bin;
	    margin_max = margin;
	  }
	}
      
      // found the best margin!
      if (margin_pos != size_type(-1)) {
	inside_oracle.features_bin[margin_pos] -= inside_forest.features_bin[margin_pos];
	
	push_back(inside_oracle.features_bin[margin_pos]);
      }
    }
  }
  
};

struct MarginViolationAll : public MarginViolation
{
  // multiple max-violation node margin
  
  void encode(const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle)
  {
    if (! forest.is_valid() || ! oracle.is_valid()) return;

    bool weights_empty = weights.empty();
    if (! weights_empty) {
      int count_non_zero = 0;
      for (feature_type::id_type id = 0; id != weights.size() && ! count_non_zero; ++ id)
	count_non_zero += (weights[id] != 0.0);
      
      weights_empty = (! count_non_zero);
    }

    if (weights_empty) {
      weight_type weight_forest;
      weight_type weight_oracle;
      
      feature_set_type features_forest;
      feature_set_type features_oracle;
      
      cicada::viterbi(oracle,
		      features_oracle,
		      weight_oracle,
		      traversal(),
		      cicada::operation::weight_function<weight_type >(weights));
      cicada::viterbi(forest,
		      features_forest,
		      weight_forest,
		      traversal(),
		      cicada::operation::weight_function<weight_type >(weights));
      
      features_oracle -= features_forest;
      
      push_back(features_oracle);
    } else {
      // First, take maximum from oracle wrt weights
      inside_oracle(weights, oracle);
      
      // Second, take maximum from forest wrt weights
      inside_forest(weights, forest);
      
      // Third, compute all the margins
      
      const size_type bin_max = utils::bithack::min(inside_forest.weights_bin.size(), inside_oracle.weights_bin.size());
      
      if (! bin_max) return;
      
      for (size_type bin = 0; bin != bin_max; ++ bin)
	if (! inside_forest.features_bin[bin].empty() && ! inside_oracle.features_bin[bin].empty()) {
	  const double margin = (cicada::semiring::log(inside_forest.weights_bin[bin])
				 - cicada::semiring::log(inside_oracle.weights_bin[bin]));
	  
	  if (margin > 0.0) {
	    inside_oracle.features_bin[bin] -= inside_forest.features_bin[bin];
	    
	    push_back(inside_oracle.features_bin[bin]);
	  }
	}
    }
  }
};

#endif
