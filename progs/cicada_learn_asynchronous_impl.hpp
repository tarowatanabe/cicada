//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_ASYNCHRONOUS_IMPL__HPP__
#define __CICADA_LEARN_ASYNCHRONOUS_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>

#include <sstream>
#include <vector>

#include "cicada/semiring.hpp"
#include "cicada/eval.hpp"

#include "cicada/kbest.hpp"
#include "cicada/operation/traversal.hpp"
#include "cicada/operation/functional.hpp"
#include "cicada/optimize_qp.hpp"
#include "cicada/optimize.hpp"
#include "cicada/feature/scorer.hpp"
#include "cicada/prune_beam.hpp"
#include "cicada/apply_exact.hpp"
#include "cicada/apply_cube_prune.hpp"
#include "cicada/viterbi.hpp"
#include "cicada/expected_ngram.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/base64.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"
#include "utils/config.hpp"
#include "utils/mathop.hpp"
#include "utils/indexed_trie.hpp"
#include "utils/getline.hpp"
#include "utils/mulvector2.hpp"
#include "utils/chunk_vector.hpp"

#include "cicada_learn_online_regularize_impl.hpp"
#include "cicada_learn_online_rate_impl.hpp"

#include <boost/tokenizer.hpp>

#include <codec/lz4.hpp>

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type  score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;

typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type> > feature_function_ptr_set_type;
typedef feature_function_ptr_set_type function_document_type;

typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;
typedef hypergraph_set_type hypergraph_document_type;

struct LearnBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;
  
};

struct LearnXBLEUBase : public LearnBase
{
  typedef cicada::semiring::Log<double> weight_type;
  
  static weight_type brevity_penalty(const double x)
  {
    typedef cicada::semiring::traits<weight_type> traits_type;

    // return (std::exp(x) - 1) / (1.0 + std::exp(1000.0 * x)) + 1.0;
      
    return ((traits_type::exp(x) - traits_type::one()) / (traits_type::one() + traits_type::exp(1000.0 * x))) + traits_type::one();
  }
    
  static weight_type derivative_brevity_penalty(const double x)
  {
    typedef cicada::semiring::traits<weight_type> traits_type;
    
    const weight_type expx     = traits_type::exp(x);
    const weight_type expxm1   = expx - traits_type::one();
    const weight_type exp1000x = traits_type::exp(1000.0 * x);
    const weight_type p1exp1000x = traits_type::one() + exp1000x;
      
    return (expx / p1exp1000x) - ((expxm1 * weight_type(1000.0) * exp1000x) / (p1exp1000x * p1exp1000x));
      
    //return expx / (1.0 + exp1000x) - boost::math::expm1(x) * (1000.0 * exp1000x) / ((1.0 + exp1000x) * (1.0 + exp1000x))
  }
    
  static weight_type clip_count(const weight_type& x, const weight_type& clip)
  {
    typedef cicada::semiring::traits<weight_type> traits_type;
    
    //return (x - clip) / (1.0 + std::exp(1000.0 * (x - clip))) + clip;
    return (weight_type(x - clip) / (traits_type::one() + traits_type::exp(1000.0 * (x - clip)))) + weight_type(clip);
  }
  
  static weight_type derivative_clip_count(const weight_type& x, const weight_type& clip)
  {
    typedef cicada::semiring::traits<weight_type> traits_type;
    
    const weight_type exp1000xmc = traits_type::exp(1000.0 * (x - clip));
    const weight_type p1exp1000xmc = exp1000xmc + traits_type::one();
    
    return (traits_type::one() / p1exp1000xmc) - ((weight_type(x - clip) * weight_type(1000.0) * exp1000xmc) / (p1exp1000xmc * p1exp1000xmc));
    
    //return 1.0 / (1.0 + exp1000x) - (x - clip) * (1000.0 * exp1000x) / ((1.0 + exp1000x) * (1.0 + exp1000x));
  }
  
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
  typedef std::vector<gradient_type, std::allocator<gradient_type> >       gradients_type;

  typedef cicada::FeatureVector<double, std::allocator<double> > gradient_xbleu_type;
  
  typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
  
  typedef cicada::Symbol       word_type;
  typedef cicada::SymbolVector ngram_type;
  
  typedef utils::indexed_trie<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> > index_set_type;
  typedef utils::simple_vector<index_set_type::id_type, std::allocator<index_set_type::id_type> > id_set_type;
  typedef std::vector<id_set_type, std::allocator<id_set_type> > id_map_type;
  
  struct Count
  {
    weight_type c;
    weight_type mu_prime;
    
    Count() : c(), mu_prime() {}
  };
  typedef Count count_type;
  typedef std::vector<count_type, std::allocator<count_type> > count_set_type;
  typedef std::vector<ngram_type, std::allocator<ngram_type> > ngram_set_type;
  
  struct CollectCounts
  {
    CollectCounts(index_set_type& __index,
		  ngram_set_type& __ngrams,
		  count_set_type& __counts,
		  id_map_type& __ids)
      : index(__index), ngrams(__ngrams), counts(__counts), ids(__ids) {}
      
    template <typename Edge, typename Weight, typename Counts>
    void operator()(const Edge& edge, const Weight& weight, Counts& __counts)
    {
	
    }
      
    template <typename Edge, typename Weight, typename Counts, typename Iterator>
    void operator()(const Edge& edge, const Weight& weight, Counts& __counts, Iterator first, Iterator last)
    {
      if (first == last) return;
	
      index_set_type::id_type id = index.root();
      for (Iterator iter = first; iter != last; ++ iter)
	id = index.push(id, *iter);
	
      if (id >= ngrams.size())
	ngrams.resize(id + 1);
      if (id >= counts.size())
	counts.resize(id + 1);
	
      counts[id].c += weight;
	
      if (ngrams[id].empty())
	ngrams[id] = ngram_type(first, last);
	
      ids[edge.id].push_back(id);
    }
      
    index_set_type& index;
    ngram_set_type& ngrams;
    count_set_type& counts;
    id_map_type& ids;
  };

  typedef cicada::semiring::Tuple<weight_type> ngram_weight_type;
  typedef cicada::semiring::Expectation<weight_type, ngram_weight_type> bleu_weight_type;
  typedef std::vector<bleu_weight_type, std::allocator<bleu_weight_type> > bleu_weights_type;
  
  struct bleu_function
  {
    typedef bleu_weight_type value_type;
      
    bleu_function(const ngram_set_type& __ngrams,
		  const count_set_type& __counts,
		  const id_map_type& __ids,
		  const weight_set_type& __weights,
		  const double& __scale)
      : ngrams(__ngrams), counts(__counts), ids(__ids),
	weights(__weights), scale(__scale) {}
    
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      const double margin = cicada::dot_product(edge.features, weights);
      const weight_type weight = cicada::semiring::traits<weight_type>::exp(margin * scale);
	
      value_type bleu(weight, ngram_weight_type(order * 2, weight_type()));
	
      id_set_type::const_iterator iter_end = ids[edge.id].end();
      for (id_set_type::const_iterator iter = ids[edge.id].begin(); iter != iter_end; ++ iter) {
	const int n = ngrams[*iter].size();
	const int index = (n - 1) << 1;
	  
	bleu.r[index] += weight;
	bleu.r[index + 1] += counts[*iter].mu_prime * weight;
      }
	
      return bleu;
    }
      
    const ngram_set_type&  ngrams;
    const count_set_type&  counts;
    const id_map_type&     ids;
    const weight_set_type& weights;
    const double           scale;
  };
    
  struct bleu_gradient_function
  {
    struct value_type
    {
      value_type(const hypergraph_type::edge_type& __edge)
	: edge(__edge) {}
	
      friend
      value_type operator*(value_type x, const bleu_weight_type& weight)
      {
	x.inside_outside = weight;
	return x;
      }
	
      bleu_weight_type inside_outside;
      const hypergraph_type::edge_type& edge;
    };
      
    bleu_gradient_function() {}
      
    value_type operator()(const hypergraph_type::edge_type& edge) const
    {
      return value_type(edge);
    }
  };
    
  struct bleu_gradient_type
  {
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > accumulated_type;
    typedef std::vector<accumulated_type, std::allocator<accumulated_type> > accumulated_set_type;

    struct value_type
    {
      value_type& operator+=(const bleu_gradient_function::value_type& x)
      {
	const double margin = cicada::dot_product(x.edge.features, impl.weights) * impl.scale_var;
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(margin * impl.scale_const);
	
	bleu_weight_type bleu(weight, ngram_weight_type(order * 2, weight_type()));
	  
	id_set_type::const_iterator iter_end = impl.ids[x.edge.id].end();
	for (id_set_type::const_iterator iter = impl.ids[x.edge.id].begin(); iter != iter_end; ++ iter) {
	  const int n = impl.ngrams[*iter].size();
	  const int index = (n - 1) << 1;
	    
	  bleu.r[index] += weight;
	  bleu.r[index + 1] += impl.counts[*iter].mu_prime * weight;
	}
	  
	bleu *= x.inside_outside;
	  
	// accumulate gradients....
	for (int n = 1; n <= order; ++ n) 
	  if (impl.matched[n] > weight_type()) {
	    const int index = (n - 1) << 1;
	    const weight_type scale_matched = bleu.r[index + 1] - bleu.p * impl.matched[n];
	    const weight_type scale_hypo    = bleu.r[index]     - bleu.p * impl.hypo[n];
	      
	    feature_set_type::const_iterator fiter_end = x.edge.features.end();
	    for (feature_set_type::const_iterator fiter = x.edge.features.begin(); fiter != fiter_end; ++ fiter)
	      if (fiter->second != 0.0) {
		const weight_type value(fiter->second * impl.scale_const);
		  
		impl.dM[n][fiter->first] += value * scale_matched;
		impl.dH[n][fiter->first] += value * scale_hypo;
	      }
	  }
	  
	return *this;
      }
	
      value_type(bleu_gradient_type& __impl) : impl(__impl) {}
	
      bleu_gradient_type& impl;
    };
      
    value_type operator[](size_t id) { return value_type(*this); }
      
    bleu_gradient_type(const ngram_set_type& __ngrams,
		       const count_set_type& __counts,
		       const id_map_type& __ids,
		       const weights_type& __matched,
		       const weights_type& __hypo,
		       const weight_set_type& __weights,
		       const double& __scale_var,
		       const double& __scale_const) 
      : ngrams(__ngrams), counts(__counts), ids(__ids),
	matched(__matched), hypo(__hypo),
	weights(__weights), scale_var(__scale_var), scale_const(__scale_const),
	dM(order + 1),
	dH(order + 1) {}
    
    const ngram_set_type&  ngrams;
    const count_set_type&  counts;
    const id_map_type&     ids;
      
    const weights_type&    matched;
    const weights_type&    hypo;
      
    const weight_set_type& weights;
    const double           scale_var;
    const double           scale_const;
      
    accumulated_set_type dM;
    accumulated_set_type dH;
  };
  
  typedef cicada::semiring::Expectation<weight_type, weight_type> entropy_weight_type;
  
  struct entropy_function
  {
    typedef entropy_weight_type value_type;
      
    entropy_function(const weight_set_type& __weights, const double& __scale) : weights(__weights), scale(__scale) {}
      
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      const double value = cicada::dot_product(edge.features, weights) * scale;
      const weight_type weight = cicada::semiring::traits<weight_type>::exp(value);
	
      return value_type(weight, weight * weight_type(value));
    }
      
    const weight_set_type& weights;
    const double scale;
  };

  struct entropy_gradient_function
  {
    struct value_type
    {
      value_type(const feature_set_type& __features, const weight_set_type& __weights, const double& __scale_var, const double& __scale_const)
	: features(__features), weights(__weights), scale_var(__scale_var), scale_const(__scale_const) {}
      
      friend
      value_type operator*(value_type x, const entropy_weight_type& weight)
      {
	x.inside_outside = weight;
	return x;
      }
	
      entropy_weight_type inside_outside;
	
      const feature_set_type& features;
      const weight_set_type& weights;
      const double scale_var;
      const double scale_const;
    };
      
    entropy_gradient_function(const weight_set_type& __weights, const double& __scale_var, const double& __scale_const)
      : weights(__weights), scale_var(__scale_var), scale_const(__scale_const) {}
      
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      return value_type(edge.features, weights, scale_var, scale_const);
    }
    
    const weight_set_type& weights;
    const double scale_var;
    const double scale_const;
  };

  struct entropy_gradient_type
  {
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > accumulated_type;

    struct proxy_type
    {
      proxy_type(accumulated_type& __dZ, accumulated_type& __dR) : dZ(__dZ), dR(__dR) {}
	
      proxy_type& operator+=(const entropy_gradient_function::value_type& x) 
      {
	const double value = cicada::dot_product(x.features, x.weights) * x.scale_var;
	const double log_p_e = value * x.scale_const;
	const weight_type p_e = cicada::semiring::traits<weight_type>::exp(log_p_e);
	
	// dZ += \lnabla p_e * x.inside_outside.p;
	// dR += (1 + \log p_e) * \nalba p_e * x.inside_outside.p + \lnabla p_e * x.inside_outside.r;
	  
	feature_set_type::const_iterator fiter_end = x.features.end();
	for (feature_set_type::const_iterator fiter = x.features.begin(); fiter != fiter_end; ++ fiter) 
	  if (fiter->second != 0.0) {
	    const weight_type value(fiter->second * x.scale_const);
	      
	    dZ[fiter->first] += value * p_e * x.inside_outside.p;
	    dR[fiter->first] += (weight_type(1.0 + log_p_e) * value * p_e * x.inside_outside.p + value * p_e * x.inside_outside.r);
	  }
	
	return *this;
      }
	
      accumulated_type& dZ;
      accumulated_type& dR;
    };
      
    typedef proxy_type value_type;
      
    proxy_type operator[](size_t id) { return proxy_type(dZ, dR); }
      
    accumulated_type dZ;
    accumulated_type dR;
  };

  typedef std::vector<entropy_weight_type, std::allocator<entropy_weight_type> > entropy_weights_type;

  LearnXBLEUBase() { clear(); }

  void clear()
  {
    counts_matched.reserve(order + 1);
    counts_hypo.reserve(order + 1);
    
    gradients_matched.reserve(order + 1);
    gradients_hypo.reserve(order + 1);

    counts_matched.resize(order + 1);
    counts_hypo.resize(order + 1);
    
    gradients_matched.resize(order + 1);
    gradients_hypo.resize(order + 1);
    gradients_entropy.clear();
    
    std::fill(counts_matched.begin(), counts_matched.end(), weight_type());
    std::fill(counts_hypo.begin(), counts_hypo.end(), weight_type());
    counts_reference = weight_type();
    counts_entropy   = weight_type();
    norm_entropy = 0;
    
    for (int n = 1; n <= order; ++ n) {
      gradients_matched[n].clear();
      gradients_hypo[n].clear();
    }
  }
  
  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer, const double& scale_var, const double& scale_const)
  {
    if (! forest.is_valid()) return;
    
    const cicada::eval::BleuScorer* bleu = dynamic_cast<const cicada::eval::BleuScorer*>(scorer.get());
    if (! bleu)
      throw std::runtime_error("we do not have bleu scorer...");
    
    // first, collect expected ngrams
    index.clear();
    counts.clear();
    ngrams.clear();
    ids.clear();
    ids.resize(forest.edges.size());
    
    cicada::expected_ngram(forest,
			   cicada::operation::weight_scaled_function<weight_type>(weights, scale_var * scale_const),
			   CollectCounts(index, ngrams, counts, ids),
			   index,
			   order);
    
    // second, commpute clipped ngram counts (\mu')
    hypo.clear();
    matched.clear();
    hypo.resize(order + 1);
    matched.resize(order + 1);
    
    for (size_type i = 0; i != ngrams.size(); ++ i) 
      if (! ngrams[i].empty()) {
	const size_type    order = ngrams[i].size();
	const weight_type& count = counts[i].c;
	const weight_type  clip = bleu->find(ngrams[i]);
	
	counts[i].mu_prime = derivative_clip_count(count, clip);
	    
	// collect counts for further inside/outside
	matched[order] += counts[i].c * counts[i].mu_prime;
	hypo[order]    += counts[i].c;
	
	// collect global counts
	counts_matched[order] += clip_count(count, clip);
	counts_hypo[order]    += counts[i].c;
      }
    
    counts_reference += bleu->reference_length(hypo[1]);
    
    // third, collect feature expectation, \hat{m} - m and \hat{h} - h
    bleu_inside.clear();
    bleu_inside.resize(forest.nodes.size(), bleu_weight_type());
    
    bleu_gradient_type bleu_gradient(ngrams, counts, ids,
				     matched, hypo,
				     weights, scale_var, scale_const);
    
    cicada::inside_outside(forest,
			   bleu_inside,
			   bleu_gradient,
			   bleu_function(ngrams, counts, ids, weights, scale_var * scale_const),
			   bleu_gradient_function());
    
    for (int n = 1; n <= order; ++ n) {
      const weight_type& Z = bleu_inside.back().p;
      const bleu_gradient_type::accumulated_set_type& dM = bleu_gradient.dM;
      const bleu_gradient_type::accumulated_set_type& dH = bleu_gradient.dH;
      
      bleu_gradient_type::accumulated_type::const_iterator miter_end = dM[n].end();
      for (bleu_gradient_type::accumulated_type::const_iterator miter = dM[n].begin(); miter != miter_end; ++ miter)
	gradients_matched[n][miter->first] += miter->second / Z;
      
      bleu_gradient_type::accumulated_type::const_iterator hiter_end = dH[n].end();
      for (bleu_gradient_type::accumulated_type::const_iterator hiter = dH[n].begin(); hiter != hiter_end; ++ hiter)
	gradients_hypo[n][hiter->first] += hiter->second / Z;
    }
    
    // forth, compute entorpy...
    entropy_inside.clear();
    entropy_inside.resize(forest.nodes.size(), entropy_weight_type());
    
    entropy_gradient_type entropy_gradient;
    
    cicada::inside_outside(forest,
			   entropy_inside,
			   entropy_gradient,
			   entropy_function(weights, scale_var * scale_const),
			   entropy_gradient_function(weights, scale_var, scale_const));
    
    const weight_type& Z = entropy_inside.back().p;
    const weight_type& R = entropy_inside.back().r;
    
    const weight_type entropy = weight_type(cicada::semiring::log(Z)) - (R / Z);
    
    counts_entropy += entropy;
    ++ norm_entropy;
    
    const entropy_gradient_type::accumulated_type& dZ = entropy_gradient.dZ;
    const entropy_gradient_type::accumulated_type& dR = entropy_gradient.dR;
    
    // compute...
    // \frac{\nabla Z}{Z} - \frac{Z \nabla \bar{r} - \bar{r} \nabla Z}{Z^2}
    
    entropy_gradient_type::accumulated_type::const_iterator ziter_end = dZ.end();
    for (entropy_gradient_type::accumulated_type::const_iterator ziter = dZ.begin(); ziter != ziter_end; ++ ziter)
      gradients_entropy[ziter->first] += ziter->second * ((cicada::semiring::traits<weight_type>::one() / Z) + R / (Z * Z));
    
    entropy_gradient_type::accumulated_type::const_iterator riter_end = dR.end();
    for (entropy_gradient_type::accumulated_type::const_iterator riter = dR.begin(); riter != riter_end; ++ riter)
      gradients_entropy[riter->first] -= riter->second / Z;
  }
  
  std::pair<double, bool> encode(gradient_xbleu_type& g)
  {
    g.clear();
    
    // check if we have any matching...
    if (counts_reference <= weight_type() || counts_hypo[1] <= weight_type() || counts_matched[1] <= weight_type())
      return std::make_pair(0.0, false);
    
    // smoothing...
    double smoothing = 1e-40;
    for (int n = 1; n <= order; ++ n) {
      if (counts_hypo[n] > weight_type() && counts_matched[n] <= weight_type())
	counts_matched[n] = std::min(weight_type(smoothing), counts_hypo[n]);
      
      smoothing *= 0.1;
    }
    
    // compute P
    double P = 0.0;
    for (int n = 1; n <= order; ++ n)
      if (counts_hypo[n] > weight_type())
	P += (1.0 / order) * (cicada::semiring::log(counts_matched[n]) - cicada::semiring::log(counts_hypo[n]));
    
    if (! std::isfinite(P))
      return std::make_pair(0.0, false);
    
    // compute C and B
    const weight_type C      = counts_reference / counts_hypo[1];
    const double      minusC = 1.0 - C;
    const weight_type B      = brevity_penalty(minusC);
    
    // for computing g...
    const weight_type exp_P = cicada::semiring::traits<weight_type>::exp(P);
    const weight_type C_dC  = C * derivative_brevity_penalty(minusC);
    
    const weight_type objective_bleu = exp_P * B;
    const weight_type factor_entropy = 1.0 / norm_entropy;
    const weight_type entropy = counts_entropy * factor_entropy;
    const weight_type factor_order = 1.0 / order;

    const double objective = - objective_bleu - temperature * entropy;
    
    // entropy...
    if (temperature != 0.0) {
      const weight_type factor_temp(- temperature);
      
      gradient_type::const_iterator eiter_end = gradients_entropy.end();
      for (gradient_type::const_iterator eiter = gradients_entropy.begin(); eiter != eiter_end; ++ eiter) {
	double& grad = g[eiter->first];
	grad = factor_temp * factor_entropy * eiter->second;

	if (! std::isfinite(grad))
	  return std::make_pair(objective, false);
      }
    }
    
    // we will collect minus gradient for minimizing negative-xBLEU
    for (int n = 1; n <= order; ++ n) 
      if (counts_hypo[n] > weight_type()) {
	const weight_type factor_matched = - (exp_P * B * factor_order) / counts_matched[n];
	const weight_type factor_hypo    = - (exp_P * B * factor_order) / counts_hypo[n];
	
	gradient_type::const_iterator miter_end = gradients_matched[n].end();
	for (gradient_type::const_iterator miter = gradients_matched[n].begin(); miter != miter_end; ++ miter) {
	  double& grad = g[miter->first];
	  grad += factor_matched * miter->second;
	  
	  if (! std::isfinite(grad))
	    return std::make_pair(objective, false);
	}
	
	gradient_type::const_iterator hiter_end = gradients_hypo[n].end();
	for (gradient_type::const_iterator hiter = gradients_hypo[n].begin(); hiter != hiter_end; ++ hiter) {
	  double& grad = g[hiter->first];
	  grad -= factor_hypo * hiter->second;

	  if (! std::isfinite(grad))
	    return std::make_pair(objective, false);
	}
      }
    
    if (counts_hypo[1] > weight_type()) {
      const weight_type factor_hypo = - (exp_P * C_dC) / counts_hypo[1];
      
      gradient_type::const_iterator hiter_end = gradients_hypo[1].end();
      for (gradient_type::const_iterator hiter = gradients_hypo[1].begin(); hiter != hiter_end; ++ hiter) {
	double& grad = g[hiter->first];
	grad += factor_hypo * hiter->second;
	
	if (! std::isfinite(grad))
	  return std::make_pair(objective, false);
      }
    }
    
    return std::make_pair(objective, ! g.empty());
  }
  
  // required for learn()
  weights_type counts_matched;
  weights_type counts_hypo;
  weight_type  counts_reference;
  weight_type  counts_entropy;
  int          norm_entropy;
  
  gradients_type gradients_matched;
  gradients_type gradients_hypo;
  gradient_type  gradients_entropy;
  
  // local variables
  index_set_type index;
  ngram_set_type ngrams;
  count_set_type counts;
  id_map_type    ids;
  
  bleu_weights_type bleu_inside;
  weights_type      matched;
  weights_type      hypo;
  entropy_weights_type entropy_inside;
};

struct LearnXBLEU : public LearnXBLEUBase
{
  LearnXBLEU(Regularize& __regularizer,
	     Rate& __rate)
    : regularizer(__regularizer),
      rate(__rate) {}
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }

  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer)
  {
    LearnXBLEUBase::encode(id, weights, forest, oracle, scorer, regularizer.scale(), scale);
  }

  void update(weight_set_type& weights, const feature_set_type& updates)
  {
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    feature_set_type::const_iterator giter_end = updates.end();
    for (feature_set_type::const_iterator giter = updates.begin(); giter != giter_end; ++ giter)
      regularizer.update(weights, giter->first, giter->second, rate(giter->first, giter->second));
    
    regularizer.postprocess(weights, eta);
  }
  
  double learn(weight_set_type& weights, feature_set_type& updates)
  {
    updates.clear();
    
    // compute gradient...
    const std::pair<double, bool> objective = LearnXBLEUBase::encode(g);
    
    if (! objective.second) {
      clear();
      return objective.first;
    }
    
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    gradient_xbleu_type::const_iterator giter_end = g.end();
    for (gradient_xbleu_type::const_iterator giter = g.begin(); giter != giter_end; ++ giter) {
      const double amount = static_cast<double>(giter->second);
      
      regularizer.update(weights, giter->first, amount, rate(giter->first, amount));
      
      updates[giter->first] = amount;
    }
    
    regularizer.postprocess(weights, eta);
    
    clear();
    
    return objective.first;
  }
  
  gradient_xbleu_type g;
  
  Regularize& regularizer;
  Rate&       rate;
};


// logistic regression base...
struct LearnSoftmaxBase : public LearnBase
{
  typedef cicada::semiring::Log<double> weight_type;
  typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
  typedef cicada::FeatureVector<double, std::allocator<double> > expectation_type;
  
  LearnSoftmaxBase() { clear(); }
  
  struct weight_function
  {
    typedef weight_type value_type;

    weight_function(const weight_set_type& __weights, const double& __scale) : weights(__weights), scale(__scale) {}
      
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e
      return cicada::semiring::traits<value_type>::exp(cicada::dot_product(edge.features, weights) * scale);
    }
      
    const weight_set_type& weights;
    const double scale;
  };
    
  struct feature_function
  {
    struct value_type
    {
      value_type(const feature_set_type& __features,
		 const weight_set_type& __weights,
		 const double& __scale_var,
		 const double& __scale_const)
	: features(__features), weights(__weights), scale_var(__scale_var), scale_const(__scale_const) {}
	
      friend
      value_type operator*(value_type x, const weight_type& weight)
      {
	x.inside_outside = weight;
	return x;
      }
	
      weight_type inside_outside;
      const feature_set_type& features;
      const weight_set_type&  weights;
      const double            scale_var;
      const double            scale_const;
    };
      
    feature_function(const weight_set_type& __weights, const double& __scale_var, const double& __scale_const)
      : weights(__weights), scale_var(__scale_var), scale_const(__scale_const) {}
      
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      return value_type(edge.features, weights, scale_var, scale_const);
    }
      
    const weight_set_type& weights;
    const double           scale_var;
    const double           scale_const;
  };
    
  struct gradients_type
  {
    struct value_type
    {
      value_type(gradient_type& __gradient) : gradient(__gradient) {}
	
      value_type& operator+=(const feature_function::value_type& x)
      {
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(x.features, x.weights) * x.scale_var * x.scale_const) * x.inside_outside;
	  
	feature_set_type::const_iterator fiter_end = x.features.end();
	for (feature_set_type::const_iterator fiter = x.features.begin(); fiter != fiter_end; ++ fiter)
	  gradient[fiter->first] += weight_type(fiter->second * x.scale_const) * weight;
	  
	return *this;
      }
	
      gradient_type& gradient;
    };
      
    value_type operator[](size_t pos)
    {
      return value_type(gradient);
    }
      
    void clear() { gradient.clear(); }
      
    gradient_type gradient;
  };
  
  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer, const double& scale_var, const double& scale_const)
  {
    gradients_forest.clear();
    gradients_oracle.clear();
    
    inside_forest.clear();
    inside_oracle.clear();
    
    inside_forest.resize(forest.nodes.size());
    inside_oracle.resize(oracle.nodes.size());
    
    cicada::inside_outside(forest, inside_forest, gradients_forest,
			   weight_function(weights, scale_var * scale_const),
			   feature_function(weights, scale_var, scale_const));
    cicada::inside_outside(oracle, inside_oracle, gradients_oracle,
			   weight_function(weights, scale_var * scale_const),
			   feature_function(weights, scale_var, scale_const));
    
    const gradient_type& gradient_forest = gradients_oracle.gradient;
    const gradient_type& gradient_oracle = gradients_oracle.gradient;
    
    weight_type& Z_forest = inside_forest.back();
    weight_type& Z_oracle = inside_oracle.back();
    
    gradient_type::const_iterator fiter_end = gradient_forest.end();
    for (gradient_type::const_iterator fiter = gradient_forest.begin(); fiter != fiter_end; ++ fiter)
      gradient[fiter->first] += fiter->second / Z_forest;
    
    gradient_type::const_iterator oiter_end = gradient_oracle.end();
    for (gradient_type::const_iterator oiter = gradient_oracle.begin(); oiter != oiter_end; ++ oiter)
      gradient[oiter->first] -= oiter->second / Z_oracle;
    
    objective -= cicada::semiring::log(Z_oracle) - cicada::semiring::log(Z_forest);
    ++ samples;
  }
  
  void clear()
  {
    gradient.clear();
    objective = 0.0;
    samples = 0;
  }
  
  gradients_type gradients_forest;
  gradients_type gradients_oracle;
  weights_type   inside_forest;
  weights_type   inside_oracle;

  gradient_type  gradient;
  double objective;
  int samples;
};

// SoftmaxL2 learner
struct LearnSoftmax : public LearnSoftmaxBase
{
  LearnSoftmax(Regularize& __regularizer,
	       Rate& __rate)
    : regularizer(__regularizer),
      rate(__rate) {}
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }
  
  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer)
  {
    LearnSoftmaxBase::encode(id, weights, forest, oracle, scorer, regularizer.scale(), scale);
  }
  
  void update(weight_set_type& weights, const feature_set_type& updates)
  {
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    feature_set_type::const_iterator giter_end = updates.end();
    for (feature_set_type::const_iterator giter = updates.begin(); giter != giter_end; ++ giter)
      regularizer.update(weights, giter->first, giter->second, rate(giter->first, giter->second));
    
    regularizer.postprocess(weights, eta);
  }

  double learn(weight_set_type& weights, feature_set_type& updates)
  {
    updates.clear();

    if (! samples) {
      clear();
      return 0.0;
    }
    
    const double eta = rate();

    const size_type k = samples;
    const double k_norm = 1.0 / k;
    
    const double objective_normalized = objective * k_norm;
    
    regularizer.preprocess(weights, eta);
    
    gradient_type::const_iterator giter_end = gradient.end();
    for (gradient_type::const_iterator giter = gradient.begin(); giter != giter_end; ++ giter) {
      const double amount = static_cast<double>(giter->second) * k_norm;
      
      regularizer.update(weights, giter->first, amount, rate(giter->first, amount));
      
      // updates!
      updates[giter->first] = amount;
    }
    
    regularizer.postprocess(weights, eta);
    
    clear();
    
    return objective_normalized;
  }

  Regularize& regularizer;
  Rate&       rate;
};

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
  
  virtual ~Margin() {}

  virtual void encode(const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle) = 0;
  
  void clear() { deltas.clear(); }
  
  delta_set_type deltas;
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
    
    // already achieved optimum...
    if (weight_oracle >= weight_forest) return;
    
    features_oracle -= features_forest;
    
    deltas.push_back(features_oracle.begin(), features_oracle.end());
  }
};

struct MarginViolation : public Margin
{
  // single max-violation node margin
  typedef cicada::semiring::Tropical<double>  weight_type;

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
      deltas.push_back(inside_oracle.features_bin[margin_pos].begin(), inside_oracle.features_bin[margin_pos].end());
    }
  }
  
};

struct MarginViolationAll : public MarginViolation
{
  // multiple max-violation node margin
  
  void encode(const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle)
  {
    if (! forest.is_valid() || ! oracle.is_valid()) return;

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
	  deltas.push_back(inside_oracle.features_bin[bin].begin(), inside_oracle.features_bin[bin].end());
	}
      }
  }
};

struct LearnMargin : public LearnBase
{
  LearnMargin(Margin& __margin)
    : margin(__margin) {}

  void clear()
  {
    margin.clear();
  }
  
  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer)
  {
    margin.encode(weights, forest, oracle);
  }
  
  Margin&     margin;
};

struct LearnHinge : public LearnMargin
{
  LearnHinge(Margin& __margin,
	     Regularize& __regularizer,
	     Rate& __rate)
    : LearnMargin(__margin),
      regularizer(__regularizer),
      rate(__rate) {}
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }
  
  void update(weight_set_type& weights, const feature_set_type& updates)
  {
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    // udpate...
    feature_set_type::const_iterator fiter_end = updates.end();
    for (feature_set_type::const_iterator fiter = updates.begin(); fiter != fiter_end; ++ fiter)
      regularizer.update(weights, fiter->first, fiter->second, rate(fiter->first, fiter->second));
    
    regularizer.postprocess(weights, eta);
  }
  
  typedef std::vector<bool, std::allocator<bool> > suffered_set_type;
  suffered_set_type suffered;

  double learn(weight_set_type& weights, feature_set_type& updates)
  {
    updates.clear();
    
    if (margin.deltas.empty()) return 0.0;
  
    size_type k = 0;
    suffered.clear();
    suffered.resize(margin.deltas.size(), false);
    
    for (size_type i = 0; i != margin.deltas.size(); ++ i) {
      const double loss = 1.0 - cicada::dot_product(margin.deltas[i].begin(), margin.deltas[i].end(), weights, 0.0) * regularizer.scale();
      const bool suffer_loss = loss > 0.0;
      
      suffered[i] = suffer_loss;
      k += suffer_loss;
    }
    
    if (k) {
      //const double k_norm = 1.0 / (features.size());
      const double k_norm = 1.0 / k; // it is wrong, but works quite well in practice
      
      // udpate...
      for (size_type i = 0; i != margin.deltas.size(); ++ i)
	if (suffered[i]) {
	  Margin::delta_set_type::const_reference::const_iterator fiter_end = margin.deltas[i].end();
	  for (Margin::delta_set_type::const_reference::const_iterator fiter = margin.deltas[i].begin(); fiter != fiter_end; ++ fiter)
	    updates[fiter->first] -= k_norm * fiter->second;
	}
      
      if (! updates.empty())
	update(weights, updates);
    }
    
    margin.clear();
    
    return 0.0;
  }

  Regularize& regularizer;
  Rate&       rate;
};

// optimized-Pegasos learner
struct LearnOHinge : public LearnMargin
{
  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;
  typedef std::vector<int, std::allocator<int> >          index_type;

  typedef Margin::delta_set_type delta_set_type;

  struct HMatrix
  {
    HMatrix(const delta_set_type& __deltas, const index_type& __index) : deltas(__deltas), index(__index) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(deltas[index[i]].begin(), deltas[index[i]].end(), deltas[index[j]].begin(), deltas[index[j]].end(), 0.0);
    }
    
    const delta_set_type& deltas;
    const index_type& index;
  };
  
  struct MMatrix
  {
    MMatrix(const delta_set_type& __deltas, const index_type& __index) : deltas(__deltas), index(__index) {}
    
    template <typename __W>
    void operator()(__W& w, const alpha_type& alpha) const
    {
      const size_type model_size = index.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  delta_set_type::const_reference::const_iterator fiter_end = deltas[index[i]].end();
	  for (delta_set_type::const_reference::const_iterator fiter = deltas[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename __W>
    double operator()(const __W& w, const size_t& i) const
    {
      double dot = 0.0;
      delta_set_type::const_reference::const_iterator fiter_end = deltas[index[i]].end();
      for (delta_set_type::const_reference::const_iterator fiter = deltas[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename __W>
    void operator()(__W& w, const double& update, const size_t& i) const
    {
      delta_set_type::const_reference::const_iterator fiter_end = deltas[index[i]].end();
      for (delta_set_type::const_reference::const_iterator fiter = deltas[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const delta_set_type& deltas;
    const index_type& index;
  };
  

  LearnOHinge(Margin& __margin,
	      Regularize& __regularizer,
	      Rate& __rate)
    : LearnMargin(__margin),
      tolerance(0.1),
      regularizer(__regularizer),
      rate(__rate) {}
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }
  
  void update(weight_set_type& weights, const feature_set_type& updates)
  {
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    // udpate...
    feature_set_type::const_iterator fiter_end = updates.end();
    for (feature_set_type::const_iterator fiter = updates.begin(); fiter != fiter_end; ++ fiter)
      regularizer.update(weights, fiter->first, fiter->second, rate(fiter->first, fiter->second));
    
    regularizer.postprocess(weights, eta);
  }
  
  double learn(weight_set_type& weights, feature_set_type& updates)
  {
    updates.clear();
    
    if (margin.deltas.empty()) return 0.0;

    const double eta = rate();
    
    alpha.clear();
    f.clear();
    index.clear();
    
    for (size_type i = 0; i != margin.deltas.size(); ++ i)
      f.push_back(cicada::dot_product(margin.deltas[i].begin(), margin.deltas[i].end(), weights, 0.0));
    
    const double weight_scale_curr = regularizer.scale();
    
    regularizer.preprocess(weights, eta);
    
    double objective = 0.0;
    for (size_type i = 0; i != f.size(); ++ i) {
      const double loss = 1.0 - f[i] * weight_scale_curr;
      
      if (loss <= 0.0) continue;
      
      f[index.size()] = - (1.0 - f[i] * regularizer.scale());
      index.push_back(i);
      objective += loss;
    }
    
    objective /= f.size();
    
    // resize of f and alpha....
    f.resize(index.size());
    alpha.resize(index.size(), 0.0);

    {
      HMatrix H(margin.deltas, index);
      MMatrix M(margin.deltas, index);
      
      cicada::optimize::QPDCD()(alpha, f, H, M, eta, tolerance);
    }
    
    size_type actives = 0;
    size_type negatives = 0;
    for (size_type i = 0; i != index.size(); ++ i)
      if (alpha[i] > 0.0) {
	delta_set_type::const_reference::const_iterator fiter_end = margin.deltas[index[i]].end();
	for (delta_set_type::const_reference::const_iterator fiter = margin.deltas[index[i]].begin(); fiter != fiter_end; ++ fiter)
	  updates[fiter->first] -= alpha[i] * fiter->second / eta;
	
	++ actives;
	negatives += f[i] > 0.0;
      }
    
    feature_set_type::const_iterator uiter_end = updates.end();
    for (feature_set_type::const_iterator uiter = updates.begin(); uiter != uiter_end; ++ uiter)
      regularizer.update(weights, uiter->first, uiter->second, eta);
    
        if (debug >= 2)
	  std::cerr << "actives: " << actives << " negatives: " << negatives << " vectors: " << alpha.size() << std::endl;
	
    regularizer.postprocess(weights, eta);
    
    margin.clear();
    
    return 0.0;
  }

  double    tolerance;

  Regularize& regularizer;
  Rate&       rate;

  alpha_type    alpha;
  f_type        f;
  index_type    index;
};

struct YieldSentence
{
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;

  template <typename Line>
  score_ptr_type operator()(const Line& x, const scorer_ptr_type& scorer) const
  {
    return scorer->score(x->yield(cicada::operation::sentence_traversal()));
  }
};

struct YieldAlignment
{
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;

  template <typename Line>
  score_ptr_type operator()(const Line& x, const scorer_ptr_type& scorer) const
  {
    std::ostringstream os;
    os << x->yield(cicada::operation::alignment_traversal());
    
    return scorer->score(sentence_type(os.str()));
  }
};

struct YieldDependency
{
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;
  
  template <typename Line>
  score_ptr_type operator()(const Line& x, const scorer_ptr_type& scorer) const
  {
    std::ostringstream os;
    os << x->yield(cicada::operation::dependency_traversal());
    
    return scorer->score(sentence_type(os.str()));
  }
};


struct ViterbiSentence
{
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;

  typedef cicada::semiring::Logprob<double>  weight_type;

  sentence_type sentence;

  score_ptr_type operator()(const weight_set_type& weights,
			    const hypergraph_type& graph,
			    const scorer_ptr_type& scorer)
  {
    weight_type weight;
    
    cicada::viterbi(graph,
		    sentence,
		    weight,
		    cicada::operation::sentence_traversal(),
		    cicada::operation::weight_function<weight_type >(weights));
    
    return scorer->score(sentence);
  }

  score_ptr_type operator()(const hypergraph_type& graph,
			    const scorer_ptr_type& scorer,
			    const feature_type& feature,
			    const double& scale)
  {
    weight_type weight;
    
    cicada::viterbi(graph,
		    sentence,
		    weight,
		    cicada::operation::sentence_traversal(),
		    cicada::operation::single_scaled_function<weight_type >(feature, scale));
    
    return scorer->score(sentence);
  }
};

struct ViterbiAlignment
{
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;
  
  typedef cicada::semiring::Logprob<double>  weight_type;
  typedef cicada::Alignment alignment_type;
  
  alignment_type alignment;
 
  score_ptr_type operator()(const weight_set_type& weights,
			    const hypergraph_type& graph,
			    const scorer_ptr_type& scorer)
  {
    weight_type weight;
    
    cicada::viterbi(graph,
		    alignment,
		    weight,
		    cicada::operation::alignment_traversal(),
		    cicada::operation::weight_function<weight_type >(weights));

    std::ostringstream os;
    os << alignment;
    
    return scorer->score(sentence_type(os.str()));
  }
 
  score_ptr_type operator()(const hypergraph_type& graph,
			    const scorer_ptr_type& scorer,
			    const feature_type& feature,
			    const double& scale)
  {
    weight_type weight;
    
    cicada::viterbi(graph,
		    alignment,
		    weight,
		    cicada::operation::alignment_traversal(),
		    cicada::operation::single_scaled_function<weight_type >(feature, scale));
    
    std::ostringstream os;
    os << alignment;
    
    return scorer->score(sentence_type(os.str()));
  }
};

struct ViterbiDependency
{
  typedef scorer_document_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_document_type::score_ptr_type  score_ptr_type;

  typedef cicada::semiring::Logprob<double>  weight_type;
  typedef cicada::Dependency dependency_type;

  dependency_type dependency;

  score_ptr_type operator()(const weight_set_type& weights,
			    const hypergraph_type& graph,
			    const scorer_ptr_type& scorer)
  {
    weight_type weight;
    
    cicada::viterbi(graph,
		    dependency,
		    weight,
		    cicada::operation::dependency_traversal(),
		    cicada::operation::weight_function<weight_type >(weights));
    
    std::ostringstream os;
    os << dependency;
    
    return scorer->score(sentence_type(os.str()));
    
  }
  
  score_ptr_type operator()(const hypergraph_type& graph,
			    const scorer_ptr_type& scorer,
			    const feature_type& feature,
			    const double& scale)
  {
    typedef cicada::semiring::Logprob<double>  weight_type;
    
    weight_type weight;
    
    cicada::viterbi(graph,
		    dependency,
		    weight,
		    cicada::operation::dependency_traversal(),
		    cicada::operation::single_scaled_function<weight_type >(feature, scale));
    
    std::ostringstream os;
    os << dependency;
    
    return scorer->score(sentence_type(os.str()));
  }
};

template <typename Viterbi>
struct OracleForest
{
  OracleForest(const int cube_size) : cube_size_(cube_size) {}
  
  Viterbi __viterbi;
  const int cube_size_;

  template <typename Generator>
  std::pair<score_ptr_type, score_ptr_type>
  operator()(const weight_set_type& weights, 
	     const hypergraph_document_type& forests,
	     const hypergraph_document_type& oracles,
	     const scorer_document_type& scorers,
	     Generator& generator)
  {
    score_ptr_type score_1best;
    score_ptr_type score_oracle;
    
    for (size_t id = 0; id != forests.size(); ++ id)
      if (forests[id].is_valid()) {
	if (score_1best)
	  *score_1best += *__viterbi(weights, forests[id], scorers[id]);
	else
	  score_1best = __viterbi(weights, forests[id], scorers[id]);
      }

    for (size_t id = 0; id != oracles.size(); ++ id)
      if (oracles[id].is_valid()) {
	if (score_oracle)
	  *score_oracle += *__viterbi(weights, oracles[id], scorers[id]);
	else
	  score_oracle = __viterbi(weights, oracles[id], scorers[id]);
      }

    return std::make_pair(score_1best, score_oracle);
  }
  
  template <typename Generator>
  std::pair<score_ptr_type, score_ptr_type>
  operator()(const weight_set_type& weights, 
	     const hypergraph_document_type& forests,
	     const scorer_document_type& scorers,
	     const function_document_type& functions,
	     hypergraph_document_type& oracles,
	     Generator& generator)
  {
    typedef std::vector<size_t, std::allocator<size_t> > id_set_type;
    
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);
    
    score_ptr_type score_1best;
    
    id_set_type ids;
    boost::random_number_generator<boost::mt19937> gen(generator);
    
    weight_set_type::feature_type feature_scorer;
    for (size_t id = 0; id != functions.size(); ++ id)
      if (functions[id] && forests[id].is_valid()) {
	feature_scorer = functions[id]->feature_name();
	ids.push_back(id);
	
	if (score_1best)
	  *score_1best += *__viterbi(weights, forests[id], scorers[id]);
	else
	  score_1best = __viterbi(weights, forests[id], scorers[id]);
      }
    
    score_ptr_type score_best;
    score_ptr_type score_curr;
    score_ptr_type score_next;
    
    score_ptr_set_type scores_best(forests.size());
    score_ptr_set_type scores_curr(forests.size());
    score_ptr_set_type scores_next(forests.size());
    
    hypergraph_document_type oracles_best(forests.size());
    hypergraph_document_type oracles_curr(forests.size());
    hypergraph_document_type oracles_next(forests.size());
    
    double objective_best = - std::numeric_limits<double>::infinity();
    double objective_curr = - std::numeric_limits<double>::infinity();
    double objective_next = - std::numeric_limits<double>::infinity();

    hypergraph_type pruned;
    
    // 10 iterations will be fine...
    for (int i = 0; i < 10; ++ i) {
      typedef cicada::semiring::Tropical<double> weight_type;
      
      for (id_set_type::const_iterator iiter = ids.begin(); iiter != ids.end(); ++ iiter) {
	const size_t id = *iiter;
	
	score_ptr_type score_removed = (score_next ? score_next->clone() : score_ptr_type());
	if (score_removed && scores_curr[id])
	  *score_removed -= *scores_curr[id];
	
	dynamic_cast<cicada::feature::Scorer*>(functions[id].get())->assign(score_removed);
	
	model_type model(functions[id]);
	
	if (cube_size_ <= 0)
	  cicada::apply_exact(model, forests[id], oracles_next[id]);
	else
	  cicada::apply_cube_prune(model, forests[id], oracles_next[id], cicada::operation::single_scaled_function<weight_type>(feature_scorer, score_factor), cube_size_);
	
	// compute pruning...
	cicada::prune_beam(oracles_next[id],
			   pruned,
			   cicada::operation::single_scaled_function<weight_type >(feature_scorer, score_factor),
			   scorer_beam);
	oracles_next[id].swap(pruned);
	
	scores_next[id] = __viterbi(oracles_next[id], scorers[id], feature_scorer, score_factor);
	
	score_ptr_type score_sample = score_removed;
	if (score_sample)
	  *score_sample += *scores_next[id];
	else
	  score_sample = scores_next[id]->clone();
	
	const double objective_sample = score_sample->score() * score_factor;
	
	if (objective_sample > objective_next || ! scores_curr[id]) {
	  objective_next = objective_sample;
	  score_next = score_sample;
	} else {
	  oracles_next[id].swap(oracles_curr[id]);
	  scores_next[id].swap(scores_curr[id]);
	}
      }

      if (objective_next > objective_best) {
	score_best = score_next->clone();
	for (size_t id = 0; id != scores_next.size(); ++ id)
	  if (scores_next[id])
	    scores_best[id] = scores_next[id]->clone();
	
	oracles_best   = oracles_next;
	objective_best = objective_next;
      }
      
      if (objective_next <= objective_curr) break;
      
      score_curr = score_next->clone();
      for (size_t id = 0; id != scores_next.size(); ++ id)
	if (scores_next[id])
	  scores_curr[id] = scores_next[id]->clone();
      
      oracles_curr   = oracles_next;
      objective_curr = objective_next;
      
      std::random_shuffle(ids.begin(), ids.end(), gen);
    }
    
    oracles.swap(oracles_best);
    
    // remove scorer feature...
    for (size_t id = 0; id != oracles.size(); ++ id) 
      if (oracles[id].is_valid()) {
	hypergraph_type::edge_set_type::iterator eiter_end = oracles[id].edges.end();
	for (hypergraph_type::edge_set_type::iterator eiter = oracles[id].edges.begin(); eiter != eiter_end; ++ eiter)
	  eiter->features.erase(feature_scorer);
      }
    
    return std::make_pair(score_1best, score_best);
  }
};

inline
void read_refset(const path_type& refset_path,
		 scorer_document_type& scorers,
		 function_document_type& functions,
		 const size_t shard_rank=0,
		 const size_t shard_size=1)
{
  typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;
  
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  sentence_document_type document;
  
  parser_type parser;
  id_sentence_type id_sentence;

  utils::compress_istream is(refset_path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iter_type iter(is);
  iter_type iter_end;
  
  while (iter != iter_end) {
    id_sentence.second.clear();
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
      if (iter != iter_end)
	throw std::runtime_error("refset parsing failed");
    
    const size_t id = id_sentence.first;
    
    if (shard_size && (id % shard_size != shard_rank)) continue;
    
    const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
    
    if (id_rank >= scorers.size())
      scorers.resize(id_rank + 1);
    
    if (id_rank >= document.size())
      document.resize(id_rank + 1);
    
    if (! scorers[id_rank])
      scorers[id_rank] = scorers.create();
    
    scorers[id_rank]->insert(id_sentence.second);
    document[id_rank].push_back(id_sentence.second);
  }

  static const hypergraph_type __hypergraph;
  static const lattice_type __lattice;
  static const span_set_type __spans;
  static const ngram_count_set_type __ngram_counts;
  
  functions.reserve(document.size());
  functions.resize(document.size());
  for (size_t seg = 0; seg != document.size(); ++ seg) {
    if (document[seg].empty())
      throw std::runtime_error("no reference at segment: " + utils::lexical_cast<std::string>(seg));
    
    functions[seg] = feature_function_type::create(scorers.parameter());
    
    if (! dynamic_cast<const cicada::feature::Scorer*>(functions[seg].get()))
      throw std::runtime_error("this is not a scorer feature");
    
    functions[seg]->assign((seg / shard_size) + shard_rank, __hypergraph, __lattice, __spans, document[seg], __ngram_counts);
  }
}

class Event
{
public:
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
public:
  Event() : buffer(), size(0) {}
  Event(const Event& x) : buffer(x.buffer), size(x.size) {}
  Event(const std::string& data) { encode(data); }
  Event& operator=(const std::string& data)
  {
    encode(data);
    return *this;
  }
  
  Event& operator=(const Event& x)
  {
    buffer = x.buffer;
    size   = x.size;
    return *this;
  }

  operator std::string() const { return decode(); }

  bool empty() const { return buffer.empty(); }
  void swap(Event& x)
  {
    buffer.swap(x.buffer);
    std::swap(size, x.size);
  }
  
  void encode(const std::string& data)
  {
    buffer.clear();
    buffer.reserve(data.size() >> 2);
    
    boost::iostreams::filtering_ostream os;
    os.push(codec::lz4_compressor());
    os.push(boost::iostreams::back_inserter(buffer));
    os.write(data.c_str(), data.size());
    os.reset();
    
    buffer_type(buffer).swap(buffer);
    
    size = data.size();
  }
  
  std::string decode() const
  {
    buffer_type output(size);
    
    boost::iostreams::filtering_istream is;
    is.push(codec::lz4_decompressor());
    is.push(boost::iostreams::array_source(&(*buffer.begin()), buffer.size()));
    is.read(&(*output.begin()), size);
    
    return std::string(output.begin(), output.end());
  }
  
private:
  typedef std::vector<char, std::allocator<char> > buffer_type;

private:
  buffer_type buffer;
  size_type   size;
};

namespace std
{
  inline
  void swap(Event& x, Event& y)
  {
    x.swap(y);
  }
};

typedef Event event_type;
typedef std::vector<event_type, std::allocator<event_type> > event_set_type;

inline
void read_events(const path_type& input_path,
		 event_set_type& events,
		 const bool directory_mode,
		 const bool id_mode,
		 const size_t shard_rank=0,
		 const size_t shard_size=1)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;

  if (directory_mode) {
    if (! boost::filesystem::is_directory(input_path))
      throw std::runtime_error("input is not directory! " + input_path.string());

    boost::spirit::qi::uint_parser<size_t, 10, 1, -1> id_parser;
    
    size_t id;
    std::string line;
    for (size_t i = 0; /**/; ++ i)
      if (shard_size <= 0 || i % shard_size == shard_rank) {
	const path_type path = input_path / (utils::lexical_cast<std::string>(i) + ".gz");
	
	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	utils::getline(is, line);
	
	if (line.empty()) continue;
	
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (id != i)
	  throw std::runtime_error("id doest not match!");
	
	const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	
	if (id_rank >= events.size())
	  events.resize(id_rank + 1);
	
	events[id_rank] = line;
      }
  } else if (id_mode) {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    boost::spirit::qi::uint_parser<size_t, 10, 1, -1> id_parser;
    
    size_t id;
    std::string line;
    while (utils::getline(is, line)) 
      if (! line.empty()) {
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (shard_size == 0 || id % shard_size == shard_rank) {
	  const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	  
	  if (id_rank >= events.size())
	    events.resize(id_rank + 1);
	  
	  events[id_rank] = line;
	}
      }
  } else {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    std::string line;
    for (size_t id = 0; utils::getline(is, line); ++ id) 
      if (shard_size == 0 || id % shard_size == shard_rank) {
	if (! line.empty())
	  events.push_back(utils::lexical_cast<std::string>(id) + " ||| " + line);
	else
	  events.push_back(std::string());
      }
  }
}


inline
path_type add_suffix(const path_type& path, const std::string& suffix)
{
  bool has_suffix_gz  = false;
  bool has_suffix_bz2 = false;
  
  path_type path_added = path;
  
  if (path.extension() == ".gz") {
    path_added = path.parent_path() / path.stem();
    has_suffix_gz = true;
  } else if (path.extension() == ".bz2") {
    path_added = path.parent_path() / path.stem();
    has_suffix_bz2 = true;
  }
  
  path_added = path_added.string() + suffix;
  
  if (has_suffix_gz)
    path_added = path_added.string() + ".gz";
  else if (has_suffix_bz2)
    path_added = path_added.string() + ".bz2";
  
  return path_added;
}

#endif
