//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_BLOCK_IMPL__HPP__
#define __CICADA_LEARN_BLOCK_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>

#include <sstream>
#include <vector>

#include "cicada_kbest_impl.hpp"
#include "cicada_mert_kbest_impl.hpp"

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
#include "cicada/viterbi.hpp"
#include "cicada/expected_ngram.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/sgi_hash_set.hpp"
#include "utils/base64.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"
#include "utils/config.hpp"
#include "utils/mathop.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_trie.hpp"


#include <boost/tokenizer.hpp>

#ifdef HAVE_SNAPPY
#include <snappy.h>
#endif

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

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


struct LearnXBLEU : public LearnBase
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
  
  typedef cicada::FeatureVectorUnordered<weight_type, std::allocator<weight_type> > gradient_type;
  typedef std::vector<gradient_type, std::allocator<gradient_type> >                gradients_type;
  
  typedef cicada::FeatureVectorUnordered<double, std::allocator<double> > expectation_type;
  
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
		  const weight_set_type& __weights)
      : ngrams(__ngrams), counts(__counts), ids(__ids),
	weights(__weights) {}
      
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
    typedef cicada::FeatureVectorUnordered<weight_type, std::allocator<weight_type> > accumulated_type;
    typedef std::vector<accumulated_type, std::allocator<accumulated_type> > accumulated_set_type;

    struct value_type
    {
      value_type& operator+=(const bleu_gradient_function::value_type& x)
      {
	const double margin = cicada::dot_product(x.edge.features, impl.weights);
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(margin * scale);
	
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
		const weight_type value(fiter->second * scale);
		  
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
		       const weight_set_type& __weights) 
      : ngrams(__ngrams), counts(__counts), ids(__ids),
	matched(__matched), hypo(__hypo),
	weights(__weights), 
	dM(order + 1),
	dH(order + 1) {}
    
    const ngram_set_type&  ngrams;
    const count_set_type&  counts;
    const id_map_type&     ids;
      
    const weights_type&    matched;
    const weights_type&    hypo;
      
    const weight_set_type& weights;
      
    accumulated_set_type dM;
    accumulated_set_type dH;
  };
  
  LearnXBLEU() { clear(); }

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
    
    std::fill(counts_matched.begin(), counts_matched.end(), weight_type());
    std::fill(counts_hypo.begin(), counts_hypo.end(), weight_type());
    counts_reference = weight_type();
    
    for (int n = 1; n <= order; ++ n) {
      gradients_matched[n].clear();
      gradients_hypo[n].clear();
    }
  }
  
  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer)
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
			   cicada::operation::weight_scaled_function<weight_type>(weights, scale),
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
				     weights);
    
    cicada::inside_outside(forest,
			   bleu_inside,
			   bleu_gradient,
			   bleu_function(ngrams, counts, ids, weights),
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
  }
  
  double encode(expectation_type& g)
  {
    g.clear();
    
    // smoothing...
    double smoothing = 1e-40;
    for (int n = 1; n <= order; ++ n) {
      if (counts_hypo[n] > weight_type() && counts_matched[n] <= weight_type())
	counts_matched[n] = smoothing;
      smoothing *= 0.1;
    }
    
    // compute P
    double P = 0.0;
    for (int n = 1; n <= order; ++ n)
      if (counts_hypo[n] > weight_type())
	P += (1.0 / order) * (cicada::semiring::log(counts_matched[n]) - cicada::semiring::log(counts_hypo[n]));
    
    // compute C and B
    const double C = counts_reference / counts_hypo[1];
    const double B = brevity_penalty(1.0 - C);
    
    // for computing g...
    const double exp_P = utils::mathop::exp(P);
    const double C_dC  = C * derivative_brevity_penalty(1.0 - C);

    const double objective = exp_P * B;
    
    // we will collect minus gradient for minimizing negative-xBLEU
    
    for (int n = 1; n <= order; ++ n) 
      if (counts_hypo[n] > weight_type()) {
	const double factor_matched = - (exp_P * B / order) / counts_matched[n];
	const double factor_hypo    = - (exp_P * B / order) / counts_hypo[n];
	
	gradient_type::const_iterator miter_end = gradients_matched[n].end();
	for (gradient_type::const_iterator miter = gradients_matched[n].begin(); miter != miter_end; ++ miter)
	  g[miter->first] += factor_matched * miter->second;
	
	gradient_type::const_iterator hiter_end = gradients_hypo[n].end();
	for (gradient_type::const_iterator hiter = gradients_hypo[n].begin(); hiter != hiter_end; ++ hiter)
	  g[hiter->first] -= factor_hypo * hiter->second;
      }
    
    if (counts_hypo[1] > weight_type()) {
      const double factor_hypo = - (exp_P * C_dC) / counts_hypo[1];
      
      gradient_type::const_iterator hiter_end = gradients_hypo[1].end();
      for (gradient_type::const_iterator hiter = gradients_hypo[1].begin(); hiter != hiter_end; ++ hiter)
	g[hiter->first] += factor_hypo * hiter->second;
    }
    
    return objective;
  }
  
  // required for learn()
  weights_type counts_matched;
  weights_type counts_hypo;
  weight_type  counts_reference;
  
  gradients_type gradients_matched;
  gradients_type gradients_hypo;
  
  // local variables
  index_set_type index;
  ngram_set_type ngrams;
  count_set_type counts;
  id_map_type    ids;
  
  bleu_weights_type bleu_inside;
  weights_type      matched;
  weights_type      hypo;
};

struct LearnXBLEUL2 : public LearnXBLEU
{
  LearnXBLEUL2(const size_type __instances) : instances(__instances), epoch(0), lambda(C), weight_scale(1.0), weight_norm(0.0) {}
  
  void initialize(weight_set_type& weights)
  {
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  
  void finalize(weight_set_type& weights)
  {
    weights *= weight_scale;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  
  double learn(weight_set_type& weights)
  {
    if (counts_reference <= weight_type()) {
      clear();
      return 0.0;
    }
    
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + batch_size - 1) / batch_size;
    const double eta = eta0 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    rescale(weights, 1.0 - eta * lambda);
    
    // compute gradient...
    const double objective = LearnXBLEU::encode(g);
    
    double a_norm = 0.0;
    double pred = 0.0;
    expectation_type::const_iterator giter_end = g.end();
    for (expectation_type::const_iterator giter = g.begin(); giter != giter_end; ++ giter) {
      // we will update "minus" value...
      
      double& x = weights[giter->first];
      const double alpha = - static_cast<double>(giter->second) * eta;
      
      a_norm += alpha * alpha;
      pred += 2.0 * x * alpha;
      
      //weight_norm += 2.0 * x * alpha * weight_scale + alpha * alpha;
      x += alpha / weight_scale;
    }
    
    // avoid numerical instability...
    weight_norm += a_norm + pred * weight_scale;
    
    if (weight_norm > 1.0 / lambda || project_weight)
      rescale(weights, std::sqrt(1.0 / (lambda * weight_norm)));
    
    if (weight_scale < 0.001 || 1000 < weight_scale)
      finalize(weights);
    
    clear();
    
    return objective;
  }
  
  void rescale(weight_set_type& weights, const double scaling)
  {
    weight_norm *= scaling * scaling;
    if (scaling != 0.0)
      weight_scale *= scaling;
    else {
      weight_scale = 1.0;
      std::fill(weights.begin(), weights.end(), 0.0);
    }
  }

  expectation_type g;
  
  size_type instances;
  
  size_type epoch;
  double    lambda;
  
  double weight_scale;
  double weight_norm;
};

struct LearnXBLEUL1 : public LearnXBLEU
{
  typedef cicada::WeightVector<double> penalty_set_type;
  
  LearnXBLEUL1(const size_type __instances) : instances(__instances), epoch(0), lambda(C), penalties(), penalty(0.0) {}

  void initialize(weight_set_type& weights)
  {

  }

  void finalize(weight_set_type& weights)
  {
    
  }

  double learn(weight_set_type& weights)
  {
    if (counts_reference <= weight_type()) {
      clear();
      return 0.0;
    }
    
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + batch_size - 1) / batch_size;
    const double eta = eta0 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    penalty += eta * lambda;
    
    // compute gradient...
    const double objective = LearnXBLEU::encode(g);
    
    expectation_type::const_iterator giter_end = g.end();
    for (expectation_type::const_iterator giter = g.begin(); giter != giter_end; ++ giter) {
      // we will update "minus" value...
      
      double& x = weights[giter->first];
      x += - static_cast<double>(giter->second) * eta;
      
      // apply penalty
      apply(x, penalties[giter->first], penalty);
    }
    
    clear();
    
    return objective;
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
  
  expectation_type g;
  
  size_type instances;
  
  size_type epoch;
  double    lambda;
  
  penalty_set_type penalties;
  double penalty;
};


// logistic regression base...
struct LearnLR : public LearnBase
{
  typedef cicada::semiring::Log<double> weight_type;
  typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
  typedef cicada::FeatureVectorUnordered<weight_type, std::allocator<weight_type> > gradient_type;
  typedef cicada::FeatureVectorUnordered<double, std::allocator<double> > expectation_type;
  
  LearnLR() { clear(); }
  
  struct weight_function
  {
    typedef weight_type value_type;

    weight_function(const weight_set_type& __weights) : weights(__weights) {}
      
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e
      return cicada::semiring::traits<value_type>::exp(cicada::dot_product(edge.features, weights));
    }
      
    const weight_set_type& weights;
  };
    
  struct feature_function
  {
    struct value_type
    {
      value_type(const feature_set_type& __features,
		 const weight_set_type& __weights)
	: features(__features), weights(__weights) {}
	
      friend
      value_type operator*(value_type x, const weight_type& weight)
      {
	x.inside_outside = weight;
	return x;
      }
	
      weight_type inside_outside;
      const feature_set_type& features;
      const weight_set_type&  weights;
    };
      
    feature_function(const weight_set_type& __weights) : weights(__weights) {}
      
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      return value_type(edge.features, weights);
    }
      
    const weight_set_type& weights;
  };
    
  struct gradients_type
  {
    struct value_type
    {
      value_type(gradient_type& __gradient) : gradient(__gradient) {}
	
      value_type& operator+=(const feature_function::value_type& x)
      {
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(x.features, x.weights)) * x.inside_outside;
	  
	feature_set_type::const_iterator fiter_end = x.features.end();
	for (feature_set_type::const_iterator fiter = x.features.begin(); fiter != fiter_end; ++ fiter)
	  gradient[fiter->first] += weight_type(fiter->second) * weight;
	  
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
  
  void encode(const size_type id, const weight_set_type& weights, const hypergraph_type& forest, const hypergraph_type& oracle, const scorer_ptr_type& scorer)
  {
    gradients_forest.clear();
    gradients_oracle.clear();
    
    inside_forest.clear();
    inside_oracle.clear();
    
    inside_forest.resize(forest.nodes.size());
    inside_oracle.resize(oracle.nodes.size());
    
    cicada::inside_outside(forest, inside_forest, gradients_forest, weight_function(weights), feature_function(weights));
    cicada::inside_outside(oracle, inside_oracle, gradients_oracle, weight_function(weights), feature_function(weights));
    
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

// SGDL2 learner
struct LearnSGDL2 : public LearnLR
{
  LearnSGDL2(const size_type __instances) : instances(__instances), epoch(0), lambda(C), weight_scale(1.0), weight_norm(0.0) {}
  
  
  void initialize(weight_set_type& weights)
  {
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }
  
  void finalize(weight_set_type& weights)
  {
    weights *= weight_scale;
    
    weight_scale = 1.0;
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
  }

  double learn(weight_set_type& weights)
  {
    if (! samples) {
      clear();
      return 0.0;
    }
    
    const size_type k = samples;
    const double k_norm = 1.0 / k;
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + batch_size - 1) / batch_size;
    const double eta = eta0 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;

    const double objective_normalized = objective * k_norm;
    
    // do we really need this...?
    rescale(weights, 1.0 - eta * lambda);
    
    double a_norm = 0.0;
    double pred = 0.0;
    gradient_type::const_iterator giter_end = gradient.end();
    for (gradient_type::const_iterator giter = gradient.begin(); giter != giter_end; ++ giter) {
      // we will update "minus" value...
      
      double& x = weights[giter->first];
      const double alpha = - static_cast<double>(giter->second) * eta * k_norm;
      
      a_norm += alpha * alpha;
      pred += 2.0 * x * alpha;
      
      //weight_norm += 2.0 * x * alpha * weight_scale + alpha * alpha;
      x += alpha / weight_scale;
    }
    
    // avoid numerical instability...
    weight_norm += a_norm + pred * weight_scale;
    
    if (weight_norm > 1.0 / lambda || project_weight)
      rescale(weights, std::sqrt(1.0 / (lambda * weight_norm)));
    
    if (weight_scale < 0.001 || 1000 < weight_scale)
      finalize(weights);
    
    clear();
    
    return objective_normalized;
  }

  void rescale(weight_set_type& weights, const double scaling)
  {
    weight_norm *= scaling * scaling;
    if (scaling != 0.0)
      weight_scale *= scaling;
    else {
      weight_scale = 1.0;
      std::fill(weights.begin(), weights.end(), 0.0);
    }
  }

  size_type instances;
  
  size_type epoch;
  double    lambda;
  
  double weight_scale;
  double weight_norm;
};


// SGDL1 learner
struct LearnSGDL1 : public LearnLR
{
  typedef cicada::WeightVector<double> penalty_set_type;
  
  // maximize...
  // L_w =  \sum \log p(y | x) - C |w|
  // 
  
  LearnSGDL1(const size_type __instances) : instances(__instances), epoch(0), lambda(C), penalties(), penalty(0.0) {}
  
  
  void initialize(weight_set_type& weights)
  {

  }

  void finalize(weight_set_type& weights)
  {
    
  }
  
  double learn(weight_set_type& weights)
  {
    if (! samples) {
      clear();
      return 0.0;
    }
    
    const size_type k = samples;
    const double k_norm = 1.0 / k;
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + batch_size - 1) / batch_size;
    const double eta = eta0 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    penalty += eta * lambda * k_norm;

    const double objective_normalized = objective * k_norm;
    
    gradient_type::const_iterator giter_end = gradient.end();
    for (gradient_type::const_iterator giter = gradient.begin(); giter != giter_end; ++ giter) {
      double& x = weights[giter->first];
      
      // update weight ... we will update "minus" value
      x += - static_cast<double>(giter->second) * eta * k_norm;
      
      // apply penalty
      apply(x, penalties[giter->first], penalty);
    }
    
    clear();
    
    return objective_normalized;
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

  size_type instances;
  
  size_type epoch;
  double    lambda;
  
  penalty_set_type penalties;
  double penalty;
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
  Viterbi __viterbi;

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
	
	cicada::apply_exact(model, forests[id], oracles_next[id]);
	
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
  Event() {}
  Event(const Event& x) : buffer(x.buffer) {}
  Event(const std::string& data) { encode(data); }
  Event& operator=(const std::string& data)
  {
    encode(data);
    return *this;
  }
  
  Event& operator=(const Event& x)
  {
    buffer = x.buffer;
    return *this;
  }

  operator std::string() const { return decode(); }

  bool empty() const { return buffer.empty(); }
  void swap(Event& x) { buffer.swap(x.buffer); }
  
  void encode(const std::string& data)
  {
#ifdef HAVE_SNAPPY
    const size_t max_length = snappy::MaxCompressedLength(data.size());
    buffer.reserve(max_length);
    buffer.resize(max_length);
    size_t length = 0;
    snappy::RawCompress(data.c_str(), data.size(), &(*buffer.begin()), &length);
    buffer.resize(length);
    buffer_type(buffer).swap(buffer);
#else
    buffer.clear();
    
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(boost::iostreams::back_insert_device<buffer_type>(buffer));
    os.write(data.c_str(), data.size());
    os.reset();
#endif
  }
  
  std::string decode() const
  {
#ifdef HAVE_SNAPPY
    std::string output;
    snappy::Uncompress(&(*buffer.begin()), buffer.size(), &output);
    return output;
#else
    std::string output;
    
    boost::iostreams::filtering_istream is;
    is.push(boost::iostreams::zlib_decompressor());
    is.push(boost::iostreams::array_source(&(*buffer.begin()), buffer.size()));

    char buf[1024];
    
    do {
      is.read(buf, 1024);
      std::copy(buf, buf + is.gcount(), std::back_inserter(output));
    } while (is);
    
    is.reset();

    return output;
#endif
  }
  
private:
  typedef std::vector<char, std::allocator<char> > buffer_type;

private:
  buffer_type buffer;
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
	std::getline(is, line);
	
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
    while (std::getline(is, line)) 
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
    for (size_t id = 0; std::getline(is, line); ++ id) 
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
