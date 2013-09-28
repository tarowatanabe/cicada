//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_ONLINE_KBEST_IMPL__HPP__
#define __CICADA_LEARN_ONLINE_KBEST_IMPL__HPP__ 1

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

#include "utils/unordered_set.hpp"
#include "utils/base64.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"
#include "utils/config.hpp"
#include "utils/mathop.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/getline.hpp"

#include "cicada_learn_online_regularize_impl.hpp"
#include "cicada_learn_online_rate_impl.hpp"

#include <boost/tokenizer.hpp>

#include <codec/lz4.hpp>

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

struct LearnBase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef hypothesis_type::feature_value_type feature_value_type;
  
  struct SampleSet
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    typedef std::vector<size_type, std::allocator<size_type> > offsets_type;

    struct Sample
    {
      typedef features_type::const_iterator const_iterator;

      Sample(const_iterator __first, const_iterator __last) : first(__first), last(__last) {}

      const_iterator begin() const { return first; }
      const_iterator end() const { return last; }
      size_type size() const { return last - first; }
      bool emtpy() const { return first == last; }
      
      const_iterator first;
      const_iterator last;
    };

    typedef Sample sample_type;
    typedef sample_type value_type;
    
    SampleSet() : features(), offsets() { offsets.push_back(0); }
    
    void clear()
    {
      features.clear();
      offsets.clear();
      offsets.push_back(0);
    }
    
    template <typename Iterator>
    void insert(Iterator first, Iterator last)
    {
      features.insert(features.end(), first, last);
      offsets.push_back(features.size());
    }
    
    sample_type operator[](size_type pos) const
    {
      return sample_type(features.begin() + offsets[pos], features.begin() + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }

    void swap(SampleSet& x)
    {
      features.swap(x.features);
      offsets.swap(x.offsets);
    }

    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    features_type features;
    offsets_type  offsets;
  };
  
  typedef SampleSet sample_set_type;
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
  
  typedef std::vector<double, std::allocator<double> > margins_type;

  typedef std::deque<sample_set_type, std::allocator<sample_set_type> >       sample_map_type;
  typedef std::deque<score_ptr_set_type, std::allocator<score_ptr_set_type> > score_ptr_map_type;

  LearnXBLEUBase() { clear(); }

  void clear()
  {
    features.clear();
    bleus.clear();
    
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
    counts_entropy   = weight_type();
    
    for (int n = 1; n <= order; ++ n) {
      gradients_matched[n].clear();
      gradients_hypo[n].clear();
    }
    
    gradients_reference.clear();
    gradients_entropy.clear();
  }
  
  template <typename Iterator>
  boost::fusion::tuple<double, double, double> gradient(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return boost::fusion::tuple<double, double, double>(0.0, 0.0, 0.0);
  }
  
  // we need features + bleu stats...
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles)
  {
    features.push_back(sample_set_type());
    bleus.push_back(score_ptr_set_type());
    
    for (size_type k = 0; k != kbests.size(); ++ k) {
      features.back().insert(kbests[k].features.begin(), kbests[k].features.end());
      bleus.back().push_back(kbests[k].score);
    }
  }


  void encode(const sample_set_type& features,
	      const score_ptr_set_type& bleus,
	      const weight_set_type& weights,
	      const double& weight_scale)
  {
    if (features.size() != bleus.size())
      throw std::runtime_error("# of training data do not match");

    // we will collect gradients from kbests
    weight_type Z;
    weight_type Z_reference;
    weight_type Z_entropy;
    weight_type dR;

    margins.clear();
    hypo.clear();
    matched.clear();
    expectation.clear();
    
    hypo.resize(order + 1);
    matched.resize(order + 1);
    
    // first pass... compute margin and Z
    for (size_type k = 0; k != features.size(); ++ k) {
      const double margin  = cicada::dot_product(weights, features[k].begin(), features[k].end(), 0.0) * weight_scale;
      
      margins.push_back(margin);
      Z += cicada::semiring::traits<weight_type>::exp(margin * scale);
    }
    
    // second pass.. compute sums
    for (size_type k = 0; k != features.size(); ++ k) {
      const double& margin = margins[k];
      const weight_type prob = cicada::semiring::traits<weight_type>::exp(margin * scale) / Z;
      
      const cicada::eval::Bleu* bleu = dynamic_cast<const cicada::eval::Bleu*>(bleus[k].get());
      if (! bleu)
	throw std::runtime_error("no bleu statistics?");
      
      // collect scaled bleu stats
      for (size_t n = 1; n <= static_cast<size_t>(order); ++ n) {
	if (n - 1 < bleu->ngrams_hypothesis.size())
	  hypo[n] += prob * bleu->ngrams_hypothesis[n - 1];
	if (n - 1 < bleu->ngrams_matched.size())
	  matched[n] += prob * bleu->ngrams_matched[n - 1];
      }
      
      // collect reference length
      Z_reference += prob * bleu->length_reference;
      
      // collect entropy...
      Z_entropy -= prob * cicada::semiring::log(prob);
      
      // collect expectation
      sample_set_type::value_type::const_iterator fiter_end = features[k].end();
      for (sample_set_type::value_type::const_iterator fiter = features[k].begin(); fiter != fiter_end; ++ fiter)
	expectation[fiter->first] += prob * weight_type(fiter->second * scale);
      
      dR += weight_type(1.0 + cicada::semiring::log(prob)) * prob;
    }
    
    // accumulate
    std::transform(hypo.begin(), hypo.end(), counts_hypo.begin(), counts_hypo.begin(), std::plus<weight_type>());
    std::transform(matched.begin(), matched.end(), counts_matched.begin(), counts_matched.begin(), std::plus<weight_type>());
    
    counts_reference += Z_reference;
    counts_entropy   += Z_entropy;
    
    gradient_type::const_iterator eiter_end = expectation.end();
    for (gradient_type::const_iterator eiter = expectation.begin(); eiter != eiter_end; ++ eiter) {
      // collect bleus...
      for (int n = 1; n <= order; ++ n) {
	gradients_hypo[n][eiter->first] -= eiter->second * hypo[n];
	gradients_matched[n][eiter->first] -= eiter->second * matched[n];
      }
      
      // reference lengths
      gradients_reference[eiter->first] -= eiter->second * Z_reference;
      
      // entropy gradient...
      gradients_entropy[eiter->first] -= - dR * eiter->second;
    }
    
    // third pass, collect gradients...
    for (size_type k = 0; k != features.size(); ++ k) {
      const double& margin = margins[k];
      const weight_type prob = cicada::semiring::traits<weight_type>::exp(margin * scale) / Z;
      const cicada::eval::Bleu* bleu = dynamic_cast<const cicada::eval::Bleu*>(bleus[k].get());
      
      sample_set_type::value_type::const_iterator fiter_end = features[k].end();
      for (sample_set_type::value_type::const_iterator fiter = features[k].begin(); fiter != fiter_end; ++ fiter) {
	const weight_type value(fiter->second * scale);
	
	// bleu statistics
	for (size_t n = 1; n <= static_cast<size_t>(order); ++ n) {
	  weight_type& grad_hypo    = gradients_hypo[n][fiter->first];
	  weight_type& grad_matched = gradients_matched[n][fiter->first];
	  
	  if (n - 1 < bleu->ngrams_hypothesis.size())
	    grad_hypo += value * prob * bleu->ngrams_hypothesis[n - 1];
	  
	  if (n - 1 < bleu->ngrams_matched.size())
	    grad_matched += value * prob * bleu->ngrams_matched[n - 1];
	}
	
	// reference lengths
	gradients_reference[fiter->first] += value * prob * bleu->length_reference;
	
	// entropy: we will collect minus values!
	gradients_entropy[fiter->first] += - weight_type(1.0 + cicada::semiring::log(prob)) * prob * value;
      }
    }
  }

  std::pair<double, bool> encode(gradient_xbleu_type& g)
  {
    g.clear();
    
    // smoothing...
    {
      double smoothing = 1e-40;
      for (int n = 1; n <= order; ++ n) {
	if (counts_hypo[n] > weight_type() && counts_matched[n] <= weight_type())
	  counts_matched[n] = std::min(weight_type(smoothing), counts_hypo[n]);
	
	smoothing *= 0.1;
      }
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
    
    // xBLEU...
    const weight_type objective_bleu = exp_P * B;
    const double      factor_instance = 1.0 / features.size();
    const weight_type entropy = counts_entropy * factor_instance;
    const weight_type factor_order = 1.0 / order;
    
    const double objective = - objective_bleu - temperature * entropy;
    
    // entropy...
    if (temperature != 0.0) {
      gradient_type::const_iterator eiter_end = gradients_entropy.end();
      for (gradient_type::const_iterator eiter = gradients_entropy.begin(); eiter != eiter_end; ++ eiter) {
	double& grad = g[eiter->first];
	grad = - temperature * factor_instance * eiter->second;

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
      const weight_type factor_ref  = - (exp_P * C_dC) / counts_reference;
      const weight_type factor_hypo = - (exp_P * C_dC) / counts_hypo[1];
      
      gradient_type::const_iterator riter_end = gradients_reference.end();
      for (gradient_type::const_iterator riter = gradients_reference.begin(); riter != riter_end; ++ riter) {
	double& grad = g[riter->first];
	grad -= factor_ref * riter->second;
	
	if (! std::isfinite(grad))
	  return std::make_pair(objective, false);
      }
      
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
  
  sample_map_type    features;
  score_ptr_map_type bleus;
  
  // required for learn()
  weights_type counts_matched;
  weights_type counts_hypo;
  weight_type  counts_reference;
  weight_type  counts_entropy;
  
  gradients_type gradients_matched;
  gradients_type gradients_hypo;
  gradient_type  gradients_reference;
  gradient_type  gradients_entropy;
  
  // local variables
  margins_type  margins;
  weights_type  matched;
  weights_type  hypo;
  gradient_type expectation;
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
  
  double learn(weight_set_type& weights)
  {
    if (features.empty()) {
      clear();
      return 0.0;
    }
    
    // collect gradient
    for (size_type k = 0; k != features.size(); ++ k)
      LearnXBLEUBase::encode(features[k], bleus[k], weights, regularizer.scale());

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
    }
    
    regularizer.postprocess(weights, eta);
    
    clear();
    
    return objective.first;
  }

  gradient_xbleu_type g;
  
  Regularize& regularizer;
  Rate&       rate;
};

struct LearnExpectedLoss : public LearnBase
{
  // lossfunction based on expected loss

  typedef std::vector<double, std::allocator<double> > margin_set_type;
  typedef std::vector<double, std::allocator<double> > loss_set_type;

  typedef std::deque<loss_set_type, std::allocator<loss_set_type> > loss_map_type;

  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::semiring::traits<weight_type> traits_type;
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > expectation_type;

  LearnExpectedLoss(Regularize& __regularizer,
		    Rate& __rate)
    : regularizer(__regularizer),
      rate(__rate) {}
  
  void clear()
  {
    features.clear();
    losses.clear();
  }
  
  template <typename Iterator>
  boost::fusion::tuple<double, double, double> gradient(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return boost::fusion::tuple<double, double, double>(0.0, 0.0, 0.0);
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles)
  {
    if (kbests.empty()) return;
    
    losses.resize(losses.size() + 1);
    losses.back().reserve(kbests.size());
    
    hypothesis_set_type::const_iterator kiter_end = kbests.end();
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
      features.insert(kiter->features.begin(), kiter->features.end());
      losses.back().push_back(kiter->loss);
    }
  }
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }

  double learn(weight_set_type& weights)
  {
    if (losses.empty() || features.empty()) {
      features.clear();
      losses.clear();
      
      return 0.0;
    }
    
    expectations.clear();

    weight_type objective;
    
    size_t pos = 0;
    for (size_t i = 0; i != losses.size(); ++ i) 
      if (! losses[i].empty()) {
	weight_type Z;
	expectations_Z.clear();
	margins.clear();
	
	const size_t pos_local = pos;
	for (size_t j = 0; j != losses[i].size(); ++ j, ++ pos) {
	  margins.push_back(cicada::dot_product(weights, features[pos].begin(), features[pos].end(), 0.0) * regularizer.scale() * scale);
	  Z += traits_type::exp(margins.back());
	}
	
	weight_type scaling_sum;
	weight_type objective_local;
	
	for (size_t j = 0, p = pos_local; j != losses[i].size(); ++ j, ++ p) {
	  const weight_type weight = traits_type::exp(margins[j]) / Z;
	  const weight_type scaling = weight_type(scale * losses[i][j]) * weight;
	  
	  scaling_sum += scaling;
	  objective_local += weight_type(losses[i][j]) * weight;
	  
	  sample_set_type::value_type::const_iterator fiter_end = features[p].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[p].begin(); fiter != fiter_end; ++ fiter) {
	    expectations[fiter->first] += weight_type(fiter->second) * scaling;
	    expectations_Z[fiter->first] += weight_type(fiter->second) * weight;
	  }
	}

	expectation_type::const_iterator eiter_end = expectations_Z.end();
	for (expectation_type::const_iterator eiter = expectations_Z.begin(); eiter != eiter_end; ++ eiter)
	  expectations[eiter->first] -= weight_type(eiter->second) * scaling_sum;
	
	objective += objective_local;
      }
    
    const size_type k = losses.size();
    const double k_norm = 1.0 / k;
    
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    // update by expectations...
    expectation_type::const_iterator eiter_end = expectations.end();
    for (expectation_type::const_iterator eiter = expectations.begin(); eiter != eiter_end; ++ eiter) {
      const double amount = static_cast<double>(eiter->second) * k_norm;

      regularizer.update(weights, eiter->first, amount, rate(eiter->first, amount));
    }
    
    regularizer.postprocess(weights, eta);
    
    features.clear();
    losses.clear();
    
    return objective * k_norm;
  }

  margin_set_type    margins;
  sample_set_type    features;
  expectation_type   expectations;
  expectation_type   expectations_Z;

  loss_map_type losses;
  
  Regularize& regularizer;
  Rate&       rate;
};

struct LearnOExpectedLoss : public LearnBase
{
  // lossfunction based on expected loss

  typedef std::vector<double, std::allocator<double> > margin_set_type;
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  
  typedef std::deque<loss_set_type, std::allocator<loss_set_type> > loss_map_type;
  
  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::semiring::traits<weight_type> traits_type;
  typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > expectation_type;

  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;
  
  struct HMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    HMatrix(const sample_set_type& __features) : features(__features) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[i].begin(), features[i].end(), features[j].begin(), features[j].end(), 0.0);
    }
    
    const sample_set_type& features;
  };
  
  struct MMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    MMatrix(const sample_set_type& __features) : features(__features) {}
    
    template <typename __W>
    void operator()(__W& w, const alpha_type& alpha) const
    {
      const size_type model_size = features.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename __W>
    double operator()(const __W& w, const size_t& i) const
    {
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename __W>
    void operator()(__W& w, const double& update, const size_t& i) const
    {
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const sample_set_type& features;
  };

  LearnOExpectedLoss(Regularize& __regularizer,
		     Rate& __rate)
    : tolerance(0.1),
      regularizer(__regularizer),
      rate(__rate) {}
  
  void clear()
  {
    features.clear();
    losses.clear();
  }
  
  template <typename Iterator>
  boost::fusion::tuple<double, double, double> gradient(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return boost::fusion::tuple<double, double, double>(0.0, 0.0, 0.0);
  }
  
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles)
  {
    if (kbests.empty()) return;
    
    losses.resize(losses.size() + 1);
    losses.back().reserve(kbests.size());
    
    hypothesis_set_type::const_iterator kiter_end = kbests.end();
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
      features.insert(kiter->features.begin(), kiter->features.end());
      losses.back().push_back(kiter->loss);
    }
  }
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }

  sample_set_type::features_type feats;
  sample_set_type features_optimize;
  alpha_type      alpha;
  f_type          f;

  feature_set_type updates;
  
  double learn(weight_set_type& weights)
  {
    if (losses.empty() || features.empty()) {
      features.clear();
      losses.clear();
      
      return 0.0;
    }
    
    features_optimize.clear();
    f.clear();
    alpha.clear();
    
    weight_type objective;
    
    size_t pos = 0;
    for (size_t i = 0; i != losses.size(); ++ i) 
      if (! losses[i].empty()) {
	weight_type Z;
	expectations.clear();
	expectations_Z.clear();
	margins.clear();
	
	const size_t pos_local = pos;
	for (size_t j = 0; j != losses[i].size(); ++ j, ++ pos) {
	  margins.push_back(cicada::dot_product(weights, features[pos].begin(), features[pos].end(), 0.0) * regularizer.scale() * scale);
	  Z += traits_type::exp(margins.back());
	}
	
	weight_type scaling_sum;
	weight_type objective_local;
	
	for (size_t j = 0, p = pos_local; j != losses[i].size(); ++ j, ++ p) {
	  const weight_type weight = traits_type::exp(margins[j]) / Z;
	  const weight_type scaling = weight_type(scale * losses[i][j]) * weight;
	  
	  scaling_sum += scaling;
	  objective_local += weight_type(losses[i][j]) * weight;
	  
	  sample_set_type::value_type::const_iterator fiter_end = features[p].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[p].begin(); fiter != fiter_end; ++ fiter) {
	    expectations[fiter->first] += weight_type(fiter->second) * scaling;
	    expectations_Z[fiter->first] += weight_type(fiter->second) * weight;
	  }
	}
	
	{
	  expectation_type::const_iterator eiter_end = expectations_Z.end();
	  for (expectation_type::const_iterator eiter = expectations_Z.begin(); eiter != eiter_end; ++ eiter)
	    expectations[eiter->first] -= weight_type(eiter->second) * scaling_sum;
	}
	
	// finaly, encode into features_optimize
	feats.clear();
	expectation_type::const_iterator eiter_end = expectations.end();
	for (expectation_type::const_iterator eiter = expectations.begin(); eiter != eiter_end; ++ eiter)
	  feats.push_back(std::make_pair(eiter->first, - double(eiter->second)));
	
	f.push_back(objective_local);
	features_optimize.insert(feats.begin(), feats.end());
	objective += objective_local;
      }

    // perform optimizatin
    
    const size_type k = losses.size();
    const double k_norm = 1.0 / k;

    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    for (size_t i = 0; i != f.size(); ++ i)
      f[i] = - (f[i] - cicada::dot_product(features_optimize[i].begin(), features_optimize[i].end(), weights, 0.0) * regularizer.scale());
    
    alpha.resize(f.size(), 0.0);
    
    {
      HMatrix H(features_optimize);
      MMatrix M(features_optimize);
      
      cicada::optimize::QPDCD()(alpha, f, H, M, eta, tolerance);
    }
    
    // update by expectations...
    updates.clear();
    size_t actives = 0;
    size_t negatives = 0;
    for (size_t i = 0; i != alpha.size(); ++ i)
      if (alpha[i] > 0.0) {
	sample_set_type::value_type::const_iterator fiter_end = features_optimize[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features_optimize[i].begin(); fiter != fiter_end; ++ fiter)
	  updates[fiter->first]  -= alpha[i] * fiter->second / eta;
	
	++ actives;
	negatives += f[i] > 0.0;
      }
    
    // updated by the merged feature-value
    // we use a constant eta...
    feature_set_type::const_iterator uiter_end = updates.end();
    for (feature_set_type::const_iterator uiter = updates.begin(); uiter != uiter_end; ++ uiter)
      regularizer.update(weights, uiter->first, uiter->second, eta);
    
    if (debug >= 2)
      std::cerr << "actives: " << actives << " negatives: " << negatives << " vectors: " << alpha.size() << std::endl;

    regularizer.postprocess(weights, eta);
    
    features.clear();
    losses.clear();
    
    return objective * k_norm;
  }
  
  margin_set_type    margins;
  sample_set_type    features;
  expectation_type   expectations;
  expectation_type   expectations_Z;
  
  loss_map_type losses;

  double tolerance;
  
  Regularize& regularizer;
  Rate&       rate;
};

struct LearnOnlineMargin : public LearnBase
{
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  
  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur3<size_t>
  {
    typedef utils::hashmurmur3<size_t> hasher_type;

    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
  typedef utils::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> >::type sentence_unique_type;
  
  void clear()
  {
    features.clear();
    losses.clear();

    history_features.clear();
    history_losses.clear();
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    
    if (kbests.empty() || oracles.empty()) return;
    
    sentences.clear();
    for (size_t o = 0; o != oracles.size(); ++ o)
      sentences.insert(oracles[o].sentence);

    features_type feats;
    
    for (size_t o = 0; o != oracles.size(); ++ o)
      for (size_t k = 0; k != kbests.size(); ++ k) {
	const hypothesis_type& oracle = oracles[o];
	const hypothesis_type& kbest  = kbests[k];
	
	if (sentences.find(kbest.sentence) != sentences.end()) continue;
	
	hypothesis_type::feature_set_type::const_iterator oiter = oracle.features.begin();
	hypothesis_type::feature_set_type::const_iterator oiter_end = oracle.features.end();
	
	hypothesis_type::feature_set_type::const_iterator kiter = kbest.features.begin();
	hypothesis_type::feature_set_type::const_iterator kiter_end = kbest.features.end();
	
	feats.clear();
	
	while (oiter != oiter_end && kiter != kiter_end) {
	  if (oiter->first < kiter->first) {
	    feats.push_back(*oiter);
	    ++ oiter;
	  } else if (kiter->first < oiter->first) {
	    feats.push_back(feature_value_type(kiter->first, - kiter->second));
	    ++ kiter;
	  } else {
	    const double value = oiter->second - kiter->second;
	    if (value != 0.0)
	      feats.push_back(feature_value_type(kiter->first, value));
	    ++ oiter;
	    ++ kiter;
	  }
	}
	
	for (/**/; oiter != oiter_end; ++ oiter)
	  feats.push_back(*oiter);
	
	for (/**/; kiter != kiter_end; ++ kiter)
	  feats.push_back(feature_value_type(kiter->first, - kiter->second));
	
	if (feats.empty()) continue;
	
	const double loss = kbest.loss - oracle.loss;
	
	if (loss <= 0.0) continue;
	
	losses.push_back(loss_rank ? 1.0 : loss);
	features.insert(feats.begin(), feats.end());
	
	// for history...
	history_losses.push_back(loss_rank ? 1.0 : loss);
	history_features.insert(feats.begin(), feats.end());
      }
  }

  template <typename Iterator>
  boost::fusion::tuple<double, double, double> gradient(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    static const double inf = std::numeric_limits<double>::infinity();
    
    double grad_pos = 0.0;
    double grad_neg = 0.0;
    for (size_t i = 0; i != history_features.size(); ++ i) {
      const double margin      = cicada::dot_product(weights,      history_features[i].begin(), history_features[i].end(), 0.0);
      const double margin_prev = cicada::dot_product(weights_prev, history_features[i].begin(), history_features[i].end(), 0.0);
      const double loss = history_losses[i];
      
      const double bi_pos = margin_prev - margin;
      const double ci_pos = loss - margin_prev;
      const double ki_pos = (bi_pos != 0.0 ? - ci_pos / bi_pos : - inf);
      
      const double bi_neg = margin_prev + margin;
      const double ci_neg = loss - margin_prev;
      const double ki_neg = (bi_neg != 0.0 ? - ci_neg / bi_neg : - inf);
      
      if (ki_pos > 0) {
	*iter = std::make_pair(ki_pos, bi_pos);
	++ iter;
      }
      
      if (ki_neg > 0) {
	*iter = std::make_pair(- ki_neg, bi_neg);
	++ iter;
      }
      
      grad_pos += bi_pos * ((bi_pos < 0.0 && ki_pos > 0.0) || (bi_pos > 0.0 && ki_pos <= 0.0));
      grad_neg += bi_neg * ((bi_neg < 0.0 && ki_neg > 0.0) || (bi_neg > 0.0 && ki_neg <= 0.0));
    }
    
    return boost::fusion::tuple<double, double, double>(grad_pos, grad_neg, history_losses.size());
  }


  sample_set_type features;
  loss_set_type   losses;
  
  sample_set_type history_features;
  loss_set_type   history_losses;
  
  sentence_unique_type sentences;
};

// Pegasos learner
struct LearnHinge : public LearnOnlineMargin
{
  
  LearnHinge(Regularize& __regularizer,
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
  
  typedef std::vector<bool, std::allocator<bool> > suffered_set_type;
  suffered_set_type suffered;

  feature_set_type updates;

  double learn(weight_set_type& weights)
  {
    if (features.empty()) return 0.0;
    
    size_type k = 0;
    suffered.clear();
    suffered.resize(features.size(), false);
    
    for (size_t i = 0; i != features.size(); ++ i) {
      const double loss = losses[i] - cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0) * regularizer.scale();
      const bool suffer_loss = loss > 0.0;
      suffered[i] = suffer_loss;
      k += suffer_loss;
    }
    
    if (! k) {
      // anyway, clear features!
      features.clear();
      losses.clear();
      
      return 0.0;
    }
    
    //const double k_norm = 1.0 / (features.size());
    const double k_norm = 1.0 / k; // it is wrong, but works quite well in practice
    
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    updates.clear();
    for (size_t i = 0; i != features.size(); ++ i) 
      if (suffered[i]) {
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter)
	  updates[fiter->first] -= k_norm * fiter->second;
      }
    
    // udpate...
    feature_set_type::const_iterator fiter_end = updates.end();
    for (feature_set_type::const_iterator fiter = updates.begin(); fiter != fiter_end; ++ fiter)
      regularizer.update(weights, fiter->first, fiter->second, rate(fiter->first, fiter->second));
    
    regularizer.postprocess(weights, eta);
    
    features.clear();
    losses.clear();
    
    return 0.0;
  }
  
  Regularize& regularizer;
  Rate&       rate;
};

// optimized-Pegasos learner
struct LearnOHinge : public LearnOnlineMargin
{
  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;
  typedef std::vector<int, std::allocator<int> >          index_type;

  struct HMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    HMatrix(const sample_set_type& __features, const index_type& __index) : features(__features), index(__index) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[index[i]].begin(), features[index[i]].end(), features[index[j]].begin(), features[index[j]].end(), 0.0);
    }
    
    const sample_set_type& features;
    const index_type& index;
  };
  
  struct MMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    MMatrix(const sample_set_type& __features, const index_type& __index) : features(__features), index(__index) {}
    
    template <typename __W>
    void operator()(__W& w, const alpha_type& alpha) const
    {
      const size_type model_size = index.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename __W>
    double operator()(const __W& w, const size_t& i) const
    {
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
      for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename __W>
    void operator()(__W& w, const double& update, const size_t& i) const
    {
      sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
      for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const sample_set_type& features;
    const index_type& index;
  };

  LearnOHinge(Regularize& __regularizer,
	      Rate& __rate)
    : tolerance(0.1),
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

  feature_set_type updates;

  double learn(weight_set_type& weights)
  {
    if (features.empty()) return 0.0;
    
    const double eta = rate();
    
    alpha.clear();
    f.clear();
    index.clear();
    
    for (size_t i = 0; i != losses.size(); ++ i)
      f.push_back(cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0));
    
    const double weight_scale_curr = regularizer.scale();

    regularizer.preprocess(weights, eta);
    
    double objective = 0.0;
    for (size_t i = 0; i != f.size(); ++ i) {
      const double loss = losses[i] - f[i] * weight_scale_curr;
      
      if (loss <= 0.0) continue;
      
      f[index.size()] = - (losses[i] - f[i] * regularizer.scale());
      index.push_back(i);
      objective += loss;
    }
    objective /= losses.size();
    
    f.resize(index.size());
    alpha.resize(index.size(), 0.0);
    
    {
      HMatrix H(features, index);
      MMatrix M(features, index);
      
      cicada::optimize::QPDCD()(alpha, f, H, M, eta, tolerance);
    }
    
    updates.clear();
    size_t actives = 0;
    size_t negatives = 0;
    for (size_t i = 0; i != index.size(); ++ i)
      if (alpha[i] > 0.0) {
	sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
	for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter)
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
    
    features.clear();
    losses.clear();
    
    return 0.0;
  }
  
  double    tolerance;

  Regularize& regularizer;
  Rate&       rate;
  
  alpha_type    alpha;
  f_type        f;
  index_type    index;
};

struct LearnPA : public LearnOnlineMargin
{
  LearnPA(const double& __lambda) : lambda(__lambda)
  {
    if (__lambda <= 0.0)
      throw std::runtime_error("lambda is <= 0?");
  }

  void initialize(weight_set_type& weights) {}
  
  void finalize(weight_set_type& weights) {}
  
  double learn(weight_set_type& weights)
  {
    double objective = 0.0;

    const double constant = 1.0 / (lambda * losses.size());
    
    for (size_t i = 0; i != losses.size(); ++ i) {
      const double margin = cicada::dot_product(weights, features[i].begin(), features[i].end(), 0.0);
      const double loss   = losses[i];
      
      const double suffered = loss - margin;

      if (suffered <= 0.0) continue;
      
      // PA-I
      const double variance = cicada::dot_product(features[i].begin(), features[i].end(), features[i].begin(), features[i].end(), 0.0);
      const double alpha = std::min(suffered / variance, constant);
      
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter)
	weights[fiter->first] += alpha * fiter->second;
      
      objective += suffered;
    }

    objective /= losses.size();
    
    features.clear();
    losses.clear();
    
    return objective;
  }
  
  double lambda;
};

struct LearnCW : public LearnOnlineMargin
{
  LearnCW(const double& __lambda) : lambda(__lambda)
  {
    if (__lambda <= 0.0)
      throw std::runtime_error("lambda is <= 0?");
  }

  void initialize(weight_set_type& weights) {}
  
  void finalize(weight_set_type& weights) {}
  
  double learn(weight_set_type& weights)
  {
    double objective = 0.0;

    covariances.allocate(1.0);
    
    for (size_t i = 0; i != losses.size(); ++ i) {
      const double margin = cicada::dot_product(weights, features[i].begin(), features[i].end(), 0.0);
      const double loss   = losses[i];
      
      const double suffered = loss - margin;
      
      if (suffered <= 0.0) continue;

      const double variance = cicada::dot_product(features[i].begin(), features[i].end(), covariances, features[i].begin(), features[i].end(), 0.0);
      
      const double theta = 1.0 + 2.0 * lambda * (margin - loss);
      const double alpha = ((- theta + std::sqrt(theta * theta - 8.0 * lambda * (margin - loss - lambda * variance))) / (4.0 * lambda * variance));
      const double beta  = (2.0 * alpha * lambda) / (1.0 + 2.0 * alpha * lambda * variance);

      if (alpha > 1e-12 && beta > 0.0) {
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	  const double var = covariances[fiter->first];
	  
	  weights[fiter->first]     += alpha * fiter->second * var;
	  covariances[fiter->first] -= beta * (var * var) * (fiter->second * fiter->second);
	}
      }
      
      objective += suffered;
    }

    objective /= losses.size();

    features.clear();
    losses.clear();
    
    return objective;
  }
  
  weight_set_type covariances;
  double lambda;
};

struct LearnAROW : public LearnOnlineMargin
{
  LearnAROW(const double& __lambda) : lambda(__lambda)
  {
    if (__lambda <= 0.0)
      throw std::runtime_error("lambda is <= 0?");
  }

  void initialize(weight_set_type& weights) {}
  
  void finalize(weight_set_type& weights) {}
  
  double learn(weight_set_type& weights)
  {
    double objective = 0.0;

    covariances.allocate(1.0);
    
    for (size_t i = 0; i != losses.size(); ++ i) {
      const double margin = cicada::dot_product(weights, features[i].begin(), features[i].end(), 0.0);
      const double loss   = losses[i];
      
      const double suffered = loss - margin;

      if (suffered <= 0.0) continue;

      const double variance = cicada::dot_product(features[i].begin(), features[i].end(), covariances, features[i].begin(), features[i].end(), 0.0);
      
      const double beta = 1.0 / (variance + lambda * losses.size());
      const double alpha = std::max(0.0, (loss - margin) * beta);
      
      if (alpha > 1e-12) {
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	  const double var = covariances[fiter->first];
	  
	  weights[fiter->first]     += alpha * fiter->second * var;
	  covariances[fiter->first] -= beta * (var * var) * (fiter->second * fiter->second);
	}
      }
      
      objective += suffered;
    }
    
    objective /= losses.size();

    features.clear();
    losses.clear();
    
    return objective;
  }
  
  weight_set_type covariances;
  double lambda;
};

struct LearnNHERD : public LearnOnlineMargin
{
  LearnNHERD(const double& __lambda) : lambda(1.0 / __lambda)
  {
    if (__lambda <= 0.0)
      throw std::runtime_error("lambda is <= 0?");
  }

  void initialize(weight_set_type& weights) {}
  
  void finalize(weight_set_type& weights) {}
  
  double learn(weight_set_type& weights)
  {
    double objective = 0.0;

    covariances.allocate(1.0);
    
    for (size_t i = 0; i != losses.size(); ++ i) {
      const double margin = cicada::dot_product(weights, features[i].begin(), features[i].end(), 0.0);
      const double loss   = losses[i];
      
      const double suffered = loss - margin;
      
      if (suffered <= 0.0) continue;
      
      const double variance = cicada::dot_product(features[i].begin(), features[i].end(), covariances, features[i].begin(), features[i].end(), 0.0);
      const double alpha = std::max(0.0, (loss - margin) / (variance + 1.0 / lambda));
      
      if (alpha > 1e-12) {
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	  const double var = covariances[fiter->first];
	
	  weights[fiter->first]     += alpha * fiter->second * var;
	  //covariances[fiter->first]  = 1.0 / ((1.0 / var) + (2.0 * lambda + lambda * lambda * variance) * fiter->second * fiter->second);
	  covariances[fiter->first]  = var / (1.0 + var * (2.0 * lambda + lambda * lambda * variance) * fiter->second * fiter->second);
	}
      }
      
      objective += suffered;
    }

    objective /= losses.size();

    features.clear();
    losses.clear();
    
    return objective;    
  }
  
  weight_set_type covariances;
  double lambda;
};


// MIRA learner
// We will run a qp solver and determine the alpha, then, translate this into w
struct LearnMIRA : public LearnOnlineMargin
{
  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;
  typedef std::vector<int, std::allocator<int> >          index_type;

  struct HMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    HMatrix(const sample_set_type& __features, const index_type& __index) : features(__features), index(__index) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[index[i]].begin(), features[index[i]].end(), features[index[j]].begin(), features[index[j]].end(), 0.0);
    }
    
    const sample_set_type& features;
    const index_type& index;
  };
  
  struct MMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    MMatrix(const sample_set_type& __features, const index_type& __index) : features(__features), index(__index) {}
    
    template <typename __W>
    void operator()(__W& w, const alpha_type& alpha) const
    {
      const size_type model_size = index.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename __W>
    double operator()(const __W& w, const size_t& i) const
    {
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
      for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename __W>
    void operator()(__W& w, const double& update, const size_t& i) const
    {
      sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
      for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const sample_set_type& features;
    const index_type& index;
  };

  LearnMIRA(const double& __lambda) : tolerance(0.1), lambda(__lambda)
  {
    if (__lambda <= 0.0)
      throw std::runtime_error("lambda is <= 0?");
  }
  
  void initialize(weight_set_type& weights) {}
  
  void finalize(weight_set_type& weights) {}
  
  double learn(weight_set_type& weights)
  {
    if (features.empty()) return 0.0;
    
    alpha.clear();
    index.clear();
    f.clear();
    
    double objective = 0.0;
    for (size_t i = 0; i != losses.size(); ++ i) {
      const double loss = losses[i] - cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0);
      
      if (loss <= 0.0) continue;
      
      f.push_back(- loss);
      index.push_back(i);
      
      objective += loss;
    }
    
    objective /= losses.size();
    
    alpha.resize(index.size(), 0.0);
    
    {
      HMatrix H(features, index);
      MMatrix M(features, index);
      
      cicada::optimize::QPDCD()(alpha, f, H, M, 1.0 / lambda, tolerance);
    }
    
    for (size_t i = 0; i != alpha.size(); ++ i)
      if (alpha[i] > 0.0) {
	// update: weights[fiter->first] += alpha[i] * fiter->second;
	
	sample_set_type::value_type::const_iterator fiter_end = features[index[i]].end();
	for (sample_set_type::value_type::const_iterator fiter = features[index[i]].begin(); fiter != fiter_end; ++ fiter)
	  weights[fiter->first] += alpha[i] * fiter->second;
      }
    
    features.clear();
    losses.clear();
    
    return objective;
  }
  
  double tolerance;
  double lambda;

  alpha_type    alpha;
  f_type        f;
  index_type    index;
};

// logistic regression base...
struct LearnSoftmaxBase : public LearnBase
{
  typedef cicada::semiring::Log<double> weight_type;

  typedef std::vector<double, std::allocator<double> > loss_set_type;
  
  struct sample_pair_type
  {
    sample_pair_type() : kbests(), oracles(), loss_kbests(), loss_oracles() {}
    
    sample_set_type kbests;
    sample_set_type oracles;

    loss_set_type loss_kbests;
    loss_set_type loss_oracles;

    void encode(const hypothesis_set_type& __kbests, const hypothesis_set_type& __oracles)
    {
      if (__kbests.empty() || __oracles.empty()) return;
      
      hypothesis_set_type::const_iterator kiter_end = __kbests.end();
      for (hypothesis_set_type::const_iterator kiter = __kbests.begin(); kiter != kiter_end; ++ kiter) {
	kbests.insert(kiter->features.begin(), kiter->features.end());
	loss_kbests.push_back(kiter->loss);
      }
      
      hypothesis_set_type::const_iterator oiter_end = __oracles.end();
      for (hypothesis_set_type::const_iterator oiter = __oracles.begin(); oiter != oiter_end; ++ oiter) {
	oracles.insert(oiter->features.begin(), oiter->features.end());
	loss_oracles.push_back(oiter->loss);
      }
    }

    typedef std::vector<double, std::allocator<double> > margin_set_type;
    
    margin_set_type __margins;
    
    template <typename Expectations>
    double encode(const weight_set_type& weights, Expectations& expectations, const double scale) const
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
      
      weight_type Z_oracle;
      weight_type Z_kbest;
      
      const double cost_factor = (softmax_margin ? 1.0 : 0.0);

      margin_set_type& margins = const_cast<margin_set_type&>(__margins);

      margins.clear();
      
      for (size_type o = 0; o != oracles.size(); ++ o) {
	margins.push_back(cicada::dot_product(weights, oracles[o].begin(), oracles[o].end(), 0.0) * scale + cost_factor * loss_oracles[o]);
	Z_oracle += traits_type::exp(margins.back());
      }
      
      for (size_type k = 0; k != kbests.size(); ++ k) {
	margins.push_back(cicada::dot_product(weights, kbests[k].begin(), kbests[k].end(), 0.0) * scale + cost_factor * loss_kbests[k]);
	Z_kbest += traits_type::exp(margins.back());
      }

      margin_set_type::const_iterator miter = margins.begin();
      
      for (size_type o = 0; o != oracles.size(); ++ o, ++ miter) {
	const weight_type weight = traits_type::exp(*miter) / Z_oracle;
	
	sample_set_type::value_type::const_iterator fiter_end = oracles[o].end();
	for (sample_set_type::value_type::const_iterator fiter = oracles[o].begin(); fiter != fiter_end; ++ fiter)
	  expectations[fiter->first] -= weight_type(fiter->second) * weight;
      }
      
      for (size_type k = 0; k != kbests.size(); ++ k, ++ miter) {
	const weight_type weight = traits_type::exp(*miter) / Z_kbest;	
	
	sample_set_type::value_type::const_iterator fiter_end = kbests[k].end();
	for (sample_set_type::value_type::const_iterator fiter = kbests[k].begin(); fiter != fiter_end; ++ fiter)
	  expectations[fiter->first] += weight_type(fiter->second) * weight;
      }
      
      return log(Z_oracle) - log(Z_kbest);
    }
  };
};


// Softmax learner
struct LearnSoftmax : public LearnSoftmaxBase
{
  typedef utils::chunk_vector<sample_pair_type, 4096 / sizeof(sample_pair_type), std::allocator<sample_pair_type> > sample_pair_set_type;
    
  LearnSoftmax(Regularize& __regularizer,
	       Rate& __rate)
    : regularizer(__regularizer),
      rate(__rate) {}
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles)
  {
    if (kbests.empty() || oracles.empty()) return;
    
    samples.push_back(sample_pair_type());
    samples.back().encode(kbests, oracles);
  }
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }

  template <typename Iterator>
  boost::fusion::tuple<double, double, double> gradient(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return boost::fusion::tuple<double, double, double>(0.0, 0.0, 0.0);
  }

  void clear() { samples.clear(); }

  double learn(weight_set_type& weights)
  {
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > expectation_type;
    
    if (samples.empty()) return 0.0;
    
    const size_type k = samples.size();
    const double k_norm = 1.0 / k;
    
    const double eta = rate();
    
    regularizer.preprocess(weights, eta);
    
    expectation_type expectations;
    
    // update... by eta / k
    double objective = 0.0;
    sample_pair_set_type::const_iterator siter_end = samples.end();
    for (sample_pair_set_type::const_iterator siter = samples.begin(); siter != siter_end; ++ siter)
      objective += siter->encode(weights, expectations, regularizer.scale());
    objective /= samples.size();
    
    // update by expectations...
    expectation_type::const_iterator eiter_end = expectations.end();
    for (expectation_type::const_iterator eiter = expectations.begin(); eiter != eiter_end; ++ eiter) {
      const double amount = static_cast<double>(eiter->second) * k_norm;
      
      regularizer.update(weights, eiter->first, amount, rate(eiter->first, amount));
    }
    
    regularizer.postprocess(weights, eta);
    
    // clear current training events..
    samples.clear();
    
    return objective;
  }
  
  sample_pair_set_type samples;
  
  Regularize& regularizer;
  Rate&       rate;
};

// Softmax learner
struct LearnOSoftmax : public LearnSoftmaxBase
{
  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;
  
  struct HMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    HMatrix(const sample_set_type& __features) : features(__features) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[i].begin(), features[i].end(), features[j].begin(), features[j].end(), 0.0);
    }
    
    const sample_set_type& features;
  };
  
  struct MMatrix
  {
    typedef LearnBase::sample_set_type sample_set_type;
    
    MMatrix(const sample_set_type& __features) : features(__features) {}
    
    template <typename __W>
    void operator()(__W& w, const alpha_type& alpha) const
    {
      const size_type model_size = features.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename __W>
    double operator()(const __W& w, const size_t& i) const
    {
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename __W>
    void operator()(__W& w, const double& update, const size_t& i) const
    {
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const sample_set_type& features;
  };

  typedef utils::chunk_vector<sample_pair_type, 4096 / sizeof(sample_pair_type), std::allocator<sample_pair_type> > sample_pair_set_type;
    
  LearnOSoftmax(Regularize& __regularizer,
		Rate& __rate)
    : tolerance(0.1),
      regularizer(__regularizer),
      rate(__rate) {}
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles)
  {
    if (kbests.empty() || oracles.empty()) return;
    
    samples.push_back(sample_pair_type());
    samples.back().encode(kbests, oracles);
  }
  
  void initialize(weight_set_type& weights)
  {
    regularizer.initialize(weights);
  }
  
  void finalize(weight_set_type& weights)
  {
    regularizer.finalize(weights);
  }

  template <typename Iterator>
  boost::fusion::tuple<double, double, double> gradient(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return boost::fusion::tuple<double, double, double>(0.0, 0.0, 0.0);
  }

  void clear() { samples.clear(); }

  sample_set_type::features_type feats;
  sample_set_type features;
  alpha_type      alpha;
  f_type          f;
  
  feature_set_type updates;

  double learn(weight_set_type& weights)
  {
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > expectation_type;
    
    if (samples.empty()) return 0.0;
    
    const double eta = rate();
    
    expectation_type expectations;
    
    features.clear();
    f.clear();
    alpha.clear();

    // use LBFGS?
    // we will minimize ||x - x'|| + loss...
    //
    
    // update... by eta / k
    double objective = 0.0;
    sample_pair_set_type::const_iterator siter_end = samples.end();
    for (sample_pair_set_type::const_iterator siter = samples.begin(); siter != siter_end; ++ siter) {
      expectations.clear();
      f.push_back(siter->encode(weights, expectations, regularizer.scale()));
      
      objective += f.back();
      
      feats.clear();
      expectation_type::const_iterator eiter_end = expectations.end();
      for (expectation_type::const_iterator eiter = expectations.begin(); eiter != eiter_end; ++ eiter)
	feats.push_back(std::make_pair(eiter->first, - double(eiter->second)));
      
      features.insert(feats.begin(), feats.end());
    }
    objective /= samples.size();
    
    // perform rescaling here!
    regularizer.preprocess(weights, eta);
    
    for (size_t i = 0; i != samples.size(); ++ i)
      f[i] = -(- f[i] - cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0) * regularizer.scale());
    
    alpha.resize(f.size(), 0.0);
    
    {
      HMatrix H(features);
      MMatrix M(features);
      
      cicada::optimize::QPDCD()(alpha, f, H, M, eta, tolerance);
    }
    
    // update by expectations...
    updates.clear();
    size_t actives = 0;
    size_t negatives = 0;
    for (size_t i = 0; i != alpha.size(); ++ i)
      if (alpha[i] > 0.0) {
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter)
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
    
    // clear current training events..
    samples.clear();
    
    return objective;
  }
  
  sample_pair_set_type samples;

  double tolerance;

  Regularize& regularizer;
  Rate&       rate;
};

struct KBestSentence
{
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::sentence_feature_traversal   traversal_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  
  typedef cicada::operation::kbest_sentence_filter_unique filter_unique_type;
  typedef cicada::operation::kbest_sentence_filter        filter_type;
  
  typedef cicada::KBest<traversal_type, function_type, filter_unique_type> derivation_unique_set_type;
  typedef cicada::KBest<traversal_type, function_type, filter_type>        derivation_set_type;

  typedef traversal_type::value_type derivation_type;
  
  void operator()(operation_set_type& operations, const std::string& input, const weight_set_type& weights, hypothesis_set_type& kbests)
  {    
    kbests.clear();
    
    if (input.empty()) return;
    
    // assign weights
    operations.assign(weights);
    
    // perform translation
    operations(input);
    
    // generate kbests...
    const hypergraph_type& graph = operations.get_data().hypergraph;
    
    if (kbest_diverse_mode) {
      derivation_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_type());

      derivation_set_type::const_iterator diter_end = derivations.end();
      for (derivation_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter)
	kbests.push_back(hypothesis_type(boost::get<0>(diter->second).begin(), boost::get<0>(diter->second).end(),
					 boost::get<1>(diter->second).begin(), boost::get<1>(diter->second).end()));
    } else {
      derivation_unique_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_unique_type(graph));
      
      derivation_unique_set_type::const_iterator diter_end = derivations.end();
      for (derivation_unique_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter)
	kbests.push_back(hypothesis_type(boost::get<0>(diter->second).begin(), boost::get<0>(diter->second).end(),
					 boost::get<1>(diter->second).begin(), boost::get<1>(diter->second).end()));
    }
  }
};

struct KBestAlignment
{
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::alignment_feature_traversal  traversal_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  
  typedef cicada::operation::kbest_alignment_filter_unique filter_unique_type;
  typedef cicada::operation::kbest_alignment_filter        filter_type;
  
  typedef cicada::KBest<traversal_type, function_type, filter_unique_type> derivation_unique_set_type;
  typedef cicada::KBest<traversal_type, function_type, filter_type>        derivation_set_type;

  typedef traversal_type::value_type derivation_type;
  
  void operator()(operation_set_type& operations, const std::string& input, const weight_set_type& weights, hypothesis_set_type& kbests)
  {    
    kbests.clear();
    
    if (input.empty()) return;
    
    // assign weights
    operations.assign(weights);
    
    // perform translation
    operations(input);
    
    // generate kbests...
    const hypergraph_type& graph = operations.get_data().hypergraph;
    
    if (kbest_diverse_mode) {
      derivation_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_type());
      
      sentence_type   sentence;

      derivation_set_type::const_iterator diter_end = derivations.end();
      for (derivation_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
	std::ostringstream os;
	os << boost::get<0>(diter->second);
	sentence.assign(os.str());
	
	kbests.push_back(hypothesis_type(sentence.begin(), sentence.end(),
					 boost::get<1>(diter->second).begin(), boost::get<1>(diter->second).end()));
      }
    } else {
      derivation_unique_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_unique_type(graph));
      
      sentence_type   sentence;

      derivation_unique_set_type::const_iterator diter_end = derivations.end();
      for (derivation_unique_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
	std::ostringstream os;
	os << boost::get<0>(diter->second);
	sentence.assign(os.str());
	
	kbests.push_back(hypothesis_type(sentence.begin(), sentence.end(),
					 boost::get<1>(diter->second).begin(), boost::get<1>(diter->second).end()));
      }
    }
  }
};

struct KBestDependency
{
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::dependency_feature_traversal traversal_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  
  typedef cicada::operation::kbest_dependency_filter_unique filter_unique_type;
  typedef cicada::operation::kbest_dependency_filter        filter_type;
  
  typedef cicada::KBest<traversal_type, function_type, filter_unique_type> derivation_unique_set_type;
  typedef cicada::KBest<traversal_type, function_type, filter_type>        derivation_set_type;

  typedef traversal_type::value_type derivation_type;
  
  void operator()(operation_set_type& operations, const std::string& input, const weight_set_type& weights, hypothesis_set_type& kbests)
  {    
    kbests.clear();
    
    if (input.empty()) return;
    
    // assign weights
    operations.assign(weights);
    
    // perform translation
    operations(input);
    
    // generate kbests...
    const hypergraph_type& graph = operations.get_data().hypergraph;

    if (kbest_diverse_mode) {
      derivation_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_type());
      
      sentence_type   sentence;

      derivation_set_type::const_iterator diter_end = derivations.end();
      for (derivation_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
	std::ostringstream os;
	os << boost::get<0>(diter->second);
	sentence.assign(os.str());
	
	kbests.push_back(hypothesis_type(sentence.begin(), sentence.end(),
					 boost::get<1>(diter->second).begin(), boost::get<1>(diter->second).end()));
      }
    } else {
      derivation_unique_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_unique_type(graph));

      sentence_type   sentence;

      derivation_unique_set_type::const_iterator diter_end = derivations.end();
      for (derivation_unique_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
	std::ostringstream os;
	os << boost::get<0>(diter->second);
	sentence.assign(os.str());
	
	kbests.push_back(hypothesis_type(sentence.begin(), sentence.end(),
					 boost::get<1>(diter->second).begin(), boost::get<1>(diter->second).end()));
      }
    }
  }
};

struct Oracle
{
  typedef std::vector<const hypothesis_type*, std::allocator<const hypothesis_type*> > oracle_set_type;
  typedef std::vector<oracle_set_type, std::allocator<oracle_set_type> > oracle_map_type;

  template <typename Generator>
  std::pair<score_ptr_type, score_ptr_type>
  operator()(const hypothesis_map_type& kbests, const hypothesis_map_type& kbests_oracle, const scorer_document_type& scorers, hypothesis_map_type& oracles, Generator& generator)
  {
    score_ptr_type score_1best;
    score_ptr_type score_oracle;
    
    // initialization...
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	hypothesis_set_type::const_iterator hiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(*hiter);
	  
	  if (! hyp.score)
	    hyp.score = scorers[id]->score(sentence_type(hyp.sentence.begin(), hyp.sentence.end()));
	}
	
	if (! score_1best)
	  score_1best = kbests[id].front().score->clone();
	else
	  *score_1best += *(kbests[id].front().score);
      }
    
    if (! score_1best)
      throw std::runtime_error("no evaluation score?");
    
    // assign loss
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	score_ptr_type score_segment = score_1best->clone();
	*score_segment -= *kbests[id].front().score;
	
	hypothesis_set_type::const_iterator hiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(*hiter);
	  
	  *score_segment += *hyp.score;
	  
	  hyp.loss = score_segment->loss();
	  
	  *score_segment -= *hyp.score;
	}
      }
    
    const size_t kbests_size = utils::bithack::min(kbests_oracle.size(), kbests.size());
    for (size_t id = 0; id != kbests_size; ++ id) 
      if (! kbests_oracle[id].empty() && ! kbests[id].empty()) {
	score_ptr_type score_segment = score_1best->clone();
	*score_segment -= *kbests[id].front().score;
	
	hypothesis_set_type::const_iterator hiter_end = kbests_oracle[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests_oracle[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(*hiter);
	  
	  if (! hyp.score)
	    hyp.score = scorers[id]->score(sentence_type(hyp.sentence.begin(), hyp.sentence.end()));
	  
	  *score_segment += *hyp.score;
	  
	  hyp.loss = score_segment->loss();
	  
	  *score_segment -= *hyp.score;
	}
	
	if (! score_oracle)
	  score_oracle = kbests_oracle[id].front().score->clone();
	else
	  *score_oracle += *(kbests_oracle[id].front().score);
      }

    if (! score_oracle)
      throw std::runtime_error("no oracle evaluation score?");

    // copy!
    oracles = kbests_oracle;

    return std::make_pair(score_1best, score_oracle);
  }
    
  template <typename Generator>
  std::pair<score_ptr_type, score_ptr_type>
  operator()(const hypothesis_map_type& kbests, const scorer_document_type& scorers, hypothesis_map_type& oracles, Generator& generator)
  {
    typedef std::vector<size_t, std::allocator<size_t> > id_set_type;

    score_ptr_type score_1best;    
    
    id_set_type ids;
    boost::random_number_generator<boost::mt19937> gen(generator);
    
    // initialization...
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	ids.push_back(id);
	
	hypothesis_set_type::const_iterator hiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(*hiter);
	  
	  if (! hyp.score)
	    hyp.score = scorers[id]->score(sentence_type(hyp.sentence.begin(), hyp.sentence.end()));
	}
	
	if (! score_1best)
	  score_1best = kbests[id].front().score->clone();
	else
	  *score_1best += *(kbests[id].front().score);
      }
    
    if (! score_1best)
      throw std::runtime_error("no evaluation score?");
    
    // assign loss
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	score_ptr_type score_segment = score_1best->clone();
	*score_segment -= *kbests[id].front().score;
	
	hypothesis_set_type::const_iterator hiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(*hiter);
	  
	  *score_segment += *hyp.score;
	  
	  hyp.loss = score_segment->loss();

	  *score_segment -= *hyp.score;
	}
      }
    
    score_ptr_type score_best;
    score_ptr_type score_curr;
    score_ptr_type score_next;
    
    oracle_map_type oracles_best(kbests.size());
    oracle_map_type oracles_curr(kbests.size());
    oracle_map_type oracles_next(kbests.size());
    
    double objective_best = - std::numeric_limits<double>::infinity();
    double objective_curr = - std::numeric_limits<double>::infinity();
    double objective_next = - std::numeric_limits<double>::infinity();
    
    // 
    // 10 iteration will be fine
    //
    for (int i = 0; i < 10; ++ i) {
      
      for (id_set_type::const_iterator iiter = ids.begin(); iiter != ids.end(); ++ iiter) {
	const size_t id = *iiter;
	
	score_ptr_type score_removed = (score_next ? score_next->clone() : score_ptr_type());
	
	if (score_removed && ! oracles_curr[id].empty())
	  *score_removed -= *(oracles_curr[id].front()->score);
	
	oracles_next[id].clear();
	
	hypothesis_set_type::const_iterator hiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  score_ptr_type score_sample;
	  
	  if (score_removed) {
	    score_sample = score_removed->clone();
	    *score_sample += *(hiter->score);
	  } else
	    score_sample = hiter->score->clone();
	  
	  const double objective_sample = score_sample->reward();
	  
	  if (objective_sample > objective_next || oracles_next[id].empty()) {
	    oracles_next[id].clear();
	    oracles_next[id].push_back(&(*hiter));
	    
	    objective_next = objective_sample;
	    score_next     = score_sample;
	  } else if (objective_sample == objective_next)
	    oracles_next[id].push_back(&(*hiter));
	}
      }
      
      if (objective_next > objective_best) {
	score_best     = score_next->clone();
	objective_best = objective_next;
	oracles_best   = oracles_next;
      }
      
      if (objective_next <= objective_curr) break;
      
      score_curr     = score_next->clone();
      objective_curr = objective_next;
      oracles_curr   = oracles_next;
      
      std::random_shuffle(ids.begin(), ids.end(), gen);
    }
    
    oracles.clear();
    oracles.resize(kbests.size());
    for (size_t id = 0; id != kbests.size(); ++ id)
      for (size_t i = 0; i != oracles_best[id].size(); ++ i)
	oracles[id].push_back(*oracles_best[id][i]);

    return std::make_pair(score_1best, score_best);
  }
};


inline
void read_refset(const path_type& refset_path,
		 scorer_document_type& scorers,
		 const size_t shard_rank=0,
		 const size_t shard_size=1)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
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
    
    if (! scorers[id_rank])
      scorers[id_rank] = scorers.create();
    
    scorers[id_rank]->insert(id_sentence.second);
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
		 const size_t shard_rank,
		 const size_t shard_size)
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
