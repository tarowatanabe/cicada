//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_BLOCK_IMPL__HPP__
#define __CICADA_LEARN_BLOCK_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <sstream>
#include <vector>

#include "cicada_kbest_impl.hpp"

#include "cicada/semiring.hpp"
#include "cicada/eval.hpp"

#include "cicada/kbest.hpp"
#include "cicada/operation/traversal.hpp"
#include "cicada/operation/functional.hpp"
#include "cicada/optimize_qp.hpp"

#include "utils/sgi_hash_set.hpp"
#include "utils/base64.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"

#include <boost/tokenizer.hpp>

#include "lbfgs.h"
#include "liblinear/linear.h"

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


//
// SVM for structured output learning  
//

struct LearnSVM : public LearnBase
{
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  typedef std::vector<double, std::allocator<double> > alpha_set_type;
  typedef std::vector<double, std::allocator<double> > bound_map_type;
  typedef std::vector<double, std::allocator<double> > f_set_type;
  
  typedef utils::chunk_vector<sample_set_type, 4096 / sizeof(sample_set_type), std::allocator<sample_set_type> > sample_map_type;
  typedef utils::chunk_vector<loss_set_type, 4096 / sizeof(loss_set_type), std::allocator<loss_set_type> >       loss_map_type;
  typedef utils::chunk_vector<alpha_set_type, 4096 / sizeof(alpha_set_type), std::allocator<alpha_set_type> >    alpha_map_type;

  typedef std::pair<size_type, size_type> pos_pair_type;
  typedef std::vector<pos_pair_type, std::allocator<pos_pair_type> > pos_pair_set_type;
  
  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#endif
  
  
  struct HMatrix
  {
    HMatrix(const pos_pair_set_type& __positions,
	    const sample_map_type& __features)
      : positions(__positions), features(__features) {}

    double operator()(int i, int j) const
    {
      const pos_pair_type& pos_i = positions[i];
      const pos_pair_type& pos_j = positions[j];
      
      return cicada::dot_product(features[pos_i.first][pos_i.second].begin(), features[pos_i.first][pos_i.second].end(),
				 features[pos_j.first][pos_j.second].begin(), features[pos_j.first][pos_j.second].end(),
				 0.0);
    }
    
    const pos_pair_set_type& positions;
    const sample_map_type&   features;
  };
  
  struct MMatrix
  {
    MMatrix(const pos_pair_set_type& __positions,
	    const sample_map_type&   __features)
      : positions(__positions), features(__features) {}
    
    template <typename W>
    void operator()(W& w, const alpha_set_type& alpha) const
    {
      alpha_set_type::const_iterator aiter = alpha.begin();
      
      const size_type model_size = features.size();
      for (size_type i = 0; i != model_size; ++ i) {
	const size_type features_size = features[i].size();
	
	for (size_type j = 0; j != features_size; ++ j, ++ aiter)
	  if (*aiter > 0.0) {
	    sample_set_type::value_type::const_iterator fiter_end = features[i][j].end();
	    for (sample_set_type::value_type::const_iterator fiter = features[i][j].begin(); fiter != fiter_end; ++ fiter) 
	      w[fiter->first] += (*aiter) * fiter->second;
	  }
      }
    }
    
    template <typename W>
    double operator()(const W& w, const size_t& i) const
    {
      const pos_pair_type& pos_i = positions[i];
      
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[pos_i.first][pos_i.second].end();
      for (sample_set_type::value_type::const_iterator fiter = features[pos_i.first][pos_i.second].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename W>
    void operator()(W& w, const double& update, const size_t& i) const
    {
      const pos_pair_type& pos_i = positions[i];
      
      sample_set_type::value_type::const_iterator fiter_end = features[pos_i.first][pos_i.second].end();
      for (sample_set_type::value_type::const_iterator fiter = features[pos_i.first][pos_i.second].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const pos_pair_set_type& positions;
    const sample_map_type&   features;
  };
  
  LearnSVM(const size_type __instances) : tolerance(0.1), lambda(C) {}

  void clear()
  {
    if (features.empty())
      features.resize(1);
    if (losses.empty())
      losses.resize(1);
    if (alphas.empty())
      alphas.resize(1);
    if (bounds.empty())
      bounds.resize(1, - std::numeric_limits<double>::infinity());
    
    features.front().clear();
    losses.front().clear();
    alphas.front().clear();
    bounds.front() = - std::numeric_limits<double>::infinity();
  }
  
  std::ostream& encode(std::ostream& os)
  {
    if (features.empty()) return os;
    
    for (size_type id = 1; id != features.size(); ++ id)
      if (! features[id].empty()) {
	
	if (features[id].size() != losses[id].size())
	  throw std::runtime_error("losses size differ");
	
	for (size_type i = 0; i != features[id].size(); ++ i) {
	  utils::encode_base64(losses[id][i], std::ostream_iterator<char>(os));
	  
	  sample_set_type::value_type::const_iterator fiter_end = features[id][i].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[id][i].begin(); fiter != fiter_end; ++ fiter) {
	    os << ' ' << fiter->first << ' ';
	    utils::encode_base64(fiter->second, std::ostream_iterator<char>(os));
	  }
	  
	  os << '\n';
	}
      }
    
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;

    if (features.empty())
      features.resize(1);
    if (losses.empty())
      losses.resize(1);
    if (alphas.empty())
      alphas.resize(1);
    
    std::string line;
    features_type feats;
    
    while (std::getline(is, line)) {
      feats.clear();
      
      const utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
      
      tokenizer_type::iterator iter     = tokenizer.begin();
      tokenizer_type::iterator iter_end = tokenizer.end();
      
      if (iter == iter_end) continue;
      
      const utils::piece loss_str = *iter;
      ++ iter;
      
      if (iter == iter_end) continue;
      
      while (iter != iter_end) {
	const utils::piece feature = *iter;
	++ iter;
	
	if (iter == iter_end) break;
	
	const utils::piece value = *iter;
	++ iter;
	
	feats.push_back(feature_value_type(feature, utils::decode_base64<double>(value)));
      }
      
      if (feats.empty()) continue;

      std::sort(feats.begin(), feats.end());
      
      features.front().insert(feats.begin(), feats.end());
      alphas.front().push_back(0.0);
      losses.front().push_back(utils::decode_base64<double>(loss_str));
    }
    
    return is;
  }

  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    
    if (kbests.empty() || oracles.empty()) return;

    const size_type id_pos = id + 1;
    
    if (id_pos >= features.size())
      features.resize(id_pos + 1);
    if (id_pos >= losses.size())
      losses.resize(id_pos + 1);
    if (id_pos >= alphas.size())
      alphas.resize(id_pos + 1);
    if (id_pos >= bounds.size())
      bounds.resize(id_pos + 1, - std::numeric_limits<double>::infinity());
        
    sentences.clear();
    for (size_t o = 0; o != oracles.size(); ++ o)
      sentences.insert(oracles[o].sentence);
    
    const double error_factor = (error_metric ? - 1.0 : 1.0);
    
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

	const double loss = (oracle.score->score() - kbest.score->score()) * error_factor;
	
	if (loss <= 0.0) continue;
	
	features[id_pos].insert(feats.begin(), feats.end());
	alphas[id_pos].push_back(0.0);
	losses[id_pos].push_back(loss_rank ? 1.0 : loss);
      }
  }
  
  void initialize(weight_set_type& weights)
  {
    
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }

  void clear_history() {}

  struct multiplies_min
  {
    multiplies_min(const double& __factor, const double& __clip) : factor(__factor), clip(__clip) {}
    
    double operator()(const double& x) const
    {
      return std::min(x * factor, clip);
    }
    
    double factor;
    double clip;
  };

  
  
  double learn(weight_set_type& weights)
  {
    positions.clear();
    alpha.clear();
    f.clear();
    
    double alpha_sum = 0.0;
    for (size_type i = 0; i != features.size(); ++ i) {
      
      if (features[i].size() != alphas[i].size())
	throw std::runtime_error("alpha size differ");
      if (features[i].size() != losses[i].size())
	throw std::runtime_error("loss size differ");
      
      for (size_type j = 0; j != features[i].size(); ++ j) {
	positions.push_back(std::make_pair(i, j));
	f.push_back(- losses[i][j]);
	alpha.push_back(alphas[i][j]);
	alpha_sum += alphas[i][j];
      }
    }

    if (alpha_sum != 0.0) {
      const double factor = 1.0 / (alpha_sum * lambda * positions.size());
      std::transform(alpha.begin(), alpha.end(), alpha.begin(), std::bind2nd(std::multiplies<double>(), factor));
      
      // hack for DCD: since we have no summation constraint, we do not have to normalize by positions.size()
      //const double factor = 1.0 / (alpha_sum * lambda);
      //const double clip = 1.0 / (lambda * positions.size());
      //std::transform(alpha.begin(), alpha.end(), alpha.begin(), multiplies_min(factor, clip));
    }
    
    cicada::optimize::QPSMO solver;
    
    HMatrix H(positions, features);
    MMatrix M(positions, features);
    
    solver(alpha, f, H, M, 1.0 / (lambda * positions.size()), tolerance);
    
    weights.clear();
    alpha_set_type::const_iterator aiter = alpha.begin();
    for (size_type i = 0; i != features.size(); ++ i) {
      if (features[i].size() != alphas[i].size())
	throw std::runtime_error("alpha size differ");
      if (features[i].size() != losses[i].size())
	throw std::runtime_error("loss size differ");
      
      for (size_type j = 0; j != features[i].size(); ++ j, ++ aiter) {
	if (*aiter > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[i][j].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[i][j].begin(); fiter != fiter_end; ++ fiter)
	    weights[fiter->first] += (*aiter) * fiter->second;
	}
	
	alphas[i][j] = *aiter;
      }
    }
    
#if 0
    // re-initialize upper bound...
    std::fill(bounds.begin(), bounds.end(), - std::numeric_limits<double>::infinity());
    
    for (size_type i = 0; i != features.size(); ++ i) {
      double& bound = bounds[i];
      
      for (size_type j = 0; j != features[i].size(); ++ j, ++ aiter)
	bound = std::max(bound, losses[i][j] - cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0));
    }
#endif
    
    return 0.0;
  }
  
  double tolerance;
  double lambda;
  
  sample_map_type features;
  loss_map_type   losses;
  alpha_map_type  alphas;
  bound_map_type  bounds;
  
  sentence_unique_type sentences;
  
  pos_pair_set_type positions;
  alpha_set_type alpha;
  f_set_type     f;
};

// Pegasos learner
struct LearnPegasos : public LearnBase
{
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  
  
  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#endif
  
  LearnPegasos(const size_type __instances) : instances(__instances), epoch(0), lambda(C), weight_scale(1.0), weight_norm(0.0) {}
  
  void clear()
  {
    features.clear();
    losses.clear();
  }

  std::ostream& encode(std::ostream& os)
  {
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    return is;
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    
    if (kbests.empty() || oracles.empty()) return;
    
    sentences.clear();
    for (size_t o = 0; o != oracles.size(); ++ o)
      sentences.insert(oracles[o].sentence);

    const double error_factor = (error_metric ? - 1.0 : 1.0);

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
	
	const double loss = (oracle.score->score() - kbest.score->score()) * error_factor;
	
	if (loss <= 0.0) continue;
	
	losses.push_back(loss_rank ? 1.0 : loss);
	features.insert(feats.begin(), feats.end());
      }
  }
  
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

  void clear_history() {}

  double learn(weight_set_type& weights)
  {
    if (features.empty()) return 0.0;
    
    const size_type k = features.size();
    const double k_norm = 1.0 / k;
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + block_size - 1) / block_size;
    const double eta = 0.2 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    rescale(weights, 1.0 - eta * lambda);
    // udpate...
    
    double a_norm = 0.0;
    double pred = 0.0;
    for (size_t i = 0; i != features.size(); ++ i) {
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	double& x = weights[fiter->first];
	const double alpha = eta * k_norm * fiter->second;
	
	a_norm += alpha * alpha;
	pred += 2.0 * x * alpha;
	
	//weight_norm += 2.0 * x * alpha * weight_scale + alpha * alpha;
	x += alpha / weight_scale;
      }
    }
    
    // avoid numerical instability...
    weight_norm += a_norm + pred * weight_scale;
    
    if (weight_norm > 1.0 / lambda)
      rescale(weights, std::sqrt(1.0 / (lambda * weight_norm)));
    
    if (weight_scale < 0.001 || 1000 < weight_scale)
      finalize(weights);
    
    features.clear();
    losses.clear();
    
    return 0.0;
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
  double    weight_scale;
  double    weight_norm;
  

  sample_set_type features;
  loss_set_type   losses;
  
  sentence_unique_type sentences;
};

// optimized-Pegasos learner
struct LearnOPegasos : public LearnBase
{
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  
  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#endif

  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;

  template <typename FeatureSet>
  struct HMatrix
  {
    HMatrix(const FeatureSet& __features) : features(__features) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[i].begin(), features[i].end(), features[j].begin(), features[j].end(), 0.0);
    }
    
    const FeatureSet& features;
  };
  
  template <typename FeatureSet>
  struct MMatrix
  {
    MMatrix(const FeatureSet& __features) : features(__features) {}
    
    template <typename W>
    void operator()(W& w, const alpha_type& alpha) const
    {
      const size_type model_size = features.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename W>
    double operator()(const W& w, const size_t& i) const
    {
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename W>
    void operator()(W& w, const double& update, const size_t& i) const
    {
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const FeatureSet& features;
  };
  
  LearnOPegasos(const size_type __instances) : instances(__instances), epoch(0), tolerance(0.1), lambda(C), weight_scale(1.0), weight_norm(0.0) {}
  
  void clear()
  {
    features.clear();
    losses.clear();
  }

  std::ostream& encode(std::ostream& os)
  {
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    return is;
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    
    if (kbests.empty() || oracles.empty()) return;
    
    sentences.clear();
    for (size_t o = 0; o != oracles.size(); ++ o)
      sentences.insert(oracles[o].sentence);

    const double error_factor = (error_metric ? - 1.0 : 1.0);

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
	
	const double loss = (oracle.score->score() - kbest.score->score()) * error_factor;
	
	if (loss <= 0.0) continue;
	
	losses.push_back(loss_rank ? 1.0 : loss);
	features.insert(feats.begin(), feats.end());
      }
  }
  
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

  void clear_history() {}

  double learn(weight_set_type& weights)
  {
    if (features.empty()) return 0.0;
    
    const size_type k = features.size();
    const double k_norm = 1.0 / k;
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + block_size - 1) / block_size;
    const double eta = 0.2 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    rescale(weights, 1.0 - eta * lambda);
    // udpate...
    
    alpha.clear();
    f.clear();
    
    alpha.reserve(losses.size());
    f.reserve(losses.size());
    
    alpha.resize(losses.size(), 0.0);
    f.resize(losses.size(), 0.0);
    
    double objective = 0.0;
    for (size_t i = 0; i != losses.size(); ++ i) {
      f[i] = - (losses[i] - cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0) * weight_scale);
      objective -= f[i];
    }
    objective /= losses.size();
    
    cicada::optimize::QPDCD solver;
    
    HMatrix<sample_set_type> H(features);
    MMatrix<sample_set_type> M(features);
    
    solver(alpha, f, H, M, eta, tolerance);
    
    double a_norm = 0.0;
    double pred = 0.0;
    for (size_t i = 0; i != losses.size(); ++ i)
      if (alpha[i] > 0.0) {
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	  double& x = weights[fiter->first];
	  const double a = alpha[i] * fiter->second;
	  
	  a_norm += a * a;
	  pred += 2.0 * x * a;
	  
	  //weight_norm += 2.0 * x * a * weight_scale + a * a;
	  x += a / weight_scale;
	}
      }
    
    // avoid numerical instability...
    weight_norm += a_norm + pred * weight_scale;
    
    if (weight_norm > 1.0 / lambda)
      rescale(weights, std::sqrt(1.0 / (lambda * weight_norm)));
    
    if (weight_scale < 0.001 || 1000 < weight_scale)
      finalize(weights);
    
    features.clear();
    losses.clear();
    
    return 0.0;
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
  double    tolerance;
  double    lambda;
  double    weight_scale;
  double    weight_norm;
  

  sample_set_type features;
  loss_set_type   losses;
  
  sentence_unique_type sentences;
  
  alpha_type    alpha;
  f_type        f;
};

// MIRA learner
// We will run a qp solver and determine the alpha, then, translate this into w
struct LearnMIRA : public LearnBase
{
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  
  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#endif

  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    f_type;

  template <typename FeatureSet>
  struct HMatrix
  {
    HMatrix(const FeatureSet& __features) : features(__features) {}
    
    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[i].begin(), features[i].end(), features[j].begin(), features[j].end(), 0.0);
    }
    
    const FeatureSet& features;
  };
  
  template <typename FeatureSet>
  struct MMatrix
  {
    MMatrix(const FeatureSet& __features) : features(__features) {}
    
    template <typename W>
    void operator()(W& w, const alpha_type& alpha) const
    {
      const size_type model_size = features.size();
      
      for (size_type i = 0; i != model_size; ++ i)
	if (alpha[i] > 0.0) {
	  sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	  for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	    w[fiter->first] += alpha[i] * fiter->second;
	}
    }
    
    template <typename W>
    double operator()(const W& w, const size_t& i) const
    {
      double dot = 0.0;
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	dot += w[fiter->first] * fiter->second;
      return dot;
    }
    
    template <typename W>
    void operator()(W& w, const double& update, const size_t& i) const
    {
      sample_set_type::value_type::const_iterator fiter_end = features[i].end();
      for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) 
	w[fiter->first] += update * fiter->second;
    }
    
    const FeatureSet& features;
  };
  
  LearnMIRA(const size_type __instances) : tolerance(0.1), lambda(C) {}
  
  void clear()
  {
    features.clear();
    losses.clear();
  }
  
  std::ostream& encode(std::ostream& os)
  {
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    return is;
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    
    if (kbests.empty() || oracles.empty()) return;
    
    sentences.clear();
    for (size_t o = 0; o != oracles.size(); ++ o)
      sentences.insert(oracles[o].sentence);

    const double error_factor = (error_metric ? - 1.0 : 1.0);

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
	
	const double loss = (oracle.score->score() - kbest.score->score()) * error_factor;
	
	if (loss <= 0.0) continue;
	
	losses.push_back(loss_rank ? 1.0 : loss);
	features.insert(feats.begin(), feats.end());
      }
  }
  
  void initialize(weight_set_type& weights)
  {
    
    
  }
  
  void finalize(weight_set_type& weights)
  {
    
  }

  void clear_history() {}
  
  double learn(weight_set_type& weights)
  {

    if (features.empty()) return 0.0;
    
    alpha.clear();
    f.clear();
    
    alpha.reserve(losses.size());
    f.reserve(losses.size());
    
    alpha.resize(losses.size(), 0.0);
    f.resize(losses.size(), 0.0);
    
    double objective = 0.0;
    for (size_t i = 0; i != losses.size(); ++ i) {
      f[i] = - (losses[i] - cicada::dot_product(features[i].begin(), features[i].end(), weights, 0.0));
      objective -= f[i];
    }
    
    objective /= losses.size();
    
    cicada::optimize::QPDCD solver;
    
    HMatrix<sample_set_type> H(features);
    MMatrix<sample_set_type> M(features);
    
    solver(alpha, f, H, M, 1.0 / (lambda * losses.size()), tolerance);
    
    for (size_t i = 0; i != losses.size(); ++ i)
      if (alpha[i] > 0.0) {
	// update: weights[fiter->first] += alpha[i] * fiter->second;
	
	sample_set_type::value_type::const_iterator fiter_end = features[i].end();
	for (sample_set_type::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter)
	  weights[fiter->first] += alpha[i] * fiter->second;
      }
    
    features.clear();
    losses.clear();
    
    return objective;
  }
  
  double tolerance;
  double lambda;

  sample_set_type features;
  loss_set_type   losses;
  
  sentence_unique_type sentences;
  
  alpha_type    alpha;
  f_type        f;
};

// logistic regression base...
struct LearnLR : public LearnBase
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

    void encode(const hypothesis_set_type& __kbests, const hypothesis_set_type& __oracles, const bool error_metric)
    {
      if (__kbests.empty() || __oracles.empty()) return;
      
      const double error_factor = (error_metric ? 1.0 : - 1.0);
      
      hypothesis_set_type::const_iterator kiter_end = __kbests.end();
      for (hypothesis_set_type::const_iterator kiter = __kbests.begin(); kiter != kiter_end; ++ kiter) {
	kbests.insert(kiter->features.begin(), kiter->features.end());
	loss_kbests.push_back(kiter->score->score() * error_factor);
      }
      
      hypothesis_set_type::const_iterator oiter_end = __oracles.end();
      for (hypothesis_set_type::const_iterator oiter = __oracles.begin(); oiter != oiter_end; ++ oiter) {
	oracles.insert(oiter->features.begin(), oiter->features.end());
	loss_oracles.push_back(oiter->score->score() * error_factor);
      }
    }
    
    template <typename Expectations>
    double encode(const weight_set_type& weights, Expectations& expectations, const double scale) const
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
      
      weight_type Z_oracle;
      weight_type Z_kbest;
      
      const double cost_factor = (softmax_margin ? 1.0 : 0.0);
      
      for (size_type o = 0; o != oracles.size(); ++ o)
	Z_oracle += traits_type::exp(cicada::dot_product(weights, oracles[o].begin(), oracles[o].end(), 0.0) * scale + cost_factor * loss_oracles[o]);
      
      for (size_type k = 0; k != kbests.size(); ++ k)
	Z_kbest += traits_type::exp(cicada::dot_product(weights, kbests[k].begin(), kbests[k].end(), 0.0) * scale + cost_factor * loss_kbests[k]);
      
      for (size_type o = 0; o != oracles.size(); ++ o) {
	const weight_type weight = traits_type::exp(cicada::dot_product(weights, oracles[o].begin(), oracles[o].end(), 0.0) * scale + cost_factor * loss_oracles[o]) / Z_oracle;
	
	sample_set_type::value_type::const_iterator fiter_end = oracles[o].end();
	for (sample_set_type::value_type::const_iterator fiter = oracles[o].begin(); fiter != fiter_end; ++ fiter)
	  expectations[fiter->first] -= weight_type(fiter->second) * weight;
      }
      
      for (size_type k = 0; k != kbests.size(); ++ k) {
	const weight_type weight = traits_type::exp(cicada::dot_product(weights, kbests[k].begin(), kbests[k].end(), 0.0) * scale + cost_factor * loss_kbests[k]) / Z_kbest;	
	
	sample_set_type::value_type::const_iterator fiter_end = kbests[k].end();
	for (sample_set_type::value_type::const_iterator fiter = kbests[k].begin(); fiter != fiter_end; ++ fiter)
	  expectations[fiter->first] += weight_type(fiter->second) * weight;
      }
      
      return log(Z_oracle) - log(Z_kbest);
    }
  };
};

// SGDL1 learner
struct LearnSGDL1 : public LearnLR
{
  typedef utils::chunk_vector<sample_pair_type, 4096 / sizeof(sample_pair_type), std::allocator<sample_pair_type> > sample_pair_set_type;
  
  typedef cicada::WeightVector<double> penalty_set_type;
  
  // maximize...
  // L_w =  \sum \log p(y | x) - C |w|
  // 
  
  LearnSGDL1(const size_type __instances) : instances(__instances), epoch(0), lambda(C), penalties(), penalty(0.0) {}
  
  void clear()
  {
    samples.clear();
  }
  
  std::ostream& encode(std::ostream& os)
  {
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    return is;
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    if (kbests.empty() || oracles.empty()) return;
    
    samples.push_back(sample_pair_type());
    samples.back().encode(kbests, oracles, error_metric);
  }
  
  void initialize(weight_set_type& weights)
  {

  }

  void finalize(weight_set_type& weights)
  {
    
  }
  
  void clear_history() {}

  double learn(weight_set_type& weights)
  {
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > expectation_type;

    if (samples.empty()) return 0.0;
    
    //
    // C = lambda * N
    //
    
    const size_type k = samples.size();
    const double k_norm = 1.0 / k;
    //const double eta = 1.0 / (lambda * (epoch + 2)); // this is an eta from pegasos
    const size_type num_samples = (instances + block_size - 1) / block_size;
    const double eta = 0.2 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    penalty += eta * lambda * k_norm;
    
    expectation_type expectations;
    
    double objective = 0.0;
    sample_pair_set_type::const_iterator siter_end = samples.end();
    for (sample_pair_set_type::const_iterator siter = samples.begin(); siter != siter_end; ++ siter)
      objective += siter->encode(weights, expectations, 1.0);
    objective /= samples.size();
    
    // update by expectations...
    expectation_type::const_iterator eiter_end = expectations.end();
    for (expectation_type::const_iterator eiter = expectations.begin(); eiter != eiter_end; ++ eiter) {
      double& x = weights[eiter->first];
      
      // update weight ... we will update "minus" value
      x += - static_cast<double>(eiter->second) * eta * k_norm;
      
      // apply penalty
      apply(x, penalties[eiter->first], penalty);
    }
    
    samples.clear();
    
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
  
  sample_pair_set_type samples;

  size_type instances;
  
  size_type epoch;
  double    lambda;
  
  penalty_set_type penalties;
  double penalty;
};

// SGDL2 learner
struct LearnSGDL2 : public LearnLR
{
  typedef utils::chunk_vector<sample_pair_type, 4096 / sizeof(sample_pair_type), std::allocator<sample_pair_type> > sample_pair_set_type;
    
  LearnSGDL2(const size_type __instances) : instances(__instances), epoch(0), lambda(C), weight_scale(1.0), weight_norm(0.0) {}
  
  void clear()
  {
    samples.clear();
  }
  
  std::ostream& encode(std::ostream& os)
  {
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    return is;
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    if (kbests.empty() || oracles.empty()) return;
    
    samples.push_back(sample_pair_type());
    samples.back().encode(kbests, oracles, error_metric);
  }
  
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

  void clear_history() {}

  double learn(weight_set_type& weights)
  {
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > expectation_type;
    
    if (samples.empty()) return 0.0;
    
    const size_type k = samples.size();
    const double k_norm = 1.0 / k;
    //const double eta = 1.0 / (lambda * (epoch + 2));  // this is an eta from pegasos
    const size_type num_samples = (instances + block_size - 1) / block_size;
    const double eta = 0.2 * std::pow(0.85, double(epoch) / num_samples); // eta from SGD-L1
    ++ epoch;
    
    rescale(weights, 1.0 - eta * lambda);
    
    expectation_type expectations;
    
    // update... by eta / k
    double objective = 0.0;
    sample_pair_set_type::const_iterator siter_end = samples.end();
    for (sample_pair_set_type::const_iterator siter = samples.begin(); siter != siter_end; ++ siter)
      objective += siter->encode(weights, expectations, weight_scale);
    objective /= samples.size();
    
    // update by expectations...
    double a_norm = 0.0;
    double pred = 0.0;
    expectation_type::const_iterator eiter_end = expectations.end();
    for (expectation_type::const_iterator eiter = expectations.begin(); eiter != eiter_end; ++ eiter) {
      // we will update "minus" value...
      
      double& x = weights[eiter->first];
      const double alpha = - static_cast<double>(eiter->second) * eta * k_norm;
      
      a_norm += alpha * alpha;
      pred += 2.0 * x * alpha;
      
      //weight_norm += 2.0 * x * alpha * weight_scale + alpha * alpha;
      x += alpha / weight_scale;
    }
    
    // avoid numerical instability...
    weight_norm += a_norm + pred * weight_scale;
    
    if (weight_norm > 1.0 / lambda)
      rescale(weights, std::sqrt(1.0 / (lambda * weight_norm)));
    
    if (weight_scale < 0.001 || 1000 < weight_scale)
      finalize(weights);
    
    // clear current training events..
    samples.clear();
    
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
  
  sample_pair_set_type samples;

  size_type instances;
  
  size_type epoch;
  double    lambda;
  
  double weight_scale;
  double weight_norm;
};

// LBFGS learner
struct LearnLBFGS : public LearnLR
{
  typedef utils::chunk_vector<sample_pair_type, 4096 / sizeof(sample_pair_type), std::allocator<sample_pair_type> > sample_pair_set_type;
  typedef utils::chunk_vector<sample_pair_set_type, 4096 / sizeof(sample_pair_set_type), std::allocator<sample_pair_set_type> > sample_pair_map_type;
  
  typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;

  LearnLBFGS(const size_type __instances) {}
  
  void clear()
  {
    samples_other.clear();
  }
  
  std::ostream& encode(std::ostream& os)
  {
    for (size_t id = 0; id != samples.size(); ++ id) 
      if (! samples[id].empty()) {
	sample_pair_set_type::const_iterator siter_end = samples[id].end();
	for (sample_pair_set_type::const_iterator siter = samples[id].begin(); siter != siter_end; ++ siter) {
	  const sample_pair_type& sample = *siter;
	  
	  if (sample.oracles.empty() || sample.kbests.empty()) continue;
	  
	  for (size_type o = 0; o != sample.oracles.size(); ++ o) {
	    os << "oracle: ";
	    utils::encode_base64(sample.loss_oracles[o], std::ostream_iterator<char>(os));
	    
	    sample_set_type::value_type::const_iterator fiter_end = sample.oracles[o].end();
	    for (sample_set_type::value_type::const_iterator fiter = sample.oracles[o].begin(); fiter != fiter_end; ++ fiter) {
	      os << ' ' << fiter->first << ' ';
	      utils::encode_base64(fiter->second, std::ostream_iterator<char>(os));
	    }
	    os << '\n';
	  }
	  
	  for (size_type k = 0; k != sample.kbests.size(); ++ k) {
	    os << "kbest: ";
	    utils::encode_base64(sample.loss_kbests[k], std::ostream_iterator<char>(os));
	    
	    sample_set_type::value_type::const_iterator fiter_end = sample.kbests[k].end();
	    for (sample_set_type::value_type::const_iterator fiter = sample.kbests[k].begin(); fiter != fiter_end; ++ fiter) {
	      os << ' ' << fiter->first << ' ';
	      utils::encode_base64(fiter->second, std::ostream_iterator<char>(os));
	    }
	    os << '\n';
	  }
	}
      }
    
    return os;
  }

  std::istream& decode(std::istream& is)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > feature_set_type;
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
    
    std::string mode = "kbest:";

    sample_set_type* psample;
    loss_set_type*   ploss;
    
    std::string line;
    feature_set_type features;
    while (std::getline(is, line)) {
      features.clear();
      
      const utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
      
      tokenizer_type::iterator iter     = tokenizer.begin();
      tokenizer_type::iterator iter_end = tokenizer.end();
      
      if (iter == iter_end) continue;
      
      const utils::piece mode_curr = *iter;
      ++ iter;
      
      if (iter == iter_end) continue;
      
      if (mode_curr != mode) {
	if (mode_curr == "oracle:") {
	  samples_other.push_back(sample_pair_type());
	  psample = &(samples_other.back().oracles);
	  ploss = &(samples_other.back().loss_oracles);
	} else {
	  psample = &(samples_other.back().kbests);
	  ploss = &(samples_other.back().loss_kbests);
	}
	
	mode = mode_curr;
      }

      const utils::piece loss_str = *iter;
      ++ iter;
      
      while (iter != iter_end) {
	const utils::piece feature = *iter;
	++ iter;
	
	if (iter == iter_end) break;
	
	const utils::piece value = *iter;
	++ iter;
	
	features.push_back(feature_value_type(feature, utils::decode_base64<double>(value)));
      }
      
      if (features.empty()) continue;
      
      psample->insert(features.begin(), features.end());
      ploss->push_back(utils::decode_base64<double>(loss_str));
    }
    
    return is;
  }
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    if (kbests.empty() || oracles.empty()) return;
    
    if (id >= samples.size())
      samples.resize(id + 1);
    
    samples[id].push_back(sample_pair_type());
    
    samples[id].back().encode(kbests, oracles, error_metric);
  }
  
  void initialize(weight_set_type& weights)
  {

  }

  void finalize(weight_set_type& weights)
  {
    
  }

  void clear_history() {}

  double learn(weight_set_type& __weights)
  {
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    
    if (regularize_l1) {
      param.orthantwise_c = C;
      param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    } else
      param.orthantwise_c = 0.0;
    
    double objective = 0.0;
    
    weights = __weights;
    weights.allocate();
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, LearnLBFGS::evaluate, 0, this, &param);
    
    __weights = weights;
    
    return objective;
  }
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    LearnLBFGS& optimizer = *((LearnLBFGS*) instance);
    
    expectation_type& expectations            = optimizer.expectations;
    weight_set_type& weights                  = optimizer.weights;
    const sample_pair_map_type& samples       = optimizer.samples;
    const sample_pair_set_type& samples_other = optimizer.samples_other;
    
    expectations.clear();
    expectations.allocate();

    double objective = 0.0;
    size_t instances = 0;

    for (size_t i = 0; i != samples_other.size(); ++ i) {
      const sample_pair_type& sample = samples_other[i];
      
      const double margin = sample.encode(weights, expectations, 1.0);
      objective -= margin;
      ++ instances;
    }
    
    for (size_t id = 0; id != samples.size(); ++ id) 
      if (! samples[id].empty()) {
	sample_pair_set_type::const_iterator siter_end = samples[id].end();
	for (sample_pair_set_type::const_iterator siter = samples[id].begin(); siter != siter_end; ++ siter) {
	  const sample_pair_type& sample = *siter;
	  
	  const double margin = sample.encode(weights, expectations, 1.0);
	  objective -= margin;
	  ++ instances;
	}
      }
    
    std::copy(expectations.begin(), expectations.begin() + n, g);
    
    objective /= instances;
    std::transform(g, g + n, g, std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    
    // L2...
    if (regularize_l2) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i) {
	g[i] += C * x[i];
	norm += x[i] * x[i];
      }
      objective += 0.5 * C * norm;
    }
    
    return objective;
  }
  
  expectation_type     expectations;
  weight_set_type      weights;
  sample_pair_map_type samples;
  sample_pair_set_type samples_other;
};

// linear learner
struct LearnLinear
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef struct model        model_type;
  typedef struct parameter    parameter_type;
  typedef struct problem      problem_type;
  typedef struct feature_node feature_node_type;
  
  typedef size_t offset_type;
  
  typedef std::vector<feature_node_type*, std::allocator<feature_node_type*> > feature_node_map_type;
  typedef std::vector<int, std::allocator<int> > label_set_type;

  //
  // typedef for unique sentences
  //
  struct hash_sentence : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    size_t operator()(const hypothesis_type::sentence_type& x) const
    {
      return hasher_type()(x.begin(), x.end(), 0);
    }
  };
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> > sentence_unique_type;
#endif
  
  static void print_string_stderr(const char *s)
  {
    std::cerr << s << std::flush;
  }

  static void print_string_none(const char *s)
  {
    
  }

  struct Encoder
  {
    typedef std::vector<offset_type, std::allocator<offset_type> > offset_set_type;
    typedef std::vector<feature_node_type, std::allocator<feature_node_type> > feature_node_set_type;

    Encoder() : offsets(), features() {}

    offset_set_type       offsets;
    feature_node_set_type features;
    
    template <typename Iterator>
    void encode(Iterator first, Iterator last)
    {
      feature_node_type feature;
      
      if (first == last) return;
      
      offsets.push_back(features.size());
      
      for (/**/; first != last; ++ first) {
	feature.index = first->first.id() + 1;
	feature.value = first->second;
	
	features.push_back(feature);
      }

      if (offsets.back() == features.size())
	offsets.pop_back();
      else {
	// termination...
	feature.index = -1;
	feature.value = 0.0;
	features.push_back(feature);
      }
    }

    void encode(const hypothesis_set_type& kbests,
		const hypothesis_set_type& oracles)
    {
      feature_node_type feature;
      
      sentence_unique_type sentences;
      for (size_t o = 0; o != oracles.size(); ++ o)
	sentences.insert(oracles[o].sentence);
      
      for (size_t o = 0; o != oracles.size(); ++ o)
	for (size_t k = 0; k != kbests.size(); ++ k) {
	  const hypothesis_type& oracle = oracles[o];
	  const hypothesis_type& kbest  = kbests[k];
	  
	  if (sentences.find(kbest.sentence) != sentences.end()) continue;
	  
	  offsets.push_back(features.size());
	  
	  hypothesis_type::feature_set_type::const_iterator oiter = oracle.features.begin();
	  hypothesis_type::feature_set_type::const_iterator oiter_end = oracle.features.end();
	
	  hypothesis_type::feature_set_type::const_iterator kiter = kbest.features.begin();
	  hypothesis_type::feature_set_type::const_iterator kiter_end = kbest.features.end();
	
	  while (oiter != oiter_end && kiter != kiter_end) {
	    if (oiter->first < kiter->first) {
	      feature.index = oiter->first.id() + 1;
	      feature.value = oiter->second;
	      features.push_back(feature);
	      ++ oiter;
	    } else if (kiter->first < oiter->first) {
	      feature.index = kiter->first.id() + 1;
	      feature.value = - kiter->second;
	      features.push_back(feature);
	      ++ kiter;
	    } else {
	      feature.index = oiter->first.id() + 1;
	      feature.value = oiter->second - kiter->second;
	      if (feature.value != 0.0)
		features.push_back(feature);
	      ++ oiter;
	      ++ kiter;
	    }
	  }
	    
	  for (/**/; oiter != oiter_end; ++ oiter) {
	    feature.index = oiter->first.id() + 1;
	    feature.value = oiter->second;
	    features.push_back(feature);
	  }
	  
	  for (/**/; kiter != kiter_end; ++ kiter) {
	    feature.index = kiter->first.id() + 1;
	    feature.value = - kiter->second;
	    features.push_back(feature);
	  }
	  
	  if (offsets.back() == features.size())
	    offsets.pop_back();
	  else {
	    // termination...
	    feature.index = -1;
	    feature.value = 0.0;
	    features.push_back(feature);
	  }
	}
    }
    
    void clear()
    {
      offsets.clear();
      features.clear();
    }
    
    void shrink()
    {
      offset_set_type(offsets).swap(offsets);
      feature_node_set_type(features).swap(features);
    }
  };
  typedef Encoder encoder_type;
  typedef utils::chunk_vector<encoder_type, 4096 / sizeof(encoder_type), std::allocator<encoder_type> > encoder_set_type; 

  LearnLinear(const size_type __instances) {}
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool error_metric=false)
  {
    if (id >= encoders.size())
      encoders.resize(id + 1);
    
    encoders[id].encode(kbests, oracles);
  }
  
  void clear()
  {
    encoder_other.clear();
  }
  
  std::ostream& encode(std::ostream& os)
  {
    for (size_type id = 0; id != encoders.size(); ++ id)
      for (size_type pos = 0; pos != encoders[id].offsets.size(); ++ pos) {
	encoder_type::feature_node_set_type::const_iterator fiter     = encoders[id].features.begin() + encoders[id].offsets[pos];
	encoder_type::feature_node_set_type::const_iterator fiter_end = (pos + 1 < encoders[id].offsets.size()
									 ? encoders[id].features.begin() + encoders[id].offsets[pos + 1] - 1
									 : encoders[id].features.end() - 1);
	for (/**/; fiter != fiter_end; ++ fiter) {
	  os << weight_set_type::feature_type(fiter->index - 1) << ' ';
	  utils::encode_base64(fiter->value, std::ostream_iterator<char>(os));
	  os << ' ';
	}
	os << '\n';
      }
    return os;
  }
  
  std::istream& decode(std::istream& is)
  {
    typedef hypothesis_type::feature_value_type feature_value_type;
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > feature_set_type;
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;

    std::string line;
    feature_set_type features;
    while (std::getline(is, line)) {
      features.clear();
      
      const utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
      
      tokenizer_type::iterator iter     = tokenizer.begin();
      tokenizer_type::iterator iter_end = tokenizer.end();
      
      while (iter != iter_end) {
	const utils::piece feature = *iter;
	++ iter;
	
	if (iter == iter_end) break;
	
	const utils::piece value = *iter;
	++ iter;
	
	features.push_back(feature_value_type(feature, utils::decode_base64<double>(value)));
      }

      if (features.empty()) continue;
      
      std::sort(features.begin(), features.end());
      encoder_other.encode(features.begin(), features.end());
    }
    return is;
  }
  
  void initialize(weight_set_type& weights)
  {

  }

  void finalize(weight_set_type& weights)
  {
    
  }

  void clear_history() {}

  double learn(weight_set_type& weights)
  {
    size_type data_size = encoder_other.offsets.size();
    for (size_type id = 0; id != encoders.size(); ++ id)
      data_size += encoders[id].offsets.size();
    
    label_set_type        labels(data_size, 1);
    feature_node_map_type features;
    features.reserve(data_size);

    for (size_type pos = 0; pos != encoder_other.offsets.size(); ++ pos)
      features.push_back(const_cast<feature_node_type*>(&(*encoder_other.features.begin())) + encoder_other.offsets[pos]);
    
    for (size_type id = 0; id != encoders.size(); ++ id)
      for (size_type pos = 0; pos != encoders[id].offsets.size(); ++ pos)
	features.push_back(const_cast<feature_node_type*>(&(*encoders[id].features.begin())) + encoders[id].offsets[pos]);
    
    problem_type problem;
    problem.l = labels.size();
    problem.n = feature_type::allocated();
    problem.y = &(*labels.begin());
    problem.x = &(*features.begin());
    problem.bias = -1;
    
    parameter_type parameter;
    parameter.solver_type = linear_solver;
    parameter.eps = eps;
    parameter.C = 1.0 / (C * labels.size()); // renormalize!
    parameter.nr_weight    = 0;
    parameter.weight_label = 0;
    parameter.weight       = 0;
    
    if (parameter.eps == std::numeric_limits<double>::infinity()) {
      if (parameter.solver_type == L2R_LR || parameter.solver_type == L2R_L2LOSS_SVC)
	parameter.eps = 0.01;
      else if (parameter.solver_type == L2R_L2LOSS_SVC_DUAL || parameter.solver_type == L2R_L1LOSS_SVC_DUAL || parameter.solver_type == MCSVM_CS || parameter.solver_type == L2R_LR_DUAL)
	parameter.eps = 0.1;
      else if (parameter.solver_type == L1R_L2LOSS_SVC || parameter.solver_type == L1R_LR)
	parameter.eps = 0.01;
    }
    
    if (debug >= 2)
      set_print_string_function(print_string_stderr);
    else
      set_print_string_function(print_string_none);
    
    const char* error_message = check_parameter(&problem, &parameter);
    if (error_message)
      throw std::runtime_error(std::string("error: ") + error_message);
    
    static const char* names[] = {"L2R_LR", "L2R_L2LOSS_SVC_DUAL", "L2R_L2LOSS_SVC", "L2R_L1LOSS_SVC_DUAL", "MCSVM_CS",
				  "L1R_L2LOSS_SVC", "L1R_LR", "L2R_LR_DUAL"};
    
    if (debug)
      std::cerr << "solver: " << names[parameter.solver_type] << std::endl;
    
    const model_type* model = train(&problem, &parameter);
    
    const double objective = model->objective * C;
    
    // it is an optimization...
    weights.clear();
    for (int j = 0; j != model->nr_feature; ++ j)
      weights[weight_set_type::feature_type(j)] = model->w[j];
    
    free_and_destroy_model(const_cast<model_type**>(&model));
    
    return objective;
  }
  
  encoder_set_type encoders;
  encoder_type     encoder_other;
};


struct KBestSentence
{
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::sentence_feature_traversal   traversal_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  typedef cicada::operation::kbest_sentence_filter_unique filter_type;
  
  typedef cicada::KBest<traversal_type, function_type, filter_type> derivation_set_type;

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
    
    derivation_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_type(graph));
    
    derivation_type derivation;
    weight_type     weight;
    
    for (int k = 0; k != kbest_size && derivations(k, derivation, weight); ++ k)
      kbests.push_back(hypothesis_type(boost::get<0>(derivation).begin(), boost::get<0>(derivation).end(),
				       boost::get<1>(derivation).begin(), boost::get<1>(derivation).end()));
  }
};

struct KBestAlignment
{
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::alignment_feature_traversal  traversal_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  typedef cicada::operation::kbest_alignment_filter_unique filter_type;
  
  typedef cicada::KBest<traversal_type, function_type, filter_type> derivation_set_type;

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
    
    derivation_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_type(graph));
    
    sentence_type   sentence;
    derivation_type derivation;
    weight_type     weight;
    
    for (int k = 0; k != kbest_size && derivations(k, derivation, weight); ++ k) {
      std::ostringstream os;
      os << boost::get<0>(derivation);
      sentence.assign(os.str());
      
      kbests.push_back(hypothesis_type(sentence.begin(), sentence.end(),
				       boost::get<1>(derivation).begin(), boost::get<1>(derivation).end()));
    }
  }
};

struct KBestDependency
{
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::dependency_feature_traversal traversal_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  typedef cicada::operation::kbest_dependency_filter_unique filter_type;
  
  typedef cicada::KBest<traversal_type, function_type, filter_type> derivation_set_type;

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
    
    derivation_set_type derivations(graph, kbest_size, traversal_type(), function_type(weights), filter_type(graph));
    
    sentence_type   sentence;
    derivation_type derivation;
    weight_type     weight;
    
    for (int k = 0; k != kbest_size && derivations(k, derivation, weight); ++ k) {
      std::ostringstream os;
      os << boost::get<0>(derivation);
      sentence.assign(os.str());
      
      kbests.push_back(hypothesis_type(sentence.begin(), sentence.end(),
				       boost::get<1>(derivation).begin(), boost::get<1>(derivation).end()));
    }
  }
};


struct Oracle
{
  typedef std::vector<const hypothesis_type*, std::allocator<const hypothesis_type*> > oracle_set_type;
  typedef std::vector<oracle_set_type, std::allocator<oracle_set_type> > oracle_map_type;

  template <typename Generator>
  std::pair<score_ptr_type, score_ptr_type>
  operator()(const hypothesis_map_type& kbests, const scorer_document_type& scorers, hypothesis_map_type& oracles, Generator& generator)
  {
    typedef std::vector<size_t, std::allocator<size_t> > id_set_type;

    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);

    score_ptr_type score_1best;
        
    score_ptr_type score_prev;
    score_ptr_type score_best;
    score_ptr_type score_curr;
    
    oracle_map_type oracles_prev(kbests.size());
    oracle_map_type oracles_best(kbests.size());
    oracle_map_type oracles_curr(kbests.size());

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
	
	if (! score_prev)
	  score_prev = kbests[id].front().score->clone();
	else
	  *score_prev += *(kbests[id].front().score);
	
	oracles_prev[id].push_back(&kbests[id].front());
      }
    
    if (score_prev)
      score_best = score_prev->clone();
    else
      throw std::runtime_error("no scores?");
    
    score_1best = score_best->clone();
    
    oracles_best = oracles_prev;
    
    double objective_prev = score_prev->score() * score_factor;
    double objective_best = objective_prev;
    double objective_curr = objective_prev;
    
    // 
    // 10 iteration will be fine
    //
    for (int i = 0; i < 10; ++ i) {
      score_curr     = score_prev->clone();
      objective_curr = objective_prev;
      oracles_curr   = oracles_prev;

      for (id_set_type::const_iterator iiter = ids.begin(); iiter != ids.end(); ++ iiter) {
	const size_t id = *iiter;
	
	score_ptr_type score = score_curr->clone();
	*score -= *(oracles_curr[id].front()->score);
	  
	hypothesis_set_type::const_iterator hiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  score_ptr_type score_sample = score->clone();
	  *score_sample += *(hiter->score);
	    
	  const double objective_sample = score_sample->score() * score_factor;
	    
	  if (objective_sample > objective_curr) {
	    oracles_curr[id].clear();
	    oracles_curr[id].push_back(&(*hiter));
	      
	    objective_curr = objective_sample;
	    score_curr     = score_sample;
	  } else if (objective_sample == objective_curr)
	    oracles_curr[id].push_back(&(*hiter));
	}
      }
      
      if (objective_curr > objective_best) {
	score_best     = score_curr->clone();
	objective_best = objective_curr;
	oracles_best   = oracles_curr;
      }
      
      if (objective_curr <= objective_prev) break;
      
      score_prev     = score_curr->clone();
      objective_prev = objective_curr;
      oracles_prev   = oracles_curr;
      
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
		 const size_t shard_rank,
		 const size_t shard_size)
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

inline
void read_samples(const path_type& input_path,
		  sample_set_type& samples,
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
	std::getline(is, line);
	
	if (line.empty()) continue;
	
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (id != i)
	  throw std::runtime_error("id doest not match!");
	
	const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	
	if (id_rank >= samples.size())
	  samples.resize(id_rank + 1);
	
	samples[id_rank] = line;
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
	  
	  if (id_rank >= samples.size())
	    samples.resize(id_rank + 1);
	  
	  samples[id_rank] = line;
	}
      }
  } else {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    std::string line;
    for (size_t id = 0; std::getline(is, line); ++ id) 
      if (shard_size == 0 || id % shard_size == shard_rank) {
	if (! line.empty())
	  samples.push_back(utils::lexical_cast<std::string>(id) + " ||| " + line);
	else
	  samples.push_back(std::string());
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
