//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_BLOCK_IMPL__HPP__
#define __CICADA_LEARN_BLOCK_IMPL__HPP__ 1

#include <sstream>
#include <vector>

#include "cicada_kbest_impl.hpp"

#include "cicada/semiring.hpp"
#include "cicada/eval.hpp"

#include "cicada/kbest.hpp"
#include "cicada/operation/traversal.hpp"
#include "cicada/operation/functional.hpp"

#include "utils/sgi_hash_set.hpp"

#include "lbfgs.h"
#include "liblinear/linear.h"

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

// LBFGS learner
struct LearnLBFGS
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef hypothesis_type::feature_set_type feature_set_type;
  typedef std::vector<feature_set_type, std::allocator<feature_set_type> > sample_type;

  struct sample_pair_type
  {
    sample_pair_type() : kbests(), oracles() {}
    
    sample_type kbests;
    sample_type oracles;
  };
  
  typedef std::vector<sample_pair_type, std::allocator<sample_pair_type> > sample_pair_set_type;
  typedef utils::chunk_vector<sample_pair_set_type, 4096 / sizeof(sample_pair_set_type), std::allocator<sample_pair_set_type> > sample_pair_map_type;
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool merge=false)
  {
    if (id >= samples.size())
      samples.resize(id + 1);
    
    if (! merge)
      samples[id].clear();
    
    samples[id].push_back(sample_pair_type());
    
    samples[id].back().kbests.reserve(kbests.size());
    samples[id].back().oracles.reserve(oracles.size());
    
    hypothesis_set_type::const_iterator kiter_end = kbests.end();
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter)
      samples[id].back().kbests.push_back(kiter->features);
    
    hypothesis_set_type::const_iterator oiter_end = oracles.end();
    for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter)
      samples[id].back().oracles.push_back(oiter->features);
  }
  
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

  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    LearnLBFGS& optimizer = *((LearnLBFGS*) instance);
    
    expectation_type& expectations = optimizer.expectations;
    weight_set_type& weights       = optimizer.weights;
    sample_pair_map_type& samples  = optimizer.samples;
    
    expectations.clear();
    expectations.allocate();

    double objective = 0.0;
    size_t instances = 0;
    
    for (size_t id = 0; id != samples.size(); ++ id) 
      if (! samples[id].empty()) {
	sample_pair_set_type::const_iterator siter_end = samples[id].end();
	for (sample_pair_set_type::const_iterator siter = samples[id].begin(); siter != siter_end; ++ siter) {
	  const sample_pair_type& sample = *siter;
	  
	  weight_type Z_oracle;
	  weight_type Z_kbest;
	  
	  sample_type::const_iterator oiter_end = sample.oracles.end();
	  for (sample_type::const_iterator oiter = sample.oracles.begin(); oiter != oiter_end; ++ oiter)
	    Z_oracle += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->begin(), oiter->end(), 0.0));
	  
	  sample_type::const_iterator kiter_end = sample.kbests.end();
	  for (sample_type::const_iterator kiter = sample.kbests.begin(); kiter != kiter_end; ++ kiter)
	    Z_kbest += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->begin(), kiter->end(), 0.0));
	  
	  for (sample_type::const_iterator oiter = sample.oracles.begin(); oiter != oiter_end; ++ oiter) {
	    const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->begin(), oiter->end(), 0.0)) / Z_oracle;
	    
	    hypothesis_type::feature_set_type::const_iterator fiter_end = oiter->end();
	    for (hypothesis_type::feature_set_type::const_iterator fiter = oiter->begin(); fiter != fiter_end; ++ fiter)
	      expectations[fiter->first] -= weight_type(fiter->second) * weight;
	  }
	  
	  for (sample_type::const_iterator kiter = sample.kbests.begin(); kiter != kiter_end; ++ kiter) {
	    const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->begin(), kiter->end(), 0.0)) / Z_kbest;
	    
	    hypothesis_type::feature_set_type::const_iterator fiter_end = kiter->end();
	    for (hypothesis_type::feature_set_type::const_iterator fiter = kiter->begin(); fiter != fiter_end; ++ fiter)
	      expectations[fiter->first] += weight_type(fiter->second) * weight;
	  }
	  
	  const double margin = log(Z_oracle) - log(Z_kbest);
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
	  
	  // termination...
	  feature.index = -1;
	  feature.value = 0.0;
	  features.push_back(feature);
	}

      shrink();
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
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool merge=false)
  {
    if (id >= encoders.size())
      encoders.resize(id + 1);
    
    if (! merge)
      encoders[id].clear();
    
    encoders[id].encode(kbests, oracles);
  }

  double learn(weight_set_type& weights)
  {
    size_type data_size = 0;
    for (size_type id = 0; id != encoders.size(); ++ id)
      data_size += encoders[id].offsets.size();
    
    label_set_type        labels(data_size, 1);
    feature_node_map_type features;
    features.reserve(data_size);
    
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

  void operator()(const hypothesis_map_type& kbests, const scorer_document_type& scorers, hypothesis_map_type& oracles)
  {
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);
        
    score_ptr_type score_prev;
    score_ptr_type score_best;
    score_ptr_type score_curr;
    
    oracle_map_type oracles_prev(kbests.size());
    oracle_map_type oracles_best(kbests.size());
    oracle_map_type oracles_curr(kbests.size());
    
    // initialization...
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
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
      
      for (size_t id = 0; id != kbests.size(); ++ id) 
	if (! kbests[id].empty()) {
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
    }
    
    oracles.clear();
    oracles.resize(kbests.size());
    for (size_t id = 0; id != kbests.size(); ++ id)
      for (size_t i = 0; i != oracles_best[id].size(); ++ i)
	oracles[id].push_back(*oracles_best[id][i]);
  }
};


#endif
