//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// PYP-IHMM!
//

//
// @InProceedings{vangael-vlachos-ghahramani:2009:EMNLP,
//   author    = {Van Gael, Jurgen  and  Vlachos, Andreas  and  Ghahramani, Zoubin},
//   title     = {The infinite {HMM} for unsupervised {PoS} tagging},
//   booktitle = {Proceedings of the 2009 Conference on Empirical Methods in Natural Language Processing},
//   month     = {August},
//   year      = {2009},
//   address   = {Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {678--687},
//   url       = {http://www.aclweb.org/anthology/D/D09/D09-1071}
// }
//

// Iterate:
//   1. compute cutoffs, and potentially break sticks, further
//   2. distribute and run DP and fill-in the new assingment
//         we can asynchronously decrement and increment, since the model is completely independent!
//   3. compute model

// TODO: use id_type(0) as BOS and do not stick-break for BOS!
// make an adjustment for induced POSs

#include <map>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>

#include "utils/unordered_map.hpp"
#include "utils/vector2.hpp"
#include "utils/chunk_vector.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/restaurant.hpp"
#include "utils/restaurant_vector.hpp"
#include "utils/stick_break.hpp"
#include "utils/pyp_parameter.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/dense_hash_set.hpp"
#include "utils/lockfree_list_queue.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

typedef cicada::Symbol    word_type;
typedef cicada::Sentence  sentence_type;
typedef cicada::Vocab     vocab_type;

struct PYP
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef uint32_t  id_type;
};

struct PYPPOS
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  typedef PYP::id_type         id_type;
  
  typedef utils::restaurant_vector< > table_transition_type;
  typedef utils::restaurant<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type > > table_emission_type;
			    
  typedef utils::chunk_vector<table_transition_type, 4096/sizeof(table_transition_type), std::allocator<table_transition_type> > transition_type;
  typedef utils::chunk_vector<table_emission_type, 4096/sizeof(table_emission_type), std::allocator<table_emission_type> >       emission_type;
  
  typedef utils::stick_break< > beta_type;

  typedef utils::pyp_parameter parameter_type;
  
  typedef utils::vector2<double, std::allocator<double > > cache_transition_type;
  typedef utils::vector2<double, std::allocator<double > > cache_emission_type;
  
  PYPPOS(const double __h,
	 const size_type classes,
	 const parameter_type& __emission,
	 const parameter_type& __transition)
    : h(__h),
      h_counts(0),
      phi0(__emission.discount, __emission.strength),
      phi(classes + 1, table_emission_type(__emission.discount,
					   __emission.strength)),
      base0(1.0 / classes),
      counts0(0),
      beta(__transition.discount, __transition.strength),
      pi0(__transition.discount, __transition.strength),
      pi(classes + 1, table_transition_type(__transition.discount,
					    __transition.strength)),
      emission0(__emission),
      emission(__emission),
      transition0(__transition),
      transition(__transition) {}
  
  
  template <typename Sampler>
  void increment(const id_type prev, const id_type next, const word_type& word, Sampler& sampler, const double temperature=1.0)
  {
    if (! next)
      throw std::runtime_error("invalid state");

    // prev = [0, # of states]
    // next = [1, # of states]
    
    // emission
    if (next >= phi.size())
      phi.resize(next + 1, table_emission_type(emission.discount, emission.strength));
    
    if (phi[next].increment(word, phi0.prob(word, h), sampler, temperature))
      if (phi0.increment(word, h, sampler, temperature))
	++ h_counts;
    
    // transition... we need to consider BOS...
    if (prev >= pi.size())
      pi.resize(prev + 1, table_transition_type(transition.discount, transition.strength));
    
    if (pi[prev].increment(next, pi0.prob(next - 1, base0), sampler, temperature))
      if (pi0.increment(next - 1, base0, sampler, temperature))
	++ counts0;
  }
  
  template <typename Sampler>
  void decrement(const id_type prev, const id_type next, const word_type& word, Sampler& sampler)
  {
    if (! next)
      throw std::runtime_error("invalid state");

    // prev = [0, # of states]
    // next = [1, # of states]
    
    // emission
    if (phi[next].decrement(word, sampler))
      if (phi0.decrement(word, sampler))
	-- h_counts;
    
    // transition
    if (pi[prev].decrement(next - 1, sampler))
      if (pi0.decrement(next - 1, sampler))
	-- counts0;
  }

  template <typename Iterator>
  void initialize_cache(Iterator first, Iterator last)
  {
    const size_type K = beta.size() + 1;

    caches_transition.clear();
    caches_emission.clear();
    
    caches_transition.resize(K, K);
    caches_emission.resize(K, word_type::allocated());
    
    for (id_type prev = 0; prev != K; ++ prev)
      for (id_type next = 1; next != K; ++ next)
	caches_transition(prev, next) = prob_transition(prev, next);
    
    for (/**/; first != last; ++ first) {
      if (first->id() >= caches_emission.size2())
	throw std::runtime_error("wrong word...?");
      
      for (id_type state = 1; state != K; ++ state)
	caches_emission(state, first->id()) = prob_emission(state, *first);
    }
  }
  
  double cache_emission(const id_type next, const word_type& word) const
  {
    if (! next || next >= caches_emission.size1())
      throw std::runtime_error("invalid state");
    
    return caches_emission(next, word.id());
  }
  
  double cache_transition(const id_type prev, const id_type next) const
  {
    if (! next || prev >= caches_transition.size1() || next >= caches_transition.size2())
      throw std::runtime_error("invalid state");
    
    return caches_transition(prev, next);
  }
  
  double cache_transition(const id_type next) const
  {
    if (! next || next - 1 >= beta.size())
      throw std::runtime_error("invalid state");
    
    return beta[next - 1];
  }
  
  double prob_emission(const id_type next, const word_type& word) const
  {
    if (! next)
      throw std::runtime_error("invalid state");
    
    const double p0 = phi0.prob(word, h);
    
    return (next < phi.size() ? phi[next].prob(word, p0) : p0);
  }
  
  double prob_transition(const id_type prev, const id_type next) const
  {
    if (! next || next - 1 >= beta.size())
      throw std::runtime_error("invalid state");
    
    const double p0 = beta[next - 1];
    
    return (prev < pi.size() ? pi[prev].prob(next - 1, p0) : p0);
  }

  double prob_transition(const id_type next) const
  {
    if (! next || next - 1 >= beta.size())
      throw std::runtime_error("invalid state");
    
    return beta[next - 1];
  }
  
  double log_likelihood() const
  {
    double logprob = std::log(h) * h_counts + std::log(base0) * counts0;
    
    logprob += phi0.log_likelihood() + transition0.log_likelihood();
    logprob += pi0.log_likelihood()  + emission0.log_likelihood();
    
    logprob += emission.log_likelihood();
    logprob += transition.log_likelihood();
    
    // we will start from 1!
    if (! phi.empty())
      for (size_type i = 1; i != phi.size(); ++ i)
	logprob += phi[i].log_likelihood();
    
    for (size_type i = 0; i != pi.size(); ++ i)
      logprob += pi[i].log_likelihood();
    
    return logprob;
  }

  struct greater_customer
  {
    greater_customer(const table_transition_type& __pi0) : pi0(__pi0) {}

    bool operator()(const size_type& x, const size_type& y) const
    {
      return pi0.size_customer(x) > pi0.size_customer(y);
    }
    
    const table_transition_type& pi0;
  };
    
  template <typename Mapping>
  void permute(Mapping& mapping)
  {
    // we will sort id by the counts...
    // zero-based mapping
    mapping.clear();
    for (size_type i = 0; i != beta.size(); ++ i)
      mapping.push_back(i);
    
    std::sort(mapping.begin(), mapping.end(), greater_customer(pi0));
    
    for (size_type i = 0; i != beta.size(); ++ i)
      std::cerr << "i=" << i << " origin: " << mapping[i] << " customer: " << pi0.size_customer(mapping[i]) << std::endl;
    
    // re-map ids....
    // actually, the mapping data will be used to re-map the training data...
    
    // re-map for transition...
    pi0.permute(mapping);
    for (size_type i = 0; i != pi.size(); ++ i)
      pi[i].permute(mapping);
    
    // + 1 for BOS
    transition_type pi_new(pi0.size() + 1, table_transition_type(transition.discount, transition.strength));
    
    pi_new.front().swap(pi.front()); // for BOS...
    for (size_type i = 1; i != pi_new.size(); ++ i)
      if (mapping[i - 1] + 1 < pi.size()) {
	std::cerr << "pi i=" << i << " from: " << mapping[i - 1] + 1 << std::endl;
	pi_new[i].swap(pi[mapping[i - 1] + 1]);
      } else
	std::cerr << "exceed... pi i=" << i << " from: " << mapping[i - 1] + 1 << std::endl;
    
    pi.swap(pi_new);
    
    // re-map for emission...
    emission_type phi_new(pi0.size() + 1, table_emission_type(emission.discount, emission.strength));
    
    // skip BOS...
    for (size_type i = 1; i != phi_new.size(); ++ i)
      if (mapping[i - 1] + 1 < phi.size()) {
	std::cerr << "phi i=" << i << " from: " << mapping[i - 1] + 1 << std::endl;
	phi_new[i].swap(phi[mapping[i - 1] + 1]);
      } else
	std::cerr << "exceed... phi i=" << i << " from: " << mapping[i - 1] + 1 << std::endl;
    
    phi.swap(phi_new);
  }

  template <typename Sampler>
  void initialize(Sampler& sampler, const id_type classes, const int num_loop = 2, const int num_iterations = 8)
  {
    sample_parameters(sampler, num_loop, num_iterations);
    
    beta.sample_parameters(classes, sampler);
    
    base0 = 1.0 / beta.size();
  }
  
  template <typename Sampler>
  void sample_sticks(const double& cutoff_min, Sampler& sampler)
  {
    double pi_max = - std::numeric_limits<double>::infinity();
    for (id_type prev = 0; prev != beta.size() + 1; ++ prev) {
      double pi_min = std::numeric_limits<double>::infinity();
      for (id_type next = 1; next != beta.size() + 1; ++ next)
	pi_min = std::min(pi_min, prob_transition(prev, next));
      
      pi_max = std::max(pi_max, pi_min);
    }
    
    while (pi_max > cutoff_min) {
      beta.increment(sampler);
      base0 = 1.0 / beta.size();
      
      const size_type K = beta.size() + 1;
      for (size_type k = 0; k != K; ++ k)
	pi_max = std::min(pi_max, prob_transition(k, K - 1));
    }
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    for (int iter = 0; iter != num_loop; ++ iter) {
      emission0.strength = sample_strength(&phi0, &phi0 + 1, sampler, emission0);
      emission0.discount = sample_discount(&phi0, &phi0 + 1, sampler, emission0);
      
      // skip BOS part...
      if (! phi.empty()) {
	emission.strength = sample_strength(phi.begin() + 1, phi.end(), sampler, emission);
	emission.discount = sample_discount(phi.begin() + 1, phi.end(), sampler, emission);
      }
      
      transition0.strength = sample_strength(&pi0, &pi0 + 1, sampler, transition0);
      transition0.discount = sample_discount(&pi0, &pi0 + 1, sampler, transition0);
      
      transition.strength = sample_strength(pi.begin(), pi.end(), sampler, transition);
      transition.discount = sample_discount(pi.begin(), pi.end(), sampler, transition);      
    }
    
    phi0.discount() = emission0.discount;
    phi0.strength() = emission0.strength;
    
    for (size_type i = 0; i != phi.size(); ++ i) {
      phi[i].discount() = emission.discount;
      phi[i].strength() = emission.strength;
    }
    
    pi0.discount() = transition0.discount;
    pi0.strength() = transition0.strength;
    
    beta.discount() = transition0.discount;
    beta.strength() = transition0.strength;
    
    for (size_type i = 0; i != pi.size(); ++ i) {
      pi[i].discount() = transition.discount;
      pi[i].strength() = transition.strength;
    }
    
    // the transition base... base0
    // + 1 for allowing infinity...
    // correct this beta's strength/discount sampling
    if (! pi0.empty()) {
      // sample beta from pi0 and base0
      // or, do we directly use the estimates in pi0???
      
      std::vector<double, std::allocator<double> > counts(pi0.size() + 1);
      for (id_type state = 0; state != pi0.size(); ++ state)
	counts[state] = pi0.size_customer(state) - pi0.size_table(state) * beta.discount();
      counts.back() = beta.strength() + counts0 * beta.discount();
      
      beta.sample_parameters(counts.begin(), counts.end(), sampler);
    }
  }

  template <typename Iterator, typename Sampler>
  double sample_strength(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double x = 0.0;
    double y = 0.0;
    
    for (/**/; first != last; ++ first) {
      x += first->sample_log_x(sampler, param.discount, param.strength);
      y += first->sample_y(sampler, param.discount, param.strength);
    }
    
    return sampler.gamma(param.strength_shape + y, param.strength_rate - x);
  }
  
  template <typename Iterator, typename Sampler>
  double sample_discount(Iterator first, Iterator last, Sampler& sampler, const parameter_type& param) const
  {
    double y = 0.0;
    double z = 0.0;
    
    for (/**/; first != last; ++ first) {
      y += first->sample_y_inv(sampler, param.discount, param.strength);
      z += first->sample_z_inv(sampler, param.discount, param.strength);
    }
    
    return sampler.beta(param.discount_alpha + y, param.discount_beta + z);
  }
  
  double                h;
  size_type             h_counts;
  table_emission_type   phi0;
  emission_type         phi;
  
  double                base0;
  size_type             counts0;
  beta_type             beta;
  table_transition_type pi0;
  transition_type       pi;
  
  parameter_type emission0;
  parameter_type emission;
  parameter_type transition0;
  parameter_type transition;

  cache_transition_type   caches_transition;
  cache_emission_type     caches_emission;
};


struct PYPGraph
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  typedef PYP::id_type         id_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > derivation_type;
  typedef std::vector<double, std::allocator<double> >   cutoff_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  typedef utils::vector2<double, std::allocator<double> > alpha_type;
  
  typedef utils::vector2<double, std::allocator<double> >    phi_type; // emission
  typedef utils::vector2<double, std::allocator<double> >    pi_type;  // transition
  
  typedef std::vector<prob_type, std::allocator<prob_type> > prob_set_type;  

  void initialize(const sentence_type& sentence, const PYPPOS& model)
  {
    const size_type T = sentence.size() + 1;
    const size_type K = model.beta.size() + 1;
    
    alpha.clear();
    alpha.reserve(T, K);
    alpha.resize(T, K);
    alpha(0, 0) = 1.0;
    
    pi.clear();
    pi.resize(K, K);
    for (id_type prev = 0; prev != K; ++ prev)
      for (id_type next = 1; next != K; ++ next)
	pi(prev, next) = model.cache_transition(prev, next);
    
    phi.clear();
    phi.resize(T, K);
    for (size_type t = 1; t != T; ++ t)
      for (id_type state = 1; state != K; ++ state)
	phi(t, state) = model.cache_emission(state, sentence[t - 1]);
    
    phi(0, 0) = 1.0;
  }
    
  logprob_type forward(const sentence_type& sentence, const PYPPOS& model, const cutoff_type& cutoff)
  {
    initialize(sentence, model);

    const size_type T = alpha.size1();
    const size_type K = alpha.size2();

    logprob_type logsum = cicada::semiring::traits<logprob_type>::one();
    for (size_type t = 1; t != T; ++ t) {
      for (id_type prev = 0; prev != K; ++ prev)
	for (id_type next = 1; next != K; ++ next)
	  if (pi(prev, next) > cutoff[t])
	    alpha(t, next) += alpha(t - 1, prev) * pi(prev, next) * phi(t, next);
      
      double scale = std::accumulate(alpha.begin(t), alpha.end(t), 0.0);
      scale = (scale == 0.0 ? 1.0 : scale);
      if (scale != 1.0)
	std::transform(alpha.begin(t), alpha.end(t), alpha.begin(t), std::bind2nd(std::multiplies<double>(), 1.0 / scale));
      
      logsum *= scale;
    }
    
    return logsum;
  }
  
  template <typename Sampler>
  logprob_type backward(Sampler& sampler, derivation_type& derivation, const double temperature)
  {
    const size_type T = alpha.size1();
    const size_type K = alpha.size2();
    
    probs.clear();
    for (id_type state = 1; state != K; ++ state)
      probs.push_back(alpha(T - 1, state));
    
    prob_set_type::const_iterator piter = sampler.draw(probs.begin(), probs.end(), temperature);
    
    logprob_type logprob = cicada::semiring::traits<logprob_type>::one();
    id_type state = (piter - probs.begin()) + 1;
    derivation[T - 1] = state;
    
    for (size_type t = T - 1; t > 1; -- t) {
      probs.clear();
      for (id_type prev = 0; prev != K; ++ prev)
	probs.push_back(alpha(t - 1, prev) * pi(prev, state) * phi(t, state));
      
      prob_set_type::const_iterator piter = sampler.draw(probs.begin(), probs.end(), temperature);
      
      logprob *= *piter / alpha(t - 1, state);
      
      state = piter - probs.begin();
      
      derivation[t - 1] = state;
    }

    if (derivation.front())
      throw std::runtime_error("wrong initial derivation");
    
    for (size_type t = 1; t != T; ++ t)
      if (! derivation[t])
	throw std::runtime_error("wrong derivation");
    
    return logprob;
  }
  
  logprob_type score(const sentence_type& sentence, const derivation_type& derivation, const PYPPOS& model)
  {
    logprob_type logprob = cicada::semiring::traits<logprob_type>::one();
    
    for (size_type t = 1; t != derivation.size(); ++ t)
      logprob *= model.prob_transition(derivation[t - 1], derivation[t]) * model.prob_emission(derivation[t], sentence[t - 1]);
    
    return logprob;
  }
  
  template <typename Sampler>
  void increment(const sentence_type& sentence, const derivation_type& derivation, PYPPOS& model, Sampler& sampler, const double temperature)
  {
    for (size_type t = 1; t != derivation.size(); ++ t)
      model.increment(derivation[t - 1],
		      derivation[t],
		      sentence[t - 1],
		      sampler,
		      temperature);
    
  }
  
  template <typename Sampler>
  void decrement(const sentence_type& sentence, const derivation_type& derivation, PYPPOS& model, Sampler& sampler)
  {
    for (size_type t = 1; t != derivation.size(); ++ t)
      model.decrement(derivation[t - 1],
		      derivation[t],
		      sentence[t - 1],
		      sampler);
  }
  
  
  phi_type    phi;
  pi_type     pi;
  alpha_type  alpha;
  
  prob_set_type    probs;
};

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef utils::sampler<boost::mt19937> sampler_type;

typedef PYP::size_type       size_type;
typedef PYP::difference_type difference_type;

typedef PYPGraph::derivation_type derivation_type;
typedef PYPGraph::cutoff_type     cutoff_type;

typedef std::vector<derivation_type, std::allocator<derivation_type> > derivation_set_type;
typedef std::vector<cutoff_type, std::allocator<cutoff_type> > cutoff_set_type;
typedef std::vector<size_type, std::allocator<size_type> > position_set_type;
typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
typedef std::vector<size_type, std::allocator<size_type> > mapping_type;
typedef utils::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_set_type;

struct TaskBeam
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;

  TaskBeam(queue_type& __queue,
	   const sentence_set_type& __training,
	   cutoff_set_type& __cutoffs,
	   derivation_set_type& __derivations,
	   const PYPPOS& __model,
	   sampler_type& __sampler,
	   double& __cutoff_min)
    : queue(__queue),
      training(__training),
      cutoffs(__cutoffs),
      derivations(__derivations),
      model(__model),
      sampler(__sampler),
      cutoff_min(__cutoff_min)
  {}
  
  void operator()()
  {
    std::vector<double, std::allocator<double> > probs;
    size_type pos;
    
    for (;;) {
      queue.pop(pos);
      
      if (pos == size_type(-1)) break;

      const sentence_type& sentence = training[pos];
      derivation_type& derivation = derivations[pos];
      cutoff_type& cutoff = cutoffs[pos];
      
      if (derivation.empty()) {
	const size_type K = model.beta.size() + 1;
	
	derivation.reserve(training[pos].size() + 1);
	derivation.push_back(0);
	
	for (size_type i = 0; i != sentence.size(); ++ i) {
	  probs.clear();
	  for (size_type state = 1; state != K; ++ state)
	    probs.push_back(model.cache_transition(state) * model.cache_emission(state, sentence[i]));
	  
	  derivation.push_back((sampler.draw(probs.begin(), probs.end()) - probs.begin()) + 1);
	}
      }
      
      cutoff.clear();
      cutoff.resize(sentence.size() + 1, 0.0);
      
      const size_type T = cutoff.size();
      
      // we compute threshold based on pi
      for (size_type t = 1; t != T; ++ t) {
	cutoff[t] = sampler.uniform(0.0, model.cache_transition(derivation[t - 1], derivation[t]));
	cutoff_min = std::min(cutoff_min, cutoff[t]);
      }
    }
  }
  
  queue_type& queue;
  
  const sentence_set_type& training;
  cutoff_set_type& cutoffs;
  derivation_set_type& derivations;
  
  const PYPPOS& model;
  sampler_type  sampler;
  double& cutoff_min;
};

struct TaskPermute
{
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;
  
  TaskPermute(queue_type& __queue,
	      derivation_set_type& __derivations,
	      const mapping_type& __mapping)
    : queue(__queue),
      derivations(__derivations)
  {
    mapping.clear();
    mapping.resize(__mapping.size() + 1, 0);
    for (size_type i = 1; i != mapping.size(); ++ i)
      mapping[__mapping[i - 1] + 1] = i;
  }
  
  void operator()()
  {
    size_type pos;
    
    for (;;) {
      queue.pop(pos);
      
      if (pos == size_type(-1)) break;
      
      for (size_type t = 0; t != derivations[pos].size(); ++ t)
	derivations[pos][t] = mapping[derivations[pos][t]];
    }
  }
  
  queue_type& queue;
  derivation_set_type& derivations;
  mapping_type mapping;
};

struct Task
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;
  
  Task(queue_type& __mapper,
       queue_type& __reducer,
       const sentence_set_type& __training,
       const cutoff_set_type& __cutoffs,
       derivation_set_type& __derivations,
       const PYPPOS& __model,
       sampler_type& __sampler)
    : mapper(__mapper),
      reducer(__reducer),
      training(__training),
      cutoffs(__cutoffs),
      derivations(__derivations),
      model(__model),
      sampler(__sampler) {}
  
  void operator()()
  {
    std::vector<double, std::allocator<double> > probs;
    size_type pos;

    derivation_type derivation_prev;
    
    for (;;) {
      mapper.pop(pos);
      
      if (pos == size_type(-1)) break;
      
      // forward
      graph.forward(training[pos], model, cutoffs[pos]);
      
      // backward
      graph.backward(sampler, derivations[pos], temperature);
      
      reducer.push(pos);
    }
  }
  
  queue_type& mapper;
  queue_type& reducer;
  
  const sentence_set_type& training;
  const cutoff_set_type& cutoffs;
  derivation_set_type& derivations;
  
  const PYPPOS& model;
  sampler_type  sampler;
  double temperature;
  
  PYPGraph graph;
};

struct less_size
{
  less_size(const sentence_set_type& __training) : training(__training) {}
  
  bool operator()(const size_type& x, const size_type& y) const
  {
    return training[x].size() < training[y].size();
  }

  const sentence_set_type& training;
};


path_set_type train_files;
path_set_type test_files;
path_type     output_file;

int classes = 16;

int samples = 30;
int baby_steps = 0;
int anneal_steps = 0;
int resample_rate = 1;
int resample_iterations = 2;
bool slice_sampling = false;

double emission_discount = 0.1;
double emission_strength = 1;

double emission_discount_prior_alpha = 1.0;
double emission_discount_prior_beta  = 1.0;
double emission_strength_prior_shape = 1.0;
double emission_strength_prior_rate  = 1.0;

double transition_discount = 0.1;
double transition_strength = 10;

double transition_discount_prior_alpha = 1.0;
double transition_discount_prior_beta  = 1.0;
double transition_strength_prior_shape = 1.0;
double transition_strength_prior_rate  = 1.0;


int threads = 1;
int debug = 0;

void options(int argc, char** argv);
void read_data(const path_set_type& paths, sentence_set_type& sentences, word_set_type& words);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    if (samples < 0)
      throw std::runtime_error("# of samples must be positive");
            
    if (resample_rate <= 0)
      throw std::runtime_error("resample rate must be >= 1");

    if (classes <= 0)
      throw std::runtime_error("zero/negative initial class size");
    
    if (train_files.empty())
      throw std::runtime_error("no training data?");
    
    if (! slice_sampling && (emission_strength < 0.0 || transition_strength < 0.0))
      throw std::runtime_error("negative strength w/o slice sampling is not supported!");
    
    sentence_set_type training;
    word_set_type     words;
    words.set_empty_key(word_type());
    
    read_data(train_files, training, words);

    if (training.empty())
      throw std::runtime_error("no training data?");
    
    derivation_set_type derivations(training.size());
    derivation_set_type derivations_prev(training.size());
    cutoff_set_type     cutoffs(training.size());
    position_set_type   positions(training.size());
    mapping_type mapping;
    
    position_set_type positions_mapped;
    position_set_type positions_reduced;

    for (size_t i = 0; i != training.size(); ++ i)
      positions[i] = i;
    
    sampler_type sampler;
    
    PYPPOS model(1.0 / words.size(),
		 classes,
		 PYPPOS::parameter_type(emission_discount,
					emission_strength,
					emission_discount_prior_alpha,
					emission_discount_prior_beta,
					emission_strength_prior_shape,
					emission_strength_prior_rate),
		 
		 PYPPOS::parameter_type(transition_discount,
					transition_strength,
					transition_discount_prior_alpha,
					transition_discount_prior_beta,
					transition_strength_prior_shape,
					transition_strength_prior_rate));
    
    model.initialize(sampler, classes, resample_iterations);

    if (debug >= 2) {
      std::cerr << "emission-base discount=" << model.emission0.discount << " strength=" << model.emission0.strength << std::endl
		<< "emission      discount=" << model.emission.discount << " strength=" << model.emission.strength << std::endl
		<< "transition-base discount=" << model.transition0.discount << " strength=" << model.transition0.strength << std::endl
		<< "transition      discount=" << model.transition.discount << " strength=" << model.transition.strength << std::endl;
      
      std::cerr << "beta:";
      for (size_t state = 0; state != model.beta.size(); ++ state)
	std::cerr << ' ' << model.beta[state];
      std::cerr <<  " remain: " << model.beta[model.beta.size()] << std::endl;
    }
    
    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);
    
    size_t baby_iter = 0;
    const size_t baby_last = utils::bithack::branch(baby_steps > 0, baby_steps, 0);
    
    bool sampling = false;
    int sample_iter = 0;

    PYPGraph graph;
    
    Task::queue_type queue_mapper;
    Task::queue_type queue_reducer;
    
    std::vector<Task, std::allocator<Task> > tasks(threads, Task(queue_mapper,
								 queue_reducer,
								 training,
								 cutoffs,
								 derivations,
								 model,
								 sampler));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    // then, learn!
    for (size_t iter = 0; sample_iter != samples; ++ iter, sample_iter += sampling) {
      
      double temperature = 1.0;
      bool anneal_finished = true;
      if (anneal_iter != anneal_last) {
	anneal_finished = false;
	temperature = double(anneal_last - anneal_iter) + 1;
	
	++ anneal_iter;
	
	if (debug >= 2)
	  std::cerr << "temperature: " << temperature << std::endl;
      }
      
      bool baby_finished = true;
      if (baby_iter != baby_last) {
	++ baby_iter;
	baby_finished = false;
      }
      
      sampling = anneal_finished && baby_finished;
      
      if (debug) {
	if (sampling)
	  std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
	else
	  std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
      }
      
      // assign temperature...
      for (size_type i = 0; i != tasks.size(); ++ i)
	tasks[i].temperature = temperature;
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      
      std::random_shuffle(positions.begin(), positions.end(), gen);
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_size(training));
      
      // sample beam...
      if (debug >= 3)
	std::cerr << "sample beam" << std::endl;
      
      model.initialize_cache(words.begin(), words.end());
      
      cutoff_type cutoff_min(threads, std::numeric_limits<double>::infinity());
      
      TaskBeam::queue_type queue_beam;
      boost::thread_group workers_beam;
      for (int i = 0; i != threads; ++ i)
	workers_beam.add_thread(new boost::thread(TaskBeam(queue_beam,
							   training,
							   cutoffs,
							   derivations,
							   model,
							   sampler,
							   cutoff_min[i])));
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
	const size_type& pos = *piter;
	
	derivations_prev[pos] = derivations[pos];

	queue_beam.push(pos);
      }
      
      for (int i = 0; i != threads; ++ i)
	queue_beam.push(size_type(-1));
      
      workers_beam.join_all();
      
      if (debug >= 3)
	std::cerr << "sample sticks" << std::endl;
      
      // sample sticks..
      model.sample_sticks(*std::min_element(cutoff_min.begin(), cutoff_min.end()), sampler);

      if (debug >= 2)
	std::cerr << "# of sticks: " << model.beta.size() << std::endl;
      
      // sample derivations...
      if (debug >= 3)
	std::cerr << "sample derivations" << std::endl;
      
      model.initialize_cache(words.begin(), words.end());
      
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter)
	queue_mapper.push(*piter);
      
      for (size_type reduced = 0; reduced != positions.size(); ++ reduced) {
	size_type pos = 0;
	queue_reducer.pop(pos);

	if (derivations[pos].size() != cutoffs[pos].size())
	  throw std::runtime_error("derivation and cutoff size differ");
	
	if (derivations[pos].size() != training[pos].size() + 1)
	  throw std::runtime_error("derivation and setnence size differ");
	
	if (! derivations_prev[pos].empty())
	  graph.decrement(training[pos], derivations_prev[pos], model, sampler);
	
	graph.increment(training[pos], derivations[pos], model, sampler, temperature);
      }

      // permute..
      if (debug >= 3)
	std::cerr << "permute" << std::endl;
      
      model.permute(mapping);
      
      TaskPermute::queue_type queue_permute;
      boost::thread_group workers_permute;
      for (int i = 0; i != threads; ++ i)
	workers_permute.add_thread(new boost::thread(TaskPermute(queue_permute,
								 derivations,
								 mapping)));
      
      
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter)
	queue_permute.push(*piter);
      
      for (int i = 0; i != threads; ++ i)
	queue_permute.push(size_type(-1));
      
      workers_permute.join_all();
      
      // sample other parameters...
      if (debug >= 3)
	std::cerr << "sample parameters" << std::endl;
	    
      model.sample_parameters(sampler, resample_iterations);
      
      if (debug >= 2) {
	std::cerr << "emission-base discount=" << model.emission0.discount << " strength=" << model.emission0.strength << std::endl
		  << "emission      discount=" << model.emission.discount << " strength=" << model.emission.strength << std::endl
		  << "transition-base discount=" << model.transition0.discount << " strength=" << model.transition0.strength << std::endl
		  << "transition      discount=" << model.transition.discount << " strength=" << model.transition.strength << std::endl;
	
	std::cerr << "beta:";
	for (size_t state = 0; state != model.beta.size(); ++ state)
	  std::cerr << ' ' << model.beta[state];
	std::cerr <<  " remain: " << model.beta[model.beta.size()] << std::endl;
      }
      
      if (debug)
	std::cerr << "log-likelihood: " << model.log_likelihood() << std::endl
		  << "classes: " << model.pi0.size() << std::endl;
    }
    
    
    for (int i = 0; i != threads; ++ i)
      queue_mapper.push(size_type(-1));
    
    workers.join_all();
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void read_data(const path_set_type& paths, sentence_set_type& sentences, word_set_type& vocab)
{
  sentences.clear();
    
  sentence_type sentence;
  for (path_set_type::const_iterator fiter = paths.begin(); fiter != paths.end(); ++ fiter) { 
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (is >> sentence) {
      sentences.push_back(sentence);
      
      vocab.insert(sentence.begin(), sentence.end());
    }
  }
}
void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("train", po::value<path_set_type>(&train_files)->multitoken(), "train file(s)")
    ("test",  po::value<path_set_type>(&test_files)->multitoken(),  "test file(s)")
    ("output", po::value<path_type>(&output_file), "output file")

    ("classes",     po::value<int>(&classes)->default_value(classes),         "# of initial classes")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    
    ("emission-discount",       po::value<double>(&emission_discount)->default_value(emission_discount),                         "discount ~ Beta(alpha,beta)")
    ("emission-discount-alpha", po::value<double>(&emission_discount_prior_alpha)->default_value(emission_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("emission-discount-beta",  po::value<double>(&emission_discount_prior_beta)->default_value(emission_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("emission-strength",       po::value<double>(&emission_strength)->default_value(emission_strength),                         "strength ~ Gamma(shape,rate)")
    ("emission-strength-shape", po::value<double>(&emission_strength_prior_shape)->default_value(emission_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("emission-strength-rate",  po::value<double>(&emission_strength_prior_rate)->default_value(emission_strength_prior_rate),   "strength ~ Gamma(shape,rate)")

    ("transition-discount",       po::value<double>(&transition_discount)->default_value(transition_discount),                         "discount ~ Beta(alpha,beta)")
    ("transition-discount-alpha", po::value<double>(&transition_discount_prior_alpha)->default_value(transition_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("transition-discount-beta",  po::value<double>(&transition_discount_prior_beta)->default_value(transition_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("transition-strength",       po::value<double>(&transition_strength)->default_value(transition_strength),                         "strength ~ Gamma(shape,rate)")
    ("transition-strength-shape", po::value<double>(&transition_strength_prior_shape)->default_value(transition_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("transition-strength-rate",  po::value<double>(&transition_strength_prior_rate)->default_value(transition_strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}

