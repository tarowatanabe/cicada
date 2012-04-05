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

#include <map>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>

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

  PYPPOS(const double __h,
	 const size_type classes,
	 const parameter_type& __emission,
	 const parameter_type& __transition)
    : h(__h),
      h_counts(0),
      phi0(__emission),
      phi(classes, table_emission_type(__emission)),
      base0(1.0 / (classes + 1)), // + 1 for BOS
      counts0(0),
      beta(__transition.discount, __transition.strength),
      pi0(__transition),
      pi(classes, table_transition_type(__transition)),
      emission0(__emission),
      emission(__emission),
      transition0(__transition),
      transition(__transition) {}
  
  template <typename Sampler>
  void increment(const id_type prev, const id_type next, const word_type& word, Sampler& sampler, const double temperature=1.0)
  {
    // emission
    if (word != vocab_type::BOS) {
      if (next >= phi.size())
	phi.resize(next + 1, table_emission_type(emission.discount, emission.strength));
      
      if (phi[next].increment(word, phi0.prob(word, h), sampler, temperature))
	if (phi0.increment(word, h, sampler, temperature))
	  ++ h_counts;
    }
    
    // transition
    if (prev >= pi.size())
      pi.resize(prev + 1, table_transition_type(transition.discount, transition.strength));
    
    // break sticks further...
    while (next >= beta.size() || prev >= beta.size())
      beta.increment(sampler);
    
    if (pi[prev].increment(next, beta[next], sampler, temperature))
      if (pi0.increment(next, base0, sampler, temperature))
	++ counts0;
  }
  
  template <typename Sampler>
  void decrement(const id_type prev, const id_type next, const word_type& word, Sampler& sampler)
  {
    // emission
    if (word != vocab_type::BOS)
      if (phi[next].decrement(word, sampler))
	if (phi0.decrement(word, sampler))
	  -- h_counts;
    
    // transition
    if (pi[prev].decrement(next, sampler))
      if (pi0.decrement(next, sampler))
	-- counts0;
  }
  
  double prob_emission(const id_type next, const word_type& word) const
  {
    const double p0 = phi0.prob(word, h);
    
    return (next < phi.size() ? phi[next].prob(word, p0) : p0);
  }
  
  double prob_transition(const id_type prev, const id_type next) const
  {
    const double p0 = beta[next];
    
    return (prev < pi.size() ? pi[prev].prob(next, p0) : p0);
  }

  double prob_transition(const id_type next) const
  {
    return beta[next];
  }

  size_type classes() const
  {
    return beta.size();
  }
  
  double log_likelihood() const
  {
    double logprob = std::log(h) * h_counts + std::log(base0) * counts0;
    
    logprob += phi0.log_likelihood() + transition0.log_likelihood();
    logprob += pi0.log_likelihood()  + emission0.log_likelihood();
    
    logprob += transition.log_likelihood();
    for (size_type i = 0; i != phi.size(); ++ i)
      logprob += phi[i].log_likelihood();
    
    logprob += emission.log_likelihood();
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
    std::cerr << "permute" << std::endl;
    
    size_type states_size = utils::bithack::max(beta.size(), pi0.size());
    states_size = utils::bithack::max(states_size, pi.size());
    states_size = utils::bithack::max(states_size, phi.size());
    for (size_type i = 0; i != pi.size(); ++ i)
      states_size = utils::bithack::max(states_size, pi[i].size());
    
    
    // we will sort id by the counts...
    mapping.clear();
    for (size_type i = 0; i != states_size; ++ i) {
      mapping.push_back(i);
      
      std::cerr << "i="<< i << " customer=" << pi0.size_customer(i) << std::endl;
    }
    
    // we will always "fix" zero for bos/eos
    std::sort(mapping.begin() + 1, mapping.end(), greater_customer(pi0));
    
    // re-map ids....
    // actually, the mapping data will be used to re-map the training data...
        
    // re-map for transition...
    pi0.permute(mapping);
    
    std::cerr << "truncated pi0: " << pi0.size() << std::endl;
    
    for (size_type i = 0; i != pi.size(); ++ i) {
      std::cerr << "(before) pi i=" << i << " " << pi[i].size() << " table: " << pi[i].size_table() << " customer: " << pi[i].size_customer()<< std::endl;
      
      pi[i].permute(mapping);
      
      std::cerr << "pi i=" << i << " " << pi[i].size() << " table: " << pi[i].size_table() << " customer: " << pi[i].size_customer()<< std::endl;
    }
    
    {
      transition_type pi_new(pi0.size(), table_transition_type(transition));
      
      for (size_type i = 0; i != pi_new.size(); ++ i)
	if (mapping[i] < pi.size())
	  pi_new[i].swap(pi[mapping[i]]);
      
      pi.swap(pi_new);
    }
    
    std::cerr << "truncated pi: " << pi.size() << std::endl;

    // re-map for emission...
    {
      emission_type phi_new(pi0.size(), table_emission_type(emission));
      
      for (size_type i = 0; i != phi_new.size(); ++ i)
	if (mapping[i] < phi.size())
	  phi_new[i].swap(phi[mapping[i]]);
      
      phi.swap(phi_new);
    }

    std::cerr << "truncated phi: " << phi.size() << std::endl;
  }

  template <typename Sampler>
  void initialize(Sampler& sampler, const id_type classes, const int num_loop = 2, const int num_iterations = 8)
  {
    // + 1 including BOS
    sample_parameters(sampler, num_loop, num_iterations);
    
    base0 = 1.0 / (classes + 1);
    
    beta.sample_parameters(classes + 1, sampler);
    
#if 0
    // sample beta from pi0 and base0
    std::vector<double, std::allocator<double> > probs(classes);
    for (id_type state = 0; state != classes; ++ state)
      probs[state] = pi0.prob(state, base0);
    
    beta.assign(probs.begin(), probs.end());
#endif
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    std::cerr << "sample parameters" << std::endl;

    for (int iter = 0; iter != num_loop; ++ iter) {
      emission0.strength = sample_strength(&phi0, &phi0 + 1, sampler, emission0);
      emission0.discount = sample_discount(&phi0, &phi0 + 1, sampler, emission0);
      
      emission.strength = sample_strength(phi.begin(), phi.end(), sampler, emission);
      emission.discount = sample_discount(phi.begin(), phi.end(), sampler, emission);
      
      transition0.strength = sample_strength(&pi0, &pi0 + 1, sampler, transition0);
      transition0.discount = sample_discount(&pi0, &pi0 + 1, sampler, transition0);
      
      transition.strength = sample_strength(pi.begin(), pi.end(), sampler, transition);
      transition.discount = sample_discount(pi.begin(), pi.end(), sampler, transition);      
    }
    
    phi0.strength() = emission0.strength;
    phi0.discount() = emission0.discount;
    
    for (size_type i = 0; i != phi.size(); ++ i) {
      phi[i].strength() = emission.strength;
      phi[i].discount() = emission.discount;
    }
    
    pi0.strength() = transition0.strength;
    pi0.discount() = transition0.discount;
    
    for (size_type i = 0; i != pi.size(); ++ i) {
      pi[i].strength() = transition.strength;
      pi[i].discount() = transition.discount;
    }
    
    // the transition base... base0
    if (! pi0.empty()) {
      base0 = 1.0 / pi0.size();
      
      beta.strength() = pi0.strength();
      beta.discount() = pi0.discount();
      
      // sample beta from pi0 and base0
      std::vector<double, std::allocator<double> > counts(pi0.size() + 1);
      for (id_type state = 0; state != pi0.size(); ++ state)
	counts[state] = pi0.size_customer(state);
      counts.back() = counts0;
      
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
  
  template <typename Sampler>
  void prune(const sentence_type& sentence, const derivation_type& derivation, PYPPOS& model, Sampler& sampler, cutoff_type& cutoff)
  {
    // first, fill-in cutoff-values
    cutoff.clear();
    cutoff.resize(sentence.size() + 2, 0.0);

    const size_type T = cutoff.size();
    const size_type K = model.beta.size();
    
    // we compute threshold based on pi
    double cutoff_min = std::numeric_limits<double>::infinity();
    for (size_type t = 1; t != T - 1; ++ t) {
      cutoff[t] = sampler.uniform(0.0, model.prob_transition(derivation[t - 1], derivation[t]));
      cutoff_min = std::min(cutoff_min, cutoff[t]);
    }
    
    // second, check to see if we need more classes!
    // we check cutoff_min to see if we have enough sticks...
    double pi_max = - std::numeric_limits<double>::infinity();
    for (size_type prev = 0; prev != K; ++ prev) {
      double pi_min = std::numeric_limits<double>::infinity();
      for (size_type next = 0; next != K; ++ next)
	pi_min = std::min(pi_min, model.prob_transition(prev, next));
      
      pi_max = std::max(pi_max, pi_min);
    }
    
    while (pi_max > cutoff_min) {
      model.beta.increment(sampler);
      
      const size_type K = model.beta.size();
      for (size_type k = 0; k != K; ++ k)
	pi_max = std::min(pi_max, model.prob_transition(k, K - 1));
    }
  }

  void initialize(const sentence_type& sentence, const PYPPOS& model)
  {
    const size_type T = sentence.size() + 2;
    const size_type K = model.beta.size();

    alpha.clear();
    alpha.reserve(T, K);
    alpha.resize(T, K);
    alpha(0, 0) = 1.0;
    
    pi.clear();
    pi.resize(K, K);
    for (id_type prev = 0; prev != K; ++ prev)
      for (id_type next = 0; next != K; ++ next)
	pi(prev, next) = model.prob_transition(prev, next);
    
    phi.clear();
    phi.resize(K, T);
    for (size_type t = 1; t != T - 1; ++ t)
      for (id_type state = 1; state != K; ++ state)
	phi(state, t) = model.prob_emission(state, sentence[t - 1]);
    
    phi(0, 0)     = 1.0;
    phi(0, T - 1) = 1.0;
  }
  
  logprob_type forward(const sentence_type& sentence, const PYPPOS& model, const cutoff_type& cutoff)
  {
    initialize(sentence, model);

    const size_type T = alpha.size1();
    const size_type K = alpha.size2();

    logprob_type logsum = cicada::semiring::traits<logprob_type>::one();
    for (size_type t = 1; t != T; ++ t) {
      for (id_type prev = 0; prev != K; ++ prev)
	for (id_type next = 0; next != K; ++ next)
	  if (pi(prev, next) > cutoff[t])
	    alpha(t, next) += alpha(t - 1, prev) * pi(prev, next) * phi(next, t);
      
      double scale = std::accumulate(alpha.begin(t), alpha.end(t), 0.0);
      scale = (scale == 0.0 ? 1.0 : scale);
      if (scale != 1.0)
	std::transform(alpha.begin(t), alpha.end(t), alpha.begin(t), std::bind2nd(std::multiplies<double>(), 1.0 / scale));
      
      logsum *= scale;
    }
    
    return logsum;
  }
  
  template <typename Sampler>
  logprob_type backward(Sampler& sampler, derivation_type& derivation)
  {
    const size_type T = alpha.size1();
    const size_type K = alpha.size2();
    
    logprob_type logprob = cicada::semiring::traits<logprob_type>::one();
    
    id_type state = 0;
    for (size_type t = T - 1; t > 1; -- t) {
      probs.clear();
      for (id_type prev = 0; prev != K; ++ prev)
	probs.push_back(alpha(t - 1, prev) * pi(prev, state) * phi(state, t));
      
      prob_set_type::const_iterator piter = sampler.draw(probs.begin(), probs.end());
      
      state = (piter - probs.begin());
      
      derivation[t - 1] = state;
      
      logprob *= *piter / alpha(t - 1, state);
    }
    
    return logprob;
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

double emission_discount = 0.9;
double emission_strength = 1;

double emission_discount_prior_alpha = 1.0;
double emission_discount_prior_beta  = 1.0;
double emission_strength_prior_shape = 1.0;
double emission_strength_prior_rate  = 1.0;

double transition_discount = 0.9;
double transition_strength = 1;

double transition_discount_prior_alpha = 1.0;
double transition_discount_prior_beta  = 1.0;
double transition_strength_prior_shape = 1.0;
double transition_strength_prior_rate  = 1.0;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);
size_t read_data(const path_set_type& paths, sentence_set_type& sentences);

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
    const size_type vocab_size = read_data(train_files, training);

    std::cerr << "vocab: " << vocab_size << std::endl;
    
    if (training.empty())
      throw std::runtime_error("no training data?");
    
    derivation_set_type derivations(training.size());
    cutoff_set_type     cutoffs(training.size());
    position_set_type   positions(training.size());
    mapping_type mapping;
    mapping_type mapping_new;
    
    for (size_t i = 0; i != training.size(); ++ i)
      positions[i] = i;

    sampler_type sampler;

    PYPPOS model(1.0 / vocab_size,
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

    PYPGraph graph;

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
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      
      std::random_shuffle(positions.begin(), positions.end(), gen);
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_size(training));

      std::vector<double, std::allocator<double> > probs;
      
      for (size_t i = 0; i != positions.size(); ++ i) {
	const size_t pos = positions[i];
	
	if (training[pos].empty()) continue;
	
	if (debug >= 3)
	  std::cerr << "training=" << pos << " classes: " << model.beta.size() << std::endl;
	
	if (derivations[pos].empty()) {
	  // it is our initial condition... sample derivation...
	  const size_type K = model.beta.size();
	  
	  derivations[pos].reserve(training[pos].size() + 2);
	  derivations[pos].push_back(0);
	  
	  for (size_type i = 0; i != training[pos].size(); ++ i) {
	    
	    probs.clear();
	    for (size_type state = 1; state != K; ++ state)
	      probs.push_back(model.prob_transition(state) * model.prob_emission(state, training[pos][i]));
	    
	    derivations[pos].push_back((sampler.draw(probs.begin(), probs.end()) - probs.begin()) + 1);
	  }
	  
	  derivations[pos].push_back(0);
	  
	  if (debug >= 3)
	    for (size_type t = 0; t != derivations[pos].size(); ++ t)
	      std::cerr << "word=" << (t == 0 || t + 1 == derivations[pos].size() ? vocab_type::BOS : training[pos][t - 1])
			<< " pos=" << derivations[pos][t]
			<< std::endl;
	  
	} else {
	  // remove from the model...
	  for (size_type t = 1; t != derivations[pos].size(); ++ t)
	    model.decrement(derivations[pos][t - 1], derivations[pos][t], t + 1 == derivations[pos].size() ? vocab_type::BOS : training[pos][t - 1], sampler);
	}
	
	// pruning...
	graph.prune(training[pos], derivations[pos], model, sampler, cutoffs[pos]);
	
	const PYPGraph::logprob_type logsum = graph.forward(training[pos], model, cutoffs[pos]);
	
	const PYPGraph::logprob_type logderivation = graph.backward(sampler, derivations[pos]);

	for (size_type t = 1; t != derivations[pos].size() - 1; ++ t)
	  if (! derivations[pos][t])
	    std::cerr << "WARNING: empty state?" << std::endl;
	
	if (derivations[pos].front() || derivations[pos].back())
	  throw std::runtime_error("wrong BOS/EOS");
	
	
	if (debug >= 3) {
	  std::cerr << "sum=" << logsum << " derivation=" << logderivation << std::endl;
	  
	  for (size_type t = 0; t != derivations[pos].size(); ++ t)
	    std::cerr << "word=" << (t == 0 || t + 1 == derivations[pos].size() ? vocab_type::BOS : training[pos][t - 1])
		      << " pos=" << derivations[pos][t]
		      << std::endl;
	}
	
	// insert into the model
	for (size_type t = 1; t != derivations[pos].size(); ++ t)
	  model.increment(derivations[pos][t - 1], derivations[pos][t], t + 1 == derivations[pos].size() ? vocab_type::BOS : training[pos][t - 1], sampler, temperature);
      }
      
      model.permute(mapping);
      
      std::cerr << "mapping: ";
      std::copy(mapping.begin(), mapping.end(), std::ostream_iterator<int>(std::cerr, " "));
      std::cerr << std::endl;

      // remap training data
      mapping_new.resize(mapping.size());
      for (size_type i = 0; i != mapping.size(); ++ i)
	mapping_new[mapping[i]] = i;
      
      
      for (size_type pos = 0; pos != derivations.size(); ++ pos)
	for (size_type t = 0; t != derivations[pos].size(); ++ t) {
	  
	  if (derivations[pos][t] >= mapping_new.size())
	    throw std::runtime_error("invalid mapping in derivation?");
	  
	  derivations[pos][t] = mapping_new[derivations[pos][t]];
	  
	  if (derivations[pos][t] >= model.classes())
	    throw std::runtime_error("invalid new mapping in derivation?");
	}
      
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
		  << "classes: " << model.classes() << std::endl;
    }
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

size_t read_data(const path_set_type& paths, sentence_set_type& sentences)
{
  typedef utils::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type vocab_type;
  
  vocab_type vocab;
  vocab.set_empty_key(word_type());
  
  sentences.clear();
    
  sentence_type sentence;
  for (path_set_type::const_iterator fiter = paths.begin(); fiter != paths.end(); ++ fiter) { 
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (is >> sentence) {
      sentences.push_back(sentence);
      
      vocab.insert(sentence.begin(), sentence.end());
    }
  }
  
  return vocab.size();
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

    ("classes", po::value<int>(&classes)->default_value(classes), "# of initial classes")
    
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

