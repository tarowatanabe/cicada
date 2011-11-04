//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// kbest learner with MPI

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"
#include "cicada_kbest_impl.hpp"
#include "cicada_text_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/base64.hpp"
#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>

#include "lbfgs.h"

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type kbest_path;
path_set_type oracle_path;
path_type weights_path;
path_type output_path = "-";
path_type output_objective_path;

path_set_type refset_files;

int iteration = 100;
bool learn_lbfgs = false;
bool learn_sgd = false;
bool learn_mira = false;
bool learn_arow = false;
bool learn_cw = false;
bool learn_pegasos = false;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

bool loss_margin = false; // margin by loss, not rank-loss
bool softmax_margin = false;
bool line_search_local = false;
bool line_search_global = false;

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool unite_kbest = false;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

void read_kbest(const  scorer_document_type& scorers,
		const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles);
void read_refset(const path_set_type& file,
		 scorer_document_type& scorers);

template <typename Optimize>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights);
template <typename Optimize, typename Generator>
double optimize_online(const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       weight_set_type& weights,
		       Generator& generator);

template <typename Optimizer>
struct OptimizeOnline;
template <typename Optimizer>
struct OptimizeOnlineMargin;
struct OptimizeLBFGS;

void bcast_weights(const int rank, weight_set_type& weights);
void send_weights(const weight_set_type& weights);
void reduce_weights(weight_set_type& weights);


int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw + learn_pegasos > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,sgd,mira,arow,cw}");
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw + learn_pegasos == 0)
      learn_lbfgs = true;

    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));

    if ((line_search_local || line_search_global) && (learn_lbfgs || learn_sgd))
      throw std::runtime_error("line-search is applicable only for non-maxent based loss");

    if (kbest_path.empty())
      throw std::runtime_error("no kbest?");
    if (oracle_path.empty())
      throw std::runtime_error("no oracke kbest?");

    scorer_document_type scorers(scorer_name);
    
    if (! refset_files.empty()) {
      read_refset(refset_files, scorers);
      
      if (! unite_kbest && kbest_path.size() > 1) {
	scorer_document_type scorers_iterative(scorer_name);
	scorers_iterative.resize(scorers.size() * kbest_path.size());
	
	for (size_t i = 0; i != kbest_path.size(); ++ i)
	  std::copy(scorers.begin(), scorers.end(), scorers_iterative.begin() + scorers.size() * i);
	
	scorers.swap(scorers_iterative);
      }
    }

    
    hypothesis_map_type kbests;
    hypothesis_map_type oracles;
    
    read_kbest(scorers, kbest_path, oracle_path, kbests, oracles);
    
    weight_set_type weights;
    if (mpi_rank == 0 && ! weights_path.empty()) {
      if (! boost::filesystem::exists(weights_path))
	throw std::runtime_error("no path? " + weights_path.string());
      
      utils::compress_istream is(weights_path, 1024 * 1024);
      is >> weights;
    }
    
    // collect features...
    for (int rank = 0; rank < mpi_size; ++ rank) {
      weight_set_type weights;
      weights.allocate();
      
      for (feature_type::id_type id = 0; id != feature_type::allocated(); ++ id)
	if (! feature_type(id).empty())
	  weights[feature_type(id)] = 1.0;
      
      bcast_weights(rank, weights);
    }
    
    if (debug && mpi_rank == 0)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;
    
    weights.allocate();

    double objective = 0.0;

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    if (learn_sgd) {
      if (regularize_l1)
	objective = optimize_online<OptimizeOnline<OptimizerSGDL1> >(kbests, oracles, weights, generator);
      else
	objective = optimize_online<OptimizeOnline<OptimizerSGDL2> >(kbests, oracles, weights, generator);
    } else if (learn_mira)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerMIRA> >(kbests, oracles, weights, generator);
    else if (learn_arow)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerAROW> >(kbests, oracles, weights, generator);
    else if (learn_cw)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerCW> >(kbests, oracles, weights, generator);
    else if (learn_pegasos)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerPegasos> >(kbests, oracles, weights, generator);
    else 
      objective = optimize_batch<OptimizeLBFGS>(kbests, oracles, weights);

    if (debug && mpi_rank == 0)
      std::cerr << "objective: " << objective << std::endl;
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_path, 1024 * 1024);
      os.precision(20);
      os << weights;
      
      if (! output_objective_path.empty()) {
	utils::compress_ostream os(output_objective_path, 1024 * 1024);
	os.precision(20);
	os << objective << '\n';
      }
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  return 0;
}

enum {
  weights_tag = 1000,
  gradients_tag,
  notify_tag,
  termination_tag,
  point_size_tag,
  point_tag,
};

inline
int loop_sleep(bool found, int non_found_iter)
{
  if (! found) {
    boost::thread::yield();
    ++ non_found_iter;
  } else
    non_found_iter = 0;
    
  if (non_found_iter >= 64) {
    struct timespec tm;
    tm.tv_sec = 0;
    tm.tv_nsec = 2000001; // above 2ms
    nanosleep(&tm, NULL);
    
    non_found_iter = 0;
  }
  return non_found_iter;
}

template <typename Optimizer>
struct OptimizeOnline
{
  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  OptimizeOnline(Optimizer& __optimizer,
		 const hypothesis_map_type& __kbests,
		 const hypothesis_map_type& __oracles)
    : optimizer(__optimizer),
      kbests(__kbests),
      oracles(__oracles)
  {
    ids.reserve(kbests.size());
    ids.resize(kbests.size());
    
    for (size_t id = 0; id != ids.size(); ++ id)
      ids[id] = id;
  }
  
  template <typename Generator>
  void operator()(Generator& generator)
  {
    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);
    
    for (size_t id = 0; id != ids.size(); ++ id)
      if (! oracles[ids[id]].empty() && ! kbests[ids[id]].empty())
	operator()(oracles[ids[id]], kbests[ids[id]]);
  }

  typedef Optimizer optimizer_type;
  
  typedef typename optimizer_type::weight_type   weight_type;
  typedef typename optimizer_type::gradient_type gradient_type;    
  
  double objective(const weight_set_type& weights) const
  {
    double obj = 0.0;
    const double cost_factor = (softmax_margin ? 1.0 : 0.0);

    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! oracles[id].empty() && ! kbests[id].empty()) {
	weight_type Z_oracle;
	weight_type Z_kbest;
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter)
	  Z_oracle += function(weights, oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss);
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter)
	  Z_kbest += function(weights, kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss);
	
	obj -= log(Z_oracle) - log(Z_kbest);
      }
    
    return obj;
  }

  double normalizer() const
  {
    size_t num = 0;
    for (size_t id = 0; id != kbests.size(); ++ id) 
      num += (! oracles[id].empty() && ! kbests[id].empty());
    
    return num;
  }
  
  template <typename Iterator>
  double operator()(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return 0.0;
  }
  
  
  template <typename Iterator>
  weight_type function(Iterator first, Iterator last, const double init) const
  {
    return cicada::semiring::traits<weight_type>::exp(cicada::dot_product(optimizer.weights, first, last, 0.0) * optimizer.weight_scale + init);
  }
  
  template <typename Iterator>
  weight_type function(const weight_set_type& weights, Iterator first, Iterator last, const double init) const
  {
    return cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, first, last, 0.0) + init);
  }
  
  
  void operator()(const hypothesis_set_type& oracles,
		  const hypothesis_set_type& kbests)
  {
    const double cost_factor = (softmax_margin ? 1.0 : 0.0);
    
    weight_type Z_oracle;
    weight_type Z_kbest;
    
    hypothesis_set_type::const_iterator oiter_end = oracles.end();
    for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter)
      Z_oracle += function(oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss);
    
    hypothesis_set_type::const_iterator kiter_end = kbests.end();
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter)
      Z_kbest += function(kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss);
    
    gradient_type gradient_oracles;
    gradient_type gradient_kbests;
    
    for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter) {
      const weight_type weight = function(oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss) / Z_oracle;
      
      hypothesis_type::feature_set_type::const_iterator fiter_end = oiter->features.end();
      for (hypothesis_type::feature_set_type::const_iterator fiter = oiter->features.begin(); fiter != fiter_end; ++ fiter)
	gradient_oracles[fiter->first] += weight_type(fiter->second) * weight;
    }
    
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
      const weight_type weight = function(kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss) / Z_kbest;
      
      hypothesis_type::feature_set_type::const_iterator fiter_end = kiter->features.end();
      for (hypothesis_type::feature_set_type::const_iterator fiter = kiter->features.begin(); fiter != fiter_end; ++ fiter)
	gradient_kbests[fiter->first] += weight_type(fiter->second) * weight;
    }
    
    optimizer(gradient_oracles,
	      gradient_kbests,
	      Z_oracle,
	      Z_kbest);
  }
  
  Optimizer& optimizer;

  const hypothesis_map_type& kbests;
  const hypothesis_map_type& oracles;

  id_set_type ids;
};

template <typename Optimizer>
struct OptimizeOnlineMargin
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
  
  typedef std::vector<double, std::allocator<double> > loss_set_type;

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

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  OptimizeOnlineMargin(Optimizer& __optimizer,
		       const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles)
    : optimizer(__optimizer)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    
    features_type feats;
    sentence_unique_type  sentences;
    
    for (size_type id = 0; id != kbests.size(); ++ id)
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	sentences.clear();
	for (size_t o = 0; o != oracles[id].size(); ++ o)
	  sentences.insert(oracles[id][o].sentence);
	
	for (size_t o = 0; o != oracles[id].size(); ++ o)
	  for (size_t k = 0; k != kbests[id].size(); ++ k) {
	    const hypothesis_type& oracle = oracles[id][o];
	    const hypothesis_type& kbest  = kbests[id][k];
	    
	    // ignore oracle translations
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

	    if (loss_margin) {
	      const double loss = kbest.loss - oracle.loss;
	      
	      // checking...
	      if (loss > 0.0) {
		features.insert(feats.begin(), feats.end());
		losses.push_back(loss);
	      }
	    } else {
	      features.insert(feats.begin(), feats.end());
	      losses.push_back(1.0);
	    }
	  }
      }
    
    loss_set_type(losses).swap(losses);
    
    ids.reserve(losses.size());
    ids.resize(losses.size());
    
    for (size_t id = 0; id != ids.size(); ++ id)
      ids[id] = id;

    optimizer.instances = losses.size();
  }
  
  typedef Optimizer optimizer_type;
  
  typedef typename optimizer_type::weight_type   weight_type;
  typedef typename optimizer_type::gradient_type gradient_type;

  template <typename Generator>
  void operator()(Generator& generator)
  {
    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);
    
    for (size_t i = 0; i != ids.size(); ++ i) {
      const size_type id = ids[i];
      
      optimizer(features[id].begin(), features[id].end(), losses[id]);
    }
  }
  
  double objective(const weight_set_type& weights) const
  {
    double obj = 0.0;
    for (size_t id = 0; id != losses.size(); ++ id) {
      const double margin = cicada::dot_product(weights, features[id].begin(), features[id].end(), 0.0);
      const double loss = losses[id];
      
      obj += (loss - margin) * (loss - margin > 0.0);
    }
    return obj;
  }
  
  double normalizer() const
  {
    return losses.size();
  }
  
  template <typename Iterator>
  double operator()(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    double grad = 0.0;
    for (size_t id = 0; id != losses.size(); ++ id) {
      const double margin      = cicada::dot_product(weights,      features[id].begin(), features[id].end(), 0.0);
      const double margin_prev = cicada::dot_product(weights_prev, features[id].begin(), features[id].end(), 0.0);
      
      const double bi = margin_prev - margin;
      const double ci = losses[id]  - margin_prev;
      
      const double ki = - ci / bi;
      
      if (ki > 0) {
	*iter = std::make_pair(ki, bi);
	++ iter;
      }
      
      if ((bi < 0.0 && ki > 0.0) || (bi > 0.0 && ki <= 0.0))
	grad += bi;
    }
    return grad;
  }
  
  template <typename Iterator>
  double function(Iterator first, Iterator last)
  {
    return cicada::dot_product(optimizer.weights, first, last, 0.0) * optimizer.weight_scale;
  }
  
  void operator()(const hypothesis_set_type& oracles,
		  const hypothesis_set_type& kbests)
  {
    // compute the best among kbests
    // compute the worst among oracles
    
    // then, compute!
    
    double score_oracle =   std::numeric_limits<double>::infinity();
    double score_kbest  = - std::numeric_limits<double>::infinity();
    hypothesis_set_type::const_iterator oiter_best;
    hypothesis_set_type::const_iterator kiter_best;
    
    hypothesis_set_type::const_iterator oiter_end = oracles.end();
    for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter) {
      const double score = function(oiter->features.begin(), oiter->features.end());
      
      if (score < score_oracle) {
	score_oracle = score;
	oiter_best = oiter;
      }
    }
    
    hypothesis_set_type::const_iterator kiter_end = kbests.end();
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
      const double score = function(kiter->features.begin(), kiter->features.end());
	    
      if (score > score_kbest) {
	score_kbest = score;
	kiter_best = kiter;
      }
    }

    if (loss_margin) {
      const double loss = kiter_best->loss - oiter_best->loss;
      
      if (loss > 0.0)
	optimizer(feature_set_type(oiter_best->features.begin(), oiter_best->features.end()),
		  feature_set_type(kiter_best->features.begin(), kiter_best->features.end()),
		  loss);
    } else
      optimizer(feature_set_type(oiter_best->features.begin(), oiter_best->features.end()),
		feature_set_type(kiter_best->features.begin(), kiter_best->features.end()));
  }

  Optimizer& optimizer;
  
  sample_set_type features;
  loss_set_type   losses;
  id_set_type     ids;
};


double l2norm(const weight_set_type& weights)
{
  double norm = 0.0;
  weight_set_type::const_iterator witer_end = weights.end();
  for (weight_set_type::const_iterator witer = weights.begin(); witer != witer_end; ++ witer)
    norm += (*witer) * (*witer);
  return norm;
}

double l2norm_diff(const weight_set_type& weights1, const weight_set_type& weights2)
{
  weight_set_type::const_iterator witer1 = weights1.begin();
  weight_set_type::const_iterator witer2 = weights2.begin();
  weight_set_type::const_iterator witer1_end = weights1.end();
  weight_set_type::const_iterator witer2_end = weights2.end();
  
  double norm = 0.0;
  if (weights1.size() <= weights2.size()) {
    for (/**/; witer1 != witer1_end; ++ witer1, ++ witer2)
      norm += ((*witer1) - (*witer2)) * ((*witer1) - (*witer2));
    for (/**/; witer2 != witer2_end; ++ witer2)
      norm += (*witer2) * (*witer2);
  } else {
    for (/**/; witer2 != witer2_end; ++ witer1, ++ witer2)
      norm += ((*witer1) - (*witer2)) * ((*witer1) - (*witer2));
    for (/**/; witer1 != witer1_end; ++ witer1)
      norm += (*witer1) * (*witer1);
  }
  
  return norm;
}

template <typename Optimize, typename Generator>
double optimize_online(const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       weight_set_type& weights,
		       Generator& generator)
{
  typedef typename Optimize::optimizer_type optimizer_type;

  typedef std::pair<double, double> point_type;
  typedef std::vector<point_type, std::allocator<point_type> > point_set_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  int instances_local = 0;
  
  for (size_t id = 0; id != kbests.size(); ++ id)
    instances_local += ! kbests[id].empty() && ! oracles[id].empty();
  
  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);

  optimizer_type optimizer(instances, C);
  Optimize opt(optimizer, kbests, oracles);

  const double norm_local = opt.normalizer();
  double norm = 0.0;
  MPI::COMM_WORLD.Reduce(&norm_local, &norm, 1, MPI::DOUBLE, MPI::SUM, 0);
  
  weight_set_type weights_init = weights;
  point_set_type points;
  point_set_type points_next;
  
  optimizer.weights = weights;
  
  if (mpi_rank == 0) {
    double objective_prev = 0.0;
    double objective = 0.0;
    
    int increased = 0;
    
    for (int iter = 0; iter < iteration; ++ iter) {
      
      for (int rank = 1; rank < mpi_size; ++ rank)
	MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
      
      bcast_weights(0, optimizer.weights);
      
      weight_set_type weights_prev = optimizer.weights;
      
      optimizer.initialize();

      opt(generator);
            
      optimizer.finalize();
      
      optimizer.weights *= (optimizer.samples + 1);
      reduce_weights(optimizer.weights);
      
      objective = 0.0;
      MPI::COMM_WORLD.Reduce(&optimizer.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      objective /= norm;
      
      int samples = 0;
      int samples_local = optimizer.samples;
      MPI::COMM_WORLD.Reduce(&samples_local, &samples, 1, MPI::INT, MPI::SUM, 0);

      const int active_size = samples;
            
      samples += mpi_size;
      
      optimizer.weights *= (1.0 / samples);
      
      // optimizer.weights is the new weights.

      if (line_search_local) {
	// perform line-search between weights_prev and optimizer.weights, and update optimizer.weights
	
	if (debug >= 3)
	  std::cerr << "line-search" << std::endl;
	
	bcast_weights(0, optimizer.weights);
	
	points.clear();
	const double grad_local = opt(optimizer.weights, weights_prev, std::back_inserter(points));
	std::sort(points.begin(), points.end());
	
	double grad = 0.0;
	MPI::COMM_WORLD.Reduce(&grad_local, &grad, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	// merge points from others... we assume that we will consume in sorted order!
	for (int rank = 1; rank < mpi_size; ++ rank) {
	  boost::iostreams::filtering_istream is;
	  is.push(utils::mpi_device_source(rank, point_tag, 1024 * 1024));
	  
	  double point;
	  double b;
	  
	  points_next.clear();
	  point_set_type::const_iterator piter = points.begin();
	  point_set_type::const_iterator piter_end = points.end();
	  
	  while (is.read((char*) &point, sizeof(double)) && is.read((char*) &b, sizeof(double))) {
	    for (/**/; piter != piter_end && piter->first < point; ++ piter)
	      points_next.push_back(*piter);
	    points_next.push_back(std::make_pair(point, b));
	  }
	  
	  // final insertion...
	  points_next.insert(points_next.end(), piter, piter_end);
	  
	  points.swap(points_next);
	  points_next.clear();
	}

	if (debug >= 3)
	  std::cerr << "point size: " << points.size() << std::endl;
		
	const double norm_w      = cicada::dot_product(optimizer.weights, optimizer.weights);
	const double dot_prod    = cicada::dot_product(weights_prev, optimizer.weights);
	const double norm_w_prev = cicada::dot_product(weights_prev, weights_prev);
	
	const double a0 = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * norm;
	const double b0 = (dot_prod - norm_w_prev) * C * norm;

	//std::cerr << "a0: "  << a0 << " b0: " << b0 << std::endl;
	
	grad += b0;

	if (debug >= 3)
	  std::cerr << "gradient: " << grad << std::endl;

	double k = 0.0;
	
	if (! points.empty() && grad < 0.0) {
	  point_set_type::const_iterator piter_end = points.end();
	  for (point_set_type::const_iterator piter = points.begin(); piter != piter_end && grad < 0.0; ++ piter) {
	    const double k_new = piter->first;
	    const double grad_new = grad + std::fabs(piter->second) + a0 * (k_new - k);
	    
	    if (grad_new >= 0) {
	      // compute intersection...
	      k = k + grad * (k - k_new) / (grad_new - grad);
	      grad = grad_new;
	      break;
	    } else {
	      k = k_new;
	      grad = grad_new;
	    }
	  }
	  
	  if (debug >= 3)
	    std::cerr << "searched gradient: " << grad << " k: " << k << std::endl;
	  
	  k = std::max(k, 0.0);
	}
	
	//const double merge_ratio = k + 0.1 * (1.0 - k);
	const double merge_ratio = (k == 0.0 ? 0.1 : k);
	
	// move to optimizer.weights * merge_ratio + weights_prev * (1.0 - merge_ratio)

	weight_set_type weights_prev_saved = weights_prev;
	
	optimizer.weights *= merge_ratio;
	weights_prev *= (1.0 - merge_ratio);
	
	optimizer.weights += weights_prev;
	
	weights_prev.swap(weights_prev_saved);
      }

      if (regularize_l2)
	objective += 0.5 * C * cicada::dot_product(optimizer.weights, optimizer.weights);
      else {
	double norm = 0.0;
	for (size_t i = 0; i < optimizer.weights.size(); ++ i)
	  norm += std::fabs(optimizer.weights[i]);
	objective += C * norm;
      }
      
      increased = utils::bithack::branch(iter && (objective > objective_prev), increased + 1, 0);
      
      //const double norm_x = std::max(1.0, l2norm(optimizer.weights));
      //const double norm_d = l2norm_diff(optimizer.weights, weights_prev);
      
      const bool converged = (active_size == 0 || (iter && (std::fabs((objective - objective_prev) / objective) < 1e-7)) || increased > 10);
      
      if (debug >= 2)
	std::cerr << "objective: " << objective << " active size: " << active_size << std::endl;
      
      if (converged) break;
      
      objective_prev = objective;
    }
    
    // send termination!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, termination_tag);

    weights.swap(optimizer.weights);
    
    bcast_weights(0, weights);
    
    if (line_search_global) {
      // perform line-search between weights_init and weights, and update weights
      
      bcast_weights(0, weights_init);
      
      if (debug >= 3)
	std::cerr << "line-search" << std::endl;
      
      points.clear();
      const double grad_local = opt(weights, weights_init, std::back_inserter(points));
      std::sort(points.begin(), points.end());
      
      double grad = 0.0;
      MPI::COMM_WORLD.Reduce(&grad_local, &grad, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      // merge points from others... we assume that we will consume in sorted order!
      for (int rank = 1; rank < mpi_size; ++ rank) {
	boost::iostreams::filtering_istream is;
	is.push(utils::mpi_device_source(rank, point_tag, 1024 * 1024));
	
	double point;
	double b;
	
	points_next.clear();
	point_set_type::const_iterator piter = points.begin();
	point_set_type::const_iterator piter_end = points.end();
	  
	while (is.read((char*) &point, sizeof(double)) && is.read((char*) &b, sizeof(double))) {
	  for (/**/; piter != piter_end && piter->first < point; ++ piter)
	    points_next.push_back(*piter);
	  points_next.push_back(std::make_pair(point, b));
	}
	  
	// final insertion...
	points_next.insert(points_next.end(), piter, piter_end);
	
	points.swap(points_next);
	points_next.clear();
      }
      
      if (debug >= 3)
	std::cerr << "point size: " << points.size() << std::endl;
      
      const double norm_w      = cicada::dot_product(weights, weights);
      const double dot_prod    = cicada::dot_product(weights_init, weights);
      const double norm_w_prev = cicada::dot_product(weights_init, weights_init);
      
      const double a0 = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * norm;
      const double b0 = (dot_prod - norm_w_prev) * C * norm;
      
      grad += b0;

      if (debug >= 3)
	std::cerr << "gradient: " << grad << std::endl;

      double k = 0.0;
	
      if (! points.empty() && grad < 0.0) {
	point_set_type::const_iterator piter_end = points.end();
	for (point_set_type::const_iterator piter = points.begin(); piter != piter_end && grad < 0.0; ++ piter) {
	  const double k_new = piter->first;
	  const double grad_new = grad + std::fabs(piter->second) + a0 * (k_new - k);
	  
	  if (grad_new >= 0) {
	    // compute intersection...
	    k = k + grad * (k - k_new) / (grad_new - grad);
	    grad = grad_new;
	    break;
	  } else {
	    k = k_new;
	    grad = grad_new;
	  }
	}
	  
	if (debug >= 3)
	  std::cerr << "searched gradient: " << grad << " k: " << k << std::endl;
	
	k = std::max(k, 0.0);
      }
      
      //const double merge_ratio = k + 0.1 * (1.0 - k);
      const double merge_ratio = (k == 0.0 ? 0.1 : k);
      
      weights *= merge_ratio;
      weights_init *= (1.0 - merge_ratio);
      
      weights += weights_init;

      // re-bcast again...
      bcast_weights(0, weights);
    }
    
    const double objective_local = opt.objective(weights);
    
    objective = 0.0;
    MPI::COMM_WORLD.Reduce(&objective_local, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
    objective /= norm;
    
    if (regularize_l2)
      objective += 0.5 * C * cicada::dot_product(weights, weights);
    else {
      double norm = 0.0;
      for (size_t i = 0; i < weights.size(); ++ i)
	norm += std::fabs(weights[i]);
      objective += C * norm;
    }
    
    return objective;
    
  } else {
    enum {
      NOTIFY = 0,
      TERMINATION,
    };
    
    MPI::Prequest requests[2];
    
    requests[NOTIFY]      = MPI::COMM_WORLD.Recv_init(0, 0, MPI::INT, 0, notify_tag);
    requests[TERMINATION] = MPI::COMM_WORLD.Recv_init(0, 0, MPI::INT, 0, termination_tag);
    
    for (int i = 0; i < 2; ++ i)
      requests[i].Start();
    
    while (1) {
      if (MPI::Request::Waitany(2, requests))
	break;
      else {
	requests[NOTIFY].Start();
	
	bcast_weights(0, optimizer.weights);

	const weight_set_type weights_prev = optimizer.weights;
	
	optimizer.initialize();
	
	opt(generator);
	
	optimizer.finalize();
	
	optimizer.weights *= (optimizer.samples + 1);
	send_weights(optimizer.weights);
	
	double objective = 0.0;
	MPI::COMM_WORLD.Reduce(&optimizer.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	int samples = 0;
	int samples_local = optimizer.samples;
	MPI::COMM_WORLD.Reduce(&samples_local, &samples, 1, MPI::INT, MPI::SUM, 0);
	
	if (line_search_local) {
	  // perform line-search between weights_prev and optimizer.weights, and update optimizer.weights
	  
	  bcast_weights(0, optimizer.weights);
	  
	  points.clear();
	  const double grad_local = opt(optimizer.weights, weights_prev, std::back_inserter(points));
	  std::sort(points.begin(), points.end());

	  double grad = 0.0;
	  MPI::COMM_WORLD.Reduce(&grad_local, &grad, 1, MPI::DOUBLE, MPI::SUM, 0);
	  
	  boost::iostreams::filtering_ostream os;
	  os.push(utils::mpi_device_sink(0, point_tag, 1024 * 1024));
	  
	  point_set_type::const_iterator piter_end = points.end();
	  for (point_set_type::const_iterator piter = points.begin(); piter != piter_end; ++ piter) {
	    os.write((char*) &(piter->first), sizeof(double));
	    os.write((char*) &(piter->second), sizeof(double));
	  }
	}
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    bcast_weights(0, weights);

    if (line_search_global) {
      // perform line-search between weights_init and weights, and update weights
      
      bcast_weights(0, weights_init);
      
      points.clear();
      const double grad_local = opt(weights, weights_init, std::back_inserter(points));
      std::sort(points.begin(), points.end());
      
      double grad = 0.0;
      MPI::COMM_WORLD.Reduce(&grad_local, &grad, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      {
	boost::iostreams::filtering_ostream os;
	os.push(utils::mpi_device_sink(0, point_tag, 1024 * 1024));
	
	point_set_type::const_iterator piter_end = points.end();
	for (point_set_type::const_iterator piter = points.begin(); piter != piter_end; ++ piter) {
	  os.write((char*) &(piter->first), sizeof(double));
	  os.write((char*) &(piter->second), sizeof(double));
	}
      }
      
      // re-bcast again...
      bcast_weights(0, weights);
    }
    
    const double objective_local = opt.objective(weights);
    double objective = 0.0;
    MPI::COMM_WORLD.Reduce(&objective_local, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);

    return 0.0;
  }

}

struct OptimizeLBFGS
{

  OptimizeLBFGS(const hypothesis_map_type& __kbests,
		const hypothesis_map_type& __oracles,
		weight_set_type& __weights,
		const size_t& __instances)
    : kbests(__kbests),
      oracles(__oracles),
      weights(__weights),
      instances(__instances) {}

  double operator()()
  {
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    
    if (regularize_l1) {
      param.orthantwise_c = C;
      param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    } else
      param.orthantwise_c = 0.0;
    
    param.max_iterations = iteration;
    
    double objective = 0.0;
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    return objective;
  }

  
  struct Task
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;
    
    Task(const weight_set_type& __weights,
	 const hypothesis_map_type& __kbests,
	 const hypothesis_map_type& __oracles,
	 const size_t& __instances)
      : weights(__weights),
	kbests(__kbests),
	oracles(__oracles),
	instances(__instances)
    {}
    
    void operator()()
    {
      g.clear();
      objective = 0.0;

      expectation_type  expectations;
      
      expectations.allocate();
      expectations.clear();
      
      const size_t id_max = utils::bithack::min(kbests.size(), oracles.size());

      const double cost_factor = (softmax_margin ? 1.0 : 0.0);
      
      for (size_t id = 0; id != id_max; ++ id)
	if (! kbests[id].empty() && ! oracles[id].empty()) {
	  
	  weight_type Z_oracle;
	  weight_type Z_kbest;
	  
	  hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	  for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter)
	    Z_oracle += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss));
	
	  hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	  for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter)
	    Z_kbest += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss));
	
	  for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	    const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->features.begin(), oiter->features.end(), cost_factor * oiter->loss)) / Z_oracle;
	  
	    hypothesis_type::feature_set_type::const_iterator fiter_end = oiter->features.end();
	    for (hypothesis_type::feature_set_type::const_iterator fiter = oiter->features.begin(); fiter != fiter_end; ++ fiter)
	      expectations[fiter->first] -= weight_type(fiter->second) * weight;
	  }
	
	  for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	    const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->features.begin(), kiter->features.end(), cost_factor * kiter->loss)) / Z_kbest;
	  
	    hypothesis_type::feature_set_type::const_iterator fiter_end = kiter->features.end();
	    for (hypothesis_type::feature_set_type::const_iterator fiter = kiter->features.begin(); fiter != fiter_end; ++ fiter)
	      expectations[fiter->first] += weight_type(fiter->second) * weight;
	  }
	
	  const double margin = log(Z_oracle) - log(Z_kbest);
	  objective -= margin;
	
	  if (debug >= 3)
	    std::cerr << " margin: " << margin << std::endl;
	}
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(expectations.begin(), expectations.end(), g.begin());
      
      objective /= instances;
      std::transform(g.begin(), g.end(), g.begin(), std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    }
    
    const weight_set_type& weights;

    const hypothesis_map_type& kbests;
    const hypothesis_map_type& oracles;
    
    size_t instances;
    
    double          objective;
    weight_set_type g;
  };

  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    typedef Task                  task_type;

    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();

    OptimizeLBFGS& optimizer = *((OptimizeLBFGS*) instance);
    
    // send notification!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
    
    bcast_weights(0, optimizer.weights);
    
    task_type task(optimizer.weights, optimizer.kbests, optimizer.oracles, optimizer.instances);
    task();
    
    // collect all the objective and gradients...
    
    reduce_weights(task.g);
    
    std::fill(g, g + n, 0.0);
    std::transform(task.g.begin(), task.g.end(), g, g, std::plus<double>());
    
    double objective = 0.0;
    MPI::COMM_WORLD.Reduce(&task.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
    
    const double objective_unregularized = objective;
    
    // L2...
    if (regularize_l2) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i) {
	g[i] += C * x[i];
	norm += x[i] * x[i];
      }
      objective += 0.5 * C * norm;
    }
    
    if (debug >= 2)
      std::cerr << "objective: " << objective << " non-regularized: " << objective_unregularized << std::endl;
    
    return objective;
  }
    
  const hypothesis_map_type& kbests;
  const hypothesis_map_type& oracles;
  
  weight_set_type& weights;
  size_t instances;
};

template <typename Optimize>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const size_t id_max = utils::bithack::min(kbests.size(), oracles.size());
  
  int instances_local = 0;
  for (size_t id = 0; id != id_max; ++ id)
    instances_local += (!kbests[id].empty()) && (!oracles[id].empty());

  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  
  if (mpi_rank == 0) {
    const double objective = Optimize(kbests, oracles, weights, instances)();
    
    // send termination!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, termination_tag);
    
    return objective;
  } else {

    enum {
      NOTIFY = 0,
      TERMINATION,
    };
    
    MPI::Prequest requests[2];

    requests[NOTIFY]      = MPI::COMM_WORLD.Recv_init(0, 0, MPI::INT, 0, notify_tag);
    requests[TERMINATION] = MPI::COMM_WORLD.Recv_init(0, 0, MPI::INT, 0, termination_tag);
    
    for (int i = 0; i < 2; ++ i)
      requests[i].Start();

    while (1) {
      if (MPI::Request::Waitany(2, requests))
	break;
      else {
	typedef typename Optimize::Task task_type;

	requests[NOTIFY].Start();

	bcast_weights(0, weights);
	
	task_type task(weights, kbests, oracles, instances);
	task();
	
	send_weights(task.g);
	
	double objective = 0.0;
	MPI::COMM_WORLD.Reduce(&task.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    return 0.0;
  }
}

void unique_kbest(hypothesis_map_type& kbests)
{
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
				  std::allocator<hypothesis_type> > hypothesis_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
			std::allocator<hypothesis_type> > hypothesis_unique_type;
#endif
  
  hypothesis_unique_type uniques;
  
  for (size_t id = 0; id != kbests.size(); ++ id) 
    if (! kbests[id].empty()) {
      uniques.clear();
      uniques.insert(kbests[id].begin(), kbests[id].end());
      
      kbests[id].clear();
      hypothesis_set_type(kbests[id]).swap(kbests[id]);
      
      kbests[id].reserve(uniques.size());
      kbests[id].insert(kbests[id].end(), uniques.begin(), uniques.end());
    }
}

void read_kbest(const scorer_document_type& scorers,
		const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles)
{  
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const bool error_metric = scorers.error_metric();
  const double loss_factor = (error_metric ? 1.0 : - 1.0);
  
  parser_type parser;
  kbest_feature_type kbest;
  
  if (unite_kbest) {
    for (path_set_type::const_iterator piter = kbest_path.begin(); piter != kbest_path.end(); ++ piter) {
      if (mpi_rank == 0 && debug)
	std::cerr << "reading kbest: " << piter->string() << std::endl;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;

	if (i >= kbests.size())
	  kbests.resize(i + 1);
	
	utils::compress_istream is(path_kbest, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  const size_t& id = boost::fusion::get<0>(kbest);

	  if (id != i)
	    throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(id));
	  	  
	  kbests[i].push_back(hypothesis_type(kbest));

	  hypothesis_type& kbest = kbests[i].back();
	  
	  if (! scorers.empty()) {
	    if (i >= scorers.size())
	      throw std::runtime_error("reference positions outof index");
	    
	    kbest.score = scorers[i]->score(sentence_type(kbest.sentence.begin(), kbest.sentence.end()));
	    kbest.loss  = kbest.score->score() * loss_factor;
	  } else
	    kbest.loss = 1;
	}
      }
    }
    
    oracles.resize(kbests.size());
    
    for (path_set_type::const_iterator piter = oracle_path.begin(); piter != oracle_path.end(); ++ piter) {
      if (mpi_rank == 0 && debug)
	std::cerr << "reading oracles: " << piter->string() << std::endl;
      
      for (size_t i = mpi_rank; i < oracles.size(); i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_oracle = (*piter) / file_name;
	
	if (! boost::filesystem::exists(path_oracle)) continue;
	
	utils::compress_istream is(path_oracle, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  if (boost::fusion::get<0>(kbest) != i)
	    throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(boost::fusion::get<0>(kbest)));
	  
	  oracles[i].push_back(hypothesis_type(kbest));

	  hypothesis_type& oracle = oracles[i].back();
	  
	  if (! scorers.empty()) {
	    if (i >= scorers.size())
	      throw std::runtime_error("reference positions outof index");
	    
	    oracle.score = scorers[i]->score(sentence_type(oracle.sentence.begin(), oracle.sentence.end()));
	    oracle.loss  = oracle.score->score() * loss_factor;
	  } else
	    oracle.loss = 0;
	}
      }
    }
    
  } else {
    // synchronous reading...
    if (kbest_path.size() != oracle_path.size())
      throw std::runtime_error("# of kbests does not match");
    
    const size_type refset_size = scorers.size() / kbest_path.size();
    
    for (size_t pos = 0; pos != kbest_path.size(); ++ pos) {
      if (mpi_rank == 0 && debug)
	std::cerr << "reading kbest: " << kbest_path[pos].string() << " with " << oracle_path[pos].string() << std::endl;
      
      const size_type refset_offset = refset_size * pos;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest  = kbest_path[pos] / file_name;
	const path_type path_oracle = oracle_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;
	if (! boost::filesystem::exists(path_oracle)) continue;

	kbests.resize(kbests.size() + 1);
	oracles.resize(oracles.size() + 1);

	const size_type refset_pos = refset_offset + i;
	
	{
	  utils::compress_istream is(path_kbest, 1024 * 1024);
	  is.unsetf(std::ios::skipws);
	  
	  iter_type iter(is);
	  iter_type iter_end;
	  
	  while (iter != iter_end) {
	    boost::fusion::get<1>(kbest).clear();
	    boost::fusion::get<2>(kbest).clear();
	    
	    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	      if (iter != iter_end)
		throw std::runtime_error("kbest parsing failed");

	    if (boost::fusion::get<0>(kbest) != i)
	      throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(boost::fusion::get<0>(kbest)));
	    
	    kbests.back().push_back(hypothesis_type(kbest));

	    hypothesis_type& kbest = kbests.back().back();
	    
	    if (! scorers.empty()) {
	    kbest.score = scorers[refset_pos]->score(sentence_type(kbest.sentence.begin(), kbest.sentence.end()));
	    kbest.loss  = kbest.score->score() * loss_factor;
	  } else
	    kbest.loss = 1;
	  }
	}

	{
	  utils::compress_istream is(path_oracle, 1024 * 1024);
	  is.unsetf(std::ios::skipws);
	  
	  iter_type iter(is);
	  iter_type iter_end;
	  
	  while (iter != iter_end) {
	    boost::fusion::get<1>(kbest).clear();
	    boost::fusion::get<2>(kbest).clear();
	    
	    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	      if (iter != iter_end)
		throw std::runtime_error("kbest parsing failed");

	    if (boost::fusion::get<0>(kbest) != i)
	      throw std::runtime_error("different id: " + utils::lexical_cast<std::string>(boost::fusion::get<0>(kbest)));
	    
	    oracles.back().push_back(hypothesis_type(kbest));
	    
	    hypothesis_type& oracle = oracles.back().back();

	    if (! scorers.empty()) {
	      oracle.score = scorers[refset_pos]->score(sentence_type(oracle.sentence.begin(), oracle.sentence.end()));
	      oracle.loss  = oracle.score->score() * loss_factor;
	    } else
	      oracle.loss = 0.0;
	  }
	}
      }
    }
  }

  // uniques...
  unique_kbest(kbests);
  unique_kbest(oracles);  
}


void reduce_weights(weight_set_type& weights)
{
  typedef utils::mpi_device_source            device_type;
  typedef boost::iostreams::filtering_istream stream_type;

  typedef boost::shared_ptr<device_type> device_ptr_type;
  typedef boost::shared_ptr<stream_type> stream_ptr_type;

  typedef std::vector<device_ptr_type, std::allocator<device_ptr_type> > device_ptr_set_type;
  typedef std::vector<stream_ptr_type, std::allocator<stream_ptr_type> > stream_ptr_set_type;

  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  device_ptr_set_type device(mpi_size);
  stream_ptr_set_type stream(mpi_size);

  for (int rank = 1; rank < mpi_size; ++ rank) {
    device[rank].reset(new device_type(rank, weights_tag, 4096));
    stream[rank].reset(new stream_type());
    
    stream[rank]->push(boost::iostreams::zlib_decompressor());
    stream[rank]->push(*device[rank]);
  }

  std::string line;
  
  int non_found_iter = 0;
  while (1) {
    bool found = false;
    
    for (int rank = 1; rank < mpi_size; ++ rank)
      while (stream[rank] && device[rank] && device[rank]->test()) {
	if (std::getline(*stream[rank], line)) {
	  const utils::piece line_piece(line);
	  tokenizer_type tokenizer(line_piece);
	  
	  tokenizer_type::iterator iter = tokenizer.begin();
	  if (iter == tokenizer.end()) continue;
	  const utils::piece feature = *iter;
	  ++ iter;
	  if (iter == tokenizer.end()) continue;
	  const utils::piece value = *iter;
	  
	  weights[feature] += utils::decode_base64<double>(value);
	} else {
	  stream[rank].reset();
	  device[rank].reset();
	}
	found = true;
      }
    
    if (std::count(device.begin(), device.end(), device_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
}


void send_weights(const weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  boost::iostreams::filtering_ostream os;
  os.push(boost::iostreams::zlib_compressor());
  os.push(utils::mpi_device_sink(0, weights_tag, 4096));
  
  for (feature_type::id_type id = 0; id < weights.size(); ++ id)
    if (! feature_type(id).empty() && weights[id] != 0.0) {
      os << feature_type(id) << ' ';
      utils::encode_base64(weights[id], std::ostream_iterator<char>(os));
      os << '\n';
    }
}

void bcast_weights(const int rank, weight_set_type& weights)
{
  typedef std::vector<char, std::allocator<char> > buffer_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == rank) {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(utils::mpi_device_bcast_sink(rank, 4096));
    
    static const weight_set_type::feature_type __empty;
    
    weight_set_type::const_iterator witer_begin = weights.begin();
    weight_set_type::const_iterator witer_end = weights.end();
    
    for (weight_set_type::const_iterator witer = witer_begin; witer != witer_end; ++ witer)
      if (*witer != 0.0) {
	const weight_set_type::feature_type feature(witer - witer_begin);
	if (feature != __empty) {
	  os << feature << ' ';
	  utils::encode_base64(*witer, std::ostream_iterator<char>(os));
	  os << '\n';
	}
      }
  } else {
    weights.clear();
    weights.allocate();
    
    boost::iostreams::filtering_istream is;
    is.push(boost::iostreams::zlib_decompressor());
    is.push(utils::mpi_device_bcast_source(rank, 4096));
    
    std::string feature;
    std::string value;
    
    while ((is >> feature) && (is >> value))
      weights[feature] = utils::decode_base64<double>(value);
  }
}

void read_refset(const path_set_type& files, scorer_document_type& scorers)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;

  if (files.empty())
    throw std::runtime_error("no reference files?");
    
  scorers.clear();

  parser_type parser;
  id_sentence_type id_sentence;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->string());

    utils::compress_istream is(*fiter, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      id_sentence.second.clear();
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
	if (iter != iter_end)
	  throw std::runtime_error("refset parsing failed");
      
      const int& id = id_sentence.first;
      
      if (id >= static_cast<int>(scorers.size()))
	scorers.resize(id + 1);
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      scorers[id]->insert(id_sentence.second);
    }
  }
}

void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("kbest",   po::value<path_set_type>(&kbest_path)->multitoken(),  "kbest path")
    ("oracle",  po::value<path_set_type>(&oracle_path)->multitoken(), "oracle kbest path")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    ("weights", po::value<path_type>(&weights_path),      "initial parameter")
    ("output",  po::value<path_type>(&output_path),       "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",   po::bool_switch(&learn_lbfgs),   "batch LBFGS algorithm")
    ("learn-sgd",     po::bool_switch(&learn_sgd),     "online SGD algorithm")
    ("learn-mira",    po::bool_switch(&learn_mira),    "online MIRA algorithm")
    ("learn-arow",    po::bool_switch(&learn_arow),    "online AROW algorithm")
    ("learn-cw",      po::bool_switch(&learn_cw),      "online CW algorithm")
    ("learn-pegasos", po::bool_switch(&learn_pegasos), "online Pegasos algorithm")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C",             po::value<double>(&C)->default_value(C), "regularization constant")
    
    ("loss-margin",    po::bool_switch(&loss_margin),        "direct loss margin")
    ("softmax-margin", po::bool_switch(&softmax_margin),     "softmax margin")

    ("line-search-local",  po::bool_switch(&line_search_local),  "perform line search locally in each iteration")
    ("line-search-global", po::bool_switch(&line_search_global), "perform line search globally")
    
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")
    
    ("unite",    po::bool_switch(&unite_kbest), "unite kbest sharing the same id")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  desc_command.add(opts_command);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {

    if (mpi_rank == 0)
      std::cout << argv[0] << " [options] [operations]\n"
		<< opts_command << std::endl;

    MPI::Finalize();
    exit(0);
  }
}
