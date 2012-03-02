//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include "cicada_mert_kbest_impl.hpp"

#include "cicada/optimize_qp.hpp"
#include "cicada/optimize.hpp"
#include "cicada/semiring/envelope.hpp"

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
#include "utils/map_file.hpp"
#include "utils/tempfile.hpp"
#include "utils/mulvector2.hpp"
#include "utils/mathop.hpp"
#include "utils/unordered_set.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>

#include "lbfgs.h"
#include "lbfgs_fortran.h"
#include "lbfgs_error.hpp"

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;
typedef std::pair<score_ptr_type, score_ptr_type> score_ptr_pair_type;

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
typedef std::vector<size_t, std::allocator<size_t> > kbest_map_type;

typedef std::vector<weight_set_type, std::allocator<weight_set_type> > weights_history_type;

path_set_type kbest_path;
path_set_type oracle_path;
path_type     weights_path;
path_set_type weights_history_path;
path_type output_path = "-";
path_type output_objective_path;

path_type bound_lower_file;
path_type bound_upper_file;

path_set_type refset_files;

int iteration = 100;
bool learn_lbfgs = false;
bool learn_sgd = false;
bool learn_xbleu = false;
bool learn_mira = false;
bool learn_nherd = false;
bool learn_arow = false;
bool learn_cw = false;
bool learn_pegasos = false;
bool learn_cp = false;
bool learn_mcp = false;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;
double scale = 1.0;
double eta0 = 0.2;
int order = 4;

bool annealing_mode = false;
bool quenching_mode = false;

double temperature = 0.0;
double temperature_start = 1000;
double temperature_end = 0.001;
double temperature_rate = 0.5;

double quench_start = 0.01;
double quench_end = 100;
double quench_rate = 10;

bool loss_margin = false; // margin by loss, not rank-loss
bool softmax_margin = false;
bool line_search = false;
bool mert_search = false;
bool mert_search_local = false;
bool sample_vector = false;
bool direct_loss = false;
bool conservative_loss = false;

bool scale_fixed = false;

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool unite_kbest = false;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

void read_kbest(const scorer_document_type& scorers,
		const path_set_type& kbest_path,
		hypothesis_map_type& kbests,
		kbest_map_type&      kbest_map);
void read_kbest(const  scorer_document_type& scorers,
		const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles,
		kbest_map_type& kbest_map);
void read_refset(const path_set_type& file,
		 scorer_document_type& scorers);

template <typename Optimize>
double optimize_cp(const scorer_document_type& scorers,
		   const hypothesis_map_type& kbests,
		   const hypothesis_map_type& oracles,
		   const kbest_map_type& kbest_map,
		   const weight_set_type& bounds_lower,
		   const weight_set_type& bounds_upper,
		   const weights_history_type& weights_history,
		   weight_set_type& weights);
template <typename Optimize>
double optimize_xbleu(const hypothesis_map_type& kbests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights);
template <typename Optimize>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      const weight_set_type& bounds_lower,
		      const weight_set_type& bounds_upper,
		      weight_set_type& weights);
template <typename Optimize, typename Generator>
double optimize_online(const scorer_document_type& scorers,
		       const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       const kbest_map_type& kbest_map,
		       const weight_set_type& bounds_lower,
		       const weight_set_type& bounds_upper,
		       weight_set_type& weights,
		       Generator& generator);

double optimize_mert(const scorer_document_type& scorers,
		     const hypothesis_map_type& kbests,
		     const kbest_map_type& kbest_map,
		     const double scale_min,
		     const double scale_max,
		     const weight_set_type& weights_prev,
		     weight_set_type& weights);

template <typename Optimizer>
struct OptimizeOnline;
template <typename Optimizer>
struct OptimizeOnlineMargin;
struct OptimizeCP;
struct OptimizeMCP;
struct OptimizeLBFGS;
struct OptimizeXBLEU;

void bcast_weights(const int rank, weight_set_type& weights);
void send_weights(const weight_set_type& weights);
void reduce_weights(weight_set_type& weights);
void reduce_score(score_ptr_pair_type& score_pair);

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw + learn_pegasos + learn_cp + learn_mcp + learn_nherd + learn_xbleu > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,sgd,mira,arow,cw}");
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw + learn_pegasos + learn_cp + learn_mcp + learn_nherd + learn_xbleu == 0)
      learn_lbfgs = true;

    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));

    if (line_search && (learn_lbfgs || learn_sgd))
      throw std::runtime_error("line-search is applicable only for non-maxent based loss");
    
    if (kbest_path.empty())
      throw std::runtime_error("no kbest?");
    if (! learn_xbleu && oracle_path.empty())
      throw std::runtime_error("no oracle kbest?");
    if (learn_xbleu && refset_files.empty())
      throw std::runtime_error("no reference translations?");
    if (learn_xbleu && order <= 0)
      throw std::runtime_error("invalid ngram order");

    if (annealing_mode) {
      if (! (temperature_end < temperature_start))
	throw std::runtime_error("temperature should start higher, then decreased");
      if (temperature_rate <= 0.0 || temperature_rate >= 1.0)
	throw std::runtime_error("temperature rate should be 0.0 < rate < 1.0: " + utils::lexical_cast<std::string>(temperature_rate));
    }
    
    if (quenching_mode) {
      if (! (quench_start < quench_end))
	throw std::runtime_error("quenching should start lower, then increased");
      if (quench_rate <= 1.0)
	throw std::runtime_error("quenching rate should be > 1.0: " + utils::lexical_cast<std::string>(quench_rate)); 
    }

    if (! bound_lower_file.empty())
      if (bound_lower_file != "-" && ! boost::filesystem::exists(bound_lower_file))
	throw std::runtime_error("no lower-bound file? " + bound_lower_file.string());
    
    if (! bound_upper_file.empty())
      if (bound_upper_file != "-" && ! boost::filesystem::exists(bound_upper_file))
	throw std::runtime_error("no upper-bound file? " + bound_upper_file.string());


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
    
    if (mert_search && scorers.empty())
      throw std::runtime_error("mert search requires evaluation scores");
    if (mert_search_local && scorers.empty())
      throw std::runtime_error("mert search requires evaluation scores");
    if (sample_vector && scorers.empty())
      throw std::runtime_error("sampling requires evaluation scores");
    if (learn_mcp && scorers.empty())
      throw std::runtime_error("MCP requires evaluation scores");
    
    hypothesis_map_type kbests;
    hypothesis_map_type oracles;
    kbest_map_type      kbest_map;
    
    if (! learn_xbleu)
      read_kbest(scorers, kbest_path, oracle_path, kbests, oracles, kbest_map);
    else
      read_kbest(scorers, kbest_path, kbests, kbest_map);
    
    weight_set_type weights;
    if (mpi_rank == 0 && ! weights_path.empty()) {
      if (! boost::filesystem::exists(weights_path))
	throw std::runtime_error("no path? " + weights_path.string());
      
      utils::compress_istream is(weights_path, 1024 * 1024);
      is >> weights;
    }

    weights_history_type weights_history;
    if (mpi_rank == 0 && ! weights_history_path.empty()) {
      weights_history.reserve(weights_history_path.size());
      
      for (size_t i = 0; i != weights_history_path.size(); ++ i) {
	const path_type& path = weights_history_path[i];
	
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error("no path? " + path.string());
	
	weights_history.push_back(weight_set_type());
	
	utils::compress_istream is(path, 1024 * 1024);
	is >> weights_history.back();
	
	if (weights_history.back().empty())
	  weights_history.pop_back();
      }
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
    
    weight_set_type bounds_lower;
    weight_set_type bounds_upper;
    
    if (! bound_lower_file.empty())
      read_bounds(bound_lower_file, bounds_lower, - std::numeric_limits<double>::infinity());
    
    if (! bound_upper_file.empty())
      read_bounds(bound_upper_file, bounds_upper,   std::numeric_limits<double>::infinity());

    if (debug && mpi_rank == 0)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;
    
    weights.allocate();

    double objective = 0.0;

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    weight_set_type weights_prev = weights;
    
    if (learn_sgd) {
      if (regularize_l1)
	objective = optimize_online<OptimizeOnline<OptimizerSGDL1> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
      else
	objective = optimize_online<OptimizeOnline<OptimizerSGDL2> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
    } else if (learn_mira)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerMIRA> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
    else if (learn_arow)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerAROW> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
    else if (learn_cw)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerCW> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
    else if (learn_pegasos)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerPegasos> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
    else if (learn_nherd)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerNHERD> >(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights, generator);
    else if (learn_cp)
      objective = optimize_cp<OptimizeCP>(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights_history, weights);
    else if (learn_mcp)
      objective = optimize_cp<OptimizeMCP>(scorers, kbests, oracles, kbest_map, bounds_lower, bounds_upper, weights_history, weights);
    else if (learn_xbleu)
      objective = optimize_xbleu<OptimizeXBLEU>(kbests, scorers, weights);
    else 
      objective = optimize_batch<OptimizeLBFGS>(kbests, oracles, bounds_lower, bounds_upper, weights);

    if (debug && mpi_rank == 0)
      std::cerr << "objective: " << objective << std::endl;
    
    if (! bounds_lower.empty()) {
      const size_t weights_size = utils::bithack::min(weights.size(), bounds_lower.size());
      
      for (size_t i = 0; i != weights_size; ++ i)
	weights[i] = std::max(weights[i], bounds_lower[i]);
    }
    
    if (! bounds_upper.empty()) {
      const size_t weights_size = utils::bithack::min(weights.size(), bounds_upper.size());
      
      for (size_t i = 0; i != weights_size; ++ i)
	weights[i] = std::min(weights[i], bounds_upper[i]);
    }
    
    if (mert_search) {
      const double objective = optimize_mert(scorers, kbests, kbest_map, 0.1, 1.1, weights_prev, weights);
      
      if (debug && mpi_rank == 0)
	std::cerr << "mert objective: " << objective << std::endl;
    }
    
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
  notify_tag,
  termination_tag,
  point_tag,
  envelope_tag,
  score_tag,
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
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef hypothesis_type::feature_value_type feature_value_type;
  
  struct SampleSet
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    typedef std::vector<size_type, std::allocator<size_type> > offsets_type;

    struct Sample
    {
      typedef const feature_value_type* const_iterator;
      
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    features_type features;
    offsets_type  offsets;
  };
  
  typedef SampleSet sample_set_type;
  
  struct sample_pair_type
  {
    typedef std::vector<double, std::allocator<double> > loss_set_type;

    sample_pair_type() : features(), loss(),  offset(0) {}
    sample_pair_type(const hypothesis_set_type& kbests,
		     const hypothesis_set_type& oracles)
      : features(), loss(),  offset(0)
    {
      loss.reserve(kbests.size() + oracles.size());

      hypothesis_set_type::const_iterator oiter_end = oracles.end();
      for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter) {
	features.insert(oiter->features.begin(), oiter->features.end());
	loss.push_back(oiter->loss);
      }
      
      offset = loss.size();
      
      hypothesis_set_type::const_iterator kiter_end = kbests.end();
      for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
	features.insert(kiter->features.begin(), kiter->features.end());
	loss.push_back(kiter->loss);
      }
      
      features.shrink();
    }
    
    size_type oracle_begin() const { return 0; }
    size_type oracle_end() const { return offset; }
    
    size_type kbest_begin() const { return offset; }
    size_type kbest_end() const { return loss.size(); }
    
    size_type size() const { return loss.size(); }
    
    sample_set_type features;
    loss_set_type   loss;
    size_type offset;
  };
  
  typedef std::vector<sample_pair_type, std::allocator<sample_pair_type> > sample_pair_set_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  OptimizeOnline(Optimizer& __optimizer,
		 const hypothesis_map_type& kbests,
		 const hypothesis_map_type& oracles,
		 const weight_set_type& __bounds_lower,
		 const weight_set_type& __bounds_upper)
    : optimizer(__optimizer),
      bounds_lower(__bounds_lower),
      bounds_upper(__bounds_upper)
  {
    samples.reserve(kbests.size());
    
    const size_t id_max = utils::bithack::min(kbests.size(), oracles.size());
    
    for (size_t id = 0; id != id_max; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty())
	samples.push_back(sample_pair_type(kbests[id], oracles[id]));
    
    ids.reserve(samples.size());
    ids.resize(samples.size());
    for (size_t id = 0; id != ids.size(); ++ id)
      ids[id] = id;
  }
  
  template <typename Generator>
  void operator()(Generator& generator)
  {
    typedef std::vector<weight_type, std::allocator<weight_type> > margin_set_type;

    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);

    const double cost_factor = (softmax_margin ? 1.0 : 0.0);
    
    margin_set_type margins;
    
    gradient_type gradient_oracles;
    gradient_type gradient_kbests;
    
    id_set_type::const_iterator iiter_end = ids.end();
    for (id_set_type::const_iterator iiter = ids.begin(); iiter != iiter_end; ++ iiter) {
      const size_type id = *iiter;
      
      weight_type Z_oracle;
      weight_type Z_kbest;
      
      margins.clear();
      margins.resize(samples[id].size());
      
      for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	const typename sample_set_type::value_type features = samples[id].features[i];
	const double loss = samples[id].loss[i];
	
	margins[i] = function(features.begin(), features.end(), cost_factor * loss);
	
	Z_oracle += margins[i];
      }
      
      for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	const typename sample_set_type::value_type features = samples[id].features[i];
	const double loss = samples[id].loss[i];
	
	margins[i] = function(features.begin(), features.end(), cost_factor * loss);
	
	Z_kbest += margins[i];
      }
      
      gradient_oracles.clear();
      gradient_kbests.clear();
      
      for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	const typename sample_set_type::value_type features = samples[id].features[i];
	const weight_type weight = margins[i] / Z_oracle;
	
	typename sample_set_type::value_type::const_iterator fiter_end = features.end();
	for (typename sample_set_type::value_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
	  gradient_oracles[fiter->first] += weight_type(fiter->second) * weight;
      }
      
      for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	const typename sample_set_type::value_type features = samples[id].features[i];
	const weight_type weight = margins[i] / Z_kbest;
	
	typename sample_set_type::value_type::const_iterator fiter_end = features.end();
	for (typename sample_set_type::value_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
	  gradient_kbests[fiter->first] += weight_type(fiter->second) * weight;
      }
      
      
      if (! bounds_lower.empty() || ! bounds_upper.empty())
	optimizer(gradient_oracles,
		  gradient_kbests,
		  Z_oracle,
		  Z_kbest,
		  bounds_lower,
		  bounds_upper);
      else
	optimizer(gradient_oracles,
		  gradient_kbests,
		  Z_oracle,
		  Z_kbest);
    }
  }

  typedef Optimizer optimizer_type;
  
  typedef typename optimizer_type::weight_type   weight_type;
  typedef typename optimizer_type::gradient_type gradient_type;    
  
  double objective(const weight_set_type& weights) const
  {
    const double cost_factor = (softmax_margin ? 1.0 : 0.0);
    
    double obj = 0.0;

    for (size_t id = 0; id != samples.size(); ++ id) {
      weight_type Z_oracle;
      weight_type Z_kbest;
      
      for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	const typename sample_set_type::value_type features = samples[id].features[i];
	const double loss = samples[id].loss[i];
	
	Z_oracle += function(weights, features.begin(), features.end(), cost_factor * loss);
      }
      
      for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	const typename sample_set_type::value_type features = samples[id].features[i];
	const double loss = samples[id].loss[i];
	
	Z_kbest += function(weights, features.begin(), features.end(), cost_factor * loss);
      }
      
      obj -= log(Z_oracle) - log(Z_kbest);
    }
    
    return obj;
  }

  double normalizer() const
  {
    return samples.size();
  }
  
  template <typename Iterator>
  std::pair<double, double> operator()(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    return std::make_pair(0.0, 0.0);
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
  
  
  
  Optimizer& optimizer;

  sample_pair_set_type samples;

  const weight_set_type& bounds_lower;
  const weight_set_type& bounds_upper;

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
      typedef const feature_value_type* const_iterator;
      
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    void flush()
    {
      
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
  typedef typename utils::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> >::type sentence_unique_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  template <typename Iterator1, typename Iterator2, typename Features>
  void construct_pair(Iterator1 oiter, Iterator1 oiter_end,
		      Iterator2 kiter, Iterator2 kiter_end,
		      Features& feats)
  {
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
  }
    
  struct greater_loss
  {
    greater_loss(const loss_set_type&   __losses) : losses(__losses) {}

    bool operator()(const size_type& x, const size_type& y) const
    {
      return losses[x] > losses[y];
    }
      
    const loss_set_type&   losses;
  };
  
  OptimizeOnlineMargin(Optimizer& __optimizer,
		       const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       const weight_set_type& __bounds_lower,
		       const weight_set_type& __bounds_upper)
    : optimizer(__optimizer),
      bounds_lower(__bounds_lower),
      bounds_upper(__bounds_upper)
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    typedef std::vector<size_type, std::allocator<size_type> > pos_set_type;
    
    features_type feats;
    sentence_unique_type  sentences;
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> gen(generator);
    
    pos_set_type    positions;
    sample_set_type features_sample;
    loss_set_type   losses_sample;

    if (sample_vector) {
      for (size_type id = 0; id != kbests.size(); ++ id)
	if (! kbests[id].empty() && ! oracles[id].empty()) {
	  features_sample.clear();
	  losses_sample.clear();
	  
	  sentences.clear();
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    sentences.insert(oracles[id][o].sentence);
	  
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    for (size_t k = 0; k != kbests[id].size(); ++ k) {
	      const hypothesis_type& oracle = oracles[id][o];
	      const hypothesis_type& kbest  = kbests[id][k];
	      
	      // ignore oracle translations
	      if (sentences.find(kbest.sentence) != sentences.end()) continue;
	      
	      const double loss = kbest.loss - oracle.loss;
	      if (loss <= 0.0) continue;
	      
	      feats.clear();
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	      
	      if (feats.empty()) continue;
	      
	      features_sample.insert(feats.begin(), feats.end());
	      losses_sample.push_back(loss);
	    }
	  
	  // second, collect data from kbests onlly, which is the same as the first examples
	  const size_type sample_size = losses_sample.size();
	  const size_type sample_size_max = sample_size << 2;

	  while (losses_sample.size() < sample_size_max) {
	    const hypothesis_type& hyp1 = kbests[id][gen(kbests[id].size())];
	    const hypothesis_type& hyp2 = kbests[id][gen(kbests[id].size())];
	    
	    const hypothesis_type& kbest  = (hyp1.loss < hyp2.loss ? hyp2 : hyp1);
	    const hypothesis_type& oracle = (hyp1.loss < hyp2.loss ? hyp1 : hyp2);
	    
	    const double loss = kbest.loss - oracle.loss;
	    if (loss <= 1e-4) continue;
	    
	    feats.clear();
	    construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	    
	    if (feats.empty()) continue;
	    
	    features_sample.insert(feats.begin(), feats.end());
	    losses_sample.push_back(loss);
	  }

	  positions.clear();
	  for (size_type i = 0; i != losses_sample.size(); ++ i)
	    positions.push_back(i);

	  std::sort(positions.begin(), positions.end(), greater_loss(losses_sample));
	  
	  for (pos_set_type::const_iterator piter = positions.begin(); piter != positions.begin() + sample_size; ++ piter) {
	    features.insert(features_sample[*piter].begin(), features_sample[*piter].end());
	    losses.push_back(loss_margin ? losses_sample[*piter] : 1.0);
	  }

	  features.flush();
	}
      
    } else {
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
	      
	      feats.clear();
	      
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	    
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
	  
	  features.flush();
	}
    }
    
    features.shrink();
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
    
    if (! bounds_lower.empty() || ! bounds_upper.empty()) {
      for (size_t i = 0; i != ids.size(); ++ i) {
	const size_type id = ids[i];
	
	optimizer(features[id].begin(), features[id].end(), bounds_lower, bounds_upper, losses[id]);
      }
    } else {
      for (size_t i = 0; i != ids.size(); ++ i) {
	const size_type id = ids[i];
	
	optimizer(features[id].begin(), features[id].end(), losses[id]);
      }
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
  std::pair<double, double> operator()(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    static const double inf = std::numeric_limits<double>::infinity();

    double grad_pos = 0.0;
    double grad_neg = 0.0;
    for (size_t id = 0; id != losses.size(); ++ id) {
      const double margin      = cicada::dot_product(weights,      features[id].begin(), features[id].end(), 0.0);
      const double margin_prev = cicada::dot_product(weights_prev, features[id].begin(), features[id].end(), 0.0);
      
      const double bi_pos = margin_prev - margin;
      const double ci_pos = losses[id]  - margin_prev;
      const double ki_pos = (bi_pos != 0.0 ? - ci_pos / bi_pos : - inf);
      
      const double bi_neg = margin_prev + margin;
      const double ci_neg = losses[id]  - margin_prev;
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
    
    return std::make_pair(grad_pos, grad_neg);
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
  
  const weight_set_type& bounds_lower;
  const weight_set_type& bounds_upper;
  
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

template <typename Points>
void reduce_points(const int rank, Points& points)
{
  Points points_next;
  
  typename Points::const_iterator piter     = points.begin();
  typename Points::const_iterator piter_end = points.end();
  
  boost::iostreams::filtering_istream is;
  is.push(utils::mpi_device_source(rank, point_tag, 1024 * 1024));
  
  std::pair<double, double> point;
  
  while (is.read((char*) &point.first, sizeof(double)) && is.read((char*) &point.second, sizeof(double))) {
    for (/**/; piter != piter_end && *piter < point; ++ piter)
      points_next.push_back(*piter);
    points_next.push_back(point);
  }
  
  points_next.insert(points_next.end(), piter, piter_end);
  
  points.swap(points_next);
}

template <typename Iterator, typename Points>
void reduce_points(Iterator first, Iterator last, Points& points)
{
  Points points_next;
  
  for (/**/; first != last; ++ first) {
    typename Points::const_iterator piter     = points.begin();
    typename Points::const_iterator piter_end = points.end();
    
    boost::iostreams::filtering_istream is;
    is.push(utils::mpi_device_source(*first, point_tag, 1024 * 1024));
    
    std::pair<double, double> point;
    
    while (is.read((char*) &point.first, sizeof(double)) && is.read((char*) &point.second, sizeof(double))) {
      for (/**/; piter != piter_end && *piter < point; ++ piter)
	points_next.push_back(*piter);
      points_next.push_back(point);
    }
    
    points_next.insert(points_next.end(), piter, piter_end);
    
    points.swap(points_next);
    points_next.clear();
  }
}

template <typename Points>
void send_points(const int rank, const Points& points)
{
  boost::iostreams::filtering_ostream os;
  os.push(utils::mpi_device_sink(rank, point_tag, 1024 * 1024));
  
  typename Points::const_iterator piter_end = points.end();
  for (typename Points::const_iterator piter = points.begin(); piter != piter_end; ++ piter) {
    os.write((char*) &(piter->first), sizeof(double));
    os.write((char*) &(piter->second), sizeof(double));
  }
}


template <typename Points>
void reduce_points(Points& points)
{
  typedef std::vector<int, std::allocator<int> > rank_set_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  rank_set_type ranks;
  int merge_size = mpi_size;
  
  while (merge_size > 1 && mpi_rank < merge_size) {
    const int reduce_size = (merge_size / 2 == 0 ? 1 : merge_size / 2);
    
    if (mpi_rank < reduce_size) {
      ranks.clear();
      for (int i = reduce_size; i < merge_size; ++ i)
	if (i % reduce_size == mpi_rank)
	  ranks.push_back(i);
      
      if (ranks.empty()) continue;
      
      if (ranks.size() == 1)
	reduce_points(ranks.front(), points);
      else
	reduce_points(ranks.begin(), ranks.end(), points);
      
    } else
      send_points(mpi_rank % reduce_size, points);
    
    merge_size = reduce_size;
  }
}

template <typename Optimize, typename Generator>
double optimize_online(const scorer_document_type& scorers,
		       const hypothesis_map_type& kbests,
		       const hypothesis_map_type& oracles,
		       const kbest_map_type& kbest_map,
		       const weight_set_type& bounds_lower,
		       const weight_set_type& bounds_upper,
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
  Optimize opt(optimizer, kbests, oracles, bounds_lower, bounds_upper);

  const double norm_local = opt.normalizer();
  double norm = 0.0;
  MPI::COMM_WORLD.Reduce(&norm_local, &norm, 1, MPI::DOUBLE, MPI::SUM, 0);
  
  weight_set_type weights_init = weights;
  point_set_type points;
  
  optimizer.weights = weights;

  bcast_weights(0, optimizer.weights);
  
  if (mpi_rank == 0) {
    weight_set_type weights_min;
    double objective_min = std::numeric_limits<double>::infinity();
    
    double objective_prev = 0.0;
    double objective = 0.0;
    
    int increased = 0;
    
    for (int iter = 0; iter < iteration; ++ iter) {
      
      for (int rank = 1; rank < mpi_size; ++ rank)
	MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
      
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
      if (! bounds_lower.empty()) {
	const size_t weights_size = utils::bithack::min(weights.size(), bounds_lower.size());
	
	for (size_t i = 0; i != weights_size; ++ i)
	  weights[i] = std::max(weights[i], bounds_lower[i]);
      }
      
      if (! bounds_upper.empty()) {
	const size_t weights_size = utils::bithack::min(weights.size(), bounds_upper.size());
	
	for (size_t i = 0; i != weights_size; ++ i)
	  weights[i] = std::min(weights[i], bounds_upper[i]);
      }


      if (line_search) {
	// perform line-search between weights_prev and optimizer.weights, and update optimizer.weights
	
	if (debug >= 3)
	  std::cerr << "line-search" << std::endl;
	
	bcast_weights(0, optimizer.weights);
	
	points.clear();
	
	const std::pair<double, double> grads = opt(optimizer.weights, weights_prev, std::back_inserter(points));
	
	std::sort(points.begin(), points.end());
	
	reduce_points(points);

	if (debug >= 3)
	  std::cerr << "point size: " << points.size() << std::endl;
	
	double grad_pos = 0.0;
	double grad_neg = 0.0;
	
	MPI::COMM_WORLD.Reduce(&grads.first,  &grad_pos, 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&grads.second, &grad_neg, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	const double norm_w      = cicada::dot_product(optimizer.weights, optimizer.weights);
	const double dot_prod    = cicada::dot_product(weights_prev, optimizer.weights);
	const double norm_w_prev = cicada::dot_product(weights_prev, weights_prev);
	
	const double a0_pos = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * norm;
	const double b0_pos = (dot_prod - norm_w_prev) * C * norm;
	
	const double a0_neg = (norm_w + 2.0 * dot_prod + norm_w_prev) * C * norm;
	const double b0_neg = (- dot_prod - norm_w_prev) * C * norm;

	//std::cerr << "a0: "  << a0 << " b0: " << b0 << std::endl;
	
	grad_pos += b0_pos;
	grad_neg += b0_neg;
	
	if (debug >= 3)
	  std::cerr << "gradient: " << grad_pos << ' ' << grad_neg << std::endl;
	
	if (grad_pos < 0.0) {
	  double k = 0.0;
	  
	  point_set_type::const_iterator piter = std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0));
	  point_set_type::const_iterator piter_end = points.end();

	  for (/**/; piter != piter_end && grad_pos < 0.0; ++ piter) {
	    const double k_new = piter->first;
	    const double grad_new = grad_pos + std::fabs(piter->second) + a0_pos * (k_new - k);
	    
	    if (grad_new >= 0) {
	      // compute intersection...
	      k = k + grad_pos * (k - k_new) / (grad_new - grad_pos);
	      grad_pos = grad_new;
	      break;
	    } else {
	      k = k_new;
	      grad_pos = grad_new;
	    }
	  }

	  if (debug >= 3)
	    std::cerr << "grad: " << grad_pos << "  k: " << k << std::endl;
	  
	  if (k > 0.0) {
	    const size_t weights_size = utils::bithack::min(optimizer.weights.size(), weights_prev.size());
	    
	    for (size_t i = 0; i != weights_size; ++ i)
	      optimizer.weights[i] = k * optimizer.weights[i] + (1.0 - k) * weights_prev[i];
	    for (size_t i = weights_size; i < optimizer.weights.size(); ++ i)
	      optimizer.weights[i] = k * optimizer.weights[i];
	    for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	      optimizer.weights[i] = (1.0 - k) * weights_prev[i];
	  }
	  
	  
	} else if (grad_neg < 0.0) {
	  double k = 0.0;
	  
	  point_set_type::const_reverse_iterator piter(std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0)));
	  point_set_type::const_reverse_iterator piter_end = points.rend();
	  
	  for (/**/; piter != piter_end && grad_neg < 0.0; ++ piter) {
	    const double k_new = - piter->first;
	    const double grad_new = grad_neg + std::fabs(piter->second) + a0_neg * (k_new - k);
	    
	    if (grad_new >= 0) {
	      // compute intersection...
	      k = k + grad_neg * (k - k_new) / (grad_new - grad_neg);
	      grad_neg = grad_new;
	      break;
	    } else {
	      k = k_new;
	      grad_neg = grad_new;
	    }
	  }
	  
	  if (debug >= 3)
	    std::cerr << "grad: " << grad_neg << "  k: " << - k << std::endl;
	  
	  if (k > 0.0) {
	    const size_t weights_size = utils::bithack::min(optimizer.weights.size(), weights_prev.size());
	    
	    for (size_t i = 0; i != weights_size; ++ i)
	      optimizer.weights[i] = - k * optimizer.weights[i] + (1.0 + k) * weights_prev[i];
	    for (size_t i = weights_size; i < optimizer.weights.size(); ++ i)
	      optimizer.weights[i] = - k * optimizer.weights[i];
	    for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	      optimizer.weights[i] = (1.0 + k) * weights_prev[i];
	  }
	}
      }

      if (mert_search_local)
	optimize_mert(scorers, kbests, kbest_map, 0.01, 2.0, weights_prev, optimizer.weights);
      
      // compute objective
      bcast_weights(0, optimizer.weights);
      
      const double objective_local = opt.objective(optimizer.weights);
      
      objective = 0.0;
      MPI::COMM_WORLD.Reduce(&objective_local, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      objective /= norm;
      
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
      
      //const bool converged = (active_size == 0 || (iter && (std::fabs((objective - objective_prev) / objective) < 1e-4)) || increased > 16);
      const bool converged = (active_size == 0 || (iter && (std::fabs((objective - objective_prev) / objective) < 1e-5)));
      
      if (debug >= 2)
	std::cerr << "objective: " << objective << " active size: " << active_size << std::endl;
      
      if (objective <= objective_min) {
	objective_min = objective;
	weights_min = optimizer.weights;
      }
      
      if (converged) break;
      
      objective_prev = objective;
    }
    
    // send termination!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, termination_tag);

    weights.swap(weights_min);
    
    return objective_min;
    
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
	
	if (line_search) {
	  // perform line-search between weights_prev and optimizer.weights, and update optimizer.weights
	  
	  bcast_weights(0, optimizer.weights);
	  
	  points.clear();
	  const std::pair<double, double> grads = opt(optimizer.weights, weights_prev, std::back_inserter(points));
	  
	  std::sort(points.begin(), points.end());
	  
	  reduce_points(points);
	  
	  double grad_pos = 0.0;
	  double grad_neg = 0.0;
	  
	  MPI::COMM_WORLD.Reduce(&grads.first,  &grad_pos, 1, MPI::DOUBLE, MPI::SUM, 0);
	  MPI::COMM_WORLD.Reduce(&grads.second, &grad_neg, 1, MPI::DOUBLE, MPI::SUM, 0);
	}

	if (mert_search_local)
	  optimize_mert(scorers, kbests, kbest_map, 0.01, 2.0, weights_prev, optimizer.weights);
	
	// compute objective
	bcast_weights(0, optimizer.weights);
	
	const double objective_local = opt.objective(optimizer.weights);
	objective = 0.0;
	MPI::COMM_WORLD.Reduce(&objective_local, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();

    return 0.0;
  }

}

struct OptimizeCP
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
      typedef const feature_value_type* const_iterator;
      
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    void flush()
    {
      
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
  typedef utils::unordered_set<hypothesis_type::sentence_type, hash_sentence, std::equal_to<hypothesis_type::sentence_type>, std::allocator<hypothesis_type::sentence_type> >::type sentence_unique_type;

  template <typename Iterator1, typename Iterator2, typename Features>
  void construct_pair(Iterator1 oiter, Iterator1 oiter_end,
		      Iterator2 kiter, Iterator2 kiter_end,
		      Features& feats)
  {
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
  }
    
  struct greater_loss
  {
    greater_loss(const loss_set_type&   __losses) : losses(__losses) {}

    bool operator()(const size_type& x, const size_type& y) const
    {
      return losses[x] > losses[y];
    }
      
    const loss_set_type&   losses;
  };

  OptimizeCP(const scorer_document_type& scorers,
	     const hypothesis_map_type& kbests,
	     const hypothesis_map_type& oracles,
	     const kbest_map_type& kbest_map) 
  {
    typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
    typedef std::vector<size_type, std::allocator<size_type> > pos_set_type;
    
    features_type feats;
    sentence_unique_type  sentences;
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    boost::random_number_generator<boost::mt19937> gen(generator);
    
    pos_set_type    positions;
    sample_set_type features_sample;
    loss_set_type   losses_sample;
    
    if (sample_vector) {
      for (size_type id = 0; id != kbests.size(); ++ id)
	if (! kbests[id].empty() && ! oracles[id].empty()) {
	  features_sample.clear();
	  losses_sample.clear();
	  
	  sentences.clear();
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    sentences.insert(oracles[id][o].sentence);
	  
	  for (size_t o = 0; o != oracles[id].size(); ++ o)
	    for (size_t k = 0; k != kbests[id].size(); ++ k) {
	      const hypothesis_type& oracle = oracles[id][o];
	      const hypothesis_type& kbest  = kbests[id][k];
	      
	      // ignore oracle translations
	      if (sentences.find(kbest.sentence) != sentences.end()) continue;
	      
	      const double loss = kbest.loss - oracle.loss;
	      if (loss <= 0.0) continue;
	      
	      feats.clear();
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	      
	      if (feats.empty()) continue;
	      
	      features_sample.insert(feats.begin(), feats.end());
	      losses_sample.push_back(loss);
	    }
	  
	  // second, collect data from kbests onlly, which is the same as the first examples
	  const size_type sample_size = losses_sample.size();
	  const size_type sample_size_max = sample_size << 2;
	  
	  while (losses_sample.size() < sample_size_max) {
	    const hypothesis_type& hyp1 = kbests[id][gen(kbests[id].size())];
	    const hypothesis_type& hyp2 = kbests[id][gen(kbests[id].size())];
	    
	    const hypothesis_type& kbest  = (hyp1.loss < hyp2.loss ? hyp2 : hyp1);
	    const hypothesis_type& oracle = (hyp1.loss < hyp2.loss ? hyp1 : hyp2);
	    
	    const double loss = kbest.loss - oracle.loss;
	    if (loss <= 1e-4) continue;
	    
	    feats.clear();
	    construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	    
	    if (feats.empty()) continue;
	    
	    features_sample.insert(feats.begin(), feats.end());
	    losses_sample.push_back(loss);
	  }

	  positions.clear();
	  for (size_type i = 0; i != losses_sample.size(); ++ i)
	    positions.push_back(i);

	  std::sort(positions.begin(), positions.end(), greater_loss(losses_sample));
	  
	  for (pos_set_type::const_iterator piter = positions.begin(); piter != positions.begin() + sample_size; ++ piter) {
	    features.insert(features_sample[*piter].begin(), features_sample[*piter].end());
	    losses.push_back(loss_margin ? losses_sample[*piter] : 1.0);
	  }

	  features.flush();
	}
      
    } else {
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
	      
	      feats.clear();
	      
	      construct_pair(oracle.features.begin(), oracle.features.end(), kbest.features.begin(), kbest.features.end(), feats);
	      
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
	  
	  features.flush();
	}
    }
    
    features.shrink();
    loss_set_type(losses).swap(losses);
  }
  

  std::pair<double, score_ptr_pair_type> operator()(const weight_set_type& weights, weight_set_type& acc) const
  {
    const double factor = 1.0 / samples;
    
    double objective = 0.0;
    for (size_t id = 0; id != losses.size(); ++ id) {
      const double margin = cicada::dot_product(weights, features[id].begin(), features[id].end(), 0.0);
      const double loss = losses[id];
      
      const double suffered = loss - margin;
      
      if (suffered > 0.0) { 
	sample_set_type::value_type::const_iterator fiter_end = features[id].end();
	for (sample_set_type::value_type::const_iterator fiter = features[id].begin(); fiter != fiter_end; ++ fiter) 
	  acc[fiter->first] += factor * fiter->second;
	
	objective += factor * suffered;
      }
    }
    
    return std::make_pair(objective, score_ptr_pair_type());
  }

  std::pair<double, score_ptr_pair_type> objective(const weight_set_type& weights) const
  {
    const double factor = 1.0 / samples;
    
    double obj = 0.0;
    for (size_t id = 0; id != losses.size(); ++ id) {
      const double margin = cicada::dot_product(weights, features[id].begin(), features[id].end(), 0.0);
      const double loss = losses[id];
      
      obj += (loss - margin) * (loss - margin > 0.0);
    }
    
    return std::make_pair(obj * factor, score_ptr_pair_type());
  }

  template <typename Iterator>
  std::pair<double, double> operator()(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    static const double inf = std::numeric_limits<double>::infinity();

    double grad_pos = 0.0;
    double grad_neg = 0.0;
    for (size_t id = 0; id != losses.size(); ++ id) {
      const double margin      = cicada::dot_product(weights,      features[id].begin(), features[id].end(), 0.0);
      const double margin_prev = cicada::dot_product(weights_prev, features[id].begin(), features[id].end(), 0.0);
      
      const double bi_pos = margin_prev - margin;
      const double ci_pos = losses[id]  - margin_prev;
      const double ki_pos = (bi_pos != 0.0 ? - ci_pos / bi_pos : - inf);
      
      const double bi_neg = margin_prev + margin;
      const double ci_neg = losses[id]  - margin_prev;
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
    
    return std::make_pair(grad_pos, grad_neg);
  }
 
  double instances() const
  {
    return losses.size();
  }
  
  template <typename Features>
  struct HMatrix
  {
    HMatrix(const Features& __features) : features(__features) {}

    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[i], features[j]);
    }
    
    const Features& features;
  };
  
  template <typename Features>
  struct MMatrix
  {
    MMatrix(const Features& __features) : features(__features) {}
    
    template <typename W>
    void operator()(W& w, const std::vector<double, std::allocator<double> >& alpha) const
    {
      std::vector<double, std::allocator<double> >::const_iterator aiter = alpha.begin();
      
      for (size_type id = 0; id != features.size(); ++ id, ++ aiter)
	if (*aiter > 0.0)
	  operator()(w, *aiter, id);
    }
    
    template <typename W>
    double operator()(const W& w, const size_t& i) const
    {
      return cicada::dot_product(w, features[i]);
    }
    
    template <typename W>
    void operator()(W& w, const double& update, const size_t& i) const
    {
      for (size_t j = 0; j != features[i].size(); ++ j)
	w[j] += update * features[i][j];
    }

    const Features& features;
  };
  
  sample_set_type features;
  loss_set_type   losses;
  double samples;
};

struct OptimizeMCP
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
      typedef const feature_value_type* const_iterator;
      
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    void flush()
    {
      
    }
    
    features_type features;
    offsets_type  offsets;
  };
  
  typedef SampleSet sample_set_type;

    
  OptimizeMCP(const scorer_document_type& __scorers,
	      const hypothesis_map_type&  __kbests,
	      const hypothesis_map_type&  __oracles,
	      const kbest_map_type&       __kbest_map) 
    : scorers(__scorers),
      kbests(__kbests),
      oracles(__oracles),
      kbest_map(__kbest_map)
  {
    // transform kbests/oracles as an expanded features for faster learning?

    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const hypothesis_type& kbest = *kiter;
	  
	  features_kbest.insert(kbest.features.begin(), kbest.features.end());
	}
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const hypothesis_type& oracle = *oiter;
	  
	  features_oracle.insert(oracle.features.begin(), oracle.features.end());
	}
      }
    
    features_kbest.shrink();
    features_oracle.shrink();
  }
  
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  typedef std::vector<double, std::allocator<double> > margin_set_type;
  typedef std::vector<size_t, std::allocator<size_t> > pos_set_type;
  typedef std::vector<const hypothesis_type*, std::allocator<const hypothesis_type*> > hypothesis_ptr_set_type;
  
  margin_set_type         kbests_margin;
  margin_set_type         oracles_margin;
  hypothesis_ptr_set_type kbests_hyp;
  hypothesis_ptr_set_type oracles_hyp;
  pos_set_type            kbests_pos;
  pos_set_type            oracles_pos;
  loss_set_type           oracles_loss;
  
  std::pair<double, score_ptr_pair_type> operator()(const weight_set_type& weights, weight_set_type& acc)
  {
#if 0
    const double factor = 1.0 / samples;
    const double inf = std::numeric_limits<double>::infinity();
    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());
    
    score_ptr_pair_type score;
    double loss = 0.0;
    double margin = 0.0;
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	const hypothesis_type* ptr_kbest = 0;
	const hypothesis_type* ptr_oracle = 0;
	double margin_kbest  = - inf;
	double margin_oracle = - inf;
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const hypothesis_type& kbest = *kiter;
	  
	  //const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), kbest.loss);
	  const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), 0.0);
	  
	  if (! ptr_kbest || margin > margin_kbest) {
	    margin_kbest = margin;
	    ptr_kbest = &kbest;
	  }
	}
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const hypothesis_type& oracle = *oiter;
	  
	  if (oracle.loss > ptr_kbest->loss) continue;
	  
	  //const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), - oracle.loss);
	  const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), 0.0);
	  
	  if (! ptr_oracle || margin > margin_oracle) {
	    margin_oracle = margin;
	    ptr_oracle = &oracle;
	  }
	}
	
	if (! direct_loss && ptr_oracle->score) {
	  if (! score.first)
	    score.first = ptr_oracle->score->clone();
	  else
	    *score.first += *(ptr_oracle->score);
	}
	
	if (ptr_kbest->score) {
	  if (! score.second)
	    score.second = ptr_kbest->score->clone();
	  else
	    *score.second += *(ptr_kbest->score);
	}
	
	//margin += (margin_oracle + ptr_oracle->loss) - (margin_kbest - ptr_kbest->loss);
	//margin += (margin_oracle + ptr_oracle->loss) - margin_kbest;
	margin += margin_oracle - margin_kbest;
	
	{
	  hypothesis_type::feature_set_type::const_iterator kiter_end = ptr_kbest->features.end();
	  for (hypothesis_type::feature_set_type::const_iterator kiter = ptr_kbest->features.begin(); kiter != kiter_end; ++ kiter)
	    acc[kiter->first] -= kiter->second * factor;
	  
	  hypothesis_type::feature_set_type::const_iterator oiter_end = ptr_oracle->features.end();
	  for (hypothesis_type::feature_set_type::const_iterator oiter = ptr_oracle->features.begin(); oiter != oiter_end; ++ oiter)
	    acc[oiter->first] += oiter->second * factor;
	}
      }
    
    return std::make_pair((loss - margin) * factor, score);
#endif
#if 0
    kbests_margin.clear();
    kbests_hyp.clear();
    
    oracles_margin.clear();
    oracles_hyp.clear();

    const double factor = 1.0 / samples;
    const double inf = std::numeric_limits<double>::infinity();
    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	
	const size_type seg = kbest_map[id];
	
	if (seg >= kbests_hyp.size())
	  kbests_hyp.resize(seg + 1, 0);
	if (seg >= kbests_margin.size())
	  kbests_margin.resize(seg + 1, - inf);

	if (seg >= oracles_hyp.size())
	  oracles_hyp.resize(seg + 1, 0);
	if (seg >= oracles_margin.size())
	  oracles_margin.resize(seg + 1, - inf);

	bool updated = false;
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const hypothesis_type& kbest = *kiter;
	  
	  //const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), kbest.loss);
	  const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), 0.0);
	  
	  if (margin > kbests_margin[seg] || ! kbests_hyp[seg]) {
	    kbests_margin[seg] = margin;
	    kbests_hyp[seg] = &kbest;
	    
	    updated = true;
	  }
	}
	
	if (! updated) continue;
	
	const double kbest_loss = (kbests_hyp[seg] ? kbests_hyp[seg]->loss : inf);
	const double kbest_margin = kbests_margin[seg];
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const hypothesis_type& oracle = *oiter;
	  
	  if (oracle.loss > kbest_loss) continue;
	  
	  //const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), - oracle.loss);
	  const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), 0.0);
	  
	  if (margin > oracles_margin[seg] || ! oracles_hyp[seg]) {
	    oracles_margin[seg] = margin;
	    oracles_hyp[seg] = &oracle;
	  }
	}
      }
    
    score_ptr_pair_type score;
    double loss = 0.0;
    double margin = 0.0;
    
    for (size_type seg = 0; seg != kbests_hyp.size(); ++ seg)
      if (oracles_hyp[seg] && kbests_hyp[seg]) {

	if (! direct_loss && oracles_hyp[seg]->score) {
	  if (! score.first)
	    score.first = oracles_hyp[seg]->score->clone();
	  else
	    *score.first += *(oracles_hyp[seg]->score);
	}
	
	if (kbests_hyp[seg]->score) {
	  if (! score.second)
	    score.second = kbests_hyp[seg]->score->clone();
	  else
	    *score.second += *(kbests_hyp[seg]->score);
	}
	
	//margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - (kbests_margin[seg] - kbests_hyp[seg]->loss);
	//margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - kbests_margin[seg];
	margin += oracles_margin[seg] - kbests_margin[seg];
	
	hypothesis_type::feature_set_type::const_iterator kiter_end = kbests_hyp[seg]->features.end();
	for (hypothesis_type::feature_set_type::const_iterator kiter = kbests_hyp[seg]->features.begin(); kiter != kiter_end; ++ kiter)
	  acc[kiter->first] -= kiter->second * factor;
	
	hypothesis_type::feature_set_type::const_iterator oiter_end = oracles_hyp[seg]->features.end();
	for (hypothesis_type::feature_set_type::const_iterator oiter = oracles_hyp[seg]->features.begin(); oiter != oiter_end; ++ oiter)
	  acc[oiter->first] += oiter->second * factor;
      }
    
    return std::make_pair((loss - margin) * factor, score);
#endif
#if 1
    kbests_margin.clear();
    kbests_hyp.clear();
    kbests_pos.clear();
    
    oracles_margin.clear();
    oracles_hyp.clear();
    oracles_pos.clear();
    oracles_loss.clear();

    const double factor = 1.0 / samples;
    const double inf = std::numeric_limits<double>::infinity();
    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());

    size_t pos_kbest = 0;
    size_t pos_oracle = 0;
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	
	const size_type seg = kbest_map[id];
	
	if (seg >= kbests_hyp.size())
	  kbests_hyp.resize(seg + 1, 0);
	if (seg >= kbests_margin.size())
	  kbests_margin.resize(seg + 1, - inf);
	if (seg >= kbests_pos.size())
	  kbests_pos.resize(seg + 1);
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter, ++ pos_kbest) {
	  const hypothesis_type& kbest = *kiter;

	  //const double margin = cicada::dot_product(weights, features_kbest[pos_kbest].begin(), features_kbest[pos_kbest].end(), kbest.loss);
	  const double margin = cicada::dot_product(weights, features_kbest[pos_kbest].begin(), features_kbest[pos_kbest].end(), 0.0);
	  
	  if (! kbests_hyp[seg] || margin > kbests_margin[seg]) {
	    kbests_margin[seg] = margin;
	    kbests_hyp[seg] = &kbest;
	    kbests_pos[seg] = pos_kbest;
	  }
	}
      }
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	const size_type seg = kbest_map[id];
	
	if (seg >= oracles_hyp.size())
	  oracles_hyp.resize(seg + 1, 0);
	if (seg >= oracles_margin.size())
	  oracles_margin.resize(seg + 1, - inf);
	if (seg >= oracles_pos.size())
	  oracles_pos.resize(seg + 1);
	if (seg >= oracles_loss.size())
	  oracles_loss.resize(seg + 1, inf);
	
	const double kbest_loss = (kbests_hyp[seg] ? kbests_hyp[seg]->loss : inf);
	const double kbest_margin = kbests_margin[seg];
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter, ++ pos_oracle) {
	  const hypothesis_type& oracle = *oiter;
	  
	  if (oracle.loss > kbest_loss) continue;
	  //if (oracle.loss > oracles_loss[seg]) continue;
	  
	  //const double margin = cicada::dot_product(weights, features_oracle[pos_oracle].begin(), features_oracle[pos_oracle].end(), 0.0);
	  const double margin = cicada::dot_product(weights, features_oracle[pos_oracle].begin(), features_oracle[pos_oracle].end(), - oracle.loss);
	  
	  if (! oracles_hyp[seg] || margin > oracles_margin[seg]) {
	    oracles_margin[seg] = margin;
	    oracles_hyp[seg] = &oracle;
	    oracles_pos[seg] = pos_oracle;
	    oracles_loss[seg] = oracle.loss;
	  }
	}
      }
    
    if (oracles_margin.size() != kbests_margin.size())
      throw std::runtime_error("margin size differ");
    if (oracles_hyp.size() != kbests_hyp.size())
      throw std::runtime_error("margin size differ");
    if (oracles_pos.size() != kbests_pos.size())
      throw std::runtime_error("margin size differ");
    
    score_ptr_pair_type score;
    double loss = 0.0;
    double margin = 0.0;
    
    for (size_type seg = 0; seg != kbests_hyp.size(); ++ seg)
      if (oracles_hyp[seg] && kbests_hyp[seg]) {

	if (! direct_loss && oracles_hyp[seg]->score) {
	  if (! score.first)
	    score.first = oracles_hyp[seg]->score->clone();
	  else
	    *score.first += *(oracles_hyp[seg]->score);
	}
	
	if (kbests_hyp[seg]->score) {
	  if (! score.second)
	    score.second = kbests_hyp[seg]->score->clone();
	  else
	    *score.second += *(kbests_hyp[seg]->score);
	}
	
	//margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - (kbests_margin[seg] - kbests_hyp[seg]->loss);
	margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - kbests_margin[seg];
	//margin += oracles_margin[seg] - kbests_margin[seg];
	
	sample_set_type::value_type::const_iterator kiter_end = features_kbest[kbests_pos[seg]].end();
	for (sample_set_type::value_type::const_iterator kiter = features_kbest[kbests_pos[seg]].begin(); kiter != kiter_end; ++ kiter)
	  acc[kiter->first] -= kiter->second * factor;

	sample_set_type::value_type::const_iterator oiter_end = features_oracle[oracles_pos[seg]].end();
	for (sample_set_type::value_type::const_iterator oiter = features_oracle[oracles_pos[seg]].begin(); oiter != oiter_end; ++ oiter)
	  acc[oiter->first] += oiter->second * factor;
      }
    
    return std::make_pair((loss - margin) * factor, score);
#endif
  }

  std::pair<double, score_ptr_pair_type> objective(const weight_set_type& weights)
  {
#if 0
    const double factor = 1.0 / samples;
    const double inf = std::numeric_limits<double>::infinity();
    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());
    
    score_ptr_pair_type score;
    double loss = 0.0;
    double margin = 0.0;
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	const hypothesis_type* ptr_kbest = 0;
	const hypothesis_type* ptr_oracle = 0;
	double margin_kbest  = - inf;
	double margin_oracle = - inf;
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const hypothesis_type& kbest = *kiter;
	  
	  //const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), kbest.loss);
	  const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), 0.0);
	  
	  if (! ptr_kbest || margin > margin_kbest) {
	    margin_kbest = margin;
	    ptr_kbest = &kbest;
	  }
	}
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const hypothesis_type& oracle = *oiter;
	  
	  if (oracle.loss > ptr_kbest->loss) continue;
	  
	  //const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), - oracle.loss);
	  const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), 0.0);

	  if (! ptr_oracle || margin > margin_oracle) {
	    margin_oracle = margin;
	    ptr_oracle = &oracle;
	  }
	}

	if (! direct_loss && ptr_oracle->score) {
	  if (! score.first)
	    score.first = ptr_oracle->score->clone();
	  else
	    *score.first += *(ptr_oracle->score);
	}

	if (ptr_kbest->score) {
	  if (! score.second)
	    score.second = ptr_kbest->score->clone();
	  else
	    *score.second += *(ptr_kbest->score);
	}
	
	//margin += (margin_oracle + ptr_oracle->loss) - (margin_kbest - ptr_kbest->loss);
	//margin += (margin_oracle + ptr_oracle->loss) - margin_kbest;
	margin += margin_oracle - margin_kbest;
      }
    
    return std::make_pair((loss - margin) * factor, score);
#endif
#if 0
    kbests_margin.clear();
    kbests_hyp.clear();
    
    oracles_margin.clear();
    oracles_hyp.clear();
    
    const double factor = 1.0 / samples;
    const double inf = std::numeric_limits<double>::infinity();
    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	
	const size_type seg = kbest_map[id];
	
	if (seg >= kbests_hyp.size())
	  kbests_hyp.resize(seg + 1, 0);
	if (seg >= kbests_margin.size())
	  kbests_margin.resize(seg + 1, - inf);
	
	if (seg >= oracles_hyp.size())
	  oracles_hyp.resize(seg + 1, 0);
	if (seg >= oracles_margin.size())
	  oracles_margin.resize(seg + 1, - inf);
	
	bool updated = false;
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter) {
	  const hypothesis_type& kbest = *kiter;

	  //const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), kbest.loss);
	  const double margin = cicada::dot_product(weights, kbest.features.begin(), kbest.features.end(), 0.0);
	  
	  if (margin > kbests_margin[seg] || ! kbests_hyp[seg]) {
	    kbests_margin[seg] = margin;
	    kbests_hyp[seg] = &kbest;
	    
	    updated = true;
	  }
	}

	if (! updated) continue;
	
	const double kbest_loss = (kbests_hyp[seg] ? kbests_hyp[seg]->loss : inf);
	const double kbest_margin = kbests_margin[seg];
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	  const hypothesis_type& oracle = *oiter;

	  if (oracle.loss > kbest_loss) continue;
	  
	  //const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), - oracle.loss);
	  const double margin = cicada::dot_product(weights, oracle.features.begin(), oracle.features.end(), 0.0);
	  
	  if (margin > oracles_margin[seg] || ! oracles_hyp[seg]) {
	    oracles_margin[seg] = margin;
	    oracles_hyp[seg] = &oracle;
	  }
	}
      }
    
    score_ptr_pair_type score;
    double loss = 0.0;
    double margin = 0.0;
    
    for (size_type seg = 0; seg != kbests_hyp.size(); ++ seg)
      if (oracles_hyp[seg] && kbests_hyp[seg]) {
	
	if (! direct_loss && oracles_hyp[seg]->score) {
	  if (! score.first)
	    score.first = oracles_hyp[seg]->score->clone();
	  else
	    *score.first += *(oracles_hyp[seg]->score);
	}
	
	if (kbests_hyp[seg]->score) {
	  if (! score.second)
	    score.second = kbests_hyp[seg]->score->clone();
	  else
	    *score.second += *(kbests_hyp[seg]->score);
	}
	
	//margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - (kbests_margin[seg] - kbests_hyp[seg]->loss);
	//margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - kbests_margin[seg];
	margin += oracles_margin[seg] - kbests_margin[seg];
      }
    
    return std::make_pair((loss - margin) * factor, score);
#endif
#if 1
    kbests_margin.clear();
    kbests_hyp.clear();
    
    oracles_margin.clear();
    oracles_hyp.clear();
    oracles_loss.clear();
    
    const double factor = 1.0 / samples;
    const double inf = std::numeric_limits<double>::infinity();
    const size_t id_size = utils::bithack::min(kbests.size(), oracles.size());
    
    size_t pos_kbest = 0;
    size_t pos_oracle = 0;
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	
	const size_type seg = kbest_map[id];
	
	if (seg >= kbests_hyp.size())
	  kbests_hyp.resize(seg + 1, 0);
	if (seg >= kbests_margin.size())
	  kbests_margin.resize(seg + 1, - inf);
	
	hypothesis_set_type::const_iterator kiter_end = kbests[id].end();
	for (hypothesis_set_type::const_iterator kiter = kbests[id].begin(); kiter != kiter_end; ++ kiter, ++ pos_kbest) {
	  const hypothesis_type& kbest = *kiter;

	  //const double margin = cicada::dot_product(weights, features_kbest[pos_kbest].begin(), features_kbest[pos_kbest].end(), kbest.loss);
	  const double margin = cicada::dot_product(weights, features_kbest[pos_kbest].begin(), features_kbest[pos_kbest].end(), 0.0);
	  
	  if (! kbests_hyp[seg] || margin > kbests_margin[seg]) {
	    kbests_margin[seg] = margin;
	    kbests_hyp[seg] = &kbest;
	  }
	}
      }
    
    for (size_t id = 0; id != id_size; ++ id) 
      if (! kbests[id].empty() && ! oracles[id].empty()) {
	const size_type seg = kbest_map[id];
	
	if (seg >= oracles_hyp.size())
	  oracles_hyp.resize(seg + 1, 0);
	if (seg >= oracles_margin.size())
	  oracles_margin.resize(seg + 1, - inf);
	if (seg >= oracles_loss.size())
	  oracles_loss.resize(seg + 1, inf);
	  
	const double kbest_loss = (kbests_hyp[seg] ? kbests_hyp[seg]->loss : inf);
	const double kbest_margin = kbests_margin[seg];
	
	hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter, ++ pos_oracle) {
	  const hypothesis_type& oracle = *oiter;

	  if (oracle.loss > kbest_loss) continue;
	  //if (oracle.loss > oracles_loss[seg]) continue;
	  
	  //const double margin = cicada::dot_product(weights, features_oracle[pos_oracle].begin(), features_oracle[pos_oracle].end(), 0.0);
	  const double margin = cicada::dot_product(weights, features_oracle[pos_oracle].begin(), features_oracle[pos_oracle].end(), - oracle.loss);
	  
	  if (! oracles_hyp[seg] || margin > oracles_margin[seg]) {
	    oracles_margin[seg] = margin;
	    oracles_hyp[seg] = &oracle;
	    oracles_loss[seg] = oracle.loss;
	  }
	}
      }

    if (oracles_margin.size() != kbests_margin.size())
      throw std::runtime_error("margin size differ");
    if (oracles_hyp.size() != kbests_hyp.size())
      throw std::runtime_error("margin size differ");
    
    score_ptr_pair_type score;
    double loss = 0.0;
    double margin = 0.0;
    
    for (size_type seg = 0; seg != kbests_hyp.size(); ++ seg)
      if (oracles_hyp[seg] && kbests_hyp[seg]) {
	
	if (! direct_loss && oracles_hyp[seg]->score) {
	  if (! score.first)
	    score.first = oracles_hyp[seg]->score->clone();
	  else
	    *score.first += *(oracles_hyp[seg]->score);
	}
	
	if (kbests_hyp[seg]->score) {
	  if (! score.second)
	    score.second = kbests_hyp[seg]->score->clone();
	  else
	    *score.second += *(kbests_hyp[seg]->score);
	}
	
	//margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - (kbests_margin[seg] - kbests_hyp[seg]->loss);
	margin += (oracles_margin[seg] + oracles_hyp[seg]->loss) - kbests_margin[seg];
	//margin += oracles_margin[seg] - kbests_margin[seg];
      }
    
    return std::make_pair((loss - margin) * factor, score);
#endif
  }

  template <typename Iterator>
  std::pair<double, double> operator()(const weight_set_type& weights, const weight_set_type& weights_prev, Iterator iter) const
  {
    
    
    return std::make_pair(0.0, 0.0);
  }
 
  double instances()
  {
#if 0
    size_type samples = 0;
    for (size_t id = 0; id != kbests.size(); ++ id) 
      samples += (! kbests[id].empty()) && (! oracles[id].empty());
    return samples;
#endif

#if 1
    kbests_hyp.clear();
    
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	const size_type seg = kbest_map[id];
	
	if (seg >= kbests_hyp.size())
	  kbests_hyp.resize(seg + 1, 0);
	
	kbests_hyp[seg] = &kbests[id].front();
      }
    
    size_type samples = 0;
    for (size_type seg = 0; seg != kbests_hyp.size(); ++ seg)    
      samples += (kbests_hyp[seg] != 0);
    return samples;
#endif
  }
  
  template <typename Features>
  struct HMatrix
  {
    HMatrix(const Features& __features) : features(__features) {}

    double operator()(int i, int j) const
    {
      return cicada::dot_product(features[i], features[j]);
    }
    
    const Features& features;
  };
  
  template <typename Features>
  struct MMatrix
  {
    MMatrix(const Features& __features) : features(__features) {}
    
    template <typename W>
    void operator()(W& w, const std::vector<double, std::allocator<double> >& alpha) const
    {
      std::vector<double, std::allocator<double> >::const_iterator aiter = alpha.begin();
      
      for (size_type id = 0; id != features.size(); ++ id, ++ aiter)
	if (*aiter > 0.0)
	  operator()(w, *aiter, id);
    }
    
    template <typename W>
    double operator()(const W& w, const size_t& i) const
    {
      return cicada::dot_product(w, features[i]);
    }
    
    template <typename W>
    void operator()(W& w, const double& update, const size_t& i) const
    {
      for (size_t j = 0; j != features[i].size(); ++ j)
	w[j] += update * features[i][j];
    }

    const Features& features;
  };
  
  const scorer_document_type& scorers;
  const hypothesis_map_type&  kbests;
  const hypothesis_map_type&  oracles;
  const kbest_map_type&       kbest_map;

  sample_set_type features_kbest;
  sample_set_type features_oracle;
  
  double samples;
};

template <typename Optimize>
double optimize_cp(const scorer_document_type& scorers,
		   const hypothesis_map_type& kbests,
		   const hypothesis_map_type& oracles,
		   const kbest_map_type& kbest_map,
		   const weight_set_type& bounds_lower,
		   const weight_set_type& bounds_upper,
		   const weights_history_type& weights_history,
		   weight_set_type& weights)
{
  typedef std::deque<weight_set_type, std::allocator<weight_set_type> > weight_queue_type;
  typedef std::vector<double, std::allocator<double> > f_set_type;
  typedef std::vector<double, std::allocator<double> > alpha_set_type;

  typedef std::pair<double, double> point_type;
  typedef std::vector<point_type, std::allocator<point_type> > point_set_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const double loss_factor = (scorers.error_metric() ? 1.0 : - 1.0);
  const double loss_oracle = (scorers.error_metric() ? 0.0 : - 1.0);
  
  Optimize opt(scorers, kbests, oracles, kbest_map);
  
  // setup norm...
  const double instances_local = opt.instances();
  opt.samples = 0.0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &opt.samples, 1, MPI::DOUBLE, MPI::SUM);
  
  weight_queue_type a;
  f_set_type        f;
  alpha_set_type    alpha;

  {
    int history_size = weights_history.size();
    MPI::COMM_WORLD.Bcast(&history_size, 1, MPI::INT, 0);
    
    if (history_size) {
      weight_set_type weights;
      
      for (int i = 0; i != history_size; ++ i) {
	if (mpi_rank == 0)
	  a.push_back(weight_set_type());
	else
	  a.resize(1);
	a.back().clear();

	if (mpi_rank == 0)
	  bcast_weights(0, const_cast<weight_set_type&>(weights_history[i]));
	else
	  bcast_weights(0, weights);
	
	std::pair<double, score_ptr_pair_type> risk_local = opt(mpi_rank == 0 ? weights_history[i] : weights, a.back());
	if (mpi_rank == 0)
	  reduce_weights(a.back());
	else
	  send_weights(a.back());
	
	// reduce objective part
	double risk = 0.0;
	MPI::COMM_WORLD.Reduce(&risk_local.first, &risk, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	// reduce score part
	reduce_score(risk_local.second);
	
	risk -= (risk_local.second.first ? risk_local.second.first->score() * loss_factor : loss_oracle);
	risk += (risk_local.second.second ? risk_local.second.second->score() * loss_factor : 0.0);
	
	if (mpi_rank == 0) {
	  // b = risk + a \cdot w
	  f.push_back(- risk - cicada::dot_product(a.back(), weights_history[i]));
	  alpha.push_back(0.0);
	}
      }
    }
  }
  
  weight_set_type weights_min;
  double objective_master_min = std::numeric_limits<double>::infinity();
  
  double objective_master_prev = std::numeric_limits<double>::infinity();
  double objective_master = 0.0;
  double objective_reduced = 0.0;
    
  // synchronize weights...
  bcast_weights(0, weights);
  
  
  // keep the previous bests
  weight_set_type weights_prev = weights;
  
  // line-search points
  point_set_type points;
  
  // at most 1000 iterations
  // at least two updates
  
  int num_minimum = 0;
  for (int iter = 0; (num_minimum <= 1 || iter < iteration) && (iter < 1000); ++ iter) {
    if (mpi_rank == 0)
      a.push_back(weight_set_type());
    else
      a.resize(1);
    a.back().clear();

    std::pair<double, score_ptr_pair_type> risk_local = opt(weights, a.back());
    if (mpi_rank == 0)
      reduce_weights(a.back());
    else
      send_weights(a.back());
    
    // reduce objective part
    double risk = 0.0;
    MPI::COMM_WORLD.Reduce(&risk_local.first, &risk, 1, MPI::DOUBLE, MPI::SUM, 0);
    
    // reduce score part
    reduce_score(risk_local.second);
    
    risk -= (risk_local.second.first ? risk_local.second.first->score() * loss_factor : loss_oracle);
    risk += (risk_local.second.second ? risk_local.second.second->score() * loss_factor : 0.0);
    
    size_t active_size = 0;
    
    if (mpi_rank == 0) {
    
      // b = risk + a \cdot w
      f.push_back(- risk - cicada::dot_product(a.back(), weights));
      alpha.push_back(0.0);
      
      // peform maximization...
      
      cicada::optimize::QPDCD solver;
      
      typename Optimize::template HMatrix<weight_queue_type> H(a);
      typename Optimize::template MMatrix<weight_queue_type> M(a);
      
      objective_reduced = solver(alpha, f, H, M, 1.0 / C, 1e-5);
      objective_reduced *= - C;
      
      weights.clear();
      active_size = 0;
      alpha_set_type::iterator aiter = alpha.begin();
      for (size_t id = 0; id != a.size(); ++ id, ++ aiter)
	if (*aiter > 0.0) {
	  for (size_t j = 0; j != a[id].size(); ++ j)
	    weights[j] += (*aiter) * a[id][j];
	  
	  ++ active_size;
	}
      
      if (debug >= 3)
	std::cerr << "active size: " << active_size << std::endl;
    }
    
    if (line_search) {
      
      if (debug >= 3 && mpi_rank == 0)
	std::cerr << "line search"  << std::endl;

      bcast_weights(0, weights);
      
      points.clear();
      const std::pair<double, double> grads = opt(weights, weights_prev, std::back_inserter(points));
      
      std::sort(points.begin(), points.end());
      
      reduce_points(points);

      double grad_pos = 0.0;
      double grad_neg = 0.0;
      
      MPI::COMM_WORLD.Reduce(&grads.first,  &grad_pos, 1, MPI::DOUBLE, MPI::SUM, 0);
      MPI::COMM_WORLD.Reduce(&grads.second, &grad_neg, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      if (mpi_rank == 0) {
	
	if (debug >= 3)
	  std::cerr << "point size: " << points.size() << std::endl;
	
	const double norm_w      = cicada::dot_product(weights, weights);
	const double dot_prod    = cicada::dot_product(weights_prev, weights);
	const double norm_w_prev = cicada::dot_product(weights_prev, weights_prev);
	
	const double a0_pos = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * opt.samples;
	const double b0_pos = (dot_prod - norm_w_prev) * C * opt.samples;
	
	const double a0_neg = (norm_w + 2.0 * dot_prod + norm_w_prev) * C * opt.samples;
	const double b0_neg = (- dot_prod - norm_w_prev) * C * opt.samples;
	
	grad_pos += b0_pos;
	grad_neg += b0_neg;

	if (debug >= 3)
	  std::cerr << "gradient: " << grad_pos << ' ' << grad_neg << std::endl;
	
	if (grad_pos < 0.0) {
	  double k = 0.0;
	  
	  point_set_type::const_iterator piter = std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0));
	  point_set_type::const_iterator piter_end = points.end();
	  
	  for (/**/; piter != piter_end && grad_pos < 0.0; ++ piter) {
	    const double k_new = piter->first;
	    const double grad_new = grad_pos + std::fabs(piter->second) + a0_pos * (k_new - k);
	    
	    if (grad_new >= 0) {
	      // compute intersection...
	      k = k + grad_pos * (k - k_new) / (grad_new - grad_pos);
	      grad_pos = grad_new;
	      break;
	    } else {
	      k = k_new;
	      grad_pos = grad_new;
	    }
	  }
	  
	  if (debug >= 3)
	    std::cerr << "grad: " << grad_pos << "  k: " << k << std::endl;

	  if (k > 0.0) {
	    const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
	    
	    for (size_t i = 0; i != weights_size; ++ i)
	      weights[i] = k * weights[i] + (1.0 - k) * weights_prev[i];
	    for (size_t i = weights_size; i < weights.size(); ++ i)
	      weights[i] = k * weights[i];
	    for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	      weights[i] = (1.0 - k) * weights_prev[i];
	  }
	} else if (grad_neg < 0.0) {
	  double k = 0.0;
	  
	  point_set_type::const_reverse_iterator piter(std::lower_bound(points.begin(), points.end(), std::make_pair(0.0, 0.0)));
	  point_set_type::const_reverse_iterator piter_end = points.rend();
	  
	  for (/**/; piter != piter_end && grad_neg < 0.0; ++ piter) {
	    const double k_new = - piter->first;
	    const double grad_new = grad_neg + std::fabs(piter->second) + a0_neg * (k_new - k);
	    
	    if (grad_new >= 0) {
	      // compute intersection...
	      k = k + grad_neg * (k - k_new) / (grad_new - grad_neg);
	      grad_neg = grad_new;
	      break;
	    } else {
	      k = k_new;
	      grad_neg = grad_new;
	    }
	  }
	  
	  if (debug >= 3)
	    std::cerr << "grad: " << grad_neg << "  k: " << - k << std::endl;
	  
	  if (k > 0.0) {
	    const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
	    
	    for (size_t i = 0; i != weights_size; ++ i)
	      weights[i] = - k * weights[i] + (1.0 + k) * weights_prev[i];
	    for (size_t i = weights_size; i < weights.size(); ++ i)
	      weights[i] = - k * weights[i];
	    for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	      weights[i] = (1.0 + k) * weights_prev[i];
	  }
	}
      }
      
      // finished line-search
    }
    
    if (mert_search_local) {
      optimize_mert(scorers, kbests, kbest_map, 0.01, 2.0, weights_prev, weights);

#if 0
      typedef cicada::optimize::LineSearch line_search_type;
      
      typedef line_search_type::segment_type          segment_type;
      typedef line_search_type::segment_set_type      segment_set_type;
      typedef line_search_type::segment_document_type segment_document_type;
      
      typedef line_search_type::value_type optimum_type;

      bcast_weights(0, weights);
      
      const weight_set_type& origin = weights_prev;
      weight_set_type direction = weights;
      direction -= weights_prev;
      
      if (mpi_rank == 0) {
        typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
        
        EnvelopeKBest::line_set_type lines;
        EnvelopeKBest envelopes(origin, direction);
        
        segment_document_type segments;
        
        for (size_t id = 0; id != kbests.size(); ++ id) 
          if (! kbests[id].empty()) {
            segments.push_back(segment_set_type());
            
            //envelopes(kbests[id], oracles[id], lines);
            envelopes(kbests[id], lines);
            
            EnvelopeKBest::line_set_type::const_iterator liter_end = lines.end();
            for (EnvelopeKBest::line_set_type::const_iterator liter = lines.begin(); liter != liter_end; ++ liter)
              segments.back().push_back(std::make_pair(liter->x, liter->hypothesis->score));
          }
	
	for (int rank = 1; rank != mpi_size; ++ rank) {
          boost::iostreams::filtering_istream is;
          is.push(boost::iostreams::zlib_decompressor());
          is.push(utils::mpi_device_source(rank, envelope_tag, 4096));
	  
          std::string line;
          int id_prev = -1;
	  
          while (std::getline(is, line)) {
            const utils::piece line_piece(line);
            tokenizer_type tokenizer(line_piece);
        
            tokenizer_type::iterator iter = tokenizer.begin();
            if (iter == tokenizer.end()) continue;
            const utils::piece id_str = *iter; 
        
            ++ iter;
            if (iter == tokenizer.end() || *iter != "|||") continue;
        
            ++ iter;
            if (iter == tokenizer.end()) continue;
            const utils::piece x_str = *iter;
        
            ++ iter;
            if (iter == tokenizer.end() || *iter != "|||") continue;
        
            ++ iter;
            if (iter == tokenizer.end()) continue;
            const utils::piece score_str = *iter;
        
            const int id = utils::lexical_cast<int>(id_str);
            if (id_prev != id)
              segments.push_back(segment_set_type());
            
            segments.back().push_back(std::make_pair(utils::decode_base64<double>(x_str),
                                                     scorer_type::score_type::decode(score_str)));
            
            id_prev = id;
          }
	}
	
	line_search_type line_search;
        
        const optimum_type optimum = line_search(segments, 0.01, 2.0, scorers.error_metric());
        
        const double update = (optimum.lower + optimum.upper) * 0.5;
        
        if (update != 0.0) {
          direction *= update;
          weights = origin;
          weights += direction;
          
          if (debug >= 2)
            std::cerr << "mert update: " << update
                      << " objective: " << optimum.objective << std::endl;
        }
	
      } else {
        EnvelopeKBest::line_set_type lines;
        EnvelopeKBest envelopes(origin, direction);
        
        boost::iostreams::filtering_ostream os;
        os.push(boost::iostreams::zlib_compressor());
        os.push(utils::mpi_device_sink(0, envelope_tag, 4096));
        
        for (size_t id = 0; id != kbests.size(); ++ id) 
          if (! kbests[id].empty()) {
            //envelopes(kbests[id], oracles[id], lines);
            envelopes(kbests[id], lines);
            
            EnvelopeKBest::line_set_type::const_iterator liter_end = lines.end();
            for (EnvelopeKBest::line_set_type::const_iterator liter = lines.begin(); liter != liter_end; ++ liter) {
              const EnvelopeKBest::line_type& line = *liter;
              
              os << id << " ||| ";
              utils::encode_base64(line.x, std::ostream_iterator<char>(os));
              os << " ||| " << line.hypothesis->score->encode() << '\n';
            }
          }
      }

#endif
      
#if 0
      typedef cicada::optimize::LineSearch line_search_type;
      
      typedef line_search_type::segment_type          segment_type;
      typedef line_search_type::segment_set_type      segment_set_type;
      typedef line_search_type::segment_document_type segment_document_type;
      
      typedef line_search_type::value_type optimum_type;

      typedef EnvelopeKBest::line_set_type line_set_type;
      typedef std::deque<line_set_type, std::allocator<line_set_type> > line_map_type;
      
      bcast_weights(0, weights);
      
      const weight_set_type& origin = weights_prev;
      weight_set_type direction = weights;
      direction -= weights_prev;
      
      EnvelopeKBest envelopes(origin, direction);
      line_map_type lines;
      
      for (size_t id = 0; id != kbests.size(); ++ id)
	if (! kbests[id].empty()) {
	  
	  if (kbest_map[id] >= lines.size())
	    lines.resize(kbest_map[id] + 1);
	  
	  envelopes(kbests[id].begin(), kbests[id].end(), std::back_inserter(lines[kbest_map[id]]));
	}
      
      if (mpi_rank == 0) {
        typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	
        segment_document_type segments;
	
	for (size_t seg = 0; seg != lines.size(); ++ seg)
	  if (! lines[seg].empty()) {
	    
	    if (seg >= segments.size())
	      segments.resize(seg + 1);
	    
	    envelopes(lines[seg]);
	    
	    EnvelopeKBest::line_set_type::const_iterator liter_end = lines[seg].end();
	    for (EnvelopeKBest::line_set_type::const_iterator liter = lines[seg].begin(); liter != liter_end; ++ liter)
	      segments[seg].push_back(std::make_pair(liter->x, liter->hypothesis->score));
	  }
	
        for (int rank = 1; rank != mpi_size; ++ rank) {
          boost::iostreams::filtering_istream is;
          is.push(boost::iostreams::zlib_decompressor());
          is.push(utils::mpi_device_source(rank, envelope_tag, 4096));

          std::string line;
          while (std::getline(is, line)) {
            const utils::piece line_piece(line);
            tokenizer_type tokenizer(line_piece);
        
            tokenizer_type::iterator iter = tokenizer.begin();
            if (iter == tokenizer.end()) continue;
            const utils::piece id_str = *iter; 
        
            ++ iter;
            if (iter == tokenizer.end() || *iter != "|||") continue;
        
            ++ iter;
            if (iter == tokenizer.end()) continue;
            const utils::piece x_str = *iter;
        
            ++ iter;
            if (iter == tokenizer.end() || *iter != "|||") continue;
        
            ++ iter;
            if (iter == tokenizer.end()) continue;
            const utils::piece score_str = *iter;
	    
            const int id = utils::lexical_cast<int>(id_str);
	    
	    if (id >= segments.size())
	      segments.resize(id + 1);
	    
            segments[id].push_back(std::make_pair(utils::decode_base64<double>(x_str),
						  scorer_type::score_type::decode(score_str)));
          }
        }
	
	line_search_type line_search;
        
        const optimum_type optimum = line_search(segments, line_search_type::RegularizeL2(C, origin, direction), 0.01, 2.0, scorers.error_metric());
	
        const double update = (optimum.lower + optimum.upper) * 0.5;
        
        if (update != 0.0) {
          direction *= update;
          weights = origin;
          weights += direction;
	  
          if (debug >= 2)
            std::cerr << "mert update: " << update
                      << " objective: " << optimum.objective << std::endl;
        }
	
      } else {
        boost::iostreams::filtering_ostream os;
        os.push(boost::iostreams::zlib_compressor());
        os.push(utils::mpi_device_sink(0, envelope_tag, 4096));
	
	for (size_t seg = 0; seg != lines.size(); ++ seg)
	  if (! lines[seg].empty()) {
	    
	    envelopes(lines[seg]);
            
            EnvelopeKBest::line_set_type::const_iterator liter_end = lines[seg].end();
            for (EnvelopeKBest::line_set_type::const_iterator liter = lines[seg].begin(); liter != liter_end; ++ liter) {
              const EnvelopeKBest::line_type& line = *liter;
              
              os << seg << " ||| ";
              utils::encode_base64(line.x, std::ostream_iterator<char>(os));
              os << " ||| " << line.hypothesis->score->encode() << '\n';
            }
          }
      }
#endif
      // finished mert-search-local
    }
    
    // current weights is the master problems weights...
    bcast_weights(0, weights);
    
    objective_master = 0.0;
    
    std::pair<double, score_ptr_pair_type> objective_master_local = opt.objective(weights);
    
    // reduce objectice part
    MPI::COMM_WORLD.Reduce(&objective_master_local.first, &objective_master, 1, MPI::DOUBLE, MPI::SUM, 0);
    
    // reduce score part
    reduce_score(objective_master_local.second);
    
    objective_master -= (objective_master_local.second.first ? objective_master_local.second.first->score() * loss_factor : loss_oracle);
    objective_master += (objective_master_local.second.second ? objective_master_local.second.second->score() * loss_factor : 0.0);
    
    objective_master += 0.5 * C * cicada::dot_product(weights, weights);

    if (objective_master <= objective_master_min) {
      weights_min = weights;
      objective_master_min = objective_master;
      
      ++ num_minimum;
    }
    MPI::COMM_WORLD.Bcast(&num_minimum, 1, MPI::INT, 0);

    if (mpi_rank == 0 && debug >= 2)
      std::cerr << "objective master: " << objective_master
		<< " reduced: " << objective_reduced
		<< " actives: " << active_size << std::endl;
    
    // check termination condition...
    int terminate = (std::fabs((objective_master - objective_reduced) / objective_master) < 0.01);
    MPI::COMM_WORLD.Bcast(&terminate, 1, MPI::INT, 0);
    
    if (terminate) break;

    MPI::COMM_WORLD.Bcast(&objective_master, 1, MPI::DOUBLE, 0);
    
    // we will update proximy when better solution found

    if (iter && (line_search || mert_search_local) && objective_master - objective_master_prev > 0.001) {
#if 1
      // we will try find the best scaling between weights_prev and weights
      // we do not update proximy!
      const double suffered_loss = objective_master - objective_master_prev;
      
      double variance = 0.0;
      
      const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
      
      for (size_t i = 0; i != weights_size; ++ i)
	variance += (weights[i] - weights_prev[i]) * (weights[i] - weights_prev[i]);
      for (size_t i = weights_size; i < weights.size(); ++ i)
	variance += weights[i] * weights[i];
      for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	variance += weights_prev[i] * weights_prev[i];
      
      //const double k = std::min(suffered_loss / variance, 1.0);
      const double k = suffered_loss / variance;
      
      if (k != 1.0) {
	for (size_t i = 0; i != weights_size; ++ i)
	  weights[i] = k * weights[i] + (1.0 - k) * weights_prev[i];
	for (size_t i = weights_size; i < weights.size(); ++ i)
	  weights[i] = k * weights[i];
	for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	  weights[i] = (1.0 - k) * weights_prev[i];
      }

      if (mpi_rank == 0 && debug >= 2) 
	std::cerr << "cutting plane ratio: " << k << std::endl;
#endif
#if 0 
      const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
      
      const double k = 0.1;
      
      for (size_t i = 0; i != weights_size; ++ i)
	weights[i] = k * weights[i] + (1.0 - k) * weights_prev[i];
      for (size_t i = weights_size; i < weights.size(); ++ i)
	weights[i] = k * weights[i];
      for (size_t i = weights_size; i < weights_prev.size(); ++ i)
	weights[i] = (1.0 - k) * weights_prev[i];
#endif
      
    } else {
      weights_prev = weights;
      objective_master_prev = objective_master;
    }
  }

  weights = weights_min;
  
  return objective_master_min;
}


struct OptimizeXBLEU
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef hypothesis_type::feature_value_type feature_value_type;
  typedef utils::mulvector2<feature_value_type, std::allocator<feature_value_type> > sample_set_type;
  
  OptimizeXBLEU(const hypothesis_map_type& __kbests,
		const sample_set_type& __features_kbest,
		const scorer_document_type& __scorers,
		weight_set_type& __weights,
		const size_t& __instances,
		const double& __lambda,
		const feature_type& __feature_scale)
    : kbests(__kbests),
      features_kbest(__features_kbest),
      scorers(__scorers),
      weights(__weights),
      instances(__instances),
      lambda(__lambda),
      feature_scale(__feature_scale)
  { }
  
  const hypothesis_map_type& kbests;
  const sample_set_type&     features_kbest;
  const scorer_document_type& scorers;
  weight_set_type& weights;
  
  size_t instances;
  double lambda;
  const feature_type feature_scale;

  double objective_opt;
  weight_set_type weights_opt;

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
    
    objective_opt = std::numeric_limits<double>::infinity();
    double objective = 0.0;
    
    const int result = lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeXBLEU::evaluate, 0, this, &param);
    
    if (debug)
      std::cerr << "lbfgs: " << lbfgs_error(result) << std::endl;
    
    // copy from opt weights!
    if (result < 0)
      weights = weights_opt;
    
    return objective;
  }

  struct Task
  {
    typedef cicada::semiring::Log<double> weight_type;
    
    weight_type brevity_penalty(const double x) const
    {
      typedef cicada::semiring::traits<weight_type> traits_type;

      // return (std::exp(x) - 1) / (1.0 + std::exp(1000.0 * x)) + 1.0;
      
      return ((traits_type::exp(x) - traits_type::one()) / (traits_type::one() + traits_type::exp(1000.0 * x))) + traits_type::one();
    }
    
    weight_type derivative_brevity_penalty(const double x) const
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
       
      const weight_type expx     = traits_type::exp(x);
      const weight_type expxm1   = expx - traits_type::one();
      const weight_type exp1000x = traits_type::exp(1000.0 * x);
      const weight_type p1exp1000x = traits_type::one() + exp1000x;
      
      return (expx / p1exp1000x) - ((expxm1 * weight_type(1000.0) * exp1000x) / (p1exp1000x * p1exp1000x));
      
      //return expx / (1.0 + exp1000x) - boost::math::expm1(x) * (1000.0 * exp1000x) / ((1.0 + exp1000x) * (1.0 + exp1000x))
    }
    
    weight_type clip_count(const weight_type& x, const weight_type& clip) const
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
      
      //return (x - clip) / (1.0 + std::exp(1000.0 * (x - clip))) + clip;
      return (weight_type(x - clip) / (traits_type::one() + traits_type::exp(1000.0 * (x - clip)))) + weight_type(clip);
    }
    
    weight_type derivative_clip_count(const weight_type& x, const weight_type& clip) const
    {
      typedef cicada::semiring::traits<weight_type> traits_type;
      
      const weight_type exp1000xmc = traits_type::exp(1000.0 * (x - clip));
      const weight_type p1exp1000xmc = exp1000xmc + traits_type::one();
      
      return (traits_type::one() / p1exp1000xmc) - ((weight_type(x - clip) * weight_type(1000.0) * exp1000xmc) / (p1exp1000xmc * p1exp1000xmc));
      
      //return 1.0 / (1.0 + exp1000x) - (x - clip) * (1000.0 * exp1000x) / ((1.0 + exp1000x) * (1.0 + exp1000x));
    }

    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef std::vector<gradient_type, std::allocator<gradient_type> > gradients_type;
    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
    
    typedef std::vector<double, std::allocator<double> > ngram_counts_type;
    typedef std::vector<weight_set_type, std::allocator<weight_set_type> > feature_counts_type;
    
    typedef std::vector<double, std::allocator<double> > margins_type;

    typedef cicada::FeatureVectorUnordered<weight_type, std::allocator<weight_type> > expectation_type;
    
    Task(const hypothesis_map_type& __kbests,
	 const sample_set_type& __features_kbest,
	 const scorer_document_type& __scorers,
	 const weight_set_type& __weights,
	 const feature_type& __feature_scale)
      : kbests(__kbests),
	features_kbest(__features_kbest),
	scorers(__scorers),
	weights(__weights),
	feature_scale(__feature_scale),
	c_matched(order + 1),
	c_hypo(order + 1),
	g_matched(order + 1),
	g_hypo(order + 1),
	g_reference(),
	g_entropy(),
	r(0),
	e(0)
    { }

    const hypothesis_map_type& kbests;
    const sample_set_type& features_kbest;
    const scorer_document_type& scorers;
    const weight_set_type& weights;
    const feature_type feature_scale;
    
    ngram_counts_type   c_matched;
    ngram_counts_type   c_hypo;
    feature_counts_type g_matched;
    feature_counts_type g_hypo;
    weight_set_type     g_reference;
    weight_set_type     g_entropy;
    double r;
    double e;
    
    void operator()()
    {
      weights_type    matched(order + 1);
      weights_type    hypo(order + 1);
      expectation_type expectation;
      
      weights_type    counts_matched(order + 1);
      weights_type    counts_hypo(order + 1);
      
      gradients_type gradients_matched(order + 1);
      gradients_type gradients_hypo(order + 1);
      
      weight_type    reference;
      weight_type    entropy;
      gradient_type  gradient_reference;
      gradient_type  gradient_entropy;
      
      margins_type  margins;

      for (size_t n = 1; n != g_matched.size(); ++ n) {
	gradients_matched[n].allocate();
	gradients_hypo[n].allocate();
	
	g_matched[n].clear();
	g_hypo[n].clear();
      }
      
      gradient_reference.allocate();
      gradient_entropy.allocate();
      g_reference.clear();
      g_entropy.clear();
      
      std::fill(counts_matched.begin(), counts_matched.end(), weight_type());
      std::fill(counts_hypo.begin(), counts_hypo.end(), weight_type());
      std::fill(c_matched.begin(), c_matched.end(), 0.0);
      std::fill(c_hypo.begin(), c_hypo.end(), 0.0);
      r = 0.0;
      e = 0.0;
      
      const double scale = weights[feature_scale];
      
      size_type features_pos = 0;
      for (size_type id = 0; id != kbests.size(); ++ id) 
	if (! kbests[id].empty()) {
	  margins.clear();
	  
	  std::fill(matched.begin(), matched.end(), weight_type());
	  std::fill(hypo.begin(), hypo.end(), weight_type());
	  expectation.clear();
	  
	  weight_type Z;
	  weight_type Z_reference;
	  weight_type Z_entropy;
	  weight_type dR;
	  
	  const size_type features_start = features_pos;
	  
	  // first pass... compute margin and Z
	  for (size_type k = 0; k != kbests[id].size(); ++ k, ++ features_pos) {
	    const hypothesis_type& kbest = kbests[id][k];
	    
	    const double margin  = cicada::dot_product(weights, features_kbest[features_pos].begin(), features_kbest[features_pos].end(), 0.0);
	    
	    margins.push_back(margin);
	    Z += cicada::semiring::traits<weight_type>::exp(margin * scale);
	  }

	  // second pass... compute sums (counts, etc.)
	  for (size_type k = 0, i = features_start; k != kbests[id].size(); ++ k, ++ i) {
	    const hypothesis_type& kbest = kbests[id][k];
	    
	    const double& margin = margins[k];
	    const weight_type prob = cicada::semiring::traits<weight_type>::exp(margin * scale) / Z;
	    
	    const cicada::eval::Bleu* bleu = dynamic_cast<const cicada::eval::Bleu*>(kbest.score.get());
	    if (! bleu)
	      throw std::runtime_error("no bleu statistics?");
	    
	    // collect scaled bleu stats
	    for (int n = 1; n <= order; ++ n) {
	      if (n - 1 < bleu->ngrams_reference.size())
		hypo[n] += prob * bleu->ngrams_reference[n - 1];
	      if (n - 1 < bleu->ngrams_hypothesis.size())
		matched[n] += prob * bleu->ngrams_hypothesis[n - 1];
	    }
	    
	    // collect reference length
	    Z_reference += prob * bleu->length_reference;
	    
	    // collect entropy...
	    Z_entropy -= prob * cicada::semiring::log(prob);
	    
	    // collect expectation
	    sample_set_type::const_reference::const_iterator fiter_end = features_kbest[i].end();
	    for (sample_set_type::const_reference::const_iterator fiter = features_kbest[i].begin(); fiter != fiter_end; ++ fiter)
	      expectation[fiter->first] += prob * weight_type(fiter->second * scale);
	    
	    expectation[feature_scale] += prob * weight_type(margin);
	    
	    dR += weight_type(1.0 + cicada::semiring::log(prob)) * prob;
	  }
	  
	  // accumulate
	  std::transform(hypo.begin(), hypo.end(), counts_hypo.begin(), counts_hypo.begin(), std::plus<weight_type>());
	  std::transform(matched.begin(), matched.end(), counts_matched.begin(), counts_matched.begin(), std::plus<weight_type>());
	  
	  reference += Z_reference;
	  entropy += Z_entropy;
	  
	  expectation_type::const_iterator eiter_end = expectation.end();
	  for (expectation_type::const_iterator eiter = expectation.begin(); eiter != eiter_end; ++ eiter) {
	    // collect bleus...
	    for (int n = 1; n <= order; ++ n) {
	      gradients_hypo[n][eiter->first] -= eiter->second * hypo[n];
	      gradients_matched[n][eiter->first] -= eiter->second * matched[n];
	    }
	    
	    // reference lengths
	    gradient_reference[eiter->first] -= eiter->second * Z_reference;
	    
	    // entropy gradient...
	    gradient_entropy[eiter->first] -= - dR * expectation[eiter->first];
	  }
	  
	  // third pass, collect gradients...
	  for (size_type k = 0, i = features_start; k != kbests[id].size(); ++ k, ++ i) {
	    const hypothesis_type& kbest = kbests[id][k];
	    
	    const double& margin = margins[k];
	    const weight_type prob = cicada::semiring::traits<weight_type>::exp(margin * scale) / Z;
	    const cicada::eval::Bleu* bleu = dynamic_cast<const cicada::eval::Bleu*>(kbest.score.get());
	    
	    // collect feature expectations etc...
	    sample_set_type::const_reference::const_iterator fiter_end = features_kbest[i].end();
	    for (sample_set_type::const_reference::const_iterator fiter = features_kbest[i].begin(); fiter != fiter_end; ++ fiter) {
	      const weight_type value(fiter->second * scale);
	      
	      // bleu statistics
	      for (int n = 1; n <= order; ++ n) {
		if (n - 1 < bleu->ngrams_reference.size())
		  gradients_hypo[n][fiter->first] += value * prob * bleu->ngrams_reference[n - 1];
		
		if (n - 1 < bleu->ngrams_hypothesis.size())
		  gradients_matched[n][fiter->first] += value * prob * bleu->ngrams_hypothesis[n - 1];
	      }
	      
	      // reference lengths
	      gradient_reference[fiter->first] += value * prob * bleu->length_reference;
	      
	      // entropy: we will collect minus values!
	      gradient_entropy[fiter->first] += - weight_type(1.0 + cicada::semiring::log(prob)) * prob * value;
	    }
	    
	    const weight_type value_scale(margin);
	    
	    for (int n = 1; n <= order; ++ n) {
	      if (n - 1 < bleu->ngrams_reference.size())
		gradients_hypo[n][feature_scale] += value_scale * prob * bleu->ngrams_reference[n - 1];
	      
	      if (n - 1 < bleu->ngrams_hypothesis.size())
		gradients_matched[n][feature_scale] += value_scale * prob * bleu->ngrams_hypothesis[n - 1];
	    }
	    
	    gradient_reference[feature_scale] += value_scale * prob * bleu->length_reference;
	    gradient_entropy[feature_scale] += - weight_type(1.0 + cicada::semiring::log(prob)) * prob * value_scale;
	  }
	}
      
      // copy from weight space to double space
      std::copy(counts_matched.begin(), counts_matched.end(), c_matched.begin());
      std::copy(counts_hypo.begin(), counts_hypo.end(), c_hypo.begin());
      
      for (size_t n = 1; n != g_matched.size(); ++ n) {
	g_matched[n].allocate();
	g_hypo[n].allocate();
	
	std::copy(gradients_matched[n].begin(), gradients_matched[n].end(), g_matched[n].begin());
	std::copy(gradients_hypo[n].begin(), gradients_hypo[n].end(), g_hypo[n].begin());
      }
      
      g_reference.allocate();
      g_entropy.allocate();
      std::copy(gradient_reference.begin(), gradient_reference.end(), g_reference.begin());
      std::copy(gradient_entropy.begin(), gradient_entropy.end(), g_entropy.begin());

      r = reference;
      e = entropy;
    }
  };
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int size,
				  const lbfgsfloatval_t step)
  {
    typedef Task task_type;
    
    const int mpi_rank = MPI::COMM_WORLD.Get_rank();
    const int mpi_size = MPI::COMM_WORLD.Get_size();
    
    OptimizeXBLEU& optimizer = *((OptimizeXBLEU*) instance);
    
    // send notification!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
    
    if (debug >= 3)
      std::cerr << "weights:" << std::endl
		<< optimizer.weights << std::flush;
    
    bcast_weights(0, optimizer.weights);
    
    task_type task(optimizer.kbests, optimizer.features_kbest, optimizer.scorers, optimizer.weights, optimizer.feature_scale);
    task();
    
    {
      task_type::ngram_counts_type c_matched(order + 1, 0.0);
      task_type::ngram_counts_type c_hypo(order + 1, 0.0);
      double                       r(0.0);
      double                       e(0.0);
      
      // reduce c_* and r 
      MPI::COMM_WORLD.Reduce(&(*task.c_matched.begin()), &(*c_matched.begin()), order + 1, MPI::DOUBLE, MPI::SUM, 0);
      MPI::COMM_WORLD.Reduce(&(*task.c_hypo.begin()), &(*c_hypo.begin()), order + 1, MPI::DOUBLE, MPI::SUM, 0);
      MPI::COMM_WORLD.Reduce(&task.r, &r, 1, MPI::DOUBLE, MPI::SUM, 0);
      MPI::COMM_WORLD.Reduce(&task.e, &e, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      task.c_matched.swap(c_matched);
      task.c_hypo.swap(c_hypo);
      task.r = r;
      task.e = e;
    }
    
    
    // reduce g_*
    for (int n = 1; n <= order; ++ n) {
      // reduce matched counts..
      reduce_weights(task.g_matched[n]);
      reduce_weights(task.g_hypo[n]);
    }
    reduce_weights(task.g_reference);
    reduce_weights(task.g_entropy);
    
    // smoothing...
    {
      double smoothing = 1e-40;
      for (int n = 1; n <= order; ++ n) {
	if (task.c_hypo[n] > 0.0 && task.c_matched[n] <= 0.0)
	  task.c_matched[n] = smoothing;
	smoothing *= 0.1;
      }
    }
    
    if (debug >= 3) {
      for (int n = 1; n <= order; ++ n)
	std::cerr << "order: " << n << " M: " << task.c_matched[n] << " H: " << task.c_hypo[n] << std::endl;
      std::cerr << "r: " << task.r << std::endl;
    }

    // compute P
    double P = 0.0;
    for (int n = 1; n <= order; ++ n)
      if (task.c_hypo[n] > 0.0)
	P += (1.0 / order) * (utils::mathop::log(task.c_matched[n]) - utils::mathop::log(task.c_hypo[n]));

    // compute C and B
    const double C = task.r / task.c_hypo[1];
    const double B = task.brevity_penalty(1.0 - C);

    // for computing g...
    const double exp_P = utils::mathop::exp(P);
    const double C_dC  = C * task.derivative_brevity_penalty(1.0 - C);
    
    // xBLEU...
    const double objective_bleu = exp_P * B;
    const double entropy = task.e / optimizer.instances;
    
    // compute g..
    std::fill(g, g + size, 0.0);
    
    // entropy
    std::transform(task.g_entropy.begin(), task.g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), - temperature / optimizer.instances));
    
    for (int n = 1; n <= order; ++ n) 
      if (task.c_hypo[n] > 0.0) {
	const double factor_matched = - (exp_P * B / order) / task.c_matched[n];
	const double factor_hypo    = - (exp_P * B / order) / task.c_hypo[n];
	
	for (size_t i = 0; i != static_cast<size_t>(size); ++ i) {
	  g[i] += factor_matched * task.g_matched[n][i];
	  g[i] -= factor_hypo * task.g_hypo[n][i];
	}
      }
    
    if (task.c_hypo[1] > 0.0) {
      const double factor_ref  = - exp_P * C_dC / task.r;
      const double factor_hypo = - exp_P * C_dC / task.c_hypo[1];
      
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i) {
	g[i] -= factor_ref  * task.g_reference[i];
	g[i] += factor_hypo * task.g_hypo[1][i];
      }
    }
    
    if (debug >= 3) {
      std::cerr << "grad:" << std::endl;
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i)
	if (g[i] != 0.0 && feature_type(i) != feature_type())
	  std::cerr << feature_type(i) << ' ' << g[i] << std::endl;
    }
    
    // we need to minimize negative bleu... + regularized by average entropy...
    double objective = - objective_bleu - temperature * entropy;
    double objective_regularized = objective;
    
    if (regularize_l2) {
      double norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(size); ++ i) {
	g[i] += optimizer.lambda * x[i] * double(i != optimizer.feature_scale.id());
	norm += x[i] * x[i] * double(i != optimizer.feature_scale.id());
      }
      
      objective_regularized += 0.5 * optimizer.lambda * norm;
      objective += 0.5 * optimizer.lambda * norm;
    } else if (regularize_l1) {
      double norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(size); ++ i)
	norm += std::fabs(x[i]);
      
      objective_regularized += optimizer.lambda * norm;
    }
    
    if (scale_fixed)
      g[optimizer.feature_scale.id()] = 0.0;
    
    if (debug >= 2)
      std::cerr << "objective: " << objective_regularized
		<< " xBLEU: " << objective_bleu
		<< " BP: " << B
		<< " entropy: " << entropy
		<< " scale: " << optimizer.weights[optimizer.feature_scale]
		<< std::endl;
    
    
    // keep the best so forth...
    if (std::isfinite(objective_regularized) && objective_regularized <= optimizer.objective_opt) {
      optimizer.objective_opt = objective_regularized;
      optimizer.weights_opt = optimizer.weights;
    }
    
    return objective;
  }
  
};

struct OptimizeLBFGS
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
      typedef const feature_value_type* const_iterator;
      
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
      return sample_type(&(*features.begin()) + offsets[pos], &(*features.begin()) + offsets[pos + 1]);
    }
    
    size_type size() const { return offsets.size() - 1; }
    bool empty() const { return offsets.size() == 1; }
    
    void shrink()
    {
      features_type(features).swap(features);
      offsets_type(offsets).swap(offsets);
    }
    
    features_type features;
    offsets_type  offsets;
  };
  
  typedef SampleSet sample_set_type;
  
  struct sample_pair_type
  {
    typedef std::vector<double, std::allocator<double> > loss_set_type;

    sample_pair_type() : features(), loss(),  offset(0) {}
    sample_pair_type(const hypothesis_set_type& kbests,
		     const hypothesis_set_type& oracles)
      : features(), loss(),  offset(0)
    {
      loss.reserve(kbests.size() + oracles.size());

      hypothesis_set_type::const_iterator oiter_end = oracles.end();
      for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter) {
	features.insert(oiter->features.begin(), oiter->features.end());
	loss.push_back(oiter->loss);
      }
      
      offset = loss.size();
      
      hypothesis_set_type::const_iterator kiter_end = kbests.end();
      for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter) {
	features.insert(kiter->features.begin(), kiter->features.end());
	loss.push_back(kiter->loss);
      }
      
      features.shrink();
    }
    
    size_type oracle_begin() const { return 0; }
    size_type oracle_end() const { return offset; }
    
    size_type kbest_begin() const { return offset; }
    size_type kbest_end() const { return loss.size(); }
    
    size_type size() const { return loss.size(); }
    
    sample_set_type features;
    loss_set_type   loss;
    size_type offset;
  };
  
  typedef std::vector<sample_pair_type, std::allocator<sample_pair_type> > sample_pair_set_type;

  OptimizeLBFGS(const sample_pair_set_type& __samples,
		const weight_set_type& __bounds_lower,
		const weight_set_type& __bounds_upper,
		weight_set_type& __weights,
		const size_t& __instances)
    : samples(__samples),
      bounds_lower(__bounds_lower),
      bounds_upper(__bounds_upper),
      weights(__weights),
      instances(__instances) {}

  double operator()()
  {
    if (regularize_l2 && (! bounds_lower.empty() || ! bounds_upper.empty())) {
      typedef Task                  task_type;
      
      const int mpi_rank = MPI::COMM_WORLD.Get_rank();
      const int mpi_size = MPI::COMM_WORLD.Get_size();
      
      std::vector<double, std::allocator<double> > g(weights.size());
      std::vector<double, std::allocator<double> > l(weights.size());
      std::vector<double, std::allocator<double> > u(weights.size());
      std::vector<int, std::allocator<int> >       nbd(weights.size());

      const double inf = std::numeric_limits<double>::infinity();

      for (size_t i = 0; i != weights.size(); ++ i) {
	bool has_lower_bound = (i < bounds_lower.size() && bounds_lower[i] > - inf);
	bool has_upper_bound = (i < bounds_upper.size() && bounds_upper[i] <   inf);
	
	if (i < bounds_lower.size())
	  l[i] = bounds_lower[i];
	
	if (i < bounds_upper.size())
	  u[i] = bounds_upper[i];
	
	nbd[i] = (has_lower_bound && has_upper_bound ? 2 : (has_lower_bound ? 1 : (has_upper_bound ? 3 : 0)));
      }
      
      LBFGS lbfgs;
      lbfgs.init(weights.size(), 7, &(*l.begin()), &(*u.begin()), &(*nbd.begin()));
      
      double objective = 0.0;
      for (;;) {
	// send notification!
	for (int rank = 1; rank < mpi_size; ++ rank)
	  MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
	
	bcast_weights(0, weights);
	
	task_type task(weights, samples, instances);
	task();
	
	// collect all the objective and gradients...
	reduce_weights(task.g);
	
	std::fill(std::copy(task.g.begin(), task.g.end(), g.begin()), g.end(), 0.0);
	
	objective = 0.0;
	MPI::COMM_WORLD.Reduce(&task.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	const double objective_unregularized = objective;
	
	double norm = 0.0;
	for (size_t i = 0; i < g.size(); ++ i) {
	  g[i] += C * weights[i];
	  norm += weights[i] * weights[i];
	}
	objective += 0.5 * C * norm;
	
	if (debug >= 2)
	  std::cerr << "objective: " << objective << " non-regularized: " << objective_unregularized << std::endl;

	const int ret = lbfgs.optimize(&(*weights.begin()), &objective, &(*g.begin()));

	if (ret <= 0)
	  break;
      }
      
      return objective;
    } else {
      lbfgs_parameter_t param;
      lbfgs_parameter_init(&param);
      
      if (regularize_l1) {
	param.orthantwise_c = C;
	param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
      } else
	param.orthantwise_c = 0.0;
      
      param.max_iterations = iteration;
      
      objective_opt = std::numeric_limits<double>::infinity();
      double objective = 0.0;
      
      const int result = lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
      
      if (debug)
	std::cerr << "lbfgs: " << lbfgs_error(result) << std::endl;
      
      // copy from opt weights!
      if (result < 0)
	weights = weights_opt;
    
      return objective;
    }
  }
  
  struct Task
  {
    typedef hypothesis_type::feature_value_type feature_value_type;
    
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;
    
    Task(const weight_set_type& __weights,
	 const sample_pair_set_type& __samples,
	 const size_t& __instances)
      : weights(__weights),
	samples(__samples),
	instances(__instances)
    {}
    
    void operator()()
    {
      typedef std::vector<double, std::allocator<double> > margin_set_type;
      
      g.clear();
      objective = 0.0;
      
      expectation_type  expectations;
      
      expectations.allocate();
      expectations.clear();
      
      margin_set_type margins;

      const double cost_factor = (softmax_margin ? 1.0 : 0.0);
      
      for (size_t id = 0; id != samples.size(); ++ id) {
	
	weight_type Z_oracle;
	weight_type Z_kbest;
	
	margins.clear();
	margins.resize(samples[id].size());
	
	for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const double loss = samples[id].loss[i];
	  
	  margins[i] = cicada::dot_product(weights, features.begin(), features.end(), cost_factor * loss);
	  
	  Z_oracle += cicada::semiring::traits<weight_type>::exp(margins[i]);
	}
	
	for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const double loss = samples[id].loss[i];
	  
	  margins[i] = cicada::dot_product(weights, features.begin(), features.end(), cost_factor * loss);
	  
	  Z_kbest += cicada::semiring::traits<weight_type>::exp(margins[i]);
	}
	
	for (size_type i = samples[id].oracle_begin(); i != samples[id].oracle_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(margins[i]) / Z_oracle;
	  
	  sample_set_type::value_type::const_iterator fiter_end = features.end();
	  for (sample_set_type::value_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
	    expectations[fiter->first] -= weight_type(fiter->second) * weight;
	}
	
	for (size_type i = samples[id].kbest_begin(); i != samples[id].kbest_end(); ++ i) {
	  const sample_set_type::value_type features = samples[id].features[i];
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(margins[i]) / Z_kbest;
	  
	  sample_set_type::value_type::const_iterator fiter_end = features.end();
	  for (sample_set_type::value_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter)
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
    const sample_pair_set_type& samples;
    
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
    
    task_type task(optimizer.weights, optimizer.samples, optimizer.instances);
    task();
    
    // collect all the objective and gradients...
    
    reduce_weights(task.g);
    
    std::fill(std::copy(task.g.begin(), task.g.end(), g), g + n, 0.0);
    
    double objective = 0.0;
    MPI::COMM_WORLD.Reduce(&task.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
    
    const double objective_unregularized = objective;
    double objective_regularized = objective;
    
    // L2...
    if (regularize_l2) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i) {
	g[i] += C * x[i];
	norm += x[i] * x[i];
      }
      
      objective_regularized += 0.5 * C * norm;
      objective += 0.5 * C * norm;
    } else if (regularize_l1) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i)
	norm += std::fabs(x[i]);
      
      objective_regularized += C * norm;
    }
    
    if (debug >= 2)
      std::cerr << "objective: " << objective_regularized << " non-regularized: " << objective_unregularized << std::endl;
    
    // keep the best so forth...
    if (std::isfinite(objective_regularized) && objective_regularized <= optimizer.objective_opt) {
      optimizer.objective_opt = objective_regularized;
      optimizer.weights_opt = optimizer.weights;
    }
    
    return objective;
  }
    
  const sample_pair_set_type& samples;

  const weight_set_type& bounds_lower;
  const weight_set_type& bounds_upper;
  
  weight_set_type& weights;
  size_t instances;
  
  double objective_opt;
  weight_set_type weights_opt;
};

template <typename Optimize>
double optimize_xbleu(const hypothesis_map_type& kbests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  const feature_type feature_scale(":feature-scale:");
  
  weights[feature_scale] = scale;
  
  typename Optimize::sample_set_type features;

  int instances_local = 0;
  for (size_t id = 0; id != kbests.size(); ++ id) 
    if (! kbests[id].empty()) {
      
      for (size_type k = 0; k != kbests[id].size(); ++ k)
	features.push_back(kbests[id][k].features.begin(), kbests[id][k].features.end());
      
      ++ instances_local;
    }
  
  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  
  if (mpi_rank == 0) {
    Optimize optimizer(kbests, features, scorers, weights, instances, C, feature_scale);
    
    double objective = 0.0;
    
    if (annealing_mode) {
      for (temperature = temperature_start; temperature >= temperature_end; temperature *= temperature_rate) {
	if (debug >= 2)
	  std::cerr << "temperature: " << temperature << std::endl;
	
	objective = optimizer();
      }
    } else 
      objective = optimizer();
    
    if (quenching_mode) {
      temperature = 0.0;
      
      for (double quench = quench_start; quench <= quench_end; quench *= quench_rate) {
	if (debug >= 2)
	  std::cerr << "quench: " << quench << std::endl;
	
	weights[feature_scale] = quench;
	
	objective = optimizer();
      }
    }
    
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
	
	task_type task(kbests, features, scorers, weights, feature_scale);
	task();
	
	typename task_type::ngram_counts_type c_matched(order + 1, 0.0);
	typename task_type::ngram_counts_type c_hypo(order + 1, 0.0);
	double                       r(0.0);
	double                       e(0.0);
	
	// reduce c_* and r 
	MPI::COMM_WORLD.Reduce(&(*task.c_matched.begin()), &(*c_matched.begin()), order + 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&(*task.c_hypo.begin()), &(*c_hypo.begin()), order + 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&task.r, &r, 1, MPI::DOUBLE, MPI::SUM, 0);
	MPI::COMM_WORLD.Reduce(&task.e, &e, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	task.c_matched.swap(c_matched);
	task.c_hypo.swap(c_hypo);
	task.r = r;
	task.e = e;
	
	// reduce g_*
	for (int n = 1; n <= order; ++ n) {
	  // reduce matched counts..
	  send_weights(task.g_matched[n]);
	  send_weights(task.g_hypo[n]);
	}
	send_weights(task.g_reference);
	send_weights(task.g_entropy);
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    return 0.0;
  }
}


template <typename Optimize>
double optimize_batch(const hypothesis_map_type& kbests,
		      const hypothesis_map_type& oracles,
		      const weight_set_type& bounds_lower,
		      const weight_set_type& bounds_upper,
		      weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const size_t id_max = utils::bithack::min(kbests.size(), oracles.size());
  
  typename Optimize::sample_pair_set_type samples;
  samples.reserve(id_max);
  
  int instances_local = 0;
  for (size_t id = 0; id != id_max; ++ id) 
    if (! kbests[id].empty() && ! oracles[id].empty()) {
      samples.push_back(typename Optimize::sample_pair_type(kbests[id], oracles[id]));
      
      ++ instances_local;
    }

  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  
  if (mpi_rank == 0) {
    const double objective = Optimize(samples, bounds_lower, bounds_upper, weights, instances)();
    
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
	
	task_type task(weights, samples, instances);
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


double optimize_mert(const scorer_document_type& scorers,
		     const hypothesis_map_type& kbests,
		     const kbest_map_type& kbest_map,
		     const double scale_min,
		     const double scale_max,
		     const weight_set_type& weights_prev,
		     weight_set_type& weights)
{
  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;
  
  typedef line_search_type::value_type optimum_type;
  
  typedef EnvelopeKBest::line_set_type line_set_type;
  typedef std::deque<line_set_type, std::allocator<line_set_type> > line_map_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  if (kbest_map.empty()) return 0.0;
  
  weight_set_type origin = weights_prev;
  weight_set_type direction;
  
  const size_t weights_size = utils::bithack::min(weights.size(), weights_prev.size());
  for (size_t i = 0; i != weights_size; ++ i)
    direction[i] = weights[i] - weights_prev[i];
  for (size_t i = weights_size; i < weights.size(); ++ i)
    direction[i] = weights[i];
  for (size_t i = weights_size; i < weights_prev.size(); ++ i)
    direction[i] = - weights_prev[i];
  
  bcast_weights(0, origin);
  bcast_weights(0, direction);

  EnvelopeKBest envelopes(origin, direction);
  line_map_type lines;
  
  for (size_t id = 0; id != kbests.size(); ++ id)
    if (! kbests[id].empty()) {
      
      if (kbest_map[id] >= lines.size())
	lines.resize(kbest_map[id] + 1);
      
      envelopes(kbests[id].begin(), kbests[id].end(), std::back_inserter(lines[kbest_map[id]]));
    }
  
  // collect all the segments...
  if (mpi_rank == 0) {
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
    
    segment_document_type segments;
    
    for (size_t seg = 0; seg != lines.size(); ++ seg)
      if (! lines[seg].empty()) {
	
	if (seg >= segments.size())
	  segments.resize(seg + 1);
	
	envelopes(lines[seg]);
	
	EnvelopeKBest::line_set_type::const_iterator liter_end = lines[seg].end();
	for (EnvelopeKBest::line_set_type::const_iterator liter = lines[seg].begin(); liter != liter_end; ++ liter)
	  segments[seg].push_back(std::make_pair(liter->x, liter->hypothesis->score));
      }
    
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::zlib_decompressor());
      is.push(utils::mpi_device_source(rank, envelope_tag, 4096));

      std::string line;
      
      while (std::getline(is, line)) {
	const utils::piece line_piece(line);
	tokenizer_type tokenizer(line_piece);
	
	tokenizer_type::iterator iter = tokenizer.begin();
	if (iter == tokenizer.end()) continue;
	const utils::piece id_str = *iter; 
	
	++ iter;
	if (iter == tokenizer.end() || *iter != "|||") continue;
	
	++ iter;
	if (iter == tokenizer.end()) continue;
	const utils::piece x_str = *iter;
	
	++ iter;
	if (iter == tokenizer.end() || *iter != "|||") continue;
	
	++ iter;
	if (iter == tokenizer.end()) continue;
	const utils::piece score_str = *iter;
	
	const size_t id = utils::lexical_cast<size_t>(id_str);
	
	if (id >= segments.size())
	  segments.resize(id + 1);
	
	segments[id].push_back(std::make_pair(utils::decode_base64<double>(x_str),
					      scorer_type::score_type::decode(score_str)));
      }
    }
    

    line_search_type line_search;
    
    const optimum_type optimum = line_search(segments, scale_min, scale_max, scorers.error_metric());

    const double update = (optimum.lower + optimum.upper) * 0.5;

    if (update != 0.0) {
      const size_t weights_size = utils::bithack::min(origin.size(), direction.size());
      
      weights.clear();
      for (size_t i = 0; i != weights_size; ++ i)
	weights[i] = origin[i] + direction[i] * update;
      for (size_t i = weights_size; i < origin.size(); ++ i)
	weights[i] = origin[i];
      for (size_t i = weights_size; i < direction.size(); ++ i)
	weights[i] = direction[i] * update;
    }
    
    if (debug >= 2)
      std::cerr << "mert update: " << update << " objective: " <<  optimum.objective << std::endl;
    
    return optimum.objective;
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(utils::mpi_device_sink(0, envelope_tag, 4096));

    for (size_t seg = 0; seg != lines.size(); ++ seg)
      if (! lines[seg].empty()) {
	
	envelopes(lines[seg]);
	
	EnvelopeKBest::line_set_type::const_iterator liter_end = lines[seg].end();
	for (EnvelopeKBest::line_set_type::const_iterator liter = lines[seg].begin(); liter != liter_end; ++ liter) {
	  const EnvelopeKBest::line_type& line = *liter;
	  
	  os << seg << " ||| ";
	  utils::encode_base64(line.x, std::ostream_iterator<char>(os));
	  os << " ||| " << line.hypothesis->score->encode() << '\n';
	}
      }
    
    return 0.0;
  }
  
}

void unique_kbest(hypothesis_map_type& kbests)
{
  typedef utils::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
			       std::allocator<hypothesis_type> >::type hypothesis_unique_type;
  
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
		hypothesis_map_type& kbests,
		kbest_map_type&      kbest_map)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  const bool error_metric = scorers.error_metric();
  const double loss_factor = (error_metric ? 1.0 : - 1.0);
  
  parser_type parser;
  kbest_feature_type kbest;

  kbest_map.clear();
  
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
    
    kbest_map.reserve(kbests.size());
    kbest_map.resize(kbests.size());
    for (size_t seg = 0; seg != kbests.size(); ++ seg)
      kbest_map[seg] = seg;
    
  } else {
    
    const size_type refset_size = scorers.size() / kbest_path.size();
    
    for (size_t pos = 0; pos != kbest_path.size(); ++ pos) {
      if (mpi_rank == 0 && debug)
	std::cerr << "reading kbest: " << kbest_path[pos].string() << std::endl;
      
      const size_type refset_offset = refset_size * pos;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_kbest  = kbest_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_kbest)) break;
	
	kbests.resize(kbests.size() + 1);
	kbest_map.push_back(i);
	
	const size_type refset_pos = refset_offset + i;
	
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
    }
  }

  // uniques...
  unique_kbest(kbests);
  kbest_map_type(kbest_map).swap(kbest_map);
}
  
void read_kbest(const scorer_document_type& scorers,
		const path_set_type& kbest_path,
		const path_set_type& oracle_path,
		hypothesis_map_type& kbests,
		hypothesis_map_type& oracles,
		kbest_map_type&      kbest_map)
{  
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const bool error_metric = scorers.error_metric();
  const double loss_factor = (error_metric ? 1.0 : - 1.0);
  
  parser_type parser;
  kbest_feature_type kbest;

  kbest_map.clear();
  
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
    
    kbest_map.reserve(kbests.size());
    kbest_map.resize(kbests.size());
    for (size_t seg = 0; seg != kbests.size(); ++ seg)
      kbest_map[seg] = seg;
    
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
	kbest_map.push_back(i);

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

  kbest_map_type(kbest_map).swap(kbest_map);
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

void reduce_score(score_ptr_pair_type& score)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      is.push(utils::mpi_device_source(rank, score_tag, 256));
      
      std::string score_str1;
      std::string score_str2;
      if (is >> score_str1 && is >> score_str2) {
	
	if (score_str1 != "{}") {
	  if (! score.first)
	    score.first = scorer_type::score_type::decode(score_str1);
	  else
	    *score.first += *scorer_type::score_type::decode(score_str1);
	}
	
	if (score_str2 != "{}") {
	  if (! score.second)
	    score.second = scorer_type::score_type::decode(score_str2);
	  else
	    *score.second += *scorer_type::score_type::decode(score_str2);
	}
      }
    }
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(utils::mpi_device_sink(0, score_tag, 256));
    
    if (score.first)
      os << score.first->encode();
    else
      os << "{}";
    os << ' ';
    if (score.second)
      os << score.second->encode();
    else
      os << "{}";
  }
}

void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("kbest",           po::value<path_set_type>(&kbest_path)->multitoken(),           "kbest path")
    ("oracle",          po::value<path_set_type>(&oracle_path)->multitoken(),          "oracle kbest path")
    ("refset",          po::value<path_set_type>(&refset_files)->multitoken(),         "reference set file(s)")
    ("weights",         po::value<path_type>(&weights_path),                           "initial parameter")
    ("weights-history", po::value<path_set_type>(&weights_history_path)->multitoken(), "parameter history")
    ("output",          po::value<path_type>(&output_path),                            "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("bound-lower", po::value<path_type>(&bound_lower_file),                     "lower bounds definition for feature weights")
    ("bound-upper", po::value<path_type>(&bound_upper_file),                     "upper bounds definition for feature weights")

    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",   po::bool_switch(&learn_lbfgs),   "batch LBFGS algorithm")
    ("learn-sgd",     po::bool_switch(&learn_sgd),     "online SGD algorithm")
    ("learn-xbleu",   po::bool_switch(&learn_xbleu),   "xBLEU algorithm")
    ("learn-mira",    po::bool_switch(&learn_mira),    "online MIRA algorithm")
    ("learn-nherd",   po::bool_switch(&learn_nherd),   "online NHERD algorithm")
    ("learn-arow",    po::bool_switch(&learn_arow),    "online AROW algorithm")
    ("learn-cw",      po::bool_switch(&learn_cw),      "online CW algorithm")
    ("learn-pegasos", po::bool_switch(&learn_pegasos), "online Pegasos algorithm")
    ("learn-cp",      po::bool_switch(&learn_cp),      "cutting plane algorithm")
    ("learn-mcp",     po::bool_switch(&learn_mcp),     "MERT by cutting plane algorithm")
    
    ("regularize-l1",      po::bool_switch(&regularize_l1),      "L1-regularization")
    ("regularize-l2",      po::bool_switch(&regularize_l2),      "L2-regularization")
    
    ("C",             po::value<double>(&C)->default_value(C), "regularization constant")
    ("scale",         po::value<double>(&scale)->default_value(scale), "scaling for weight")
    ("eta0",          po::value<double>(&eta0),                        "\\eta_0 for decay")
    ("order",         po::value<int>(&order)->default_value(order),    "ngram order for xBLEU")
    
    ("annealing", po::bool_switch(&annealing_mode), "annealing")
    ("quenching", po::bool_switch(&quenching_mode), "quenching")
    
    ("temperature",       po::value<double>(&temperature)->default_value(temperature),             "temperature")
    ("temperature-start", po::value<double>(&temperature_start)->default_value(temperature_start), "start temperature for annealing")
    ("temperature-end",   po::value<double>(&temperature_end)->default_value(temperature_end),     "end temperature for annealing")
    ("temperature-rate",  po::value<double>(&temperature_rate)->default_value(temperature_rate),   "annealing rate")

    ("quench-start", po::value<double>(&quench_start)->default_value(quench_start), "start quench for annealing")
    ("quench-end",   po::value<double>(&quench_end)->default_value(quench_end),     "end quench for annealing")
    ("quench-rate",  po::value<double>(&quench_rate)->default_value(quench_rate),   "quenching rate")
    
    ("loss-margin",       po::bool_switch(&loss_margin),       "direct loss margin")
    ("softmax-margin",    po::bool_switch(&softmax_margin),    "softmax margin")
    ("line-search",       po::bool_switch(&line_search),       "perform line search in each iteration")
    ("mert-search",       po::bool_switch(&mert_search),       "perform one-dimensional mert")
    ("mert-search-local", po::bool_switch(&mert_search_local), "perform local one-dimensional mert")
    ("sample-vector",     po::bool_switch(&sample_vector),     "perform samling")
    ("direct-loss",       po::bool_switch(&direct_loss),       "compute loss by directly treating hypothesis score")
    ("conservative-loss", po::bool_switch(&conservative_loss), "conservative loss")
    
    ("scale-fixed", po::bool_switch(&scale_fixed), "fixed scaling")

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
