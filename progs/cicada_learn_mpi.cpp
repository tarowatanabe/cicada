//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// learning from hypergraphs...
//
// we assume two inputs, one for partition and the other for marginals
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"
#include "cicada_text_impl.hpp"

#include "cicada/prune.hpp"
#include "cicada/operation/functional.hpp"
#include "cicada/expected_ngram.hpp"
#include "cicada/symbol_vector.hpp"

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
#include "utils/mathop.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_trie.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>
#include <boost/math/special_functions/expm1.hpp>

#include "lbfgs.h"
#include "lbfgs_error.hpp"

typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type forest_path;
path_set_type intersected_path;
path_set_type refset_path;
path_type weights_path;
path_set_type weights_history_path;
path_type output_path = "-";
path_type output_objective_path;

path_type bound_lower_file;
path_type bound_upper_file;

int iteration = 100;

bool learn_lbfgs = false;
bool learn_sgd = false;
bool learn_xbleu = false;
bool learn_mira = false;
bool learn_nherd = false;
bool learn_arow = false;
bool learn_cw = false;
bool learn_pegasos = false;

bool regularize_l1 = false;
bool regularize_l2 = false;
bool regularize_entropy = false;
double C = 1.0;
double C2 = 1.0;
double scale = 1.0;
double eta0 = 0.2;
int order = 4;

double temperature = 0.5;
double temperature_start = 1000;
double temperature_end = 0.001;
double temperature_rate = 0.5;

bool loss_margin = false; // margin by loss, not rank-loss
bool softmax_margin = false;
bool annealing_mode = false;

// scorers
std::string scorer_name = "bleu:order=4,exact=true";
bool scorer_list = false;

bool unite_forest = false;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

template <typename Optimize>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights);
template <typename Optimize>
double optimize_xbleu(const hypergraph_set_type& forests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights);

template <typename Optimize, typename Generator>
double optimize_online(const hypergraph_set_type& graphs_forest,
		       const hypergraph_set_type& graphs_intersected,
		       weight_set_type& weights,
		       Generator& generator);

template <typename Optimizer>
struct OptimizeOnline;
template <typename Optimizer>
struct OptimizeOnlineMargin;
struct OptimizeXBLEU;
struct OptimizeLBFGS;

void read_refset(const path_set_type& files,
		 scorer_document_type& scorers);
void read_forest(const path_set_type& forest_path,
		 const scorer_document_type& scorers,
		 hypergraph_set_type& forest,
		 scorer_document_type& scorers_forest);
void read_forest(const path_set_type& forest_path,
		 const path_set_type& intersected_path,
		 hypergraph_set_type& graphs_forest,
		 hypergraph_set_type& graphs_intersected);
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
    
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw + learn_pegasos + learn_nherd + learn_xbleu > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,sgd,mira,arow,cw}");
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw + learn_pegasos + learn_nherd + learn_xbleu == 0)
      learn_lbfgs = true;

    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));
    
    if (forest_path.empty())
      throw std::runtime_error("no forest?");
    if (! learn_xbleu && intersected_path.empty())
      throw std::runtime_error("no intersected forest?");
    if (learn_xbleu && refset_path.empty())
      throw std::runtime_error("no reference translations?");
    if (learn_xbleu && order <= 0)
      throw std::runtime_error("invalid ngram order");
    
    if (! bound_lower_file.empty())
      if (bound_lower_file != "-" && ! boost::filesystem::exists(bound_lower_file))
	throw std::runtime_error("no lower-bound file? " + bound_lower_file.string());
    
    if (! bound_upper_file.empty())
      if (bound_upper_file != "-" && ! boost::filesystem::exists(bound_upper_file))
	throw std::runtime_error("no upper-bound file? " + bound_upper_file.string());

    scorer_document_type scorers(scorer_name);
    if (! refset_path.empty())
      read_refset(refset_path, scorers);
    
    hypergraph_set_type graphs_forest;
    hypergraph_set_type graphs_intersected;
    
    if (! learn_xbleu)
      read_forest(forest_path, intersected_path, graphs_forest, graphs_intersected);
    else {
      scorer_document_type scorers_forest(scorer_name);
      
      read_forest(forest_path, scorers, graphs_forest, scorers_forest);
      
      scorers_forest.swap(scorers);
    }
    
    
    weight_set_type weights;
    if (mpi_rank ==0 && ! weights_path.empty()) {
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
    
    if (learn_sgd) {
      if (regularize_l1)
	objective = optimize_online<OptimizeOnline<OptimizerSGDL1> >(graphs_forest, graphs_intersected, weights, generator);
      else
	objective = optimize_online<OptimizeOnline<OptimizerSGDL2> >(graphs_forest, graphs_intersected, weights, generator);
    } else if (learn_mira)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerMIRA> >(graphs_forest, graphs_intersected, weights, generator);
    else if (learn_arow)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerAROW> >(graphs_forest, graphs_intersected, weights, generator);
    else if (learn_cw)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerCW> >(graphs_forest, graphs_intersected, weights, generator);
    else if (learn_pegasos)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerPegasos> >(graphs_forest, graphs_intersected, weights, generator);
    else if (learn_nherd)
      objective = optimize_online<OptimizeOnlineMargin<OptimizerNHERD> >(graphs_forest, graphs_intersected, weights, generator);
    else if (learn_xbleu)
      objective = optimize_xbleu<OptimizeXBLEU>(graphs_forest, scorers, weights);
    else
      objective = optimize_batch<OptimizeLBFGS>(graphs_forest, graphs_intersected, weights);

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
  
  OptimizeOnline(Optimizer& __optimizer)
    : optimizer(__optimizer) {}
  
  typedef Optimizer optimizer_type;
  
  typedef typename optimizer_type::weight_type   weight_type;
  typedef typename optimizer_type::gradient_type gradient_type;    
  
  typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
    
  struct gradients_type
  {
    typedef gradient_type value_type;

    template <typename Index>
    gradient_type& operator[](Index)
    {
      return gradient;
    }

    void clear() { gradient.clear(); }
      
    gradient_type gradient;
  };
    
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
    typedef gradient_type value_type;

    feature_function(const weight_set_type& __weights, const double& __scale) : weights(__weights), scale(__scale) {}

    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      // p_e r_e
      gradient_type grad;
	
      const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(edge.features, weights) * scale);
	
      feature_set_type::const_iterator fiter_end = edge.features.end();
      for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	if (fiter->second != 0.0)
	  grad[fiter->first] = weight_type(fiter->second) * weight;
      
      return grad;
    }
    
    const weight_set_type& weights;
    const double scale;
  };
  
  
  void operator()(const hypergraph_type& hypergraph_intersected,
		  const hypergraph_type& hypergraph_forest)
  {
    gradients.clear();
    gradients_intersected.clear();
    
    inside.clear();
    inside_intersected.clear();
    
    inside.reserve(hypergraph_forest.nodes.size());
    inside.resize(hypergraph_forest.nodes.size(), weight_type());
    cicada::inside_outside(hypergraph_forest, inside, gradients,
			   weight_function(optimizer.weights, optimizer.weight_scale),
			   feature_function(optimizer.weights, optimizer.weight_scale));
    
    inside_intersected.reserve(hypergraph_intersected.nodes.size());
    inside_intersected.resize(hypergraph_intersected.nodes.size(), weight_type());
    cicada::inside_outside(hypergraph_intersected, inside_intersected, gradients_intersected,
			   weight_function(optimizer.weights, optimizer.weight_scale),
			   feature_function(optimizer.weights, optimizer.weight_scale));
    
    gradient_type& gradient = gradients.gradient;
    weight_type& Z = inside.back();
    
    gradient_type& gradient_intersected = gradients_intersected.gradient;
    weight_type& Z_intersected = inside_intersected.back();
    
    gradient /= Z;
    gradient_intersected /= Z_intersected;
    
    optimizer(gradients_intersected.gradient,
	      gradients.gradient,
	      Z_intersected,
	      Z);
  }
  
  Optimizer& optimizer;

  gradients_type gradients;
  gradients_type gradients_intersected;
  weights_type   inside;
  weights_type   inside_intersected;
};


template <typename Optimizer>
struct OptimizeOnlineMargin
{
  
  OptimizeOnlineMargin(Optimizer& __optimizer)
    : optimizer(__optimizer) {}
  
  typedef Optimizer optimizer_type;
  
  typedef typename optimizer_type::weight_type   weight_type;
  typedef typename optimizer_type::gradient_type gradient_type;    
  
  typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;
  
  struct Accumulated
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > accumulated_type;
    
    typedef accumulated_type value_type;
    
    accumulated_type& operator[](size_t index)
    {
      return accumulated;
    }
    
    void clear() { accumulated.clear(); }
    
    value_type accumulated;
  };
  typedef Accumulated accumulated_type;

  struct count_function
  {
    typedef cicada::semiring::Log<double> value_type;
    
    template <typename Edge>
    value_type operator()(const Edge& x) const
    {
      return cicada::semiring::traits<value_type>::exp(0.0);
    }
  };
  
  
  struct feature_count_function
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > accumulated_type;
    typedef accumulated_type value_type;
    
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      accumulated_type accumulated;
      
      feature_set_type::const_iterator fiter_end = edge.features.end();
      for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	if (fiter->second != 0.0)
	  accumulated[fiter->first] = weight_type(fiter->second);
      
      return accumulated;
    }
  };
  
  typedef std::vector<typename count_function::value_type, std::allocator<typename count_function::value_type> > count_set_type;
  
  void operator()(const hypergraph_type& hypergraph_intersected,
		  const hypergraph_type& hypergraph_forest)
  {
    typedef cicada::operation::weight_scaled_function<cicada::semiring::Tropical<double> > function_type;
    
    if (margin_kbest > 0)
      cicada::prune_kbest(hypergraph_forest, pruned_forest, function_type(optimizer.weights, optimizer.weight_scale), margin_kbest);
    else
      cicada::prune_beam(hypergraph_forest, pruned_forest, function_type(optimizer.weights, optimizer.weight_scale), margin_beam);
    
    if (margin_kbest > 0)
      cicada::prune_kbest(hypergraph_intersected, pruned_intersected, function_type(optimizer.weights, - optimizer.weight_scale), margin_kbest);
    else
      cicada::prune_beam(hypergraph_intersected, pruned_intersected, function_type(optimizer.weights, - optimizer.weight_scale), margin_beam);
    
    counts_intersected.clear();
    counts_forest.clear();
    
    counts_intersected.resize(pruned_intersected.nodes.size());
    counts_forest.resize(pruned_forest.nodes.size());
    
    accumulated_intersected.clear();
    accumulated_forest.clear();
    
    cicada::inside_outside(pruned_intersected, counts_intersected, accumulated_intersected, count_function(), feature_count_function());
    cicada::inside_outside(pruned_forest,      counts_forest,      accumulated_forest,      count_function(), feature_count_function());
    
    features_intersected.assign(accumulated_intersected.accumulated);
    features_forest.assign(accumulated_forest.accumulated);
    
    features_intersected *= (1.0 / double(counts_intersected.back()));
    features_forest      *= (1.0 / double(counts_forest.back()));
    
    
    // use the collected features...!
    optimizer(features_intersected, features_forest);
  }
  
  Optimizer& optimizer;

  hypergraph_type pruned_intersected;
  hypergraph_type pruned_forest;

  count_set_type counts_intersected;
  count_set_type counts_forest;

  accumulated_type accumulated_intersected;
  accumulated_type accumulated_forest;
  
  feature_set_type features_intersected;
  feature_set_type features_forest;
  
  double margin_beam;
  int    margin_kbest;
};

template <typename Optimize, typename Generator>
double optimize_online(const hypergraph_set_type& graphs_forest,
		       const hypergraph_set_type& graphs_intersected,
		       weight_set_type& weights,
		       Generator& generator)
{
  typedef std::vector<int, std::allocator<int> > id_set_type;
  typedef typename Optimize::optimizer_type optimizer_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  id_set_type ids(graphs_forest.size());
  int instances_local = 0;
  
  for (size_t id = 0; id != ids.size(); ++ id) {
    ids[id] = id;
    instances_local += (graphs_intersected[id].is_valid() && graphs_forest[id].is_valid());
  }
  
  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  
  optimizer_type optimizer(instances, C);
  Optimize opt(optimizer);
  
  optimizer.weights = weights;
  
  if (mpi_rank == 0) {
    double objective = 0.0;
    
    for (int iter = 0; iter < iteration; ++ iter) {
      
      for (int rank = 1; rank < mpi_size; ++ rank)
	MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, notify_tag);
      
      bcast_weights(0, optimizer.weights);
      
      optimizer.initialize();
      
      for (size_t id = 0; id != ids.size(); ++ id)
	if (graphs_intersected[ids[id]].is_valid() && graphs_forest[ids[id]].is_valid())
	  opt(graphs_intersected[ids[id]], graphs_forest[ids[id]]);
      
      optimizer.finalize();
      
      boost::random_number_generator<Generator> gen(generator);
      std::random_shuffle(ids.begin(), ids.end(), gen);
      
      optimizer.weights *= (optimizer.samples + 1);
      reduce_weights(optimizer.weights);
      
      objective = 0.0;
      MPI::COMM_WORLD.Reduce(&optimizer.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      int samples = 0;
      int samples_local = (optimizer.samples + 1);
      MPI::COMM_WORLD.Reduce(&samples_local, &samples, 1, MPI::INT, MPI::SUM, 0);
      
      optimizer.weights *= (1.0 / samples);
      
      if (debug >= 2)
	std::cerr << "objective: " << objective << std::endl;
    }
    
    // send termination!
    for (int rank = 1; rank < mpi_size; ++ rank)
      MPI::COMM_WORLD.Send(0, 0, MPI::INT, rank, termination_tag);

    weights.swap(optimizer.weights);
    
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
	
	optimizer.initialize();
	
	for (size_t id = 0; id != ids.size(); ++ id)
	  if (graphs_intersected[ids[id]].is_valid() && graphs_forest[ids[id]].is_valid())
	    opt(graphs_intersected[ids[id]], graphs_forest[ids[id]]);
	
	optimizer.finalize();
	
	boost::random_number_generator<Generator> gen(generator);
	std::random_shuffle(ids.begin(), ids.end(), gen);
	
	optimizer.weights *= (optimizer.samples + 1);
	send_weights(optimizer.weights);
	
	double objective = 0.0;
	MPI::COMM_WORLD.Reduce(&optimizer.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	int samples = 0;
	int samples_local = (optimizer.samples + 1);
	MPI::COMM_WORLD.Reduce(&samples_local, &samples, 1, MPI::INT, MPI::SUM, 0);
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    return 0.0;
  }
}

struct OptimizeXBLEU
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  OptimizeXBLEU(const hypergraph_set_type& __forests,
		const scorer_document_type& __scorers,
		weight_set_type& __weights,
		const size_t& __instances,
		const double& __lambda)
    : forests(__forests),
      scorers(__scorers),
      weights(__weights),
      instances(__instances),
      lambda(__lambda) {}
  
  const hypergraph_set_type& forests;
  const scorer_document_type& scorers;
  weight_set_type& weights;
  
  size_t instances;
  double lambda;

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
    
    if (debug >= 3)
      std::cerr << "lbfgs weights:" << std::endl
		<< weights << std::flush;

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

    
    typedef cicada::Symbol       word_type;
    typedef cicada::SymbolVector ngram_type;
    
    struct Count
    {
      weight_type c;
      weight_type mu_prime;
      
      Count() : c(), mu_prime() {}
    };
    typedef Count count_type;
    typedef std::vector<count_type, std::allocator<count_type> > count_set_type;
    typedef std::vector<ngram_type, std::allocator<ngram_type> > ngram_set_type;
    
    typedef utils::indexed_trie<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> > index_set_type;

    typedef std::vector<double, std::allocator<double> > ngram_counts_type;
    typedef std::vector<weight_set_type, std::allocator<weight_set_type> > feature_counts_type;
    
    struct CollectCounts
    {
      CollectCounts(index_set_type& __index,
		    ngram_set_type& __ngrams,
		    count_set_type& __counts)
	: index(__index), ngrams(__ngrams), counts(__counts) {}
      
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
      }
      
      index_set_type& index;
      ngram_set_type& ngrams;
      count_set_type& counts;
    };
    
    struct CollectExpectation
    {
      CollectExpectation(const index_set_type& __index,
			 const ngram_set_type& __ngrams,
			 const count_set_type& __counts,
			 const weights_type& __matched,
			 const weights_type& __hypo,
			 const weight_set_type& __weights,
			 const double& __scale,
			 gradients_type& __gradients_matched,
			 gradients_type& __gradients_hypo)
	: index(__index), ngrams(__ngrams), counts(__counts), matched(__matched), hypo(__hypo),
	  weights(__weights), scale(__scale), 
	  gradients_matched(__gradients_matched),
	  gradients_hypo(__gradients_hypo)
      {}

      template <typename Edge, typename Weight, typename Counts>
      void operator()(const Edge& edge, const Weight& weight, Counts& __counts)
      {
	for (int n = 1; n <= order; ++ n) 
	  if (hypo[n] > weight_type()) {
	    const weight_type scale_matched = weight * matched[n];
	    const weight_type scale_hypo    = weight * hypo[n];
	    
	    feature_set_type::const_iterator fiter_end = edge.features.end();
	    for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	      if (fiter->second != 0.0) {
		const weight_type value(fiter->second * scale);
		
		gradients_matched[n][fiter->first] -= value * scale_matched;
		gradients_hypo[n][fiter->first]    -= value * scale_hypo;
	      }
	  }
      }
      
      template <typename Edge, typename Weight, typename Counts, typename Iterator>
      void operator()(const Edge& edge, const Weight& weight, Counts& __counts, Iterator first, Iterator last)
      {
	if (first == last) return;
	
	index_set_type::id_type id = index.root();
	for (Iterator iter = first; iter != last; ++ iter) {
	  id = index.find(id, *iter);
	  if (index.is_root(id))
	    throw std::runtime_error("no ngram?");
	}
	
	if (ngrams[id].empty())
	  throw std::runtime_error("no ngram storage?");
	
	const size_type order = ngrams[id].size();
	
	const weight_type scale_matched = weight * counts[id].mu_prime;
	const weight_type scale_hypo    = weight;

	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0) {
	    const weight_type value(fiter->second * scale);
	    
	    gradients_matched[order][fiter->first] += value * scale_matched;
	    gradients_hypo[order][fiter->first]    += value * scale_hypo;
	  }
      }
      
      const index_set_type& index;
      const ngram_set_type& ngrams;
      const count_set_type& counts;
      const weights_type& matched;
      const weights_type& hypo;
      
      const weight_set_type& weights;
      const double& scale;
      
      gradients_type& gradients_matched;
      gradients_type& gradients_hypo;
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
	value_type(const feature_set_type& __features, const weight_set_type& __weights, const double& __scale)
	  : features(__features), weights(__weights), scale(__scale) {}
	
	friend
	value_type operator*(const value_type& x, const entropy_weight_type& weight)
	{
	  const_cast<value_type&>(x).inside_outside = weight;
	  return x;
	}
	
	entropy_weight_type inside_outside;
	
	const feature_set_type& features;
	const weight_set_type& weights;
	const double scale;
      };
      
      entropy_gradient_function(const weight_set_type& __weights, const double& __scale)
	: weights(__weights), scale(__scale) {}
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	return value_type(edge.features, weights, scale);
      }
      
      const weight_set_type& weights;
      const double scale;
    };
    

    struct entropy_gradient_type
    {
      typedef cicada::FeatureVectorUnordered<weight_type, std::allocator<weight_type> > accumulated_type;

      struct proxy_type
      {
	proxy_type(accumulated_type& __dZ, accumulated_type& __dR) : dZ(__dZ), dR(__dR) {}
	
	proxy_type& operator+=(const entropy_gradient_function::value_type& x) 
	{
	  const double value = cicada::dot_product(x.features, x.weights);
	  const double log_p_e = value * x.scale;
	  const weight_type p_e = cicada::semiring::traits<weight_type>::exp(log_p_e);

	  // dZ += \lnabla p_e * x.inside_outside.p;
	  // dR += (1 + \log p_e) * \nalba p_e * x.inside_outside.p + \lnabla p_e * x.inside_outside.r;
	  
	  feature_set_type::const_iterator fiter_end = x.features.end();
	  for (feature_set_type::const_iterator fiter = x.features.begin(); fiter != fiter_end; ++ fiter) 
	    if (fiter->second != 0.0) {
	      const weight_type value(fiter->second * x.scale);
	      
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
    

    Task(const hypergraph_set_type& __forests,
	 const scorer_document_type& __scorers,
	 const weight_set_type& __weights)
      : forests(__forests), scorers(__scorers), weights(__weights),
	c_matched(order + 1),
	c_hypo(order + 1),
	g_matched(order + 1),
	g_hypo(order + 1),
	r(0) {}
    
    void operator()()
    {
      const word_type __tmp;
      
      index_set_type index;
      ngram_set_type ngrams;
      count_set_type counts;
      weights_type   matched(order + 1);
      weights_type   hypo(order + 1);
      
      weights_type   counts_matched(order + 1);
      weights_type   counts_hypo(order + 1);
      gradients_type gradients_matched(order + 1);
      gradients_type gradients_hypo(order + 1);
      
      weight_type          entropy;
      entropy_weights_type entropy_inside;
      gradient_type        gradient_entropy;
      
      for (size_t n = 0; n != g_matched.size(); ++ n) {
	gradients_matched[n].allocate();
	gradients_hypo[n].allocate();
	
	g_matched[n].clear();
	g_hypo[n].clear();
      }

      gradient_entropy.allocate();
      g_entropy.clear();
      
      std::fill(counts_matched.begin(), counts_matched.end(), weight_type());
      std::fill(counts_hypo.begin(), counts_hypo.end(), weight_type());
      std::fill(c_matched.begin(), c_matched.end(), 0.0);
      std::fill(c_hypo.begin(), c_hypo.end(), 0.0);
      r = 0.0;
      e = 0.0;
      
      for (size_t id = 0; id != forests.size(); ++ id) {
	const hypergraph_type& forest = forests[id];
	
	if (! forest.is_valid()) continue;
	
	const cicada::eval::BleuScorer* scorer = dynamic_cast<const cicada::eval::BleuScorer*>(scorers[id].get());
	
	if (! scorer)
	  throw std::runtime_error("we do not have bleu scorer...");
	
	// here, we will implement forest xBLEU...
	
	
	// first, collect expected ngrams
	index.clear();
	counts.clear();
	ngrams.clear();
		
	cicada::expected_ngram(forest,
			       cicada::operation::weight_scaled_function<weight_type>(weights, scale),
			       CollectCounts(index, ngrams, counts),
			       index,
			       order);
	
	// second, commpute clipped ngram counts (\mu')
	std::fill(matched.begin(), matched.end(), weight_type());
	std::fill(hypo.begin(), hypo.end(), weight_type());
	
	for (size_type i = 0; i != ngrams.size(); ++ i) 
	  if (! ngrams[i].empty()) {
	    const size_type    order = ngrams[i].size();
	    const weight_type& count = counts[i].c;
	    const weight_type  clip = scorer->find(ngrams[i]);
	    
	    counts[i].mu_prime = derivative_clip_count(count, clip);
	    
	    // collect counts for further inside/outside
	    matched[order] += counts[i].c * counts[i].mu_prime;
	    hypo[order]    += counts[i].c;
	    
	    // collect global counts
	    counts_matched[order] += clip_count(count, clip);
	    counts_hypo[order]    += counts[i].c;
	  }
	
	r += scorer->reference_length(hypo[1]);
	
	if (debug >= 4)
	  for (int n = 1; n <= order; ++ n)
	    std::cerr << "order: " << n << " matched: " << matched[n] << " hypo: " << hypo[n] << std::endl;
	
	// third, collect feature expectation, \hat{m} - m and \hat{h} - h
	cicada::expected_ngram(forest,
			       cicada::operation::weight_scaled_function<weight_type>(weights, scale),
			       CollectExpectation(index, ngrams, counts, matched, hypo,
						  weights, scale,
						  gradients_matched, gradients_hypo),
			       index,
			       order);

	
	// forth, compute entorpy...
	entropy_inside.clear();
	entropy_inside.resize(forest.nodes.size(), entropy_weight_type());

	entropy_gradient_type entropy_gradient;
	
	cicada::inside_outside(forest,
			       entropy_inside,
			       entropy_gradient,
			       entropy_function(weights, scale),
			       entropy_gradient_function(weights, scale));
	
	const weight_type& Z = entropy_inside.back().p;
	const weight_type& R = entropy_inside.back().r;
	
	const weight_type entropy_segment = weight_type(cicada::semiring::log(Z)) - (R / Z);
	
	if (debug >= 4)
	  std::cerr << "entropy: " << double(entropy_segment) << std::endl;
	
	entropy += entropy_segment;
	
	const entropy_gradient_type::accumulated_type& dZ = entropy_gradient.dZ;
	const entropy_gradient_type::accumulated_type& dR = entropy_gradient.dR;

	// compute...
	// \frac{\nabla Z}{Z} - \frac{Z \nabla \bar{r} - \bar{r} \nabla Z}{Z^2}
	
	entropy_gradient_type::accumulated_type::const_iterator ziter_end = dZ.end();
	for (entropy_gradient_type::accumulated_type::const_iterator ziter = dZ.begin(); ziter != ziter_end; ++ ziter)
	  gradient_entropy[ziter->first] += ziter->second * ((cicada::semiring::traits<weight_type>::one() / Z) + R / (Z * Z));
	
	entropy_gradient_type::accumulated_type::const_iterator riter_end = dR.end();
	for (entropy_gradient_type::accumulated_type::const_iterator riter = dR.begin(); riter != riter_end; ++ riter)
	  gradient_entropy[riter->first] -= riter->second / Z;
      }
      
      std::copy(counts_matched.begin(), counts_matched.end(), c_matched.begin());
      std::copy(counts_hypo.begin(), counts_hypo.end(), c_hypo.begin());
      
      for (size_t n = 1; n != g_matched.size(); ++ n) {
	g_matched[n].allocate();
	g_hypo[n].allocate();
	
	std::copy(gradients_matched[n].begin(), gradients_matched[n].end(), g_matched[n].begin());
	std::copy(gradients_hypo[n].begin(), gradients_hypo[n].end(), g_hypo[n].begin());
      }
      
      g_entropy.allocate();
      std::copy(gradient_entropy.begin(), gradient_entropy.end(), g_entropy.begin());
      
      e = entropy;
    }
    
    const hypergraph_set_type& forests;
    const scorer_document_type& scorers;
    const weight_set_type& weights;
    
    ngram_counts_type   c_matched;
    ngram_counts_type   c_hypo;
    feature_counts_type g_matched;
    feature_counts_type g_hypo;
    weight_set_type     g_entropy;
    double r;
    double e;
  };
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int size,
				  const lbfgsfloatval_t step)
  {
    typedef Task                  task_type;
    
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
    
    task_type task(optimizer.forests, optimizer.scorers, optimizer.weights);
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
    reduce_weights(task.g_entropy);
    
    // smoothing...
#if 1
    {
      double smoothing = 1e-40;
      for (int n = 1; n <= order; ++ n) {
	if (task.c_hypo[n] > 0.0 && task.c_matched[n] <= 0.0)
	  task.c_matched[n] = smoothing;
	smoothing *= 0.1;
      }
    }
#endif

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
    
    //std::cerr << "P: " << P << " B: " << B << " C: " << C << std::endl;

        // xBLEU...
    const double objective_bleu = exp_P * B;
    const double entropy = task.e / optimizer.instances;
    
    // compute g..
    std::fill(g, g + size, 0.0);
    //std::copy(task.g.begin(), task.g.end(), g);

    // entropy
    if (regularize_entropy) {
      // 0.5 * (entropy - C2)^2
      // 
      // thus, derivative is: (entropy - C2) * \nabla entropy
      // we need to consider average of entropy...
      
      std::transform(task.g_entropy.begin(), task.g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), (entropy - C2) * (1.0 / optimizer.instances)));
    } else
      std::transform(task.g_entropy.begin(), task.g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), - temperature / optimizer.instances));
    
    for (int n = 1; n <= order; ++ n)  {
      if (task.c_hypo[n] > 0.0) {
	const double factor_matched = - (exp_P * B / order) / task.c_matched[n];
	const double factor_hypo    = - (exp_P * B / order) / task.c_hypo[n];
	
	for (size_t i = 0; i != static_cast<size_t>(size); ++ i) {
	  g[i] += factor_matched * task.g_matched[n][i];
	  g[i] -= factor_hypo * task.g_hypo[n][i];
	}
      }
    }
    
    if (task.c_hypo[1] > 0.0) {
      // I think the missed exp(P) is a bug in Rosti et al. (2011)
      const double factor = - exp_P * C_dC / task.c_hypo[1];
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i)
	g[i] += factor * task.g_hypo[1][i];
    }
    
    if (debug >= 3) {
      std::cerr << "grad:" << std::endl;
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i)
	if (g[i] != 0.0 && feature_type(i) != feature_type())
	  std::cerr << feature_type(i) << ' ' << g[i] << std::endl;
    }
        
    // we need to minimize negative bleu... + regularized by average entropy...
    double objective = - objective_bleu + (regularize_entropy ? 0.5 * (entropy - C2) * (entropy - C2) : - temperature * entropy);
    //double objective = - objective_bleu;
    
    if (regularize_l2) {
      double norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(size); ++ i) {
	g[i] += optimizer.lambda * x[i];
	norm += x[i] * x[i];
      }
      objective += 0.5 * optimizer.lambda * norm;
    }
    
    if (debug >= 2)
      std::cerr << "objective: " << objective
		<< " xBLEU: " << objective_bleu
		<< " BP: " << B
		<< " entropy: " << entropy
		<< std::endl;

#if 1
    // keep the best so forth...
    if (objective <= optimizer.objective_opt) {
      optimizer.objective_opt = objective;
      optimizer.weights_opt = optimizer.weights;
    }
#endif
    
    return objective;
  }
};
  

struct OptimizeLBFGS
{

  OptimizeLBFGS(const hypergraph_set_type& __graphs_forest,
		const hypergraph_set_type& __graphs_intersected,
		weight_set_type& __weights,
		const size_t& __instances)
    : graphs_forest(__graphs_forest),
      graphs_intersected(__graphs_intersected),
      weights(__weights),
      instances(__instances) {}

  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
  size_t instances;

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
    
    const int result = lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    if (debug)
      std::cerr << "lbfgs: " << lbfgs_error(result) << std::endl;

    return objective;
  }

  
  struct Task
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_static_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    Task(const weight_set_type& __weights,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected,
	 const size_t& __instances)
      : weights(__weights),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected),
	instances(__instances) {}
    
    struct gradients_type
    {
      typedef gradient_type value_type;

      template <typename Index>
      gradient_type& operator[](Index)
      {
	return gradient;
      }

      void clear() { gradient.clear(); }
      
      gradient_type gradient;
    };
    
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
      typedef gradient_type value_type;

      feature_function(const weight_set_type& __weights) : weights(__weights) {}
      

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e r_e
	gradient_type grad;
	
	const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(edge.features, weights));
	
	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0)
	    grad[fiter->first] = weight_type(fiter->second) * weight;
	
	return grad;
      }
      
      const weight_set_type& weights;
    };
    
    void operator()()
    {
      hypergraph_type hypergraph;
            
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;
      
      gradient_static_type  feature_expectations;

      g.clear();
      objective = 0.0;
      
      for (size_t i = 0; i != graphs_forest.size(); ++ i) {
	const hypergraph_type& hypergraph_forest      = graphs_forest[i];
	const hypergraph_type& hypergraph_intersected = graphs_intersected[i];

	if (! hypergraph_forest.is_valid() || ! hypergraph_intersected.is_valid()) continue;
	  
	gradients.clear();
	gradients_intersected.clear();
	  
	inside.clear();
	inside_intersected.clear();

	inside.reserve(hypergraph_forest.nodes.size());
	inside.resize(hypergraph_forest.nodes.size(), weight_type());
	inside_outside(hypergraph_forest, inside, gradients, weight_function(weights), feature_function(weights));
	  
	inside_intersected.reserve(hypergraph_intersected.nodes.size());
	inside_intersected.resize(hypergraph_intersected.nodes.size(), weight_type());
	inside_outside(hypergraph_intersected, inside_intersected, gradients_intersected, weight_function(weights), feature_function(weights));
	  
	gradient_type& gradient = gradients.gradient;
	weight_type& Z = inside.back();
	  
	gradient_type& gradient_intersected = gradients_intersected.gradient;
	weight_type& Z_intersected = inside_intersected.back();
	  
	gradient /= Z;
	gradient_intersected /= Z_intersected;
	  
	feature_expectations -= gradient_intersected;
	feature_expectations += gradient;
	  
	const double margin = log(Z_intersected) - log(Z);
	  
	objective -= margin;
	  
	if (debug >= 3)
	  std::cerr << "margin: " << margin << std::endl;
      }
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(feature_expectations.begin(), feature_expectations.end(), g.begin());
      
      // normalize!
      objective /= instances;
      std::transform(g.begin(), g.end(), g.begin(), std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    }
    
    const weight_set_type& weights;

    const hypergraph_set_type& graphs_forest;
    const hypergraph_set_type& graphs_intersected;
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
    
    task_type task(optimizer.weights, optimizer.graphs_forest, optimizer.graphs_intersected, optimizer.instances);
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
};

template <typename Optimize>
double optimize_xbleu(const hypergraph_set_type& forests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
    
  int instances_local = 0;
  for (size_t id = 0; id != forests.size(); ++ id)
    instances_local += forests[id].is_valid();
  
  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  
  if (mpi_rank == 0) {
    Optimize optimizer(forests, scorers, weights, instances, C);
    
    const double objective = optimizer();

    if (debug >= 3)
      std::cerr << "final weights:" << std::endl
		<< weights << std::flush;
    
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
	
	task_type task(forests, scorers, weights);
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
	
	// reduce g_*
	for (int n = 1; n <= order; ++ n) {
	  // reduce matched counts..
	  send_weights(task.g_matched[n]);
	  send_weights(task.g_hypo[n]);
	}
	send_weights(task.g_entropy);
	
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    return 0.0;
  }
}

template <typename Optimize>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  const size_t id_max = utils::bithack::min(graphs_forest.size(), graphs_intersected.size());
  
  int instances_local = 0;
  for (size_t id = 0; id != id_max; ++ id)
    instances_local += (graphs_intersected[id].is_valid() && graphs_forest[id].is_valid());
  
  int instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  
  if (mpi_rank == 0) {
    Optimize optimizer(graphs_forest, graphs_intersected, weights, instances);
    
    const double objective = optimizer();
    
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
	
	task_type task(weights, graphs_forest, graphs_intersected, instances);
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

void read_forest(const path_set_type& forest_path,
		 const scorer_document_type& scorers,
		 hypergraph_set_type& forests,
		 scorer_document_type& scorers_forest)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (unite_forest) {
    size_t id;
    
    std::string line;
    hypergraph_type graph;
    
    for (path_set_type::const_iterator piter = forest_path.begin(); piter != forest_path.end(); ++ piter) {
    
      if (mpi_rank == 0 && debug)
	std::cerr << "reading forest: " << piter->string() << std::endl;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_forest = (*piter) / file_name;
      
	if (! boost::filesystem::exists(path_forest)) break;
	
	utils::compress_istream is(path_forest);
	std::getline(is, line);
	
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end  = line.end();
	
	if (! parse_id(id, iter, end))
	  throw std::runtime_error("invalid id input: " + path_forest.string());
	if (id != i)
	  throw std::runtime_error("invalid id input: " + path_forest.string());
	
	if (id >= forests.size())
	  forests.resize(id + 1);
	
	if (! graph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path_forest.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path_forest.string());
	
	forests[id].unite(graph);
      }
    }

    if (forests.size() > scorers.size())
      throw std::runtime_error("invalid scorers");

    scorers_forest = scorers;
  } else {
    size_t id;
    
    std::string line;
    
    for (size_t pos = 0; pos != forest_path.size(); ++ pos) {
      if (mpi_rank == 0 && debug)
	std::cerr << "reading forest: " << forest_path[pos].string() << std::endl;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_forest = forest_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_forest)) break;
	
	utils::compress_istream is(path_forest);
	std::getline(is, line);
	
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end  = line.end();
	
	if (! parse_id(id, iter, end))
	  throw std::runtime_error("invalid id input: " + path_forest.string());
	if (id != i)
	  throw std::runtime_error("invalid id input: " + path_forest.string());
	
	forests.push_back(hypergraph_type());
	
	if (id >= scorers.size())
	  throw std::runtime_error("invalid scorers");
	
	scorers_forest.push_back(scorers[id]);
	
	if (! forests.back().assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path_forest.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path_forest.string());
      }
    }
  }
}


void read_forest(const path_set_type& forest_path,
		 const path_set_type& intersected_path,
		 hypergraph_set_type& graphs_forest,
		 hypergraph_set_type& graphs_intersected)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  if (unite_forest) {
    size_t id_forest;
    size_t id_intersected;
  
    std::string line;

    hypergraph_type graph;
  
    for (path_set_type::const_iterator piter = forest_path.begin(); piter != forest_path.end(); ++ piter) {
    
      if (mpi_rank == 0 && debug)
	std::cerr << "reading forest: " << piter->string() << std::endl;

      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
      
	const path_type path_forest = (*piter) / file_name;
      
	if (! boost::filesystem::exists(path_forest)) break;
      
	utils::compress_istream is(path_forest);
	std::getline(is, line);
      
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end  = line.end();
      
	if (! parse_id(id_forest, iter, end))
	  throw std::runtime_error("invalid id input: " + path_forest.string());
	if (id_forest != i)
	  throw std::runtime_error("invalid id input: " + path_forest.string());
      
	if (id_forest >= graphs_forest.size())
	  graphs_forest.resize(id_forest + 1);
      
	if (! graph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path_forest.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path_forest.string());
      
	graphs_forest[id_forest].unite(graph);
      }
    }
  
    graphs_intersected.resize(graphs_forest.size());
  
    for (path_set_type::const_iterator piter = intersected_path.begin(); piter != intersected_path.end(); ++ piter) {
    
      if (mpi_rank == 0 && debug)
	std::cerr << "reading intersected forest: " << piter->string() << std::endl;

      for (size_t i = mpi_rank; i < graphs_intersected.size(); i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
      
	const path_type path_intersected = (*piter) / file_name;
      
	if (! boost::filesystem::exists(path_intersected)) continue;
      
	utils::compress_istream is(path_intersected);
	std::getline(is, line);
      
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end = line.end();
      
	if (! parse_id(id_intersected, iter, end))
	  throw std::runtime_error("invalid id input" + path_intersected.string());
	if (id_intersected != i)
	  throw std::runtime_error("invalid id input: " + path_intersected.string());
      
	if (! graph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path_intersected.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path_intersected.string());
	
	graphs_intersected[id_intersected].unite(graph);
      }
    }
  } else {
    if (forest_path.size() != intersected_path.size())
      throw std::runtime_error("# of forest does not match");
    
    size_t id_forest;
    size_t id_intersected;
    
    std::string line;
    
    for (size_t pos = 0; pos != forest_path.size(); ++ pos) {
      
      if (mpi_rank == 0 && debug)
	std::cerr << "reading forest: " << forest_path[pos].string() << " with " << intersected_path[pos].string() << std::endl;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_forest      = forest_path[pos] / file_name;
	const path_type path_intersected = intersected_path[pos] / file_name;
	
	if (! boost::filesystem::exists(path_forest)) break;
	if (! boost::filesystem::exists(path_intersected)) continue;
	
	{
	  utils::compress_istream is(path_forest);
	  std::getline(is, line);
	  
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator end = line.end();
	  
	  if (! parse_id(id_forest, iter, end))
	    throw std::runtime_error("invalid id input: " + path_forest.string());
	  if (id_forest != i)
	    throw std::runtime_error("invalid id input: " + path_forest.string());
	  
	  graphs_forest.push_back(hypergraph_type());
	  
	  if (! graphs_forest.back().assign(iter, end))
	    throw std::runtime_error("invalid graph format" + path_forest.string());
	  if (iter != end)
	    throw std::runtime_error("invalid id ||| graph format" + path_forest.string());
	}
	
	{
	  utils::compress_istream is(path_intersected);
	  std::getline(is, line);
	  
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator end = line.end();
	  
	  if (! parse_id(id_intersected, iter, end))
	    throw std::runtime_error("invalid id input" + path_intersected.string());
	  if (id_intersected != i)
	    throw std::runtime_error("invalid id input: " + path_intersected.string());
	  
	  graphs_intersected.push_back(hypergraph_type());
	  
	  if (! graphs_intersected.back().assign(iter, end))
	    throw std::runtime_error("invalid graph format" + path_intersected.string());
	  if (iter != end)
	    throw std::runtime_error("invalid id ||| graph format" + path_intersected.string());
	}
      }
    }
  }
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


void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("forest",      po::value<path_set_type>(&forest_path)->multitoken(),      "forest path(s)")
    ("intersected", po::value<path_set_type>(&intersected_path)->multitoken(), "intersected forest path(s)")
    ("oracle",      po::value<path_set_type>(&intersected_path)->multitoken(), "oracle forest path(s) (an alias for --intersected)")
    ("refset",      po::value<path_set_type>(&refset_path)->multitoken(),      "reference translation(s)")
    ("weights",     po::value<path_type>(&weights_path),      "initial parameter")
    ("weights-history", po::value<path_set_type>(&weights_history_path)->multitoken(), "parameter history")
    ("output",      po::value<path_type>(&output_path),       "output parameter")
    
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
    
    ("regularize-l1",      po::bool_switch(&regularize_l1),      "L1-regularization")
    ("regularize-l2",      po::bool_switch(&regularize_l2),      "L2-regularization")
    ("regularize-entropy", po::bool_switch(&regularize_entropy), " entropy regularization")
    
    ("C",             po::value<double>(&C)->default_value(C),         "regularization constant")
    ("C2",            po::value<double>(&C2)->default_value(C2),       "an alternative regularization constant")
    ("scale",         po::value<double>(&scale)->default_value(scale), "scaling for weight")
    ("eta0",          po::value<double>(&eta0),                        "\\eta_0 for decay")
    ("order",         po::value<int>(&order)->default_value(order),    "ngram order for xBLEU")
    
    ("temperature",       po::value<double>(&temperature)->default_value(temperature),             "temperature")
    ("temperature-start", po::value<double>(&temperature_start)->default_value(temperature_start), "start temperature for annealing")
    ("temperature-end",   po::value<double>(&temperature_end)->default_value(temperature_end),     "end temperature for annealing")
    ("temperature-rate",  po::value<double>(&temperature_rate)->default_value(temperature_rate),   "annealing rate")

    ("annealing", po::bool_switch(&annealing_mode), "annealing")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")
    
    ("unite",    po::bool_switch(&unite_forest), "unite forest sharing the same id")
    
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
