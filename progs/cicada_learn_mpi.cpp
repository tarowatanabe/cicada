//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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

#include "cicada/prune.hpp"
#include "cicada/operation/functional.hpp"

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

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random.hpp>

#include "lbfgs.h"

typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type forest_path;
path_set_type intersected_path;
path_type weights_path;
path_type output_path = "-";
path_type output_objective_path;

int iteration = 100;
bool learn_lbfgs = false;
bool learn_sgd = false;
bool learn_mira = false;
bool learn_arow = false;
bool learn_cw = false;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

bool unite_forest = false;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

template <typename Optimize>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
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
struct OptimizeLBFGS;

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
    
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,sgd,mira,arow,cw}");
    if (int(learn_lbfgs) + learn_sgd + learn_mira + learn_arow + learn_cw == 0)
      learn_lbfgs = true;

    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("either L1 or L2 regularization");

    if (forest_path.empty())
      throw std::runtime_error("no forest?");
    if (intersected_path.empty())
      throw std::runtime_error("no intersected forest?");

    hypergraph_set_type graphs_forest;
    hypergraph_set_type graphs_intersected;
        
    read_forest(forest_path, intersected_path, graphs_forest, graphs_intersected);

    if (debug && mpi_rank == 0)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    weight_set_type weights;
    if (mpi_rank ==0 && ! weights_path.empty()) {
      if (! boost::filesystem::exists(weights_path))
	throw std::runtime_error("no path? " + weights_path.string());
      
      utils::compress_istream is(weights_path, 1024 * 1024);
      is >> weights;
    }
    
    weights.allocate();

    double objective = 0.0;

    boost::mt19937 generator;
    generator.seed(time(0) * getpid());
    
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
    else
      objective = optimize_batch<OptimizeLBFGS>(graphs_forest, graphs_intersected, weights);

    if (debug && mpi_rank == 0)
      std::cerr << "objective: " << objective << std::endl;
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_path, 1024 * 1024);
      os.precision(20);
      os << weights;
      
      if (! outpub_objective_path.empty()) {
	utils::compress_ostream os(output_objective_xpath, 1024 * 1024);
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
      cicada::prune_kbest(hypergraph_forest, pruned_forest, function_type(optimizer.weights, 1.0), margin_kbest);
    else
      cicada::prune_beam(hypergraph_forest, pruned_forest, function_type(optimizer.weights, 1.0), margin_beam);
    
    if (margin_kbest > 0)
      cicada::prune_kbest(hypergraph_intersected, pruned_intersected, function_type(optimizer.weights, - 1.0), margin_kbest);
    else
      cicada::prune_beam(hypergraph_intersected, pruned_intersected, function_type(optimizer.weights, - 1.0), margin_beam);
    
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
  for (size_t id = 0; id != ids.size(); ++ id)
    ids[id] = id;
  
  int instances = graphs_forest.size();
  if (! unite_forest) {
    instances = 0;
    const int instances_local = graphs_forest.size();
    
    MPI::COMM_WORLD.Allreduce(&instances_local, &instances, 1, MPI::INT, MPI::SUM);
  }
  
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
      
      optimizer.weights *= optimizer.samples;
      reduce_weights(optimizer.weights);
      
      objective = 0.0;
      MPI::COMM_WORLD.Reduce(&optimizer.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      int samples = 0;
      int samples_local = optimizer.samples;
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
	
	optimizer.weights *= optimizer.samples;
	send_weights(optimizer.weights);
	
	double objective = 0.0;
	MPI::COMM_WORLD.Reduce(&optimizer.objective, &objective, 1, MPI::DOUBLE, MPI::SUM, 0);
	
	int samples = 0;
	int samples_local = optimizer.samples;
	MPI::COMM_WORLD.Reduce(&samples_local, &samples, 1, MPI::INT, MPI::SUM, 0);
      }
    }
    
    if (requests[NOTIFY].Test())
      requests[NOTIFY].Cancel();
    
    return 0.0;
  }
}
  

struct OptimizeLBFGS
{

  OptimizeLBFGS(const hypergraph_set_type& __graphs_forest,
		const hypergraph_set_type& __graphs_intersected,
		weight_set_type& __weights)
    : graphs_forest(__graphs_forest),
      graphs_intersected(__graphs_intersected),
      weights(__weights) {}

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
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_static_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    Task(const weight_set_type& __weights,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected)
      : weights(__weights),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected) {}
    
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
      const int mpi_rank = MPI::COMM_WORLD.Get_rank();
      const int mpi_size = MPI::COMM_WORLD.Get_size();

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
    }
    
    const weight_set_type& weights;

    const hypergraph_set_type& graphs_forest;
    const hypergraph_set_type& graphs_intersected;
    
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
    
    task_type task(optimizer.weights, optimizer.graphs_forest, optimizer.graphs_intersected);
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
    
  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
};

template <typename Optimize>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    const double objective = Optimize(graphs_forest, graphs_intersected, weights)();
    
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
	
	task_type task(weights, graphs_forest, graphs_intersected);
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
	std::string::const_iterator end = line.end();
      
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

  // collect features...
  for (int rank = 0; rank < mpi_size; ++ rank) {
    weight_set_type weights;
    weights.allocate();
    
    for (feature_type::id_type id = 0; id != feature_type::allocated(); ++ id)
      if (! feature_type(id).empty())
	weights[feature_type(id)] = 1.0;
    
    bcast_weights(rank, weights);
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
    device[rank].reset(new device_type(rank, weights_tag, 1024 * 1024));
    stream[rank].reset(new stream_type());
    
    stream[rank]->push(boost::iostreams::gzip_decompressor());
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
  os.push(boost::iostreams::gzip_compressor());
  os.push(utils::mpi_device_sink(0, weights_tag, 1024 * 1024));
  
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
    os.push(utils::mpi_device_bcast_sink(rank, 1024));
    
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
    is.push(utils::mpi_device_bcast_source(rank, 1024));
    
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
    ("forest",      po::value<path_set_type>(&forest_path)->multitoken(),       "forest path(s)")
    ("intersected", po::value<path_set_type>(&intersected_path)->multitoken(),  "intersected forest path(s)")
    ("weights",     po::value<path_type>(&weights_path),      "initial parameter")
    ("output",      po::value<path_type>(&output_path),       "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",  po::bool_switch(&learn_lbfgs),  "batch LBFGS algorithm")
    ("learn-sgd",    po::bool_switch(&learn_sgd),    "online SGD algorithm")
    ("learn-mira",   po::bool_switch(&learn_mira),   "online MIRA algorithm")
    ("learn-arow",   po::bool_switch(&learn_arow),   "online AROW algorithm")
    ("learn-cw",     po::bool_switch(&learn_cw),     "online CW algorithm")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C"            , po::value<double>(&C),           "regularization constant")

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
