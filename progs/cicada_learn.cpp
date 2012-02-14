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
#include <numeric>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"
#include "cicada_text_impl.hpp"

#include "cicada/operation/functional.hpp"
#include "cicada/expected_ngram.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
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
#include <boost/random.hpp>

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

bool learn_sgd = false;
bool learn_lbfgs = false;
bool learn_mira = false;
bool learn_xbleu = false;

bool regularize_l1 = false;
bool regularize_l2 = false;
bool regularize_entropy = false;
double C = 1.0;
double C2 = 1.0;
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
bool scale_fixed = false;

// scorers
std::string scorer_name = "bleu:order=4,exact=true";
bool scorer_list = false;

bool unite_forest = false;

int threads = 2;

int debug = 0;

#include "cicada_learn_impl.hpp"

void options(int argc, char** argv);

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

template <typename Optimizer>
double optimize_xbleu(const hypergraph_set_type& forests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights);
template <typename Optimizer>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights);
template <typename Optimizer, typename Generator>
double optimize_online(const hypergraph_set_type& graphs_forest,
		       const hypergraph_set_type& graphs_intersected,
		       weight_set_type& weights,
		       Generator& generator);
struct OptimizeLBFGS;

template <typename Optimizer, typename Generator>
struct OptimizeOnline;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(learn_lbfgs) + learn_sgd > 1)
      throw std::runtime_error("eitehr learn-{lbfgs,sgd}");
    if (int(learn_lbfgs) + learn_sgd == 0)
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
    
    threads = utils::bithack::max(1, threads);

    scorer_document_type scorers(scorer_name);
    if (! refset_path.empty())
      read_refset(refset_path, scorers);

    hypergraph_set_type graphs_forest;
    hypergraph_set_type graphs_intersected;
    
    if (! intersected_path.empty())
      read_forest(forest_path, intersected_path, graphs_forest, graphs_intersected);
    else {
      scorer_document_type scorers_forest(scorer_name);
      
      read_forest(forest_path, scorers, graphs_forest, scorers_forest);
      
      scorers_forest.swap(scorers);
    }
    
    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    weight_set_type weights;
    if (! weights_path.empty()) {
      if (! boost::filesystem::exists(weights_path))
	throw std::runtime_error("no path? " + weights_path.string());
      
      utils::compress_istream is(weights_path, 1024 * 1024);
      is >> weights;
    }
        
    weight_set_type bounds_lower;
    weight_set_type bounds_upper;
    
    if (! bound_lower_file.empty())
      read_bounds(bound_lower_file, bounds_lower, - std::numeric_limits<double>::infinity());
    
    if (! bound_upper_file.empty())
      read_bounds(bound_upper_file, bounds_upper,   std::numeric_limits<double>::infinity());
    
    weights.allocate();
    
    double objective = 0.0;
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());

    if (learn_sgd) {
      if (regularize_l1)
	objective = optimize_online<OptimizeOnline<OptimizerSGDL1, boost::mt19937> >(graphs_forest, graphs_intersected, weights, generator);
      else
	objective = optimize_online<OptimizeOnline<OptimizerSGDL2, boost::mt19937> >(graphs_forest, graphs_intersected, weights, generator);
    } else if (learn_xbleu)
      objective = optimize_xbleu<OptimizeXBLEU>(graphs_forest, scorers, weights);
    } else
      objective = optimize_batch<OptimizeLBFGS>(graphs_forest, graphs_intersected, weights);
    
    if (debug)
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
    
    utils::compress_ostream os(output_path, 1024 * 1024);
    os.precision(20);
    os << weights;
    
    if (! output_objective_path.empty()) {
      utils::compress_ostream os(output_objective_path, 1024 * 1024);
      os.precision(20);
      os << objective << '\n';
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Optimizer, typename Generator>
struct OptimizeOnline
{
  typedef Optimizer optimizer_type;
  typedef Generator generator_type;
  typedef std::vector<optimizer_type, std::allocator<optimizer_type> > optimizer_set_type;
  
  OptimizeOnline(const hypergraph_set_type& __graphs_forest,
		 const hypergraph_set_type& __graphs_intersected,
		 weight_set_type& __weights,
		 generator_type& __generator)
    : graphs_forest(__graphs_forest),
      graphs_intersected(__graphs_intersected),
      weights(__weights),
      generator(__generator) {}
  
  struct Task
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    typedef Optimizer optimizer_type;

    typedef typename optimizer_type::weight_type   weight_type;
    typedef typename optimizer_type::gradient_type gradient_type;    
    
    Task(queue_type& __queue,
	 optimizer_type& __optimizer,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected)
      : queue(__queue),
	optimizer(__optimizer),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected)
    {}

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
    
    void operator()()
    {
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;

      optimizer.initialize();
      
      int id;
      while (1) {
	queue.pop(id);
	if (id < 0) break;
	
	gradients.clear();
	gradients_intersected.clear();
	
	inside.clear();
	inside_intersected.clear();

	if (! graphs_forest[id].is_valid() || ! graphs_intersected[id].is_valid()) continue;
	
	inside.reserve(graphs_forest[id].nodes.size());
	inside.resize(graphs_forest[id].nodes.size(), weight_type());
	
	inside_intersected.reserve(graphs_intersected[id].nodes.size());
	inside_intersected.resize(graphs_intersected[id].nodes.size(), weight_type());
	
	cicada::inside_outside(graphs_forest[id], inside, gradients,
			       weight_function(optimizer.weights, optimizer.weight_scale),
			       feature_function(optimizer.weights, optimizer.weight_scale));
	
	cicada::inside_outside(graphs_intersected[id], inside_intersected, gradients_intersected,
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
      
      optimizer.finalize();
    }
    
    queue_type&     queue;
    optimizer_type& optimizer;

    const hypergraph_set_type& graphs_forest;
    const hypergraph_set_type& graphs_intersected;
  };
  
  double operator()()
  {
    typedef Task task_type;
    typedef typename task_type::queue_type queue_type;
    typedef std::vector<int, std::allocator<int> > id_set_type;

    id_set_type ids;
    const int id_max = utils::bithack::min(graphs_forest.size(), graphs_intersected.size());
    for (int id = 0; id != id_max; ++ id)
      if (graphs_forest[id].is_valid() && graphs_intersected[id].is_valid())
	ids.push_back(id);
    
    optimizer_set_type optimizers(threads, optimizer_type(ids.size(), C));
    
    queue_type queue;
    
    weight_set_type weights_mixed;
    double objective = 0.0;
    
    for (int iter = 0; iter < iteration; ++ iter) {
      
      boost::thread_group workers;
      for (int i = 0; i < threads; ++ i)
	workers.add_thread(new boost::thread(task_type(queue, optimizers[i], graphs_forest, graphs_intersected)));
      
      for (size_t pos = 0; pos != ids.size(); ++ pos)
	queue.push(ids[pos]);
      
      for (int i = 0; i < threads; ++ i)
	queue.push(-1);
      
      boost::random_number_generator<Generator> gen(generator);
      std::random_shuffle(ids.begin(), ids.end(), gen);
      
      workers.join_all();
      
      // collect weights from optimizers and perform averaging...
      weights_mixed.clear();
      size_t samples = 0;
      
      typename optimizer_set_type::iterator oiter_end = optimizers.end();
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
	oiter->weights *= (oiter->samples + 1);
	
	weights_mixed += oiter->weights;
	samples       += (oiter->samples + 1);
      }
      
      weights_mixed *= (1.0 / samples);
      
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter)
	oiter->weights = weights_mixed;
      
      objective = 0.0;
      for (typename optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter)
	objective += oiter->objective;
      if (debug >= 2)
	std::cerr << "objective: " << objective << std::endl;
    }
    
    weights.swap(weights_mixed);
    
    return objective;
  }

  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
  generator_type& generator;
};

template <typename Optimizer, typename Generator>
double optimize_online(const hypergraph_set_type& graphs_forest,
		       const hypergraph_set_type& graphs_intersected,
		       weight_set_type& weights,
		       Generator& generator)
{
  return Optimizer(graphs_forest, graphs_intersected, weights, generator)();
}

struct OptimizeXBLEU
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  OptimizeXBLEU(const hypergraph_set_type& __forests,
		const scorer_document_type& __scorers,
		weight_set_type& __weights,
		const double& __lambda,
		const feature_type& __feature_scale)
    : forests(__forests),
      scorers(__scorers),
      weights(__weights),
      lambda(__lambda),
      feature_scale(__feature_scale) {}
  
  const hypergraph_set_type& forests;
  const scorer_document_type& scorers;
  weight_set_type& weights;
  
  double lambda;
  const feature_type& feature_scale;

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

    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef std::vector<gradient_type, std::allocator<gradient_type> > gradients_type;
    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    
    typedef cicada::Symbol       word_type;
    typedef cicada::SymbolVector ngram_type;

    typedef utils::indexed_trie<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> > index_set_type;
    typedef utils::simple_vector<index_set_type::id_type, std::allocator<index_set_type::id_type> > id_set_type;
    typedef std::vector<id_set_type, std::allocator<id_set_type> > id_map_type;

    // queue...
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    
    struct Count
    {
      weight_type c;
      weight_type mu_prime;
      
      Count() : c(), mu_prime() {}
    };
    typedef Count count_type;
    typedef std::vector<count_type, std::allocator<count_type> > count_set_type;
    typedef std::vector<ngram_type, std::allocator<ngram_type> > ngram_set_type;
    

    typedef std::vector<double, std::allocator<double> > ngram_counts_type;
    typedef std::vector<weight_set_type, std::allocator<weight_set_type> > feature_counts_type;
    
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
		    const weight_set_type& __weights,
		    const double& __scale)
	: ngrams(__ngrams), counts(__counts), ids(__ids),
	  weights(__weights), scale(__scale) {}
      
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
      const double           scale;
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
	  const weight_type weight = cicada::semiring::traits<weight_type>::exp(margin * impl.scale);
	  const weight_type value_scale(margin);
	  
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
		  const weight_type value(fiter->second * impl.scale);
		  
		  impl.dM[n][fiter->first] += value * scale_matched;
		  impl.dH[n][fiter->first] += value * scale_hypo;
		}
	      
	      impl.dM[n][impl.feature_scale] += value_scale * scale_matched;
	      impl.dH[n][impl.feature_scale] += value_scale * scale_hypo;
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
			 const weight_set_type& __weights,
			 const double& __scale,
			 const feature_type& __feature_scale) 
	: ngrams(__ngrams), counts(__counts), ids(__ids),
	  matched(__matched), hypo(__hypo),
	  weights(__weights), scale(__scale), feature_scale(__feature_scale),
	  dM(order + 1),
	  dH(order + 1) {}
      
      
      const ngram_set_type&  ngrams;
      const count_set_type&  counts;
      const id_map_type&     ids;
      
      const weights_type&    matched;
      const weights_type&    hypo;
      
      const weight_set_type& weights;
      const double           scale;
      const feature_type     feature_scale;
      
      accumulated_set_type dM;
      accumulated_set_type dH;
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
	value_type(const feature_set_type& __features, const weight_set_type& __weights, const double& __scale, const feature_type& __feature_scale)
	  : features(__features), weights(__weights), scale(__scale), feature_scale(__feature_scale) {}
	
	friend
	value_type operator*(value_type x, const entropy_weight_type& weight)
	{
	  x.inside_outside = weight;
	  return x;
	}
	
	entropy_weight_type inside_outside;
	
	const feature_set_type& features;
	const weight_set_type& weights;
	const double scale;
	const feature_type& feature_scale;
      };
      
      entropy_gradient_function(const weight_set_type& __weights, const double& __scale, const feature_type& __feature_scale)
	: weights(__weights), scale(__scale), feature_scale(__feature_scale) {}
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	return value_type(edge.features, weights, scale, feature_scale);
      }
      
      const weight_set_type& weights;
      const double scale;
      const feature_type& feature_scale;
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
	  const weight_type value_scale(value);
	  
	  // dZ += \lnabla p_e * x.inside_outside.p;
	  // dR += (1 + \log p_e) * \nalba p_e * x.inside_outside.p + \lnabla p_e * x.inside_outside.r;
	  
	  feature_set_type::const_iterator fiter_end = x.features.end();
	  for (feature_set_type::const_iterator fiter = x.features.begin(); fiter != fiter_end; ++ fiter) 
	    if (fiter->second != 0.0) {
	      const weight_type value(fiter->second * x.scale);
	      
	      dZ[fiter->first] += value * p_e * x.inside_outside.p;
	      dR[fiter->first] += (weight_type(1.0 + log_p_e) * value * p_e * x.inside_outside.p + value * p_e * x.inside_outside.r);
	    }
	  
	  dZ[x.feature_scale] += value_scale * p_e * x.inside_outside.p;
	  dR[x.feature_scale] += (weight_type(1.0 + log_p_e) * value_scale * p_e * x.inside_outside.p + value_scale * p_e * x.inside_outside.r);
	  
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
    

    Task(queue_type& __queue,
	 const hypergraph_set_type& __forests,
	 const scorer_document_type& __scorers,
	 const weight_set_type& __weights,
	 const feature_type& __feature_scale)
      : queue(__queue),
	forests(__forests), scorers(__scorers), weights(__weights), feature_scale(__feature_scale),
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
      id_map_type    ids;
      
      weights_type   matched(order + 1);
      weights_type   hypo(order + 1);
      
      weights_type   counts_matched(order + 1);
      weights_type   counts_hypo(order + 1);
      gradients_type gradients_matched(order + 1);
      gradients_type gradients_hypo(order + 1);

      bleu_weights_type bleu_inside;
      
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

      const double scale = weights[feature_scale];
      
      for (;;) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
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
	ids.clear();
	
	ids.resize(forest.edges.size());
		
	cicada::expected_ngram(forest,
			       cicada::operation::weight_scaled_function<weight_type>(weights, scale),
			       CollectCounts(index, ngrams, counts, ids),
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
	bleu_inside.clear();
	bleu_inside.resize(forest.nodes.size(), bleu_weight_type());
	
	bleu_gradient_type bleu_gradient(ngrams, counts, ids,
					 matched, hypo,
					 weights, scale, feature_scale);
	
	cicada::inside_outside(forest,
			       bleu_inside,
			       bleu_gradient,
			       bleu_function(ngrams, counts, ids, weights, scale),
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

	
	// forth, compute entorpy...
	entropy_inside.clear();
	entropy_inside.resize(forest.nodes.size(), entropy_weight_type());

	entropy_gradient_type entropy_gradient;
	
	cicada::inside_outside(forest,
			       entropy_inside,
			       entropy_gradient,
			       entropy_function(weights, scale),
			       entropy_gradient_function(weights, scale, feature_scale));
	
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
    
    queue_type& queue;
    
    const hypergraph_set_type& forests;
    const scorer_document_type& scorers;
    const weight_set_type& weights;
    const feature_type& feature_scale;
    
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
    typedef task_type::queue_type queue_type;
    
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
    
    OptimizeXBLEU& optimizer = *((OptimizeXBLEU*) instance);
    
    queue_type queue;
    task_set_type tasks(threads, task_type(queue,
					   optimizer.forests,
					   optimizer.scorers,
					   optimizer.weights,
					   optimizer.feature_scale));
    
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    size_type instances = 0;
    for (size_t id = 0; id != optimizer.forests.size(); ++ id)
      if (optimizer.forests[id].is_valid()) {
	queue.push(id);
	++ instances;
      }
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    // clear g...
    std::fill(g, g + size, 0.0);
    
    workers.join_all();

    
    task_type::ngram_counts_type c_matched(order + 1, 0.0);
    task_type::ngram_counts_type c_hypo(order + 1, 0.0);
    
    task_type::feature_counts_type g_matched(order + 1);
    task_type::feature_counts_type g_hypo(order + 1);
    weight_set_type g_entropy;

    double r(0.0);
    double e(0.0);
    
    
    for (int i = 0; i < threads; ++ i) {
      std::transform(tasks[i].c_matched.begin(), tasks[i].c_matched.end(), c_matched.begin(), c_matched.begin(), std::plus<double>());
      std::transform(tasks[i].c_hypo.begin(), tasks[i].c_hypo.end(), c_hypo.begin(), c_hypo.begin(), std::plus<double>());
      
      for (int n = 1; n <= order; ++ n) {
	g_matched[n] += tasks[i].g_matched[n];
	g_hypo[n] += tasks[i].g_hypo[n];
      }
      
      g_entropy += tasks[i].g_entropy;
      
      r += tasks[i].r;
      e += tasks[i].e;
    }

    
    // smoothing...
    {
      double smoothing = 1e-40;
      for (int n = 1; n <= order; ++ n) {
	if (c_hypo[n] > 0.0 && c_matched[n] <= 0.0)
	  c_matched[n] = smoothing;
	smoothing *= 0.1;
      }
    }
    
    // compute P
    double P = 0.0;
    for (int n = 1; n <= order; ++ n)
      if (c_hypo[n] > 0.0)
	P += (1.0 / order) * (utils::mathop::log(c_matched[n]) - utils::mathop::log(c_hypo[n]));
    
    // compute C and B
    const double C = r / c_hypo[1];
    const double B = task_type::brevity_penalty(1.0 - C);
    
    // for computing g...
    const double exp_P = utils::mathop::exp(P);
    const double C_dC  = C * task_type::derivative_brevity_penalty(1.0 - C);
    
    //std::cerr << "P: " << P << " B: " << B << " C: " << C << std::endl;

    // xBLEU...
    const double objective_bleu = exp_P * B;
    const double entropy = e / instances;
    
    // compute g..
    std::fill(g, g + size, 0.0);

    // entropy
    if (regularize_entropy) {
      // 0.5 * (entropy - C2)^2
      // 
      // thus, derivative is: (entropy - C2) * \nabla entropy
      // we need to consider average of entropy...
      
      std::transform(g_entropy.begin(), g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), (entropy - C2) / instances));
    } else
      std::transform(g_entropy.begin(), g_entropy.end(), g, std::bind2nd(std::multiplies<double>(), - temperature / instances));
    
    for (int n = 1; n <= order; ++ n) 
      if (c_hypo[n] > 0.0) {
	const double factor_matched = - (exp_P * B / order) / c_matched[n];
	const double factor_hypo    = - (exp_P * B / order) / c_hypo[n];
	
	for (size_t i = 0; i != static_cast<size_t>(size); ++ i) {
	  g[i] += factor_matched * g_matched[n][i];
	  g[i] -= factor_hypo * g_hypo[n][i];
	}
      }
    
    if (c_hypo[1] > 0.0) {
      // I think the missed exp(P) is a bug in Rosti et al. (2011)
      const double factor = - exp_P * C_dC / c_hypo[1];
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i)
	g[i] += factor * g_hypo[1][i];
    }
    
    if (debug >= 3) {
      std::cerr << "grad:" << std::endl;
      for (size_t i = 0; i != static_cast<size_t>(size); ++ i)
	if (g[i] != 0.0 && feature_type(i) != feature_type())
	  std::cerr << feature_type(i) << ' ' << g[i] << std::endl;
    }
        
    // we need to minimize negative bleu... + regularized by average entropy...
    double objective = - objective_bleu + (regularize_entropy ? 0.5 * (entropy - C2) * (entropy - C2) : - temperature * entropy);
    
    if (regularize_l2) {
      double norm = 0.0;
      for (size_t i = 0; i < static_cast<size_t>(size); ++ i) {
	g[i] += optimizer.lambda * x[i] * double(i != optimizer.feature_scale.id());
	norm += x[i] * x[i] * double(i != optimizer.feature_scale.id());
      }
      objective += 0.5 * optimizer.lambda * norm;
    }
    
    if (scale_fixed)
      g[optimizer.feature_scale.id()] = 0.0;
    
    if (debug >= 2)
      std::cerr << "objective: " << objective
		<< " xBLEU: " << objective_bleu
		<< " BP: " << B
		<< " entropy: " << entropy
		<< " scale: " << optimizer.weights[optimizer.feature_scale]
		<< std::endl;

    // keep the best so forth...
    if (objective <= optimizer.objective_opt) {
      optimizer.objective_opt = objective;
      optimizer.weights_opt = optimizer.weights;
    }
    
    return objective;
  }
};


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
    
    const int result = lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    if (debug)
      std::cerr << "lbfgs: " << lbfgs_error(result) << std::endl;
    
    return objective;
  }

  
  struct Task
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
    
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > gradient_static_type;
    
    typedef std::vector<weight_type, std::allocator<weight_type> > weights_type;

    Task(queue_type&            __queue,
	 const weight_set_type& __weights,
	 const hypergraph_set_type& __graphs_forest,
	 const hypergraph_set_type& __graphs_intersected,
	 const size_t& __instances)
      : queue(__queue),
	weights(__weights),
	graphs_forest(__graphs_forest),
	graphs_intersected(__graphs_intersected),
	instances(__instances)
    {}

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

    void operator()()
    {
      gradients_type gradients;
      gradients_type gradients_intersected;
      weights_type   inside;
      weights_type   inside_intersected;
      
      gradient_static_type  feature_expectations;

      g.clear();
      objective = 0.0;
      
      while (1) {
	int id = 0;
	queue.pop(id);
	if (id < 0) break;
	
	gradients.clear();
	gradients_intersected.clear();
	
	inside.clear();
	inside_intersected.clear();
	
	inside.reserve(graphs_forest[id].nodes.size());
	inside.resize(graphs_forest[id].nodes.size(), weight_type());
	
	inside_intersected.reserve(graphs_intersected[id].nodes.size());
	inside_intersected.resize(graphs_intersected[id].nodes.size(), weight_type());
	
	cicada::inside_outside(graphs_forest[id], inside, gradients, weight_function(weights), feature_function(weights));
	cicada::inside_outside(graphs_intersected[id], inside_intersected, gradients_intersected, weight_function(weights), feature_function(weights));
	
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
	  std::cerr << "id: " << id << " margin: " << margin << std::endl;
      }
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(feature_expectations.begin(), feature_expectations.end(), g.begin());
      
      // normalize!
      objective /= instances;
      std::transform(g.begin(), g.end(), g.begin(), std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    }

    queue_type&            queue;
    
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
    typedef task_type::queue_type queue_type;
    
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

    OptimizeLBFGS& optimizer = *((OptimizeLBFGS*) instance);
    
    const int id_max = utils::bithack::min(optimizer.graphs_forest.size(), optimizer.graphs_intersected.size());
    size_t instances = 0;
    for (int id = 0; id != id_max; ++ id)
      instances += (optimizer.graphs_forest[id].is_valid() && optimizer.graphs_intersected[id].is_valid());
    
    queue_type queue;
    
    task_set_type tasks(threads, task_type(queue, optimizer.weights, optimizer.graphs_forest, optimizer.graphs_intersected, instances));

    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    for (int id = 0; id != id_max; ++ id)
      if (optimizer.graphs_forest[id].is_valid() && optimizer.graphs_intersected[id].is_valid())
	queue.push(id);
    
    // collect all the objective and gradients...
    double objective = 0.0;
    std::fill(g, g + n, 0.0);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    workers.join_all();
    
    for (int i = 0; i < threads; ++ i) {
      objective += tasks[i].objective;
      std::transform(tasks[i].g.begin(), tasks[i].g.end(), g, g, std::plus<double>());
    }
        
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
      std::cerr << "objective: " << objective << std::endl;
    
    return objective;
  }
  
  const hypergraph_set_type& graphs_forest;
  const hypergraph_set_type& graphs_intersected;
  
  weight_set_type& weights;
};

template <typename Optimizer>
double optimize_xbleu(const hypergraph_set_type& forests,
		      const scorer_document_type& scorers,
		      weight_set_type& weights)
{
  const feature_type feature_scale(":feature-scale:");
  
  weights[feature_scale] = scale;
  
  Optimizer optimizer(forests, scorers, weights, C, feature_scale);
  
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
    
  return objective;
}


template <typename Optimizer>
double optimize_batch(const hypergraph_set_type& graphs_forest,
		      const hypergraph_set_type& graphs_intersected,
		      weight_set_type& weights)
{
  return Optimizer(graphs_forest, graphs_intersected, weights)();
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
  if (unite_forest) {
    size_t id;
    
    std::string line;
    hypergraph_type graph;
    
    for (path_set_type::const_iterator piter = forest_path.begin(); piter != forest_path.end(); ++ piter) {
    
      if (debug)
	std::cerr << "reading forest: " << piter->string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
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
      if (debug)
	std::cerr << "reading forest: " << forest_path[pos].string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
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
  if (unite_forest) {
    size_t id_forest;
    size_t id_intersected;
    
    std::string line;
    hypergraph_type graph;
    
    for (path_set_type::const_iterator piter = forest_path.begin(); piter != forest_path.end(); ++ piter) {
      if (debug)
	std::cerr << "reading forest: " << piter->string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
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
    
      if (debug)
	std::cerr << "reading intersected forest: " << piter->string() << std::endl;

      for (size_t i = 0; i < graphs_intersected.size(); ++ i) {
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
      if (debug)
	std::cerr << "reading forest: " << forest_path[pos].string() << " with " << intersected_path[pos].string() << std::endl;
      
      for (size_t i = 0; /**/; ++ i) {
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


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("forest",      po::value<path_set_type>(&forest_path)->multitoken(),      "forest path")
    ("intersected", po::value<path_set_type>(&intersected_path)->multitoken(), "intersected forest path")
    ("oracle",      po::value<path_set_type>(&intersected_path)->multitoken(), "oracle forest path(s) (an alias for --intersected)")
    ("refset",      po::value<path_set_type>(&refset_path)->multitoken(),      "reference translation(s)")
    ("weights",     po::value<path_type>(&weights_path),      "initial parameter")
    ("weights-history", po::value<path_set_type>(&weights_history_path)->multitoken(), "parameter history")
    ("output",      po::value<path_type>(&output_path),       "output parameter")
    
    ("output-objective", po::value<path_type>(&output_objective_path), "output final objective")
    
    ("bound-lower", po::value<path_type>(&bound_lower_file),                     "lower bounds definition for feature weights")
    ("bound-upper", po::value<path_type>(&bound_upper_file),                     "upper bounds definition for feature weights")

    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max # of iterations")
    
    ("learn-lbfgs",  po::bool_switch(&learn_lbfgs),  "batch LBFGS algorithm")
    ("learn-sgd",    po::bool_switch(&learn_sgd),    "online SGD algorithm")
    ("learn-xbleu",  po::bool_switch(&learn_xbleu),  "xBLEU algorithm")
    
    ("regularize-l1",      po::bool_switch(&regularize_l1),      "L1-regularization")
    ("regularize-l2",      po::bool_switch(&regularize_l2),      "L2-regularization")
    ("regularize-entropy", po::bool_switch(&regularize_entropy), " entropy regularization")
    
    ("C",             po::value<double>(&C)->default_value(C),         "regularization constant")
    ("C2",            po::value<double>(&C2)->default_value(C2),       "an alternative regularization constant")
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

    ("scale-fixed", po::bool_switch(&scale_fixed), "fixed scaling")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")
    
    ("unite",    po::bool_switch(&unite_forest), "unite forest sharing the same id")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  desc_command.add(opts_command);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options] [operations]\n"
	      << opts_command << std::endl;
    exit(0);
  }
}
