
//
// refset format:
// 0 |||  reference translatin for source sentence 0
// 0 |||  another reference
// 1 |||  reference translation for source sentence 1
//

#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <stdexcept>
#include <numeric>
#include <algorithm>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"
#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"
#include "cicada/viterbi.hpp"
#include "cicada/sentence_vector.hpp"

#include "cicada/apply.hpp"
#include "cicada/model.hpp"

#include "cicada/feature/bleu.hpp"
#include "cicada/feature/bleu_linear.hpp"
#include "cicada/parameter.hpp"

#include "cicada/eval.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"
#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/thread.hpp>

#include "lbfgs.h"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Symbol   symbol_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef cicada::Model model_type;

typedef rule_type::feature_set_type    feature_set_type;
typedef cicada::WeightVector<double>   weight_set_type;
typedef feature_set_type::feature_type feature_type;

typedef cicada::FeatureFunction feature_function_type;
typedef feature_function_type::feature_function_ptr_type feature_function_ptr_type;

typedef std::vector<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;


typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type> > feature_function_ptr_set_type;

typedef cicada::SentenceVector sentence_set_type;
typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type  score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;


path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";

std::string scorer_name = "bleu:order=4,exact=true";
int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

bool oracle_loss = false;
bool apply_exact = false;
int cube_size = 200;
bool softmax_margin = false;

int threads = 4;

int debug = 0;

void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 feature_function_ptr_set_type& features);
void read_refset(const path_set_type& file,
		 scorer_document_type& scorers,
		 sentence_document_type& sentences);
void compute_oracles(const hypergraph_set_type& graphs,
		     const feature_function_ptr_set_type& features,
		     const scorer_document_type& scorers);
double optimize(const hypergraph_set_type& graphs,
		const feature_function_ptr_set_type& features,
		weight_set_type& weights);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("you cannot use both of L1 and L2...");
    
    threads = utils::bithack::max(threads, 1);
    
    // read reference set
    scorer_document_type   scorers(scorer_name);
    sentence_document_type sentences;
    
    read_refset(refset_files, scorers, sentences);
    
    if (debug)
      std::cerr << "# of references: " << sentences.size() << std::endl;

    // read test set
    if (tstset_files.empty())
      tstset_files.push_back("-");

    if (debug)
      std::cerr << "reading hypergraphs" << std::endl;

    hypergraph_set_type       graphs(sentences.size());
    feature_function_ptr_set_type features(sentences.size());
    
    read_tstset(tstset_files, graphs, sentences, features);

    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    if (oracle_loss)
      compute_oracles(graphs, features, scorers);
    
    weight_set_type weights;
    
    const double objective = optimize(graphs, features, weights);
    
    if (debug)
      std::cerr << "objective: " << objective << std::endl;
    
    utils::compress_ostream os(output_file);
    os.precision(10);
    os << weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct OptimizeLBFGS
{

  OptimizeLBFGS(const hypergraph_set_type&           __graphs,
		const feature_function_ptr_set_type& __features,
		weight_set_type&                      __weights)
    : graphs(__graphs),
      features(__features),
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
    weights.clear();
    weights.allocate();
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, OptimizeLBFGS::evaluate, 0, this, &param);
    
    return objective;
  }

  
  struct Task
  {
    typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
        
    Task(queue_type&            __queue,
	 const hypergraph_set_type&           __graphs,
	 const feature_function_ptr_set_type& __features,
	 const weight_set_type&               __weights)
      : queue(__queue),
	graphs(__graphs),
	features(__features),
	weights(__weights) {}
    
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > gradient_type;
    typedef cicada::WeightVector<weight_type, std::allocator<weight_type> >  expectation_type;
    
    typedef std::vector<weight_type, std::allocator<weight_type> > inside_set_type;

    typedef cicada::semiring::Tropical<double> bleu_weight_type;
    typedef std::vector<bleu_weight_type, std::allocator<bleu_weight_type> > bleu_weight_set_type;
        
    struct gradient_set_type
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

    struct bleu_weight_function
    {
      typedef cicada::semiring::Tropical<double> value_type;
      
      bleu_weight_function(const weight_set_type& __weights) : weights(__weights) {}

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
      }
      
      const weight_set_type& weights;
    };

    struct weight_set_function
    {
      typedef cicada::semiring::Logprob<double> value_type;

      weight_set_function(const weight_set_type& __weights)
	: weights(__weights) {}

      const weight_set_type& weights;
      
      value_type operator()(const feature_set_type& x) const
      {
	return cicada::semiring::traits<value_type>::log(x.dot(weights));
      }
    };

    struct weight_function
    {
      typedef weight_type value_type;

      weight_function(const weight_set_type& __weights)
	: weights(__weights) {}

      const weight_set_type& weights;

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e
	return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
      }
    };

    struct weight_max_function
    {
      typedef weight_type value_type;

      weight_max_function(const weight_set_type& __weights, const bleu_weight_set_type& __bleus, const bleu_weight_type& __max_bleu)
	: weights(__weights), bleus(__bleus), max_bleu(__max_bleu) {}

      const weight_set_type&      weights;
      const bleu_weight_set_type& bleus;
      const bleu_weight_type      max_bleu;

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e
	if (log(max_bleu) - log(bleus[edge.id]) >= 1e-4)
	  return value_type();
	else
	  return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
      }
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
	
	const weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights));
	
	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0)
	    grad[fiter->first] = weight_type(fiter->second) * weight;
	
	return grad;
      }
      
      const weight_set_type& weights;
    };

    struct feature_max_function
    {
      typedef gradient_type value_type;

      feature_max_function(const weight_set_type& __weights, const bleu_weight_set_type& __bleus, const bleu_weight_type& __max_bleu)
	: weights(__weights), bleus(__bleus), max_bleu(__max_bleu) {}

      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	// p_e r_e
	if (log(max_bleu) - log(bleus[edge.id]) >= 1e-4)
	  return gradient_type();

	gradient_type grad;
	
	const weight_type weight = cicada::semiring::traits<weight_type>::log(edge.features.dot(weights));
	
	feature_set_type::const_iterator fiter_end = edge.features.end();
	for (feature_set_type::const_iterator fiter = edge.features.begin(); fiter != fiter_end; ++ fiter)
	  if (fiter->second != 0.0)
	    grad[fiter->first] = weight_type(fiter->second) * weight;
	
	return grad;
      }
      
      const weight_set_type&      weights;
      const bleu_weight_set_type& bleus;
      const bleu_weight_type      max_bleu;
    };


    void operator()()
    {
      g.clear();
      objective = 0.0;

      weight_set_type weights_reward(weights);
      weight_set_type weights_penalty(weights);
      weight_set_type weights_bleu;

      weight_set_type::feature_type feature_bleu;
      for (int i = 0; i < features.size(); ++ i)
	if (features[i]) {
	  feature_bleu = features[i]->feature_name();
	  break;
	}
      
      weights_reward[feature_bleu]  =  1.0;
      weights_penalty[feature_bleu] = -1.0;
      weights_bleu[feature_bleu]    =  1.0;
      
      hypergraph_type graph_reward;
      hypergraph_type graph_penalty;

      bleu_weight_set_type bleus_inside;
      bleu_weight_set_type bleus_inside_outside;

      inside_set_type   inside;
      gradient_set_type gradients;

      inside_set_type   inside_correct;
      gradient_set_type gradients_correct;

      expectation_type  feature_expectations;
      
      int id = 0;
      while (1) {
	queue.pop(id);
	if (id < 0) break;
	
	if (! graphs[id].is_valid()) continue;
	
	model_type model;
	model.push_back(features[id]);
	
	if (! apply_exact) {
	  cicada::apply_cube_prune(model, graphs[id], graph_reward, weight_set_function(weights_reward), cube_size);
	  cicada::apply_cube_prune(model, graphs[id], graph_penalty, weight_set_function(weights_penalty), cube_size);
	  
	  graph_reward.unite(graph_penalty);
	  graph_penalty.clear();
	}
	
	const hypergraph_type& graph = (apply_exact ? graphs[id] : graph_reward);
	
	// compute inside/outside by bleu using tropical semiring...
	bleus_inside.clear();
	bleus_inside_outside.clear();
	
	bleus_inside.resize(graph.nodes.size());
	bleus_inside_outside.resize(graph.edges.size());
	
	cicada::inside_outside(graph, bleus_inside, bleus_inside_outside, bleu_weight_function(weights_bleu), bleu_weight_function(weights_bleu));
	
	const bleu_weight_type bleu_max = *std::max_element(bleus_inside_outside.begin(), bleus_inside_outside.end());
	
	// then, inside/outside to collect potentials...
	
	inside.clear();
	inside_correct.clear();
	
	gradients.clear();
	gradients_correct.clear();
	
	inside.resize(graph.nodes.size());
	inside_correct.resize(graph.nodes.size());
	
	if (softmax_margin) {
	  cicada::inside_outside(graph, inside, gradients, weight_function(weights_penalty), feature_function(weights_penalty));
	  
	  cicada::inside_outside(graph, inside_correct, gradients_correct,
				 weight_max_function(weights_penalty, bleus_inside_outside, bleu_max),
				 feature_max_function(weights_penalty, bleus_inside_outside, bleu_max));
	} else {
	  cicada::inside_outside(graph, inside, gradients, weight_function(weights), feature_function(weights));
	  
	  cicada::inside_outside(graph, inside_correct, gradients_correct,
				 weight_max_function(weights, bleus_inside_outside, bleu_max),
				 feature_max_function(weights, bleus_inside_outside, bleu_max));
	}
	
	
	gradient_type& gradient = gradients.gradient;
	weight_type& Z = inside.back();
	
	gradient_type& gradient_correct = gradients_correct.gradient;
	weight_type& Z_correct = inside_correct.back();
	
	gradient /= Z;
	gradient_correct /= Z_correct;
	
	feature_expectations -= gradient_correct;
	feature_expectations += gradient;
	
	const double margin = log(Z_correct) - log(Z);
	
	objective -= margin;
	
	if (debug >= 3)
	  std::cerr << "id: " << id
		    << " bleu: " << log(bleu_max)
		    << " correct: " << log(Z_correct)
		    << " partition: " << log(Z)
		    << " margin: " << margin
		    << std::endl;
      }
      
      // transform feature_expectations into g...
      g.allocate();
      
      std::copy(feature_expectations.begin(), feature_expectations.end(), g.begin());

      g[feature_bleu] = 0.0;
    }
    
    queue_type&            queue;
    
    const hypergraph_set_type&           graphs;
    const feature_function_ptr_set_type& features;
    const weight_set_type&               weights;
    
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
    
    typedef boost::shared_ptr<task_type> task_ptr_type;
    typedef std::vector<task_ptr_type, std::allocator<task_ptr_type> > task_set_type;

    OptimizeLBFGS& optimizer = *((OptimizeLBFGS*) instance);
    
    weight_set_type::feature_type feature_bleu;
    for (int i = 0; i < optimizer.features.size(); ++ i)
      if (optimizer.features[i]) {
	feature_bleu = optimizer.features[i]->feature_name();
	break;
      }
    optimizer.weights[feature_bleu] = 0.0;
        
    queue_type queue(optimizer.graphs.size());
    
    task_set_type tasks(threads);
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i) {
      tasks[i].reset(new task_type(queue, optimizer.graphs, optimizer.features, optimizer.weights));
      workers.add_thread(new boost::thread(boost::ref(*tasks[i])));
    }
    
    for (int sample = 0; sample < optimizer.graphs.size(); ++ sample)
      queue.push(sample);
    
    // collect all the objective and gradients...
    double objective = 0.0;
    std::fill(g, g + n, 0.0);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    workers.join_all();
    
    for (int i = 0; i < threads; ++ i) {
      objective += tasks[i]->objective;
      std::transform(tasks[i]->g.begin(), tasks[i]->g.end(), g, g, std::plus<double>());
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
  
  const hypergraph_set_type&           graphs;
  const feature_function_ptr_set_type& features;
  weight_set_type&                     weights;
};

double optimize(const hypergraph_set_type& graphs,
		const feature_function_ptr_set_type& features,
		weight_set_type& weights)
{
  return OptimizeLBFGS(graphs, features, weights)();
}


struct TaskOracle
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  TaskOracle(queue_type&                          __queue,
	     const hypergraph_set_type&           __graphs,
	     const feature_function_ptr_set_type& __features,
	     const scorer_document_type&          __scorers,
	     score_ptr_set_type&                  __scores)
    : queue(__queue),
      graphs(__graphs),
      features(__features),
      scorers(__scorers),
      scores(__scores)
  {
    score_optimum.reset();
    
    score_ptr_set_type::const_iterator siter_end = scores.end();
    for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter) 
      if (*siter) {
	if (! score_optimum)
	  score_optimum = (*siter)->clone();
	else
	  *score_optimum += *(*siter);
      } 
  }
  
  
  struct bleu_function
  {
    typedef rule_type::feature_set_type feature_set_type;
    
    typedef cicada::semiring::Logprob<double> value_type;

    bleu_function(const weight_set_type::feature_type& __name, const double& __factor)
      : name(__name), factor(__factor) {}

    weight_set_type::feature_type name;
    double factor;
    
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      return cicada::semiring::traits<value_type>::log(edge.features[name] * factor);
    }
    
  };
  
  struct kbest_traversal
  {
    typedef sentence_type value_type;
    
    template <typename Edge, typename Iterator>
    void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
    {
      // extract target-yield, features
      
      yield.clear();
    
      rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
      for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
	if (titer->is_non_terminal()) {
	  const int pos = titer->non_terminal_index() - 1;
	  yield.insert(yield.end(), (first + pos)->begin(), (first + pos)->end());
	} else if (*titer != vocab_type::EPSILON)
	  yield.push_back(*titer);
    }
  };

  struct weight_bleu_function
  {
    typedef cicada::semiring::Logprob<double> value_type;
    
    weight_bleu_function(const weight_set_type::feature_type& __name, const double& __factor)
      : name(__name), factor(__factor) {}
    
    weight_set_type::feature_type name;
    double factor;
    
    value_type operator()(const feature_set_type& x) const
    {
      return cicada::semiring::traits<value_type>::log(x[name] * factor);
    }
  };
  
  void operator()()
  {
    // we will try maximize    
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);
    
    weight_set_type::feature_type feature_bleu;
    for (int i = 0; i < features.size(); ++ i)
      if (features[i]) {
	feature_bleu = features[i]->feature_name();
	break;
      }
    
    double objective_optimum = (score_optimum
				? score_optimum->score().first * score_factor
				: - std::numeric_limits<double>::infinity());

    hypergraph_type graph_oracle;
    
    int id = 0;
    while (1) {
      queue.pop(id);
      if (id < 0) break;
      
      if (! graphs[id].is_valid()) continue;

      score_ptr_type score_curr;
      if (score_optimum)
	score_curr = score_optimum->clone();
      
      if (scores[id])
	*score_curr -= *scores[id];
      
      cicada::feature::Bleu*       __bleu = dynamic_cast<cicada::feature::Bleu*>(features[id].get());
      cicada::feature::BleuLinear* __bleu_linear = dynamic_cast<cicada::feature::BleuLinear*>(features[id].get());
      
      if (__bleu)
	__bleu->assign(score_curr);
      else
	__bleu_linear->assign(score_curr);
      
      model_type model;
      model.push_back(features[id]);
      
      cicada::apply_cube_prune(model, graphs[id], graph_oracle, weight_bleu_function(feature_bleu, score_factor), cube_size);
      
      cicada::semiring::Logprob<double> weight;
      sentence_type sentence;
      cicada::viterbi(graph_oracle, sentence, weight, kbest_traversal(), bleu_function(feature_bleu, score_factor));
      
      score_ptr_type score_sample = scorers[id]->score(sentence);
      if (score_curr)
	*score_curr += *score_sample;
      else
	score_curr = score_sample;
      
      const double objective = score_curr->score().first * score_factor;
      
      if (objective > objective_optimum || ! scores[id]) {
	score_optimum = score_curr;
	objective_optimum = objective;
	scores[id] = score_sample;
      }
    }
  }
  
  score_ptr_type score_optimum;
  
  queue_type&                          queue;
  const hypergraph_set_type&           graphs;
  const feature_function_ptr_set_type& features;
  const scorer_document_type&          scorers;
  score_ptr_set_type&                  scores;
};

void compute_oracles(const hypergraph_set_type& graphs,
		     const feature_function_ptr_set_type& features,
		     const scorer_document_type& scorers)
{
  typedef TaskOracle            task_type;
  typedef task_type::queue_type queue_type;

  typedef boost::shared_ptr<task_type> task_ptr_type;
  typedef std::vector<task_ptr_type, std::allocator<task_ptr_type> > task_set_type;
  
  score_ptr_set_type scores(graphs.size());

  score_ptr_type score_optimum;
  double objective_optimum = - std::numeric_limits<double>::infinity();
  
  const bool error_metric = scorers.error_metric();
  const double score_factor = (error_metric ? - 1.0 : 1.0);
  
  for (int iter = 0; iter < 5; ++ iter) {
    queue_type queue(graphs.size());
    
    task_set_type tasks(threads);
    for (int i = 0; i < threads; ++ i)
      tasks[i].reset(new task_type(queue, graphs, features, scorers, scores));
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(*tasks[i])));
    
    
    for (int id = 0; id < graphs.size(); ++ id)
      queue.push(id);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    workers.join_all();
    
    score_optimum.reset();
    score_ptr_set_type::const_iterator siter_end = scores.end();
    for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter) 
      if (*siter) {
	if (! score_optimum)
	  score_optimum = (*siter)->clone();
	else
	  *score_optimum += *(*siter);
      } 
    
    const double objective = score_optimum->score().first * score_factor;
    if (debug)
      std::cerr << "oracle score: " << objective << std::endl;
    
    if (objective <= objective_optimum) break;
    
    objective_optimum = objective;
  }
  
  for (int id = 0; id < graphs.size(); ++ id)
    if (features[id]) {
      if (! scores[id])
	throw std::runtime_error("no scores?");
      
      score_ptr_type score_curr = score_optimum->clone();
      *score_curr -= *scores[id];
      
      cicada::feature::Bleu*       __bleu = dynamic_cast<cicada::feature::Bleu*>(features[id].get());
      cicada::feature::BleuLinear* __bleu_linear = dynamic_cast<cicada::feature::BleuLinear*>(features[id].get());
      
      if (__bleu)
	__bleu->assign(score_curr);
      else
	__bleu_linear->assign(score_curr);
    }
}


void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 feature_function_ptr_set_type& features)
{
  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (debug)
      std::cerr << "file: " << *titer << std::endl;
    
    if (boost::filesystem::is_directory(*titer)) {
      
      for (int i = 0; /**/; ++ i) {
	const path_type path = (*titer) / (boost::lexical_cast<std::string>(i) + ".gz");

	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	
	int id;
	std::string sep;
	hypergraph_type hypergraph;
	
	while (is >> id >> sep >> hypergraph) {
	
	  if (sep != "|||")
	    throw std::runtime_error("format error?: " + path.file_string());
	
	  if (id >= graphs.size())
	    throw std::runtime_error("tstset size exceeds refset size?" + boost::lexical_cast<std::string>(id) + ": " + path.file_string());
	
	  graphs[id].unite(hypergraph);
	}
      }
    } else {
      
      utils::compress_istream is(*titer, 1024 * 1024);
      
      int id;
      std::string sep;
      hypergraph_type hypergraph;
      
      while (is >> id >> sep >> hypergraph) {
	
	if (sep != "|||")
	  throw std::runtime_error("format error?: " + titer->file_string());
	
	if (id >= graphs.size())
	  throw std::runtime_error("tstset size exceeds refset size?" + boost::lexical_cast<std::string>(id) + ": " + titer->file_string());
	
	graphs[id].unite(hypergraph);
      }
    }
  }
  
  typedef cicada::Parameter parameter_type;
    
  for (int id = 0; id < graphs.size(); ++ id) {
    if (graphs[id].goal == hypergraph_type::invalid)
      std::cerr << "invalid graph at: " << id << std::endl;
    else {
      features[id] = feature_function_type::create(scorer_name);
      
      cicada::feature::Bleu*       __bleu = dynamic_cast<cicada::feature::Bleu*>(features[id].get());
      cicada::feature::BleuLinear* __bleu_linear = dynamic_cast<cicada::feature::BleuLinear*>(features[id].get());
      
      if (! __bleu && ! __bleu_linear)
	throw std::runtime_error("invalid bleu feature function...");

      static const cicada::Lattice       __lattice;
      static const cicada::SpanVector    __spans;
      static const cicada::NGramCountSet __ngram_counts;
      
      if (__bleu)
	__bleu->assign(id, graphs[id], __lattice, __spans, sentences[id], __ngram_counts);
      else
	__bleu_linear->assign(id, graphs[id], __lattice, __spans, sentences[id], __ngram_counts);
      
      if (apply_exact) {
	model_type model;
	model.push_back(features[id]);
	
	hypergraph_type graph_reward;
	
	cicada::apply_exact(model, graphs[id], graph_reward);
	
	graphs[id].swap(graph_reward);
      }
    }
  }
}

void read_refset(const path_set_type& files,
		 scorer_document_type& scorers,
		 sentence_document_type& sentences)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;

  if (files.empty())
    throw std::runtime_error("no reference files?");
  
  sentences.clear();

  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->file_string());

    utils::compress_istream is(*fiter, 1024 * 1024);
    
    std::string line;
    
    while (std::getline(is, line)) {
      tokenizer_type tokenizer(line);
    
      tokenizer_type::iterator iter = tokenizer.begin();
      if (iter == tokenizer.end()) continue;
    
      const int id = boost::lexical_cast<int>(*iter);
      ++ iter;
    
      if (iter == tokenizer.end()) continue;
      if (*iter != "|||") continue;
      ++ iter;

      if (id >= scorers.size())
	scorers.resize(id + 1);
      
      if (id >= sentences.size())
	sentences.resize(id + 1);
      
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      sentences[id].push_back(sentence_type(iter, tokenizer.end()));
      scorers[id]->insert(sentences[id].back());
    }
  }  
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s) (in hypergraph format)")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    
    ("iteration",          po::value<int>(&iteration),          "# of mert iteration")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "regularization via L1")
    ("regularize-l2", po::bool_switch(&regularize_l2), "regularization via L2")
    ("C"            , po::value<double>(&C),           "regularization constant")

    ("oracle-loss", po::bool_switch(&oracle_loss), "loss from oracle translations")
    ("apply-exact", po::bool_switch(&apply_exact), "exact feature applicatin w/o pruning")
    ("cube-size",   po::value<int>(&cube_size),    "cube-pruning size")

    ("softmax-margin", po::bool_switch(&softmax_margin), "softmax-margin")
    
    ("threads", po::value<int>(&threads), "# of threads")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
