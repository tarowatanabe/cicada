//
// online learning using the margin between forest
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
#include "cicada/span.hpp"

#include "cicada/apply.hpp"
#include "cicada/model.hpp"

#include "cicada/feature/bleu.hpp"
#include "cicada/feature/bleu_linear.hpp"
#include "cicada/parameter.hpp"

#include "cicada/eval.hpp"

#include "cicada_impl.hpp"
#include "cicada_learn_online_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"
#include "utils/space_separator.hpp"
#include "utils/sgi_hash_map.hpp"

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>


typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type> > feature_function_ptr_set_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::scorer_ptr_type scorer_ptr_type;
typedef scorer_type::score_ptr_type  score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;


typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file = "-";
path_type output_file = "-";

bool input_id_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_directory_mode = false;

path_type weights_file;

std::string symbol_goal         = vocab_type::S;
std::string symbol_non_terminal = vocab_type::X;

grammar_file_set_type grammar_mutable_files;
grammar_file_set_type grammar_static_files;

bool grammar_glue_straight = false;
bool grammar_glue_inverted = false;
bool grammar_insertion = false;
bool grammar_deletion = false;

feature_parameter_set_type feature_parameters;
bool feature_list = false;

op_set_type ops;
bool op_list = false;

std::string scorer_name    = "bleu:order=4,exact=false";

bool learn_regression = false;
bool learn_factored = false;

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;
double loss_scale = 100;

double loss_margin = 0.001;
double score_margin = 0.001;

int batch_size = 1;
bool reranking = false;
bool asynchronous_vectors = false;
bool mix_weights = false;
bool average_weights = false;

bool apply_exact = false;
int  cube_size = 200;

int threads = 4;

int debug = 0;



void optimize(weight_set_type& weights);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (input_lattice_mode && input_forest_mode)
      throw std::runtime_error("input can be sentence, lattice or forest");
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("you cannot use both of L1 and L2...");

    if (int(learn_regression) + learn_factored > 1)
      throw std::runtime_error("you can learn one of learn-regression|learn-factored");
    if (int(learn_regression) + learn_factored == 0)
      learn_regression = true;

    if (feature_list) {
      std::cout << cicada::FeatureFunction::lists();
      return 0;
    }

    if (op_list) {
      std::cout << OperationSet::lists();
      return 0;
    }

    srandom(time(0) * getpid());
    
    threads = utils::bithack::max(threads, 1);
    
    weight_set_type weights;
    if (boost::filesystem::exists(weights_file)) {
      utils::compress_istream is(weights_file);
      is >> weights;
    }
    
    optimize(weights);
    
    utils::compress_ostream os(output_file);
    os.precision(20);
    os << weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}



template <typename Optimizer>
struct Task
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;
  typedef Optimizer optimizer_type;

  typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;
  
  Task(queue_type& __queue,
       queue_type& __queue_send,
       queue_type& __queue_recv,
       optimizer_type& __optimizer)
    : queue(__queue), queue_send(__queue_send), queue_recv(__queue_recv),
      optimizer(__optimizer),
      batch_current(0),
      score_1best(), score(), scores(), norm(0.0) { initialize(); }
  
  queue_type&         queue;
  queue_type&         queue_send;
  queue_type&         queue_recv;

  optimizer_type&    optimizer;

  int batch_current;

  score_ptr_type     score_1best;
  score_ptr_type     score;
  score_ptr_set_type scores;
  double norm;
  std::vector<int, std::allocator<int> > norms;

  hypergraph_set_type  hypergraph_oracles;
  
  grammar_type grammar;
  model_type model;
  
  typedef cicada::semiring::Logprob<double> weight_type;
  
  struct weight_set_function
  {
    typedef cicada::semiring::Logprob<double> value_type;
    
    weight_set_function(const weight_set_type& __weights, const double& __scale)
      : weights(__weights), scale(__scale) {}
    
    const weight_set_type& weights;
    const double scale;
    
    template <typename Edge>
    value_type operator()(const Edge& x) const
    {
      return cicada::semiring::traits<value_type>::log(x.features.dot(weights) * scale);
    }
    
    value_type operator()(const feature_set_type& x) const
    {
      return cicada::semiring::traits<value_type>::log(x.dot(weights) * scale);
    }

  };
  
  struct count_function
  {
    typedef int value_type;
    
    template <typename Edge>
    value_type operator()(const Edge& x) const
    {
      return 1;
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

  struct accumulated_set_unique_type
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
  
  struct accumulated_set_type
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > accumulated_type;
    typedef accumulated_type value_type;
    
    typedef std::vector<value_type, std::allocator<value_type> > value_set_type;

    accumulated_set_type(size_t x) : accumulated(x) {}
    
    accumulated_type& operator[](size_t index)
    {
      return accumulated[index];
    }
    
    void resize(size_t x) { accumulated.resize(x); }
    
    void clear() { accumulated.clear(); }

    size_t size() const { return accumulated.size(); }
    
    value_set_type accumulated;
  };
  
  struct bleu_function
  {
    typedef cicada::semiring::Tropical<double> value_type;
    
    bleu_function(const feature_type& __feature_name, const double __scale)
      : feature_name(__feature_name), scale(__scale) {}
    
    const feature_type feature_name;
    const double scale;
    
    template <typename Edge>
    value_type operator()(const Edge& x) const
    {
      return cicada::semiring::traits<value_type>::log(x.features[feature_name] * scale);
    }
  };

  typedef std::vector<double, std::allocator<double> > label_collection_type;
  typedef std::vector<double, std::allocator<double> > margin_collection_type;
  typedef std::vector<feature_set_type, std::allocator<feature_set_type> > feature_collection_type;

  typedef std::vector<int, std::allocator<int> > count_set_type;
  typedef boost::tuple<sentence_type, feature_set_type> yield_type;
  
  void prune_hypergraph(model_type& model_bleu,
			model_type& model_sparse,
			const hypergraph_type& hypergraph,
			hypergraph_type& modified,
			yield_type& yield, 
			const weight_set_type& weights,
			const weight_set_type& weights_prune,
			const double margin)
  {
    cicada::apply_cube_prune(model_bleu, hypergraph, modified, weight_set_function(weights, 1.0), cube_size);
    
    cicada::beam_prune(modified, weight_set_scaled_function<cicada::semiring::Tropical<double> >(weights_prune, 1.0), margin);
    
    if (! model_sparse.empty()) {
      model_sparse.apply_feature(true);
      cicada::apply_exact(model_sparse, modified);
      model_sparse.apply_feature(false);
    }
    
    weight_type weight;
    
    cicada::viterbi(modified, yield, weight, kbest_traversal(), weight_set_function(weights, 1.0));
  }

  void add_support_vectors_regression(const hypergraph_type& hypergraph_reward,
				      const hypergraph_type& hypergraph_penalty,
				      const feature_type& feature_name,
				      label_collection_type& labels,
				      margin_collection_type& margins,
				      feature_collection_type& features)
  {
    count_set_type counts_reward(hypergraph_reward.nodes.size());
    count_set_type counts_penalty(hypergraph_penalty.nodes.size());
    
    accumulated_set_unique_type accumulated_reward_unique;
    accumulated_set_unique_type accumulated_penalty_unique;
    
    cicada::inside_outside(hypergraph_reward,  counts_reward,  accumulated_reward_unique,  count_function(), feature_count_function());
    cicada::inside_outside(hypergraph_penalty, counts_penalty, accumulated_penalty_unique, count_function(), feature_count_function());
    
    feature_set_type features_reward(accumulated_reward_unique.accumulated.begin(), accumulated_reward_unique.accumulated.end());
    feature_set_type features_penalty(accumulated_penalty_unique.accumulated.begin(), accumulated_penalty_unique.accumulated.end());
    
    features_reward  *= (1.0 / counts_reward.back());
    features_penalty *= (1.0 / counts_penalty.back());
    
    features.push_back(features_reward - features_penalty);
    
    const double bleu_score = features.back()[feature_name];
    
    features.back()["bias"] = 1.0;
    features.back().erase(feature_name);
    
    labels.push_back(1.0);
    margins.push_back(bleu_score * norm * loss_scale);
  }
  

  void add_support_vectors_factored(const hypergraph_type& hypergraph_reward,
				    const hypergraph_type& hypergraph_penalty,
				    const feature_type& feature_name,
				    label_collection_type& labels,
				    margin_collection_type& margins,
				    feature_collection_type& features)
  {
    typedef std::pair<int, int> span_type;
    typedef std::vector<span_type, std::allocator<span_type> > span_set_type;

    typedef std::vector<typename bleu_function::value_type, std::allocator<typename bleu_function::value_type> > bleu_set_type;

    typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;

#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<span_type, edge_set_type, utils::hashmurmur<size_t>, std::equal_to<span_type>,
      std::allocator<std::pair<const span_type, edge_set_type> > > cluster_type;
#else
    typedef sgi::hash_map<span_type, edge_set_type, utils::hashmurmur<size_t>, std::equal_to<span_type>,
      std::allocator<std::pair<const span_type, edge_set_type> > > cluster_type;
#endif


    count_set_type counts_reward(hypergraph_reward.nodes.size());
    count_set_type counts_penalty(hypergraph_penalty.nodes.size());
    
    accumulated_set_type accumulated_reward(hypergraph_reward.edges.size());
    accumulated_set_type accumulated_penalty(hypergraph_penalty.edges.size());
    
    cicada::inside_outside(hypergraph_reward,  counts_reward,  accumulated_reward,  count_function(), feature_count_function());
    cicada::inside_outside(hypergraph_penalty, counts_penalty, accumulated_penalty, count_function(), feature_count_function());
    
    bleu_set_type bleu_reward(hypergraph_reward.nodes.size());
    bleu_set_type bleu_penalty(hypergraph_penalty.nodes.size());
    
    bleu_set_type bleu_edge_reward(hypergraph_reward.edges.size());
    bleu_set_type bleu_edge_penalty(hypergraph_penalty.edges.size());
    
    cicada::inside_outside(hypergraph_reward,  bleu_reward,  bleu_edge_reward,  bleu_function(feature_name,   1.0), bleu_function(feature_name,   1.0));
    cicada::inside_outside(hypergraph_penalty, bleu_penalty, bleu_edge_penalty, bleu_function(feature_name, - 1.0), bleu_function(feature_name, - 1.0));
    
    span_set_type spans_reward(hypergraph_reward.nodes.size());
    span_set_type spans_penalty(hypergraph_penalty.nodes.size());
    
    cicada::node_span(hypergraph_reward,  spans_reward);
    cicada::node_span(hypergraph_penalty, spans_penalty);
    
    // we will split features into spans...
    cluster_type cluster_reward;
    cluster_type cluster_penalty;
    
    {
      hypergraph_type::node_set_type::const_iterator riter_end = hypergraph_reward.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator riter = hypergraph_reward.nodes.begin(); riter != riter_end; ++ riter) {
	edge_set_type& edges = cluster_reward[spans_reward[riter->id]];
	edges.insert(edges.end(), riter->edges.begin(), riter->edges.end());
      }
      
      hypergraph_type::node_set_type::const_iterator piter_end = hypergraph_penalty.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator piter = hypergraph_penalty.nodes.begin(); piter != piter_end; ++ piter) {
	edge_set_type& edges = cluster_penalty[spans_penalty[piter->id]];
	edges.insert(edges.end(), piter->edges.begin(), piter->edges.end());
      }
    }
    
    cluster_type::const_iterator riter_end = cluster_reward.end();
    for (cluster_type::const_iterator riter = cluster_reward.begin(); riter != riter_end; ++ riter) {

      // collect intersection only...!
      cluster_type::const_iterator piter = cluster_penalty.find(riter->first);
      if (piter != cluster_penalty.end()) {
	
	const edge_set_type& edges_reward  = riter->second;
	const edge_set_type& edges_penalty = piter->second;

	feature_set_type features_reward;
	feature_set_type features_penalty;
	
	double bleu_reward = 0;
	double bleu_penalty = 0;
	
	{
	  edge_set_type::const_iterator eiter_end = edges_reward.end();
	  for (edge_set_type::const_iterator eiter = edges_reward.begin(); eiter != eiter_end; ++ eiter) {
	    features_reward += feature_set_type(accumulated_reward[*eiter].begin(), accumulated_reward[*eiter].end());
	    bleu_reward += bleu_edge_reward[*eiter] * norm * loss_scale;
	  }
	}

	{
	  edge_set_type::const_iterator eiter_end = edges_penalty.end();
	  for (edge_set_type::const_iterator eiter = edges_penalty.begin(); eiter != eiter_end; ++ eiter) {
	    features_penalty += feature_set_type(accumulated_penalty[*eiter].begin(), accumulated_penalty[*eiter].end());
	    bleu_penalty -= bleu_edge_penalty[*eiter] * norm * loss_scale;
	  }
	}
	
	features_reward *= (1.0 / counts_reward.back());
	bleu_reward *= (1.0 / edges_reward.size());
      
	features_penalty *= (1.0 / counts_penalty.back());
	bleu_penalty *= (1.0 / edges_penalty.size());
      
	features.push_back(features_reward - features_penalty);
	
	features.back()["bias"] = 1.0;
	features.back().erase(feature_name);

	labels.push_back(1.0);
	margins.push_back(bleu_reward - bleu_penalty);
      } 
    }
    
  }
  
  void operator()()
  {
    scorer_ptr_type           scorer(scorer_type::create(scorer_name));
    feature_function_ptr_type feature_function(feature_function_type::create(scorer_name));

    cicada::feature::Bleu* __bleu = dynamic_cast<cicada::feature::Bleu*>(feature_function.get());
    if (! __bleu)
      throw std::runtime_error("not supported feature");
    
    model_type model_bleu;
    model_bleu.push_back(feature_function);

    weight_set_type weights_bleu;
    weights_bleu[__bleu->feature_name()] = 1.0;

    model_type model_sparse;
    for (model_type::const_iterator iter = model.begin(); iter != model.end(); ++ iter)
      if ((*iter)->sparse_feature())
	model_sparse.push_back(*iter);
    
    OperationSet operations(ops.begin(), ops.end(),
			    grammar,
			    model,
			    symbol_goal,
			    symbol_non_terminal,
			    grammar_insertion,
			    grammar_deletion,
			    true,
			    input_lattice_mode,
			    input_forest_mode,
			    true,
			    false,
			    debug);

    weight_set_type& weights = optimizer.weights;

    if (! reranking)
      operations.assign(weights);
    
    hypergraph_type hypergraph_reward;
    hypergraph_type hypergraph_penalty;
    
    yield_type  yield_viterbi;
    weight_type weight_viterbi;
    
    yield_type  yield_reward;
    weight_type weight_reward;
    
    yield_type  yield_penalty;
    weight_type weight_penalty;

    
    label_collection_type   labels;
    margin_collection_type  margins;
    feature_collection_type features;

    sentence_set_type targets;

    bool terminated_merge = false;
    bool terminated_sample = false;

    std::string buffer;
    while (! terminated_merge || ! terminated_sample) {
      
      while (! terminated_merge && queue_recv.pop_swap(buffer, true)) {
	if (buffer.empty()) {
	  terminated_merge = true;
	  boost::thread::yield();
	  break;
	} else {
	  
	  size_t id = size_t(-1);
	  int source_length;
	  decode_support_vectors(buffer, id, source_length, boost::get<0>(yield_viterbi), hypergraph_reward, hypergraph_penalty, targets);

	  if (id == size_t(-1))
	    throw std::runtime_error("invalid encoded feature vector");

	  if (id >= scores.size())
	    scores.resize(id + 1);
	  
	  if (id >= norms.size())
	    norms.resize(id + 1);
	  
	  // remove "this" score
#if 0
	  if (score && scores[id])
	    *score -= *scores[id];

	  norm += source_length;
	  norm -= norms[id];
	  norms[id] = source_length;
#endif
	  	  
	  norm *= 0.9;
	  if (score)
	    *score *= 0.9;
	  norm += source_length;
	  
	  if (score_1best && scores[id])
	    *score_1best -= *scores[id];

	  scorer->clear();
	  sentence_set_type::const_iterator titer_end = targets.end();
	  for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
	    scorer->insert(*titer);
	    
	  
	  scores[id] = scorer->score(boost::get<0>(yield_viterbi));
	  if (! score)
	    score = scores[id]->clone();
	  else
	    *score += *scores[id];
	  
	  if (! score_1best)
	    score_1best = scores[id]->clone();
	  else
	    *score_1best += *scores[id];

	  
	  if (debug >= 2) {
	    const std::pair<double, double> bleu_1best = score_1best->score();
	    const std::pair<double, double> bleu_viterbi = score->score();
	    
	    std::cerr << "bleu: " << bleu_1best.first
		      << " peanlty: " << bleu_1best.second
		      << " viterbi: " << bleu_viterbi.first
		      << " penalty: " << bleu_viterbi.second
		      << std::endl;
	  }

	  if (learn_factored)
	    add_support_vectors_factored(hypergraph_reward, hypergraph_penalty, __bleu->feature_name(), labels, margins, features);
	  else
	    add_support_vectors_regression(hypergraph_reward, hypergraph_penalty, __bleu->feature_name(), labels, margins, features);
	  
	  ++ batch_current;
	  
	  if (batch_current >= batch_size && ! labels.empty()) {
	    if (debug)
	      std::cerr << "# of support vectors: " << labels.size() << std::endl;
	    
	    optimizer(labels, margins, features);
	    
	    batch_current = 0;
	    labels.clear();
	    margins.clear();
	    features.clear();
	  }
	}
      }
      
      if (terminated_sample) continue;
      if (! queue.pop_swap(buffer, true)) {
	boost::thread::yield();
	continue;
      }
      
      if (buffer.empty()) {
	terminated_sample = true;
	queue_send.push(std::string());
	boost::thread::yield();
	continue;
      }
      
      model.apply_feature(false);
      
      operations(buffer);
      
      // operations.hypergraph contains result...
      const size_t& id = operations.id;
      const lattice_type& lattice = operations.lattice;
      const hypergraph_type& hypergraph = operations.hypergraph;
      const sentence_set_type& targets = operations.targets;
      
      if (debug)
	std::cerr << "id: " << id << std::endl;
      
      // compute source-length
      int source_length = lattice.shortest_distance();
      if (hypergraph.is_valid()) {
	// we will enumerate forest structure... and collect min-size...
	std::vector<source_length_function::value_type, std::allocator<source_length_function::value_type> > lengths(hypergraph.nodes.size());
	
	cicada::inside(hypergraph, lengths, source_length_function());
	
	source_length = - log(lengths.back());
      }
            
      // collect max-feature from hypergraph
      cicada::viterbi(hypergraph, yield_viterbi, weight_viterbi, kbest_traversal(), weight_set_function(weights, 1.0));
      
      // update scores...
      if (id >= scores.size())
	scores.resize(id + 1);

      if (id >= norms.size())
	norms.resize(id + 1);
      
#if 0
      if (score && scores[id])
        *score -= *scores[id];
      
      norm += source_length;
      norm -= norms[id];
      norms[id] = source_length;
#endif
      
      norm *= 0.9;
      if (score)
	*score *= 0.9;
      norm += source_length;
      
      // create scorers...
      scorer->clear();
      __bleu->clear();
      sentence_set_type::const_iterator titer_end = targets.end();
      for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer) {
	scorer->insert(*titer);
	__bleu->insert(source_length, *titer);
      }
      __bleu->insert(score);
      
      // compute bleu-rewarded instance
      weights[__bleu->feature_name()] =  loss_scale * norm;
      
      if (id >= hypergraph_oracles.size())
	hypergraph_oracles.resize(id + 1);
      
      hypergraph_oracles[id].unite(hypergraph);
      {
	hypergraph_type hypergraph_reward;
      
	prune_hypergraph(model_bleu, model_sparse, hypergraph_oracles[id], hypergraph_reward, yield_reward, weights, weights_bleu, loss_margin);
      
	hypergraph_oracles[id].swap(hypergraph_reward);
      }
      
      const hypergraph_type& hypergraph_reward = hypergraph_oracles[id];
      
      // compute bleu-penalty hypergraph
      weights[__bleu->feature_name()] = - loss_scale * norm;
      
      prune_hypergraph(model_bleu, model_sparse, hypergraph, hypergraph_penalty, yield_penalty, weights, weights, score_margin);
      
      // erase unused weights...
      weights.erase(__bleu->feature_name());
      
      score_ptr_type score_reward  = scorer->score(boost::get<0>(yield_reward));
      score_ptr_type score_penalty = scorer->score(boost::get<0>(yield_penalty));
      if (score) {
	*score_reward  += *score;
	*score_penalty += *score;
      }

      if (score_1best && scores[id])
	*score_1best -= *scores[id];
      
      scores[id] = scorer->score(boost::get<0>(yield_viterbi));
      if (! score)
	score = scores[id]->clone();
      else
	*score += *scores[id];

      if (! score_1best)
	score_1best = scores[id]->clone();
      else
	*score_1best += *scores[id];
      
      if (debug) {
	const std::pair<double, double> bleu_1best = score_1best->score();
	const std::pair<double, double> bleu_viterbi = score->score();
	
	std::cerr << "viterbi: " << boost::get<0>(yield_viterbi) << std::endl;
	
	std::cerr << "hypergraph density:"
		  << " oracle: " << (double(hypergraph_reward.edges.size()) / hypergraph_reward.nodes.size())
		  << " violated: " << (double(hypergraph_penalty.edges.size()) / hypergraph_penalty.nodes.size())
		  << std::endl;
	
	std::cerr << "bleu: " << bleu_1best.first
		  << " peanlty: " << bleu_1best.second
		  << " viterbi: " << bleu_viterbi.first
		  << " penalty: " << bleu_viterbi.second
		  << " oracle: " << score_reward->score().first
		  << " violated: " << score_penalty->score().first
		  << std::endl;
      }

      if (learn_factored)
	add_support_vectors_factored(hypergraph_reward, hypergraph_penalty, __bleu->feature_name(), labels, margins, features);
      else
	add_support_vectors_regression(hypergraph_reward, hypergraph_penalty, __bleu->feature_name(), labels, margins, features);
      
      ++ batch_current;
      
      if (batch_current >= batch_size && ! labels.empty()) {
	if (debug)
	  std::cerr << "# of support vectors: " << labels.size() << std::endl;
	
	optimizer(labels, margins, features);
	
	batch_current = 0;
	labels.clear();
	margins.clear();
	features.clear();
      }
      
      if (asynchronous_vectors) {
	encode_support_vectors(buffer, id, source_length, boost::get<0>(yield_viterbi), hypergraph_reward, hypergraph_penalty, targets);
	
	queue_send.push_swap(buffer);
      }
    }

    if (! labels.empty()) {
      if (debug)
	std::cerr << "# of support vectors: " << labels.size() << std::endl;
      
      optimizer(labels, margins, features);

      batch_current = 0;
      labels.clear();
      margins.clear();
      features.clear();
    }

    
    // final function call...
    optimizer.finalize();
  }
  
  void initialize()
  {
    // read grammars...
    for (grammar_file_set_type::const_iterator fiter = grammar_static_files.begin(); fiter != grammar_static_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarStatic(*fiter)));
    
    if (debug)
      std::cerr << "loaded mutable grammar: " << grammar.size() << std::endl;
    
    for (grammar_file_set_type::const_iterator fiter = grammar_mutable_files.begin(); fiter != grammar_mutable_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarMutable(*fiter)));

    if (debug)
      std::cerr << "loaded static grammar: " << grammar.size() << std::endl;
    
    if (grammar_glue_straight || grammar_glue_inverted)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										  symbol_non_terminal,
										  grammar_glue_straight,
										  grammar_glue_inverted)));

    if (debug)
      std::cerr << "grammar: " << grammar.size() << std::endl;
    
    // read features...
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));
    model.initialize();
  }  
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

void optimize(weight_set_type& weights)
{
  typedef OptimizeMIRA optimizer_type;
  typedef std::vector<optimizer_type, std::allocator<optimizer_type> > optimizer_set_type;
  
  typedef Task<optimizer_type>         task_type;
  typedef boost::shared_ptr<task_type> task_ptr_type;
  typedef std::vector<task_ptr_type, std::allocator<task_ptr_type> > task_ptr_set_type;
  
  typedef task_type::queue_type queue_type;
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;
  
  typedef std::vector<std::string, std::allocator<std::string> > sample_set_type;
  
  // read all the training data...
  sample_set_type samples;
  if (input_directory_mode) {
    std::string line;
      
    boost::filesystem::directory_iterator iter_end;
    for (boost::filesystem::directory_iterator iter(input_file); iter != iter_end; ++ iter) {
      utils::compress_istream is(*iter, 1024 * 1024);
      
      if (std::getline(is, line) && ! line.empty())
	samples.push_back(line);
    }
  } else {
    utils::compress_istream is(input_file, 1024 * 1024);
    
    size_t id = 0;

    std::string line;
    while (std::getline(is, line))
      if (! line.empty()) {
	if (! input_id_mode)
	  samples.push_back(boost::lexical_cast<std::string>(id) + " ||| " + line);
	else
	  samples.push_back(line);
	++ id;
      }
  }

  if (debug)
    std::cerr << "# of samples: " << samples.size() << std::endl;
  
  queue_type queue(samples.size() + threads);
  queue_ptr_set_type queue_reduce(threads);
  queue_ptr_set_type queue_bcast(threads);
  
  for (int i = 0; i < threads; ++ i) {
    queue_reduce[i].reset(new queue_type());
    queue_bcast[i].reset(new queue_type());
  }

  optimizer_set_type optimizers(threads, optimizer_type(weights, C, debug));
  
  task_ptr_set_type tasks(threads);
  for (int i = 0; i < threads; ++ i)
    tasks[i].reset(new task_type(queue, *queue_reduce[i], *queue_bcast[i], optimizers[i]));
  
  weight_set_type weights_mixed;
  weight_set_type weights_accumulated;
  double norm_accumulated = 0;
  
  for (int iter = 0; iter < iteration; ++ iter) {
    
    if (debug)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(*tasks[i])));
    
    sample_set_type::const_iterator siter_end = samples.end();
    for (sample_set_type::const_iterator siter = samples.begin(); siter != siter_end; ++ siter)
      queue.push(*siter);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(std::string());

    std::vector<bool, std::allocator<bool> > terminated(threads, false);
    
    std::string buffer;
    int non_found_iter = 0;
    while (1) {
      bool found = false;
      
      for (int i = 0; i < threads; ++ i)
	if (! terminated[i] && queue_reduce[i]->pop_swap(buffer, true)) {
	  if (buffer.empty())
	    terminated[i] = true;
	  else {
	    for (int j = 0; j < threads; ++ j)
	      if (j != i)
		queue_bcast[i]->push(buffer);
	  }
	  
	  found = true;
	}
      
      if (std::count(terminated.begin(), terminated.end(), true) == threads)
	break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    for (int i = 0; i < threads; ++ i)
      queue_bcast[i]->push(std::string());
    
    workers.join_all();
    
    std::random_shuffle(samples.begin(), samples.end());
    
    // merge vector...
    weights_mixed.clear();
    
    optimizer_set_type::iterator oiter_end = optimizers.end();
    for (optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
      weights_mixed       += oiter->weights;
      weights_accumulated += oiter->accumulated;
      norm_accumulated    += oiter->updated;
    }
    
    weights_mixed *= (1.0 / tasks.size());
    
    {
      // output mixed weights...
      bool has_suffix_gz  = false;
      bool has_suffix_bz2 = false;
      
      path_type path_output = output_file;
      
      if (path_output.extension() == ".gz") {
	path_output = path_output.parent_path() / path_output.stem();
	has_suffix_gz = true;
      } else if (path_output.extension() == ".bz2") {
	path_output = path_output.parent_path() / path_output.stem();
	has_suffix_bz2 = true;
      }
      
      if (has_suffix_gz)
	path_output = path_output.file_string() + '.' + boost::lexical_cast<std::string>(iter + 1) + ".gz";
      else if (has_suffix_bz2)
	path_output = path_output.file_string() + '.' + boost::lexical_cast<std::string>(iter + 1) + ".bz2";
      else
	path_output = path_output.file_string() + '.' + boost::lexical_cast<std::string>(iter + 1);

      if (average_weights) {
	weights = weights_accumulated;
	weights /= norm_accumulated;
	
	utils::compress_ostream os(path_output, 1024 * 1024);
	os.precision(20);
	os << weights;
      } else {
	utils::compress_ostream os(path_output, 1024 * 1024);
	os.precision(20);
	os << weights_mixed;
      }
    }
    
    for (optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
      if (mix_weights)
	oiter->weights = weights_mixed;
      oiter->accumulated.clear();
      oiter->updated = 1;
    }
  }
  
  weights = weights_mixed;
  
  if (average_weights) {
    weights = weights_accumulated;
    weights /= norm_accumulated;
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    // options for input/output format
    ("input-id",         po::bool_switch(&input_id_mode),         "id-prefixed input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")

    ("weights", po::value<path_type>(&weights_file), "initial weights")
    
    // grammar
    ("goal",           po::value<std::string>(&symbol_goal)->default_value(symbol_goal),                 "goal symbol")
    ("non-terminal",   po::value<std::string>(&symbol_non_terminal)->default_value(symbol_non_terminal), "default non-terminal symbol")
    ("grammar",        po::value<grammar_file_set_type >(&grammar_mutable_files)->composing(),           "grammar file(s)")
    ("grammar-static", po::value<grammar_file_set_type >(&grammar_static_files)->composing(),            "static binary grammar file(s)")
    
    // special handling
    ("grammar-glue-straight", po::bool_switch(&grammar_glue_straight), "add straight hiero glue rule")
    ("grammar-glue-inverted", po::bool_switch(&grammar_glue_inverted), "add inverted hiero glue rule")
    ("grammar-insertion",     po::bool_switch(&grammar_insertion),     "source-to-target transfer grammar")
    ("grammar-deletion",      po::bool_switch(&grammar_deletion),      "source-to-<epsilon> transfer grammar")
    
    // models...
    ("feature-function",      po::value<feature_parameter_set_type >(&feature_parameters)->composing(), "feature function(s)")
    ("feature-function-list", po::bool_switch(&feature_list),                                           "list of available feature function(s)")

    //operatins...
    ("operation",      po::value<op_set_type>(&ops)->composing(), "operations")
    ("operation-list", po::bool_switch(&op_list),                 "list of available operation(s)")
    
    // learning related..
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")

    ("learn-regression", po::bool_switch(&learn_regression), "learn by regression")
    ("learn-factored",   po::bool_switch(&learn_factored),   "learn by edge-factored linear classification")
    
    ("iteration",          po::value<int>(&iteration),          "# of mert iteration")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "regularization via L1")
    ("regularize-l2", po::bool_switch(&regularize_l2), "regularization via L2")
    ("C"            , po::value<double>(&C),           "regularization constant")

    ("loss-scale",    po::value<double>(&loss_scale)->default_value(loss_scale),     "loss scaling")
    
    ("loss-margin",   po::value<double>(&loss_margin)->default_value(loss_margin),   "loss margin for oracle forest")
    ("score-margin",  po::value<double>(&score_margin)->default_value(score_margin), "score margin for hypothesis forest")
    
    ("batch-size",           po::value<int>(&batch_size)->default_value(batch_size), "batch size")
    ("reranking",            po::bool_switch(&reranking),                            "learn by forest reranking")
    ("asynchronous-vectors", po::bool_switch(&asynchronous_vectors),                 "asynchrounsly merge support vectors")
    ("average-weights",      po::bool_switch(&average_weights),                      "average weight vectors")
    ("mix-weights",          po::bool_switch(&mix_weights),                          "mixing weight vectors at every epoch")
    
    ("apply-exact", po::bool_switch(&apply_exact), "exact feature applicatin w/o pruning")
    ("cube-size",   po::value<int>(&cube_size),    "cube-pruning size")
    
    ("threads", po::value<int>(&threads), "# of threads")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config", po::value<path_type>(), "configuration file")
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;
  
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  if (variables.count("config")) {
    utils::compress_istream is(variables["config"].as<path_type>());
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
