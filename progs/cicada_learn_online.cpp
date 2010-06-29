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

std::string algorithm_name = "perceptron";
std::string scorer_name    = "bleu:order=4,exact=false";

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;
double loss_scale = 100;

bool bundle = false;
double loss_margin = 0.001;
double score_margin = 0.001;

bool average_vectors = false;
bool mix_vectors = false;

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
  
  Task(queue_type& __queue,
       queue_type& __queue_send,
       queue_type& __queue_recv,
       optimizer_type& __optimizer)
    : queue(__queue), queue_send(__queue_send), queue_recv(__queue_recv),
      optimizer(__optimizer) { initialize(); }
  
  queue_type&         queue;
  queue_type&         queue_send;
  queue_type&         queue_recv;

  optimizer_type&    optimizer;

  score_ptr_type     score;
  score_ptr_set_type scores;
  
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
  
  struct accumulated_set_type
  {
    typedef cicada::semiring::Log<double> weight_type;
    typedef cicada::FeatureVector<weight_type, std::allocator<weight_type> > accumulated_type;
    typedef accumulated_type value_type;
    
    template <typename Index>
    accumulated_type& operator[](Index)
    {
      return accumulated;
    }
    
    void clear() { accumulated.clear(); }
    
    accumulated_type accumulated;
  };
  
  
  void operator()()
  {
    typedef boost::tuple<sentence_type, feature_set_type> yield_type;

    scorer_ptr_type           scorer(scorer_type::create(scorer_name));
    feature_function_ptr_type feature_function(feature_function_type::create(scorer_name));

    cicada::feature::Bleu* __bleu = dynamic_cast<cicada::feature::Bleu*>(feature_function.get());
    if (! __bleu)
      throw std::runtime_error("not supported feature");
    
    model_type model_bleu;
    model_bleu.push_back(feature_function);

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

    operations.assign(weights);
    
    hypergraph_type hypergraph_reward;
    hypergraph_type hypergraph_penalty;
    
    yield_type  yield_viterbi;
    weight_type weight_viterbi;
    
    yield_type  yield_reward;
    weight_type weight_reward;
    
    yield_type  yield_penalty;
    weight_type weight_penalty;

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
	  double loss;
	  decode_feature_vectors(buffer, id, loss, boost::get<0>(yield_viterbi), boost::get<1>(yield_reward), boost::get<1>(yield_penalty));

	  if (id == size_t(-1))
	    throw std::runtime_error("invalid encoded feature vector");
	  
	  optimizer(loss, boost::get<1>(yield_reward), boost::get<1>(yield_penalty));
	  
	  if (id >= scores.size())
	    scores.resize(id + 1);
	  
	  // remove "this" score
	  if (score && scores[id])
	    *score -= *scores[id];
	  
	  scores[id] = scorer->score(boost::get<0>(yield_viterbi));
	  if (! score)
	    score = scores[id]->clone();
	  else
	    *score += *scores[id];
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
      
      // remove "this" score
      if (score && scores[id])
	*score -= *scores[id];
      
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
      weights[__bleu->feature_name()] =  loss_scale * source_length * scores.size();
      
      // cube-pruning for bleu computation
      cicada::apply_cube_prune(model_bleu, hypergraph, hypergraph_reward, weight_set_function(weights, 1.0), cube_size);
      
      // prune by bleu-score
      cicada::beam_prune(hypergraph_reward, weight_set_scaled_function<cicada::semiring::Tropical<double> >(weights, 1.0), loss_margin);
      
      // apply sparce features again
      if (! model_sparse.empty()) {
	model_sparse.apply_feature(true);
	cicada::apply_exact(model_sparse, hypergraph_reward, weight_set_function(weights, 1.0), cube_size);
      }
      
      cicada::viterbi(hypergraph_reward, yield_reward, weight_reward, kbest_traversal(), weight_set_function(weights, 1.0));
      if (bundle) {
	std::vector<int, std::allocator<int> > counts(hypergraph_reward.nodes.size());
	accumulated_set_type accumulated;
	
	cicada::inside_outside(hypergraph_reward, counts, accumulated, count_function(), feature_count_function());
	
	accumulated.accumulated /= counts.back();
	boost::get<1>(yield_reward).assign(accumulated.accumulated.begin(), accumulated.accumulated.end());
      }
      
      // compute bleu-penalty hypergraph
      weights[__bleu->feature_name()] = - loss_scale * source_length * scores.size();
      
      // cube-pruning for bleu computation
      cicada::apply_cube_prune(model_bleu, hypergraph, hypergraph_penalty, weight_set_function(weights, 1.0), cube_size);
      
      // prune...
      cicada::beam_prune(hypergraph_penalty, weight_set_scaled_function<cicada::semiring::Tropical<double> >(weights, 1.0), score_margin);
      
      // apply sparce features again
      if (! model_sparse.empty()) {
	model_sparse.apply_feature(true);
	cicada::apply_exact(model_sparse, hypergraph_penalty, weight_set_function(weights, 1.0), cube_size);
      }
      
      cicada::viterbi(hypergraph_penalty, yield_penalty, weight_penalty, kbest_traversal(), weight_set_function(weights, 1.0));
      if (bundle) {
	std::vector<int, std::allocator<int> > counts(hypergraph_penalty.nodes.size());
	accumulated_set_type accumulated;
	
	cicada::inside_outside(hypergraph_penalty, counts, accumulated, count_function(), feature_count_function());
	
	accumulated.accumulated /= counts.back();
	boost::get<1>(yield_penalty).assign(accumulated.accumulated.begin(), accumulated.accumulated.end());
      }

      const double bleu_reward  = boost::get<1>(yield_reward)[__bleu->feature_name()]; 
      const double bleu_penalty = boost::get<1>(yield_penalty)[__bleu->feature_name()];
      
      const double bleu_loss =  bleu_reward - bleu_penalty;
      const double loss = loss_scale * source_length * bleu_loss * scores.size();

      scores[id] = scorer->score(boost::get<0>(yield_viterbi));
      if (! score)
	score = scores[id]->clone();
      else
	*score += *scores[id];
      
      const std::pair<double, double> bleu_viterbi = score->score();
      
      if (debug) {
	std::cerr << "viterbi:  " << boost::get<0>(yield_viterbi) << std::endl
		  << "oracle:   " << boost::get<0>(yield_reward) << std::endl
		  << "violated: " << boost::get<0>(yield_penalty) << std::endl;
	
	std::cerr << "hypergraph density:"
		  << " oracle: " << (double(hypergraph_reward.edges.size()) / hypergraph_reward.nodes.size())
		  << " violated: " << (double(hypergraph_penalty.edges.size()) / hypergraph_penalty.nodes.size())
		  << std::endl;
	
	std::cerr << "bleu: " << bleu_viterbi.first
		  << " peanlty: " << bleu_viterbi.second
		  << " loss: " << bleu_loss
		  << " oracle: " << bleu_reward
		  << " violated: " << bleu_penalty
		  << " scaled: " << loss
		  << std::endl;
      }
      
      // reset bleu scores...
      weights.erase(__bleu->feature_name());
      boost::get<1>(yield_reward).erase(__bleu->feature_name());
      boost::get<1>(yield_penalty).erase(__bleu->feature_name());
      
      optimizer(loss, boost::get<1>(yield_reward), boost::get<1>(yield_penalty));
      
      encode_feature_vectors(buffer, id, loss, boost::get<0>(yield_viterbi), boost::get<1>(yield_reward), boost::get<1>(yield_penalty));
      queue_send.push_swap(buffer);
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

  optimizer_set_type optimizers(threads, optimizer_type(C, debug));
  
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
    weights.clear();
    weights_mixed.clear();
    double norm = 0.0;
    
    optimizer_set_type::iterator oiter_end = optimizers.end();
    for (optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
      weight_set_type& weights_local = oiter->weights;
      
      // mixing weights...
      weights_mixed += weights_local;
      
      weights_local *= oiter->updated;
      weights_local -= oiter->accumulated;
      
      weights += weights_local;
      norm    += oiter->updated;
    }
    
    weights_accumulated += weights;
    norm_accumulated += norm;
    
    // mixing...
    weights_mixed *= (1.0 / tasks.size());
    
    for (optimizer_set_type::iterator oiter = optimizers.begin(); oiter != oiter_end; ++ oiter) {
      oiter->weights = weights_mixed;
      oiter->accumulated.clear();
      oiter->updated = 1;
    }
  }
  
  weights = weights_accumulated;
  weights /= norm_accumulated;
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
    ("algorithm",   po::value<std::string>(&algorithm_name)->default_value(algorithm_name),     "learning algorithm")
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    
    ("iteration",          po::value<int>(&iteration),          "# of mert iteration")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "regularization via L1")
    ("regularize-l2", po::bool_switch(&regularize_l2), "regularization via L2")
    ("C"            , po::value<double>(&C),           "regularization constant")

    ("loss-scale",    po::value<double>(&loss_scale)->default_value(loss_scale),     "loss scaling")
    
    ("bundle",        po::bool_switch(&bundle),                                      "bundle support vectors from hypergraphs")
    ("loss-margin",   po::value<double>(&loss_margin)->default_value(loss_margin),   "loss margin for oracle forest")
    ("score-margin",  po::value<double>(&score_margin)->default_value(score_margin), "score margin for hypothesis forest")
    
    ("average-vectors", po::bool_switch(&average_vectors), "average vectors")
    ("mix-vectors",     po::bool_switch(&mix_vectors),     "mixing at every epoch")
    
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
