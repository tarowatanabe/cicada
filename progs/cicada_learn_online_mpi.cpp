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

#include "utils/base64.hpp"
#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"

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
std::string scorer_name = "bleu:order=4,exact=false";

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;
double loss_margin = 0.01;
double score_margin = 0.01;
double loss_scale = 100;

bool apply_exact = false;
int cube_size = 200;

int debug = 0;

void optimize(OperationSet& operations, model_type& model, weight_set_type& weights);

void bcast_weights(const int rank, weight_set_type& weights);
void send_weights(const weight_set_type& weights);
void reduce_weights(weight_set_type& weights);
void options(int argc, char** argv);

enum {
  
  weights_tag = 1000,
  sample_tag,
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

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    options(argc, argv);
    
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
    
    // read grammars...
    grammar_type grammar;
    for (grammar_file_set_type::const_iterator fiter = grammar_static_files.begin(); fiter != grammar_static_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarStatic(*fiter)));

    if (debug && mpi_rank == 0)
      std::cerr << "loaded mutable grammar: " << grammar.size() << std::endl;
    
    for (grammar_file_set_type::const_iterator fiter = grammar_mutable_files.begin(); fiter != grammar_mutable_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarMutable(*fiter)));

    if (debug && mpi_rank == 0)
      std::cerr << "loaded static grammar: " << grammar.size() << std::endl;
    
    if (grammar_glue_straight || grammar_glue_inverted)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										  symbol_non_terminal,
										  grammar_glue_straight,
										  grammar_glue_inverted)));

    if (debug && mpi_rank == 0)
      std::cerr << "grammar: " << grammar.size() << std::endl;
    
    // read features...
    model_type model;
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));
    model.initialize();
    
    if (debug && mpi_rank == 0)
      std::cerr << "feature functions: " << model.size() << std::endl;

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
    
    weight_set_type weights;
    
    optimize(operations, model, weights);
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_file);
      os.precision(20);
      os << weights;
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  return 0;
}


struct Task
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  Task(queue_type& __queue,
       OperationSet& __operations,
       model_type& __model)
    : queue(__queue),
      operations(__operations),
      model(__model),
      weights(), weights_accumulated(), updated(1) { }

  queue_type&        queue;
  OperationSet&      operations;
  model_type&        model;
  
  weight_set_type    weights;
  weight_set_type    weights_accumulated;
  long               updated;
  score_ptr_type     score;
  score_ptr_set_type scores;
  
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

  struct diff_norm
  {
    double operator()(const double& x, const double& y) const{
      return (x - y) * (x - y);
    }
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

    weight_set_type weights_bleu;
    weights_bleu[__bleu->feature_name()] = 1.0;
    
    operations.assign(weights);
    
    hypergraph_type hypergraph_reward;
    hypergraph_type hypergraph_penalty;      
    
    std::string line;
    
    while (1) {
      queue.pop(line);
      if (line.empty()) break;

      model.apply_feature(false);
      
      operations(line);
      
      // operations.hypergraph contains result...
      const size_t& id = operations.id;
      const lattice_type& lattice = operations.lattice;
      const hypergraph_type& hypergraph = operations.hypergraph;
      const sentence_set_type& targets = operations.targets;

      std::cerr << "id: " << id << std::endl;

      // compute source-length
      int source_length = lattice.shortest_distance();
      if (hypergraph.is_valid()) {
	// we will enumerate forest structure... and collect min-size...
	std::vector<source_length_function::value_type, std::allocator<source_length_function::value_type> > lengths(hypergraph.nodes.size());
	
	cicada::inside(hypergraph, lengths, source_length_function());
	
	source_length = - log(lengths.back());
      }
      
      std::cerr << "source length: " << source_length << std::endl;

      // collect max-feature from hypergraph
      yield_type  yield_viterbi;
      weight_type weight_viterbi;
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
      yield_type  yield_reward;
      weight_type weight_reward;
      
      weights[__bleu->feature_name()] =  loss_scale * source_length;
      // cube-pruning for bleu computation
      cicada::apply_cube_prune(model_bleu, hypergraph, hypergraph_reward, weight_set_function(weights, 1.0), cube_size);
      
      // prune by bleu-score
      cicada::beam_prune<cicada::semiring::Tropical<double>, weight_set_type>(hypergraph_reward, weights_bleu, 1.0, loss_margin);
      
      // apply sparce features again
      if (! model_sparse.empty()) {
	model_sparse.apply_feature(true);
	cicada::apply_exact(model_sparse, hypergraph_reward, weight_set_function(weights, 1.0), cube_size);
      }
      
      // compute viterbi
      cicada::viterbi(hypergraph_reward, yield_reward, weight_reward, kbest_traversal(), weight_set_function(weights_bleu, 1.0));
      
      // compute bleu-penalty hypergraph
      yield_type  yield_penalty;
      weight_type weight_penalty;
      
      weights[__bleu->feature_name()] = - loss_scale * source_length;
      // cube-pruning for bleu computation
      cicada::apply_cube_prune(model_bleu, hypergraph, hypergraph_penalty, weight_set_function(weights, 1.0), cube_size);
      
      // prune...
      cicada::beam_prune<cicada::semiring::Tropical<double>, weight_set_type>(hypergraph_penalty, weights, 1.0, score_margin);
      
      // apply sparce features again
      if (! model_sparse.empty()) {
	model_sparse.apply_feature(true);
	cicada::apply_exact(model_sparse, hypergraph_penalty, weight_set_function(weights, 1.0), cube_size);
      }
      
      // compute biterbi
      cicada::viterbi(hypergraph_penalty, yield_penalty, weight_penalty, kbest_traversal(), weight_set_function(weights, 1.0));
      
      // reset bleu scores...
      weights.erase(__bleu->feature_name());
      boost::get<1>(yield_reward).erase(__bleu->feature_name());
      boost::get<1>(yield_penalty).erase(__bleu->feature_name());
      
      score_ptr_type score_reward  = scorer->score(boost::get<0>(yield_reward));
      score_ptr_type score_penalty = scorer->score(boost::get<0>(yield_penalty));
      scores[id] = scorer->score(boost::get<0>(yield_viterbi));
      
      if (score) {
	*score_reward  += *score;
	*score_penalty += *score;
      }
      
      if (! score)
	score = scores[id]->clone();
      else
	*score += *scores[id];
      
      const std::pair<double, double> bleu_viterbi = score->score();
      const std::pair<double, double> bleu_reward  = score_reward->score();
      const std::pair<double, double> bleu_penalty = score_penalty->score();
      
      if (debug)
	std::cerr << "viterbi:  " << boost::get<0>(yield_viterbi) << std::endl
		  << "oracle:   " << boost::get<0>(yield_reward) << std::endl
		  << "violated: " << boost::get<0>(yield_penalty) << std::endl
		  << "viterbi score:  " << bleu_viterbi.first << " penalty: " << bleu_viterbi.second << std::endl
		  << "oracle score:   " << bleu_reward.first  << " penalty: " << bleu_reward.second << std::endl
		  << "violated score: " << bleu_penalty.first << " penalty: " << bleu_penalty.second << std::endl;
      
      const double loss   = loss_scale * source_length * (bleu_reward.first - bleu_penalty.first);
      const double margin = boost::get<1>(yield_penalty).dot(weights) - boost::get<1>(yield_reward).dot(weights);
      const double norm   = boost::get<1>(yield_penalty).dot(boost::get<1>(yield_reward), diff_norm());
      
      const double alpha = std::max(0.0, std::min(1.0 / C, (loss - margin) / norm));
      
      if (loss - margin > 0.0) {
	std::cerr << "loss: " << loss << " margin: " << margin << " norm: " << norm << " alpha: " << alpha << std::endl;

	// update...
	feature_set_type::const_iterator riter_end = boost::get<1>(yield_reward).end();
	for (feature_set_type::const_iterator riter = boost::get<1>(yield_reward).begin(); riter != riter_end; ++ riter) {
	  weights[riter->first] += riter->second * alpha;
	  weights_accumulated[riter->first] += riter->second * alpha * updated;
	}
	
	feature_set_type::const_iterator piter_end = boost::get<1>(yield_penalty).end();
	for (feature_set_type::const_iterator piter = boost::get<1>(yield_penalty).begin(); piter != piter_end; ++ piter) {
	  weights[piter->first] -= piter->second * alpha;
	  weights_accumulated[piter->first] -= piter->second * alpha * updated;
	}
	
	++ updated;
      }
    }
  }
  
  
};

void optimize(OperationSet& operations, model_type& model, weight_set_type& weights)
{
  typedef Task                  task_type;
  typedef task_type::queue_type queue_type;

  typedef std::vector<std::string, std::allocator<std::string> > sample_set_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  sample_set_type samples;
    
  // read all the training data...
  if (mpi_rank == 0) {
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
  }
  
  queue_type queue(1);
  task_type task(queue, operations, model);

  weight_set_type weights_mixed;
  weight_set_type weights_accumulated;
  double norm_accumulated = 0;
  
  for (int iter = 0; iter < iteration; ++ iter) {

    boost::thread thread(boost::ref(task));
    
    if (mpi_rank == 0) {
      typedef utils::mpi_ostream ostream_type;
      typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
      
      std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > stream(mpi_size);
      
      for (int rank = 1; rank < mpi_size; ++ rank)
	stream[rank].reset(new ostream_type(rank, sample_tag, 4096));
      
      sample_set_type::const_iterator siter = samples.begin();
      sample_set_type::const_iterator siter_end = samples.end();
      
      int non_found_iter = 0;
      while (siter != siter_end) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size && siter != siter_end; ++ rank)
	  if (stream[rank]->test()) {
	    stream[rank]->write(*siter);
	    ++ siter;
	    
	    found = true;
	  }
	
	if (queue.empty() && siter != siter_end) {
	  queue.push(*siter);
	  ++ siter;
	  
	  found = true;
	}
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }

      while (1) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size; ++ rank) 
	  if (stream[rank] && stream[rank]->test()) {
	    if (! stream[rank]->terminated())
	      stream[rank]->terminate();
	    else
	      stream[rank].reset();
	    
	    found = true;
	  }
	
	if (std::count(stream.begin(), stream.end(), ostream_ptr_type()) == mpi_size
	    && queue.push(std::string(), true))
	  break;
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
      
      std::random_shuffle(samples.begin(), samples.end());
      
    } else {
      utils::mpi_istream is(0, sample_tag, 4096, true);
      
      std::string line;
      while (is.read(line)) {
	queue.push_swap(line);
	queue.wait_empty();
	is.ready();
      }
      queue.push(std::string());
    }
    
    thread.join();
    
    // merge weights...
    
    weights_mixed = task.weights;
    if (mpi_rank == 0)
      reduce_weights(weights_mixed);
    else
      send_weights(task.weights);
    
    task.weights *= task.updated;
    task.weights -= task.weights_accumulated;
    
    weights = task.weights;
    if (mpi_rank == 0)
      reduce_weights(weights);
    else
      send_weights(task.weights);
    
    long updated_accumulated = 0;
    MPI::COMM_WORLD.Reduce(&task.updated, &updated_accumulated, 1, MPI::LONG, MPI::SUM, 0);
    
    weights_accumulated += weights;
    norm_accumulated += updated_accumulated;
    
    weights_mixed *= (1.0 / mpi_size);
    
    bcast_weights(0, weights_mixed);
    
    task.weights = weights_mixed;
    task.weights_accumulated.clear();
    task.updated = 1;
  }
  
  weights = weights_accumulated;
  weights /= norm_accumulated;
}


void reduce_weights(weight_set_type& weights)
{
  typedef utils::mpi_device_source            device_type;
  typedef boost::iostreams::filtering_istream stream_type;

  typedef boost::shared_ptr<device_type> device_ptr_type;
  typedef boost::shared_ptr<stream_type> stream_ptr_type;

  typedef std::vector<device_ptr_type, std::allocator<device_ptr_type> > device_ptr_set_type;
  typedef std::vector<stream_ptr_type, std::allocator<stream_ptr_type> > stream_ptr_set_type;

  typedef boost::tokenizer<utils::space_separator> tokenizer_type;

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
	  tokenizer_type tokenizer(line);
	  
	  tokenizer_type::iterator iter = tokenizer.begin();
	  if (iter == tokenizer.end()) continue;
	  std::string feature = *iter;
	  ++ iter;
	  if (iter == tokenizer.end()) continue;
	  std::string value = *iter;
	  
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
    if (! feature_type(id).empty() && weights[id] != 0.0)
      os << feature_type(id) << ' ' << utils::encode_base64(weights[id]) << '\n';
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
	if (feature != __empty)
	  os << feature << ' ' << utils::encode_base64(*witer) << '\n';
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
    
    ("loss-margin",   po::value<double>(&loss_margin)->default_value(loss_margin),   "loss margin for oracle forest")
    ("score-margin",  po::value<double>(&score_margin)->default_value(score_margin), "score margin for hypothesis forest")
    ("loss-scale",    po::value<double>(&loss_scale)->default_value(loss_scale),     "loss scaling")
    
    ("apply-exact", po::bool_switch(&apply_exact), "exact feature applicatin w/o pruning")
    ("cube-size",   po::value<int>(&cube_size),    "cube-pruning size")
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
