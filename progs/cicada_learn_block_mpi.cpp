//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// New: blocked learning
//
// In each iteration, we sample a block and decode, and update "support vectors" and learn.
// For learning, we employ lbfgs or liblinear for faster training (but we will use reranking training)
// We will parallelize in each block, and perform oracle computation in parallel.

//
// The idea is taken from the cicada-learn*.sh scripts, but adapted so that:
// we will discard or keep old support vectors
// we will compute error metric wrt "block"
// run in online fashion for faster convergence on larger data
// we will keep only support vectors, not actual translations
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <deque>

#include "cicada_impl.hpp"
#include "cicada_kbest_impl.hpp"
#include "cicada_text_impl.hpp"

#include "cicada/semiring.hpp"
#include "cicada/eval.hpp"

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
#include "utils/sgi_hash_set.hpp"

#include "lbfgs.h"
#include "liblinear/linear.h"

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;
typedef std::vector<std::string, std::allocator<std::string> > sample_set_type;

path_type input_file;

bool input_id_mode = false;
bool input_bitext_mode = false;
bool input_sentence_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_span_mode = false;
bool input_alignment_mode = false;
bool input_dependency_mode = false;
bool input_directory_mode = false;

std::string symbol_goal         = vocab_type::S;

grammar_file_set_type grammar_files;
bool grammar_list = false;

grammar_file_set_type tree_grammar_files;
bool tree_grammar_list = false;

feature_parameter_set_type feature_parameters;
bool feature_list = false;

op_set_type ops;
bool op_list = false;

// learning options...
path_type refset_file;
path_type output_file;  // weights output
path_type weights_file;

// scorers
std::string scorer_name = "bleu:order=4,exact=true";
bool yield_sentence   = false;
bool yield_alignment  = false;
bool yield_dependnecy = false;

// learning parameters
int iteration = 10;
int block_size = 64;
int kbest_size = 1000;

// solver parameters
bool learn_lbfgs  = false;
bool learn_linear = false;
int linear_solver = L2R_L2LOSS_SVC_DUAL;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1e-3;
double eps = std::numeric_limits<double>::infinity();

// additional misc parameters...
bool mix_weights_mode     = false; // mix weights after each iteration
bool average_weights_mode = false; // dump averaged weights at the final iteration
bool merge_samples_mode   = false; // merge all the samples, instead of "replacing"
bool dump_weights_mode    = false; // dump current weights... for debugging purpose etc.

int debug = 0;

struct LearnLBFGS;
struct LearnLinear;

void options(int argc, char** argv);

void read_refset(const path_type& refset_path,
		 scorer_document_type& scorers,
		 const size_t shard_rank = 0,
		 const size_t shard_size = 0);
void read_samples(const path_type& input_path,
		  sample_set_type& samples,
		  const bool directory_mode,
		  const bool id_mode,
		  const size_t shard_rank = 0,
		  const size_t shard_size = 0);
template <typename Learner>
void cicada_learn(operation_set_type& operations,
		  const sample_set_type& samples,
		  const scorer_document_type& scorers,
		  weight_set_type& weights);
void synchronize();

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    options(argc, argv);
    
    if (feature_list) {
      if (mpi_rank == 0)
	std::cout << cicada::FeatureFunction::lists();
      return 0;
    }

    if (op_list) {
      if (mpi_rank == 0)
	std::cout << operation_set_type::lists();
      return 0;
    }

    if (grammar_list) {
      if (mpi_rank == 0)
	std::cout << grammar_type::lists();
      return 0;
    }

    if (tree_grammar_list) {
      if (mpi_rank == 0)
	std::cout << tree_grammar_type::lists();
      return 0;
    }


    // read grammars...
    grammar_type grammar(grammar_files.begin(), grammar_files.end());
    if (debug && mpi_rank == 0)
      std::cerr << "grammar: " << grammar.size() << std::endl;

    tree_grammar_type tree_grammar(tree_grammar_files.begin(), tree_grammar_files.end());
    if (mpi_rank == 0 && debug)
      std::cerr << "tree grammar: " << tree_grammar.size() << std::endl;
    
    // read features...
    model_type model;
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));
    model.initialize();

    if (debug && mpi_rank == 0)
      std::cerr << "feature functions: " << model.size() << std::endl;
    
    operation_set_type operations(ops.begin(), ops.end(),
				  model,
				  grammar,
				  tree_grammar,
				  symbol_goal,
				  true,
				  input_sentence_mode,
				  input_lattice_mode,
				  input_forest_mode,
				  input_span_mode,
				  input_alignment_mode,
				  input_dependency_mode,
				  input_bitext_mode,
				  true,
				  debug);
    
    // make sure to synchronize here... otherwise, badthink may happen...
    if (mpi_rank == 0 && ! operations.get_output_data().directory.empty()) {
      const path_type& directory = operations.get_output_data().directory;
      
      if (boost::filesystem::exists(directory) && ! boost::filesystem::is_directory(directory))
	utils::filesystem::remove_all(directory);
      
      boost::filesystem::create_directories(directory);
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(directory); iter != iter_end; ++ iter)
	utils::filesystem::remove_all(*iter);
    }
    
    synchronize();
    
    ::sync();
    
    // read input data
    sample_set_type samples;
    read_samples(input_file, samples, input_directory_mode, input_id_mode, mpi_rank, mpi_size);
    
    // read reference data
    scorer_document_type scorers(scorer_name);
    read_refset(refset_file, scorers, mpi_rank, mpi_size);
    
    if (scorers.size() != samples.size())
      throw std::runtime_error("training sample size and reference translation size does not match");
    
    // weights
    weight_set_type weights;
    if (! weights_file.empty()) {
      if (weights_file != "-" && ! boost::filesystem::exists(weights_file))
	throw std::runtime_error("no weights file? " + weights_file.string());
      
      if (mpi_rank == 0) {
	utils::compress_istream is(weights_file, 1024 * 1024);
	is >> weights;
      }
    }
    
    // perform learning...
    if (learn_lbfgs)
      cicada_learn<LearnLBFGS>(operations, samples, scorers, weights);
    else
      cicada_learn<LearnLinear>(operations, samples, scorers, weights);
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_file, 1024 * 1024);
      os << weights;
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
}

enum {
  sample_tag = 1000,
  result_tag,
  notify_tag,
  stat_tag,
};

inline
int loop_sleep(bool found, int non_found_iter)
{
  if (! found) {
    boost::thread::yield();
    ++ non_found_iter;
  } else
    non_found_iter = 0;
    
  if (non_found_iter >= 50) {
    struct timespec tm;
    tm.tv_sec = 0;
    tm.tv_nsec = 2000001;
    nanosleep(&tm, NULL);
      
    non_found_iter = 0;
  }
  return non_found_iter;
}


void synchronize()
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  if (mpi_rank == 0) {
    std::vector<MPI::Request, std::allocator<MPI::Request> > request_recv(mpi_size);
    std::vector<MPI::Request, std::allocator<MPI::Request> > request_send(mpi_size);
    std::vector<bool, std::allocator<bool> > terminated_recv(mpi_size, false);
    std::vector<bool, std::allocator<bool> > terminated_send(mpi_size, false);
    
    terminated_recv[0] = true;
    terminated_send[0] = true;
    for (int rank = 1; rank != mpi_size; ++ rank) {
      request_recv[rank] = MPI::COMM_WORLD.Irecv(0, 0, MPI::INT, rank, notify_tag);
      request_send[rank] = MPI::COMM_WORLD.Isend(0, 0, MPI::INT, rank, notify_tag);
    }
    
    int non_found_iter = 0;
    for (;;) {
      bool found = false;
      
      for (int rank = 1; rank != mpi_size; ++ rank)
	if (! terminated_recv[rank] && request_recv[rank].Test()) {
	  terminated_recv[rank] = true;
	  found = true;
	}
      
      for (int rank = 1; rank != mpi_size; ++ rank)
	if (! terminated_send[rank] && request_send[rank].Test()) {
	  terminated_send[rank] = true;
	  found = true;
	}
      
      if (std::count(terminated_send.begin(), terminated_send.end(), true) == mpi_size
	  && std::count(terminated_recv.begin(), terminated_recv.end(), true) == mpi_size) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  } else {
    MPI::Request request_send = MPI::COMM_WORLD.Isend(0, 0, MPI::INT, 0, notify_tag);
    MPI::Request request_recv = MPI::COMM_WORLD.Irecv(0, 0, MPI::INT, 0, notify_tag);
    
    bool terminated_send = false;
    bool terminated_recv = false;
    
    int non_found_iter = 0;
    for (;;) {
      bool found = false;
      
      if (! terminated_send && request_send.Test()) {
	terminated_send = true;
	found = true;
      }
      
      if (! terminated_recv && request_recv.Test()) {
	terminated_recv = true;
	found = true;
      }
      
      if (terminated_send && terminated_recv) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
}

void read_refset(const path_type& refset_path,
		 scorer_document_type& scorers,
		 const size_t shard_rank = 0,
		 const size_t shard_size = 0)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  parser_type parser;
  id_sentence_type id_sentence;

  utils::compress_istream is(refset_path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iter_type iter(is);
  iter_type iter_end;
  
  while (iter != iter_end) {
    id_sentence.second.clear();
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
      if (iter != iter_end)
	throw std::runtime_error("refset parsing failed");
    
    const size_t id = id_sentence.first;
    
    if (shard_size && (id % shard_size != shard_rank)) continue;
    
    const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
    
    if (id_rank >= scorers.size())
      scorers.resize(id_rank + 1);
    
    if (! scorers[id_rank])
      scorers[id_rank] = scorers.create();
    
    scorers[id_rank]->insert(id_sentence.second);
  }
}

void read_samples(const path_type& input_path,
		  sample_set_type& samples,
		  const bool directory_mode,
		  const bool id_mode,
		  const size_t shard_rank = 0,
		  const size_t shard_size = 0)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;

  if (directory_mode) {
    if (! boost::filesystem::is_directory(input_path))
      throw std::runtime_error("input is not directory! " + input_path.string());

    boost::spirit::qi::uint_parser<size_t, 10, 1, -1> id_parser;
    
    size_t id;
    std::string line;
    for (size_t i = 0; /**/; ++ i)
      if (shard_size <= 0 || i % shard_size == shard_rank) {
	const path_type path = input_path / (utils::lexical_cast<std::string>(i) + ".gz");
	
	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	std::getline(is, line);
	
	if (line.empty()) continue;
	
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (id != i)
	  throw std::runtime_error("id doest not match!");
	
	const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	
	if (id_rank >= samples.size())
	  samples.resize(id_rank + 1);
	
	samples[id_rank] = line;
      }
  } else if (id_mode) {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    boost::spirit::qi::uint_parser<size_t, 10, 1, -1> id_parser;
    
    size_t id;
    std::string line;
    while (std::getline(is, line)) 
      if (! line.empty()) {
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (shard_size == 0 || id % shard_size == shard_rank) {
	  const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	  
	  if (id_rank >= samples.size())
	    samples.resize(id_rank + 1);
	  
	  samples[id_rank] = line;
	}
      }
  } else {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    std::string line;
    for (size_t id = 0; std::getine(is, line); ++ id) 
      if (shard_size == 0 || id % shard_size == shard_rank) {
	if (! line.empty())
	  samples.push_back(utils::lexical_cast<std::string>(id) + " ||| " + line);
	else
	  samples.push_back(std::string());
      }
  }
}

// LBFGS learner
struct LearnLBFGS
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef hypothesis_type::feature_set_type feature_set_type;
  typedef std::vector<feature_set_type, std::allocator<feature_set_type> > sample_type;

  struct sample_pair_type
  {
    sample_pair_type() : kbests(), oracles() {}
    
    sample_type kbests;
    sample_type oracles;
  };
  
  typedef std::vector<sample_pair_type, std::allocator<sample_pair_type> > sample_pair_set_type;
  typedef utils::chunk_vector<sample_pair_set_type, 4096 / sizeof(sample_pair_set_type), std::allocator<sample_pair_set_type> > sample_pair_map_type;
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool merge=false)
  {
    if (id >= samples.size())
      samples.resize(id + 1);
    
    if (! merge)
      samples[id].clear();
    
    samples[id].push_back(sample_pair_type());
    
    samples[id].back().kbests.reserve(kbests.size());
    samples[id].back().oracles.reserve(oracles.size());
    
    hypothesis_set_type::const_iterator kiter_end = kbests.end();
    for (hypothesis_set_type::const_iterator kiter = kbests.begin(); kiter != kiter_end; ++ kiter)
      samples[id].back().kbests.push_back(kiter->features);
    
    hypothesis_set_type::const_iterator oiter_end = oracles.end();
    for (hypothesis_set_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter)
      samples[id].back().oracles.push_back(oiter->features);
  }
  
  double learn(weight_set_type& __weights)
  {
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);
    
    if (regularize_l1) {
      param.orthantwise_c = C;
      param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;
    } else
      param.orthantwise_c = 0.0;
    
    double objective = 0.0;
    
    weights = __weights;
    weights.allocate();
    
    lbfgs(weights.size(), &(*weights.begin()), &objective, LearnLBFGS::evaluate, 0, this, &param);
    
    __weights = weights;
    
    return objective;
  }

  typedef cicada::semiring::Log<double> weight_type;
  typedef cicada::WeightVector<weight_type, std::allocator<weight_type> > expectation_type;
  
  static lbfgsfloatval_t evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
  {
    LearnLBFGS& optimizer = *((LearnLBFGS*) instance);
    
    expectation_type& expectations = optimizer.expectations;
    weight_set_type& weights       = optimizer.weights;
    sample_pair_map_type& samples  = optimizer.samples;
    
    expectations.clear();
    expectations.allocate();

    double objective = 0.0;
    size_t instances = 0;
    
    for (size_t id = 0; id != samples.size(); ++ id) 
      if (! samples[id].empty()) {
	sample_pair_set_type::const_iterator siter_end = samples[id].end();
	for (sample_pair_set_type::const_iterator siter = samples[id].begin(); siter != siter_end; ++ siter) {
	  const sample_pair_type& sample = *siter;
	  
	  weight_type Z_oracle;
	  weight_type Z_kbest;
	  
	  sample_type::const_iterator oiter_end = sample.oracles.end();
	  for (sample_type::const_iterator oiter = sample.oracles.begin(); oiter != oiter_end; ++ oiter)
	    Z_oracle += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->begin(), oiter->end(), 0.0));
	  
	  sample_type::const_iterator kiter_end = sample.kbests.end();
	  for (sample_type::const_iterator kiter = sample.kbests.begin(); kiter != kiter_end; ++ kiter)
	    Z_kbest += cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->begin(), kiter->end(), 0.0));
	  
	  for (sample_type::const_iterator oiter = sample.oracles.begin(); oiter != oiter_end; ++ oiter) {
	    const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, oiter->begin(), oiter->end(), 0.0)) / Z_oracle;
	    
	    hypothesis_type::feature_set_type::const_iterator fiter_end = oiter->end();
	    for (hypothesis_type::feature_set_type::const_iterator fiter = oiter->begin(); fiter != fiter_end; ++ fiter)
	      expectations[fiter->first] -= weight_type(fiter->second) * weight;
	  }
	  
	  for (sample_type::const_iterator kiter = sample.kbests.begin(); kiter != kiter_end; ++ kiter) {
	    const weight_type weight = cicada::semiring::traits<weight_type>::exp(cicada::dot_product(weights, kiter->begin(), kiter->end(), 0.0)) / Z_kbest;
	    
	    hypothesis_type::feature_set_type::const_iterator fiter_end = kiter->end();
	    for (hypothesis_type::feature_set_type::const_iterator fiter = kiter->begin(); fiter != fiter_end; ++ fiter)
	      expectations[fiter->first] += weight_type(fiter->second) * weight;
	  }
	  
	  const double margin = log(Z_oracle) - log(Z_kbest);
	  objective -= margin;
	  ++ instances;
	}
      }
    
    std::copy(expectations.begin(), expectations.begin() + n, g);
    
    objective /= instances;
    std::transform(g, g + n, g, std::bind2nd(std::multiplies<double>(), 1.0 / instances));
    
    // L2...
    if (regularize_l2) {
      double norm = 0.0;
      for (int i = 0; i < n; ++ i) {
	g[i] += C * x[i];
	norm += x[i] * x[i];
      }
      objective += 0.5 * C * norm;
    }
    
    return objective;
  }
  
  expectation_type     expectations;
  weight_set_type      weights;
  sample_pair_map_type samples;
};

// linear learner
struct LearnLinear
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef struct model        model_type;
  typedef struct parameter    parameter_type;
  typedef struct problem      problem_type;
  typedef struct feature_node feature_node_type;
  
  typedef size_t offset_type;
  
  typedef std::vector<feature_node_type*, std::allocator<feature_node_type*> > feature_node_map_type;
  typedef std::vector<int, std::allocator<int> > label_set_type;

  //
  // typedef for unique sentences
  //
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
  
  static void print_string_stderr(const char *s)
  {
    std::cerr << s << std::flush;
  }

  static void print_string_none(const char *s)
  {
    
  }

  struct Encoder
  {
    typedef std::vector<offset_type, std::allocator<offset_type> > offset_set_type;
    typedef std::vector<feature_node_type, std::allocator<feature_node_type> > feature_node_set_type;

    Encoder() : offsets(), features() {}

    offset_set_type       offsets;
    feature_node_set_type features;

    void encode(const hypothesis_set_type& kbests,
		const hypothesis_set_type& oracles)
    {
      feature_node_type feature;
      
      sentence_unique_type sentences;
      for (size_t o = 0; o != oracles.size(); ++ o)
	sentences.insert(oracles[o].sentence);
      
      for (size_t o = 0; o != oracles.size(); ++ o)
	for (size_t k = 0; k != kbests.size(); ++ k) {
	  const hypothesis_type& oracle = oracles[o];
	  const hypothesis_type& kbest  = kbests[k];
	  
	  if (sentences.find(kbest.sentence) != sentences.end()) continue;
	  
	  offsets.push_back(features.size());
	  
	  hypothesis_type::feature_set_type::const_iterator oiter = oracle.features.begin();
	  hypothesis_type::feature_set_type::const_iterator oiter_end = oracle.features.end();
	
	  hypothesis_type::feature_set_type::const_iterator kiter = kbest.features.begin();
	  hypothesis_type::feature_set_type::const_iterator kiter_end = kbest.features.end();
	
	  while (oiter != oiter_end && kiter != kiter_end) {
	    if (oiter->first < kiter->first) {
	      feature.index = oiter->first.id() + 1;
	      feature.value = oiter->second;
	      features.push_back(feature);
	      ++ oiter;
	    } else if (kiter->first < oiter->first) {
	      feature.index = kiter->first.id() + 1;
	      feature.value = - kiter->second;
	      features.push_back(feature);
	      ++ kiter;
	    } else {
	      feature.index = oiter->first.id() + 1;
	      feature.value = oiter->second - kiter->second;
	      if (feature.value != 0.0)
		features.push_back(feature);
	      ++ oiter;
	      ++ kiter;
	    }
	  }
	    
	  for (/**/; oiter != oiter_end; ++ oiter) {
	    feature.index = oiter->first.id() + 1;
	    feature.value = oiter->second;
	    features.push_back(feature);
	  }
	  
	  for (/**/; kiter != kiter_end; ++ kiter) {
	    feature.index = kiter->first.id() + 1;
	    feature.value = - kiter->second;
	    features.push_back(feature);
	  }
	  
	  // termination...
	  feature.index = -1;
	  feature.value = 0.0;
	  features.push_back(feature);
	}

      shrink();
    }
    
    void clear()
    {
      offsets.clear();
      features.clear();
    }
    
    void shrink()
    {
      offset_set_type(offsets).swap(offsets);
      feature_node_set_type(features).swap(features);
    }
  };
  typedef Encoder encoder_type;
  typedef utils::chunk_vector<encoder_type, 4096 / sizeof(encoder_type), std::allocator<encoder_type> > encoder_set_type; 
  
  void encode(const size_type id, const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, const bool merge=false)
  {
    if (id >= encoders.size())
      encoders.resize(id + 1);
    
    if (! merge)
      encoders[id].clear();
    
    encoders[id].encode(kbests, oracles);
  }

  double learn(weight_set_type& weights)
  {
    size_type data_size = 0;
    for (size_type id = 0; id != encoders.size(); ++ id)
      data_size += encoders[id].size();
    
    label_set_type        labels(data_size, 1);
    feature_node_map_type features;
    features.reserve(data_size);
    
    for (size_type id = 0; id != encoders.size(); ++ id)
      for (size_type pos = 0; pos != encoders[id].offsets.size(); ++ pos)
	features.push_back(const_cast<feature_node_type*>(&(*encoders[id].features.begin())) + encoders[id].offsets[pos]);
    
    problem_type problem;
    problem.l = labels.size();
    problem.n = feature_type::allocated();
    problem.y = &(*labels.begin());
    problem.x = &(*features.begin());
    problem.bias = -1;
    
    parameter_type parameter;
    parameter.solver_type = linear_solver;
    parameter.eps = eps;
    parameter.C = 1.0 / (C * labels.size()); // renormalize!
    parameter.nr_weight    = 0;
    parameter.weight_label = 0;
    parameter.weight       = 0;
    
    if (parameter.eps == std::numeric_limits<double>::infinity()) {
      if (parameter.solver_type == L2R_LR || parameter.solver_type == L2R_L2LOSS_SVC)
	parameter.eps = 0.01;
      else if (parameter.solver_type == L2R_L2LOSS_SVC_DUAL || parameter.solver_type == L2R_L1LOSS_SVC_DUAL || parameter.solver_type == MCSVM_CS || parameter.solver_type == L2R_LR_DUAL)
	parameter.eps = 0.1;
      else if (parameter.solver_type == L1R_L2LOSS_SVC || parameter.solver_type == L1R_LR)
	parameter.eps = 0.01;
    }
    
    if (debug >= 2)
      set_print_string_function(print_string_stderr);
    else
      set_print_string_function(print_string_none);
    
    const char* error_message = check_parameter(&problem, &parameter);
    if (error_message)
      throw std::runtime_error(std::string("error: ") + error_message);
    
    static const char* names[] = {"L2R_LR", "L2R_L2LOSS_SVC_DUAL", "L2R_L2LOSS_SVC", "L2R_L1LOSS_SVC_DUAL", "MCSVM_CS",
				  "L1R_L2LOSS_SVC", "L1R_LR", "L2R_LR_DUAL"};
    
    if (debug)
      std::cerr << "solver: " << names[parameter.solver_type] << std::endl;
    
    const model_type* model = train(&problem, &parameter);
    
    const double objective = model->objective * C;
    
    // it is an optimization...
    weights.clear();
    for (int j = 0; j != model->nr_feature; ++ j)
      weights[weight_set_type::feature_type(j)] = model->w[j];
    
    free_and_destroy_model(const_cast<model_type**>(&model));
    
    return objective;
  }
  
  encoder_set_type encoders;
};

template <typename Learner>
void cicada_learn(operation_set_type& operations,
		  const sample_set_type& samples,
		  const scorer_document_type& scorers,
		  weight_set_type& weights)
{
  typedef std::vector<size_t, std::allocator<size_t> > segment_set_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  Learner         learner;
  KBestGenerator  kbest_generator;
  OracleGenerator oracle_generator;
  
  segment_set_type segments(samples.size());
  for (size_t seg = 0; seg != segments.size(); ++ seg)
    segments[seg] = seg;
  
  // random number generator...
  boost::mt19937 generator;
  generator.seed(time(0) * getpid());
  boost::random_number_generator<generator_type> gen(generator);
  
  weight_set_type weights_total;
  int             learned_total = 0;

  hypothesis_set_type kbests;
  
  segment_set_type     segments_block;
  hypothesis_map_type  kbests_block;
  hypothesis_map_type  oracles_block;
  scorer_document_type scorers_block(scorers);
  
  for (int iter = 0; iter != iteration; ++ iter) {
    // perform learning
    
    segment_set_type::const_iterator siter     = segments.begin();
    segment_set_type::const_iterator siter_end = segments.end();
    
    while (siter != siter_end) {
      segments_block.clear();
      kbests_block.clear();
      oracles_block.clear();
      scorers_block.clear();
      
      segment_set_type::const_iterator siter_last = std::min(siter + block_size, siter_end);
      for (/**/; siter != siter_last; ++ siter) {
	const size_t id = *siter;
	
	if (samples[id].empty() || ! scorers[id]) continue;
	
	kbest_generator(operations, samples[id], kbests);
	
	// no translations?
	if (kbests.empty()) continue;
	
	segments_block.push_back(id);
	kbests_block.push_back(kbests);
	scorers_block.push_back(scorers[id]);
      }
      
      // no kbests?
      if (segments_block.empty()) continue;
      
      // oracle computation
      oracle_generator(kbests_block, scorers_block, oracles_block);
      
      // encode into learner...
      for (size_t i = 0; i != kbests_block.size(); ++ i)
	learner.encode(segments_block[i], kbests_block[i], oracles_block[i], merge_samples_mode);
      
      // perform learning...
      learner.learn(weights);
      
      // keep totals...
      weights_total += weights;
      ++ learned_total;
    }
    
    // randomize..
    std::random_shuffle(segments.begin(), segments.end(), gen);
    
    if (mix_weights_mode) {
      
    }
    
    if (dump_weights_mode && mpi_rank == 0) {
      if (average_weights_mode) {
	
      } else {
	
      }
    }
  }
  
  if (average_weights_mode) {
    
  } else {
    
  }
}

struct deprecated
{
  deprecated(const boost::program_options::options_description& __desc)
    : desc(__desc) {}
  
  template <typename Tp>
  void operator()(const Tp& x) const
  {
    std::cout << desc << std::endl;
    exit(1);
  }
  
  const boost::program_options::options_description& desc;
};

void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    
    // options for input/output format
    ("input-id",         po::bool_switch(&input_id_mode),         "id-prefixed input")
    ("input-bitext",     po::bool_switch(&input_bitext_mode),     "target sentence prefixed input")
    ("input-sentence",   po::bool_switch(&input_sentence_mode),   "sentence input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-span",       po::bool_switch(&input_span_mode),       "span input")
    ("input-alignment",  po::bool_switch(&input_alignment_mode),  "alignment input")
    ("input-dependency", po::bool_switch(&input_dependency_mode), "dependency input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")
    
    // grammar
    ("goal",              po::value<std::string>(&symbol_goal)->default_value(symbol_goal),    "goal symbol")
    ("grammar",           po::value<grammar_file_set_type >(&grammar_files)->composing(),      "grammar specification(s)")
    ("grammar-list",      po::bool_switch(&grammar_list),                                      "list of available grammar specifications")
    ("tree-grammar",      po::value<grammar_file_set_type >(&tree_grammar_files)->composing(), "tree grammar specification(s)")
    ("tree-grammar-list", po::bool_switch(&tree_grammar_list),                                 "list of available grammar specifications")
    
    // models...
    ("feature-function",      po::value<feature_parameter_set_type >(&feature_parameters)->composing(), "feature function(s)")
    ("feature-function-list", po::bool_switch(&feature_list),                                           "list of available feature function(s)")
    
    //operatins...
    ("operation",      po::value<op_set_type>(&ops)->composing(), "operations")
    ("operation-list", po::bool_switch(&op_list),                 "list of available operation(s)");
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config", po::value<path_type>(), "configuration file")
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

    po::options_description opts_deprecated("deprecated options");
  opts_deprecated.add_options()
    ("non-terminal",          po::value<std::string>()->notifier(deprecated(opts_deprecated)), "see --grammar-list")
    ("grammar-static",        po::value<grammar_file_set_type >()->composing()->notifier(deprecated(opts_deprecated)), "use --grammar ")
    ("grammar-glue-straight", po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar glue:straight=[true|false],inverted=[true|false],non-terminal=[x]")
    ("grammar-glue-inverted", po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar glue:straight=[true|false],inverted=[true|false],non-terminal=[x]")
    ("grammar-insertion",     po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar insetion:non-terminal=[x]")
    ("grammar-deletion",      po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar deletion:non-terminal=[x]")
    ("tree-grammar-static",   po::value<grammar_file_set_type >()->composing()->notifier(deprecated(opts_deprecated)),  "use --tree-grammar")
    ("tree-grammar-fallback", po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --tree-grammar fallback:non-terminal=[x]");
  
  po::options_description desc_config;
  po::options_description desc_command;
  po::options_description desc_visible;

  desc_config.add(opts_config).add(opts_deprecated);
  desc_command.add(opts_config).add(opts_command).add(opts_deprecated);
  desc_visible.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  if (variables.count("config")) {
    const path_type path_config = variables["config"].as<path_type>();
    if (! boost::filesystem::exists(path_config))
      throw std::runtime_error("no config file: " + path_config.string());
	
    utils::compress_istream is(path_config);
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {

    if (mpi_rank == 0)
      std::cout << argv[0] << " [options]\n"
		<< desc_visible << std::endl;

    MPI::Finalize();
    exit(0);
  }
}

