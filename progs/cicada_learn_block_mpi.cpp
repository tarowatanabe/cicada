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
#include <sstream>
#include <cstdlib>

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
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>

#include "lbfgs.h"
#include "liblinear/linear.h"

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

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
bool yield_dependency = false;

// learning parameters
int iteration = 10;
int block_size = 64;
int kbest_size = 1000;

// solver parameters
bool learn_lbfgs  = false;
bool learn_pegasos = false;
bool learn_opegasos = false;
bool learn_mira   = false;
bool learn_sgd    = false;
bool learn_svm    = false;
bool learn_linear = false;
int linear_solver = L2R_L2LOSS_SVC_DUAL;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1e-3;
double eps = std::numeric_limits<double>::infinity();

// additional misc parameters...
bool loss_rank = false; // loss by rank
bool softmax_margin = false;
bool merge_vectors_mode  = false; // merge all the vectors from others
bool line_search_mode = false;    // perform line-search
bool dump_weights_mode   = false; // dump current weights... for debugging purpose etc.


int debug = 0;

#include "cicada_learn_block_impl.hpp"

// forward declarations...

void options(int argc, char** argv);

template <typename Learner, typename KBestGenerator, typename OracleGenerator>
void cicada_learn(operation_set_type& operations,
		  const sample_set_type& samples,
		  const scorer_document_type& scorers,
		  weight_set_type& weights);
void synchronize();

void bcast_weights(weight_set_type& weights);
void reduce_weights(weight_set_type& weights);
void reduce_score(score_ptr_type& score, const int tag);

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


    // check parameters for learning...
    // 
    //
    if (refset_file.empty() || ! boost::filesystem::exists(refset_file))
      throw std::runtime_error("no refset file? " + refset_file.string());
    
    if (output_file.empty())
      throw std::runtime_error("no output?");
    
    if (int(yield_sentence) + yield_alignment + yield_dependency > 1)
      throw std::runtime_error("you can specify either --yield-{sentence,alignment,dependency}");
    if (int(yield_sentence) + yield_alignment + yield_dependency == 0)
      yield_sentence = true;
    
    if (int(learn_lbfgs) + learn_mira + learn_sgd + learn_linear + learn_svm + learn_pegasos + learn_opegasos> 1)
      throw std::runtime_error("you can specify either --learn-{lbfgs,mira,sgd,linear,svm}");
    if (int(learn_lbfgs) + learn_mira + learn_sgd + learn_linear + learn_svm + learn_pegasos + learn_opegasos == 0)
      learn_lbfgs = true;
    
    if (int(regularize_l1) + regularize_l2 > 1)
      throw std::runtime_error("either L1 or L2 regularization");
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));

    if (block_size <= 0)
      throw std::runtime_error("block size must be possitive: " + utils::lexical_cast<std::string>(block_size));
    if (kbest_size <= 0)
      throw std::runtime_error("kbest size must be possitive: " + utils::lexical_cast<std::string>(kbest_size));
    
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
				  debug ? debug - 1 : debug);
    
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
    if (yield_sentence) {
      if (learn_lbfgs)
	cicada_learn<LearnLBFGS, KBestSentence, Oracle>(operations, samples, scorers, weights);
      else if (learn_mira)
	cicada_learn<LearnMIRA, KBestSentence, Oracle>(operations, samples, scorers, weights);
      else if (learn_sgd && regularize_l1)
	cicada_learn<LearnSGDL1, KBestSentence, Oracle>(operations, samples, scorers, weights);
      else if (learn_sgd && regularize_l2)
	cicada_learn<LearnSGDL2, KBestSentence, Oracle>(operations, samples, scorers, weights);      
      else if (learn_linear)
	cicada_learn<LearnLinear, KBestSentence, Oracle>(operations, samples, scorers, weights);
      else if (learn_svm)
	cicada_learn<LearnSVM, KBestSentence, Oracle>(operations, samples, scorers, weights);
      else if (learn_pegasos)
	cicada_learn<LearnPegasos, KBestSentence, Oracle>(operations, samples, scorers, weights);
      else
	cicada_learn<LearnOPegasos, KBestSentence, Oracle>(operations, samples, scorers, weights);
    } else if (yield_alignment) {
      if (learn_lbfgs)
	cicada_learn<LearnLBFGS, KBestAlignment, Oracle>(operations, samples, scorers, weights);
      else if (learn_mira)
	cicada_learn<LearnMIRA, KBestAlignment, Oracle>(operations, samples, scorers, weights);
      else if (learn_sgd && regularize_l1)
	cicada_learn<LearnSGDL1, KBestAlignment, Oracle>(operations, samples, scorers, weights);
      else if (learn_sgd && regularize_l2)
	cicada_learn<LearnSGDL2, KBestAlignment, Oracle>(operations, samples, scorers, weights);      
      else if (learn_linear)
	cicada_learn<LearnLinear, KBestAlignment, Oracle>(operations, samples, scorers, weights);
      else if (learn_svm)
	cicada_learn<LearnSVM, KBestAlignment, Oracle>(operations, samples, scorers, weights);
      else if (learn_pegasos)
	cicada_learn<LearnPegasos, KBestAlignment, Oracle>(operations, samples, scorers, weights);
      else
	cicada_learn<LearnOPegasos, KBestAlignment, Oracle>(operations, samples, scorers, weights);
    } else if (yield_dependency) {
      if (learn_lbfgs)
	cicada_learn<LearnLBFGS, KBestDependency, Oracle>(operations, samples, scorers, weights);
      else if (learn_mira)
	cicada_learn<LearnMIRA, KBestDependency, Oracle>(operations, samples, scorers, weights);
      else if (learn_sgd && regularize_l1)
	cicada_learn<LearnSGDL1, KBestDependency, Oracle>(operations, samples, scorers, weights);
      else if (learn_sgd && regularize_l2)
	cicada_learn<LearnSGDL2, KBestDependency, Oracle>(operations, samples, scorers, weights);      
      else if (learn_linear)
	cicada_learn<LearnLinear, KBestDependency, Oracle>(operations, samples, scorers, weights);
      else if (learn_svm)
	cicada_learn<LearnSVM, KBestDependency, Oracle>(operations, samples, scorers, weights);
      else if (learn_pegasos)
	cicada_learn<LearnPegasos, KBestDependency, Oracle>(operations, samples, scorers, weights);
      else
	cicada_learn<LearnOPegasos, KBestDependency, Oracle>(operations, samples, scorers, weights);
    } else
      throw std::runtime_error("invalid yield");
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_file, 1024 * 1024);
      os.precision(20);
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
  weights_tag = 1000,
  score_1best_tag,
  score_oracle_tag,
  notify_tag,
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


struct Dumper
{
  typedef std::pair<path_type, weight_set_type > value_type;
  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;
  
  Dumper(queue_type& __queue)
    : queue(__queue) {}
  
  void operator()()
  {
    value_type value;
    
    while (1) {
      queue.pop_swap(value);
      if (value.first.empty()) break;
      
      utils::compress_ostream os(value.first, 1024 * 1024);
      os.precision(20);
      os << value.second;
    }
  }

  queue_type& queue;
};

template <typename Learner, typename KBestGenerator, typename OracleGenerator>
void cicada_learn(operation_set_type& operations,
		  const sample_set_type& samples,
		  const scorer_document_type& scorers,
		  weight_set_type& weights)
{
  typedef std::pair<double, double> point_type;
  typedef std::vector<point_type, std::allocator<point_type> > point_set_type;
  
  typedef std::vector<size_t, std::allocator<size_t> > segment_set_type;
  
  typedef Dumper dumper_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  Learner         learner(samples.size());
  KBestGenerator  kbest_generator;
  OracleGenerator oracle_generator;
  
  segment_set_type segments(samples.size());
  for (size_t seg = 0; seg != segments.size(); ++ seg)
    segments[seg] = seg;
  
  // random number generator...
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> gen(generator);
  
  hypothesis_set_type kbests;
  
  segment_set_type     segments_block;
  hypothesis_map_type  kbests_block;
  hypothesis_map_type  oracles_block;
  scorer_document_type scorers_block(scorers);

  const bool error_metric = scorers.error_metric();
  
  dumper_type::queue_type queue_dumper;
  std::auto_ptr<boost::thread> thread_dumper(mpi_rank == 0 ? new boost::thread(dumper_type(queue_dumper)) : 0);
  
  // first, bcast weights...
  bcast_weights(weights);

  // previous weights...
  weight_set_type weights_prev;
  point_set_type points;
  point_set_type points_next;
    
  // start training
  for (int iter = 0; iter != iteration; ++ iter) {
    // perform learning

    if (debug && mpi_rank == 0)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    int updated = 0;
    score_ptr_type score_1best;
    score_ptr_type score_oracle;

    weights_prev = weights;

    learner.initialize(weights);
    
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
	
	kbest_generator(operations, samples[id], weights, kbests);
	
	// no translations?
	if (kbests.empty()) continue;
	
	segments_block.push_back(id);
	kbests_block.push_back(kbests);
	scorers_block.push_back(scorers[id]);
      }
      
      // no kbests?
      if (segments_block.empty()) continue;
      
      // oracle computation
      std::pair<score_ptr_type, score_ptr_type> scores = oracle_generator(kbests_block, scorers_block, oracles_block, generator);


      if (! score_1best)
	score_1best = scores.first;
      else
	*score_1best += *(scores.first);
      
      if (! score_oracle)
	score_oracle = scores.second;
      else
	*score_oracle += *(scores.second);
      
      if (debug >= 2)
	std::cerr << "devset block       1best: " << scores.first->score() << " oracle: " << scores.second->score() << std::endl
		  << "devset accumulated 1best: " << score_1best->score()  << " oracle: " << score_oracle->score() << std::endl;
      
      // encode into learner...
      for (size_t i = 0; i != kbests_block.size(); ++ i)
	learner.encode(segments_block[i], kbests_block[i], oracles_block[i], error_metric);
      
      // perform learning...
      learner.learn(weights);
      
      // keep totals...
      ++ updated;
    }

    if (debug) {
      reduce_score(score_1best, score_1best_tag);
      reduce_score(score_oracle, score_oracle_tag);
      
      if (mpi_rank == 0)
	std::cerr << "devset total: 1best: " << score_1best->score() << " oracle: " << score_oracle->score() << std::endl;
    }
    
    // randomize..
    std::random_shuffle(segments.begin(), segments.end(), gen);
    
    // merge vectors...
    if (merge_vectors_mode) {
      if (debug && mpi_rank == 0)
	std::cerr << "merge vectors" << std::endl;
      
      learner.clear();
      
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank == mpi_rank) {
	  boost::iostreams::filtering_ostream os;
	  os.push(boost::iostreams::zlib_compressor());
	  os.push(utils::mpi_device_bcast_sink(rank, 1024 * 1024));
	  
	  learner.encode(os);
	} else {
	  boost::iostreams::filtering_istream is;
	  is.push(boost::iostreams::zlib_decompressor());
	  is.push(utils::mpi_device_bcast_source(rank, 1024 * 1024));
	  
	  learner.decode(is);
	}
      
      if (debug && mpi_rank == 0)
	std::cerr << "learning by merged vectors" << std::endl;
      
      // perform learning...
      learner.learn(weights);
    }

    learner.finalize(weights);
    
    // mix weights
    if (debug && mpi_rank == 0)
      std::cerr << "mix weights" << std::endl;

    if (merge_vectors_mode) {
      // simply averaging, since we have distributed vectors
      reduce_weights(weights);
      
      bcast_weights(weights);
      
      weights *= 1.0 / mpi_size;
    } else {
      ++ updated; // avoid zero...
      weights *= updated;
      
      reduce_weights(weights);
      
      bcast_weights(weights);
      
      int updated_total = 0;
      MPI::COMM_WORLD.Allreduce(&updated, &updated_total, 1, MPI::INT, MPI::SUM);
      
      weights *= 1.0 / updated_total;
    }
    
    // perform line-search....
    if (line_search_mode) {
      
      points.clear();
      const std::pair<double, double> grad_norm = learner.gradient(weights, weights_prev, std::back_inserter(points));
      
      double grad = 0.0;
      MPI::COMM_WORLD.Reduce(&grad_norm.first, &grad, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      double loss_norm = 0.0;
      MPI::COMM_WORLD.Reduce(&grad_norm.second, &loss_norm, 1, MPI::DOUBLE, MPI::SUM, 0);
      
      std::sort(points.begin(), points.end());
      
      // merge points from others...
      if (mpi_rank == 0) {
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
	
	// compute here...
	const double norm_w      = cicada::dot_product(weights, weights);
	const double dot_prod    = cicada::dot_product(weights_prev, weights);
	const double norm_w_prev = cicada::dot_product(weights_prev, weights_prev);
	
	const double a0 = (norm_w - 2.0 * dot_prod + norm_w_prev) * C * loss_norm;
	const double b0 = (dot_prod - norm_w_prev) * C * loss_norm;
	
	grad += b0;
	
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
	
	  k = std::max(k, 0.0);
	}
	
	const double merge_ratio = k + 0.1 * (1.0 - k);
	
	weights *= merge_ratio;
	weights_prev *= (1.0 - merge_ratio);
	
	weights += weights_prev;
      } else {
	boost::iostreams::filtering_ostream os;
	os.push(utils::mpi_device_sink(0, point_tag, 1024 * 1024));
	
	point_set_type::const_iterator piter_end = points.end();
	for (point_set_type::const_iterator piter = points.begin(); piter != piter_end; ++ piter) {
	  os.write((char*) &(piter->first), sizeof(double));
	  os.write((char*) &(piter->second), sizeof(double));
	}
      }
      
      // receive new weights...
      bcast_weights(weights);
    }
    
    // clear history for line-search...
    learner.clear_history();
    
    // dump...
    if (dump_weights_mode && mpi_rank == 0)
      queue_dumper.push(std::make_pair(add_suffix(output_file, "." + utils::lexical_cast<std::string>(iter + 1)), weights));
  }
  
  if (mpi_rank == 0) {
    queue_dumper.push(std::make_pair(path_type(), weight_set_type()));
    thread_dumper->join();
  }
}

void reduce_score(score_ptr_type& score, const int tag)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      is.push(utils::mpi_device_source(rank, tag, 256));
      
      std::string score_str;
      is >> score_str;
      
      if (! score)
	score = scorer_type::score_type::decode(score_str);
      else
	*score += *scorer_type::score_type::decode(score_str);
    }
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(utils::mpi_device_sink(0, tag, 256));
    os << score->encode();
  }
}

void send_weights(const int rank, const weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  boost::iostreams::filtering_ostream os;
  os.push(boost::iostreams::zlib_compressor());
  os.push(utils::mpi_device_sink(rank, weights_tag, 4096));
  
  for (feature_type::id_type id = 0; id < weights.size(); ++ id)
    if (! feature_type(id).empty() && weights[id] != 0.0) {
      os << feature_type(id) << ' ';
      utils::encode_base64(weights[id], std::ostream_iterator<char>(os));
      os << '\n';
    }
}

void reduce_weights(const int rank, weight_set_type& weights)
{
  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  boost::iostreams::filtering_istream is;
  is.push(boost::iostreams::zlib_decompressor());
  is.push(utils::mpi_device_source(rank, weights_tag, 4096));
  
  std::string line;
  
  while (std::getline(is, line)) {
    const utils::piece line_piece(line);
    tokenizer_type tokenizer(line_piece);
    
    tokenizer_type::iterator iter = tokenizer.begin();
    if (iter == tokenizer.end()) continue;
    const utils::piece feature = *iter;
    ++ iter;
    if (iter == tokenizer.end()) continue;
    const utils::piece value = *iter;
    
    weights[feature] += utils::decode_base64<double>(value);
  }
}

template <typename Iterator>
void reduce_weights(Iterator first, Iterator last, weight_set_type& weights)
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

  device_ptr_set_type device;
  stream_ptr_set_type stream;
  
  for (/**/; first != last; ++ first) {
    device.push_back(device_ptr_type(new device_type(*first, weights_tag, 4096)));
    stream.push_back(stream_ptr_type(new stream_type()));
    
    stream.back()->push(boost::iostreams::zlib_decompressor());
    stream.back()->push(*device.back());
  }
  
  std::string line;
  
  int non_found_iter = 0;
  while (1) {
    bool found = false;
    
    for (size_t i = 0; i != device.size(); ++ i)
      while (stream[i] && device[i] && device[i]->test()) {
	if (std::getline(*stream[i], line)) {
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
	  stream[i].reset();
	  device[i].reset();
	}
	found = true;
      }
    
    if (std::count(device.begin(), device.end(), device_ptr_type()) == static_cast<int>(device.size())) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
}


void reduce_weights(weight_set_type& weights)
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
	reduce_weights(ranks.front(), weights);
      else
	reduce_weights(ranks.begin(), ranks.end(), weights);
      
    } else
      send_weights(mpi_rank % reduce_size, weights);
    
    merge_size = reduce_size;
  }
}

void bcast_weights(const int rank, weight_set_type& weights)
{
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

void bcast_weights(weight_set_type& weights)
{
  bcast_weights(0, weights);
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

  po::options_description opts_learn("learning options");
  opts_learn.add_options()
    ("refset",  po::value<path_type>(&refset_file),  "refset")
    ("output",  po::value<path_type>(&output_file),  "model (or weights) output")
    ("weights", po::value<path_type>(&weights_file), "initial model (or weights)")
    
    ("scorer",           po::value<std::string>(&scorer_name)->default_value(scorer_name), "evaluation scorer")
    ("yield-sentence",   po::bool_switch(&yield_sentence),                                 "sentence yield")
    ("yield-alignment",  po::bool_switch(&yield_alignment),                                "alignment yield")
    ("yield-dependency", po::bool_switch(&yield_dependency),                               "dependency yield")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration),   "learning iterations")
    ("block",     po::value<int>(&block_size)->default_value(block_size), "block (or batch, bin) size")
    ("kbest",     po::value<int>(&kbest_size)->default_value(kbest_size), "kbest size")
    
    ("learn-lbfgs",    po::bool_switch(&learn_lbfgs),    "batch LBFGS algorithm")
    ("learn-mira",     po::bool_switch(&learn_mira),     "online MIRA algorithm")
    ("learn-pegasos",  po::bool_switch(&learn_pegasos),  "online Pegasos algorithm")
    ("learn-opegasos", po::bool_switch(&learn_opegasos), "online optimized-Pegasos algorithm")
    ("learn-sgd",      po::bool_switch(&learn_sgd),      "online SGD algorithm")
    ("learn-svm",      po::bool_switch(&learn_svm),      "SVM for structured output")
    ("learn-linear",   po::bool_switch(&learn_linear),   "liblinear algorithm")
    ("solver",         po::value<int>(&linear_solver),   "liblinear solver type (default: 1)\n"
     " 0: \tL2-regularized logistic regression (primal)\n"
     " 1: \tL2-regularized L2-loss support vector classification (dual)\n"
     " 2: \tL2-regularized L2-loss support vector classification (primal)\n"
     " 3: \tL2-regularized L1-loss support vector classification (dual)\n"
     " 5: \tL1-regularized L2-loss support vector classification\n"
     " 6: \tL1-regularized logistic regression\n"
     " 7: \tL2-regularized logistic regression (dual)")
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C",             po::value<double>(&C)->default_value(C), "regularization constant")
    ("eps",           po::value<double>(&eps),                 "tolerance for liblinear")
    
    ("loss-rank",      po::bool_switch(&loss_rank),          "rank loss")
    ("softmax-margin", po::bool_switch(&softmax_margin),     "softmax margin")
    ("merge-vector",   po::bool_switch(&merge_vectors_mode), "merge vectors from others")
    ("line-search",    po::bool_switch(&line_search_mode),   "perform line search")
    ("dump-weights",   po::bool_switch(&dump_weights_mode),  "dump mode (or weights) during iterations")
    ;
    

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
  desc_command.add(opts_config).add(opts_learn).add(opts_command).add(opts_deprecated);
  desc_visible.add(opts_config).add(opts_learn).add(opts_command);
  
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

