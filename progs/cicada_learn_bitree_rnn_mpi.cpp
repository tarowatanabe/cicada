//
//  Copyright(C) 2013-2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// mini-batch kbest-based tree-rnn learner
//

//
// This implementation is motivated by (Chiang et al., 2008) and (Chiang et al., 2009)
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
#include "cicada_output_impl.hpp"

#include "cicada/eval/score.hpp"
#include "cicada/format.hpp"
#include "cicada/signature.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/tokenizer.hpp"
#include "cicada/matcher.hpp"
#include "cicada/feature/frontier_bitree_rnn.hpp"

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
#include "utils/mpi_traits.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

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
path_type output_feature;

op_set_type ops;
bool op_list = false;

bool scorer_list = false;
bool format_list = false;
bool signature_list = false;
bool stemmer_list   = false;
bool tokenizer_list = false;
bool matcher_list = false;

// learning options...
path_type refset_file;
path_type output_weights_file; // weights output
path_type output_model_file;   // model output
path_type weights_file;

// scorers
std::string scorer_name = "bleu:order=4,exact=true";

// learning parameters
int iteration = 10;
int batch_size = 8;
int kbest_size = 1000;
bool kbest_unique_mode = true;
double kbest_diversity = 0.0;

bool optimize_adagrad = false;
bool optimize_sgd = false;

double lambda = 0.0;
double eta0 = 0.1;

// additional misc parameters...
int merge_history = 0;
bool mix_none_mode = false;
bool mix_average_mode = false;
bool mix_select_mode = false;
bool dump_weights_mode   = false; // dump current weights... for debugging purpose etc.

int debug = 0;

#include "cicada_learn_bitree_rnn_impl.hpp"

// forward declarations...

void options(int argc, char** argv);
void merge_features();
void merge_statistics(const operation_set_type& operations, operation_set_type::statistics_type& statistics);

template <typename Learner>
void cicada_learn(const Learner& learner,
		  operation_set_type& operations,
		  const model_type& model,
		  const event_set_type& events,
		  const scorer_document_type& scorers,
		  weight_set_type& weights,
		  tree_rnn_type& theta);
void synchronize();

void bcast_weights(const int rank, weight_set_type& weights, tree_rnn_type& theta);
void bcast_weights(weight_set_type& weights, tree_rnn_type& theta);
void reduce_weights(weight_set_type& weights, tree_rnn_type& theta);
void reduce_score_pair(score_ptr_type& score_1best, score_ptr_type& score_oracle);

void send_weights(const int rank, const weight_set_type& weights, const tree_rnn_type& theta);
void recv_weights(const int rank, weight_set_type& weights, tree_rnn_type& theta);

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    options(argc, argv);
    
    if (feature_list || op_list || grammar_list || tree_grammar_list || scorer_list || format_list || signature_list || stemmer_list || tokenizer_list || matcher_list) {
      
      if (mpi_rank == 0) {
	if (scorer_list)
	  std::cout << cicada::eval::Scorer::lists();
	
	if (feature_list)
	  std::cout << cicada::FeatureFunction::lists();

	if (format_list)
	  std::cout << cicada::Format::lists();

	if (grammar_list)
	  std::cout << grammar_type::lists();

	if (matcher_list)
	  std::cout << cicada::Matcher::lists();
      
	if (op_list)
	  std::cout << operation_set_type::lists();

	if (signature_list)
	  std::cout << cicada::Signature::lists();
      
	if (stemmer_list)
	  std::cout << cicada::Stemmer::lists();

	if (tokenizer_list)
	  std::cout << cicada::Tokenizer::lists();
      
	if (tree_grammar_list)
	  std::cout << tree_grammar_type::lists();
      }
      
      return 0;
    }    
    
    if (boost::filesystem::exists(input_file)) {
      if (boost::filesystem::is_directory(input_file))
	input_directory_mode = true;
      else if (input_directory_mode)
	throw std::runtime_error("non directory input: " + input_file.string());
    }
    
    // check parameters for learning...
    // 
    //
    if (refset_file.empty() || ! boost::filesystem::exists(refset_file))
      throw std::runtime_error("no refset file? " + refset_file.string());
    
    if (output_weights_file.empty())
      throw std::runtime_error("no output?");

    if (output_model_file.empty())
      throw std::runtime_error("no output?");
    
    if (int(optimize_sgd) + optimize_adagrad > 1)
      throw std::runtime_error("either one of optimize-{sgd,adagrad}");
    
    if (int(optimize_sgd) + optimize_adagrad == 0)
      optimize_sgd = true;
    
    if (lambda < 0)
      throw std::runtime_error("regularization constant must be positive");

    if (eta0 < 0)
      throw std::runtime_error("learning rate constant must be positive");
    
    if (batch_size <= 0)
      throw std::runtime_error("batch size must be possitive: " + utils::lexical_cast<std::string>(batch_size));
    if (kbest_size <= 0)
      throw std::runtime_error("kbest size must be possitive: " + utils::lexical_cast<std::string>(kbest_size));
    
    if (int(mix_none_mode) + mix_average_mode + mix_select_mode > 1)
      throw std::runtime_error("you can specify only one of mix-{none,average,select}");
    
    if (int(mix_none_mode) + mix_average_mode + mix_select_mode == 0)
      mix_none_mode = true;
    
    // random number seed
    ::srandom(utils::random_seed());

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
    
    if (mpi_rank == 0 && debug)
      std::cerr << "operations: " << operations.size() << std::endl;

    // make sure to synchronize here... otherwise, badthink may happen...
    if (mpi_rank == 0 && ! operations.get_output_data().directory.empty())
      prepare_directory(operations.get_output_data().directory);
    
    synchronize();
    
    ::sync();
    
    // read input data
    event_set_type events;
    read_events(input_file, events, input_directory_mode, input_id_mode, mpi_rank, mpi_size);
    
    // read reference data
    scorer_document_type scorers(scorer_name);
    read_refset(refset_file, scorers, mpi_rank, mpi_size);
    
    if (scorers.size() != events.size())
      throw std::runtime_error("training sample size and reference translation size does not match");
    
    // weights and rnn parameters
    weight_set_type weights;
    
    if (! weights_file.empty()) {
      if (weights_file != "-" && ! boost::filesystem::exists(weights_file))
	throw std::runtime_error("no weights file? " + weights_file.string());
      
      if (mpi_rank == 0) {
	utils::compress_istream is(weights_file, 1024 * 1024);
	is >> weights;
      }
    }
        
    // check if we have correct feature!
    const cicada::feature::FrontierBiTreeRNN* tree_rnn_feature = 0;
    
    model_type::const_iterator miter_end = model.end();
    for (model_type::const_iterator miter = model.begin(); miter != miter_end; ++ miter)
      if (dynamic_cast<cicada::feature::FrontierBiTreeRNN*>(miter->get())) {
	if (tree_rnn_feature)
	  throw std::runtime_error("We do not allow multiple tree-rnn features!");
	
	tree_rnn_feature = dynamic_cast<cicada::feature::FrontierBiTreeRNN*>(miter->get());
      }
    
    if (! tree_rnn_feature)
      throw std::runtime_error("no tree-rnn feature?");
    
    tree_rnn_type& theta = tree_rnn_feature->model();
    
    // perform learning...
    if (optimize_adagrad)
      cicada_learn(LearnAdaGrad(theta, lambda, eta0), operations, model, events, scorers, weights, theta);
    else 
      cicada_learn(LearnSGD(theta, lambda, eta0), operations, model, events, scorers, weights, theta);
    
    // output model...
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_weights_file, 1024 * 1024);
      os.precision(20);
      os << weights;
      
      theta.write(output_model_file);
    }

    synchronize();
    
    ::sync();
    
    operation_set_type::statistics_type statistics;
    merge_statistics(operations, statistics);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "statistics"<< '\n'
		<< statistics;
    
    if (! output_feature.empty()) {
      merge_features();
      
      if (mpi_rank == 0) {
	utils::compress_ostream os(output_feature, 1024 * 1024);
	
	for (feature_type::id_type id = 0; id != feature_type::allocated(); ++ id)
	  if (! feature_type(id).empty())
	    os << feature_type(id) << '\n';
      }
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
  score_tag,
  notify_tag,
  update_tag,
  feature_tag,
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

void merge_features()
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    typedef utils::mpi_device_source            device_type;
    typedef boost::iostreams::filtering_istream stream_type;
    
    typedef boost::shared_ptr<device_type> device_ptr_type;
    typedef boost::shared_ptr<stream_type> stream_ptr_type;
    
    typedef std::vector<device_ptr_type, std::allocator<device_ptr_type> > device_ptr_set_type;
    typedef std::vector<stream_ptr_type, std::allocator<stream_ptr_type> > stream_ptr_set_type;
    
    device_ptr_set_type device(mpi_size);
    stream_ptr_set_type stream(mpi_size);
    
    for (int rank = 1; rank != mpi_size; ++ rank) {
      device[rank].reset(new device_type(rank, feature_tag, 1024 * 1024));
      stream[rank].reset(new stream_type());
      
      stream[rank]->push(boost::iostreams::zlib_decompressor());
      stream[rank]->push(*device[rank]);
    }
    
    std::string line;
    
    int non_found_iter = 0;
    while (1) {
      bool found = false;
      
      for (int rank = 1; rank != mpi_size; ++ rank)
	while (stream[rank] && device[rank] && device[rank]->test()) {
	  if (std::getline(*stream[rank], line)) {
	    if (! line.empty())
	      feature_type(line);
	  } else {
	    stream[rank].reset();
	    device[rank].reset();
	  }
	  
	  found = true;
	}
      
      if (std::count(device.begin(), device.end(), device_ptr_type()) == mpi_size) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }

  } else {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(utils::mpi_device_sink(0, feature_tag, 1024 * 1024));
    
    for (feature_type::id_type id = 0; id != feature_type::allocated(); ++ id)
      if (! feature_type(id).empty())
	os << feature_type(id) << '\n';
  }
}

void merge_statistics(const operation_set_type& operations,
		      operation_set_type::statistics_type& statistics)
{
  typedef operation_set_type::statistics_type statistics_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  statistics.clear();
  
  if (mpi_rank == 0) {
    typedef utils::mpi_device_source            device_type;
    typedef boost::iostreams::filtering_istream stream_type;
    
    typedef boost::shared_ptr<device_type> device_ptr_type;
    typedef boost::shared_ptr<stream_type> stream_ptr_type;
    
    typedef std::vector<device_ptr_type, std::allocator<device_ptr_type> > device_ptr_set_type;
    typedef std::vector<stream_ptr_type, std::allocator<stream_ptr_type> > stream_ptr_set_type;

    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
    
    device_ptr_set_type device(mpi_size);
    stream_ptr_set_type stream(mpi_size);
    
    for (int rank = 1; rank != mpi_size; ++ rank) {
      device[rank].reset(new device_type(rank, stat_tag, 4096));
      stream[rank].reset(new stream_type());
      
      stream[rank]->push(boost::iostreams::zlib_decompressor());
      stream[rank]->push(*device[rank]);
    }
    
    statistics = operations.get_statistics();

    std::string line;
    
    int non_found_iter = 0;
    while (1) {
      bool found = false;
      
      for (int rank = 1; rank != mpi_size; ++ rank)
	while (stream[rank] && device[rank] && device[rank]->test()) {
	  if (std::getline(*stream[rank], line)) {
	    const utils::piece line_piece(line);
	    tokenizer_type tokenizer(line_piece);
	    
	    tokenizer_type::iterator iter = tokenizer.begin();
	    if (iter == tokenizer.end()) continue;
	    const utils::piece name = *iter;
	    
	    ++ iter;
	    if (iter == tokenizer.end()) continue;
	    const utils::piece count = *iter;
	    
	    ++ iter;
	    if (iter == tokenizer.end()) continue;
	    const utils::piece node = *iter;
	
	    ++ iter;
	    if (iter == tokenizer.end()) continue;
	    const utils::piece edge = *iter;
	
	    ++ iter;
	    if (iter == tokenizer.end()) continue;
	    const utils::piece user_time = *iter;
	    
	    ++ iter;
	    if (iter == tokenizer.end()) continue;
	    const utils::piece cpu_time = *iter;

	    ++ iter;
	    if (iter == tokenizer.end()) continue;
	    const utils::piece thread_time = *iter;
	
	    statistics[name] += statistics_type::statistic_type(utils::lexical_cast<statistics_type::count_type>(count),
								utils::lexical_cast<statistics_type::count_type>(node),
								utils::lexical_cast<statistics_type::count_type>(edge),
								utils::decode_base64<statistics_type::second_type>(user_time),
								utils::decode_base64<statistics_type::second_type>(cpu_time),
								utils::decode_base64<statistics_type::second_type>(thread_time));
	  } else {
	    stream[rank].reset();
	    device[rank].reset();
	  }
	  
	  found = true;
	}
      
      if (std::count(device.begin(), device.end(), device_ptr_type()) == mpi_size) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(utils::mpi_device_sink(0, stat_tag, 4096));
    
    statistics_type::const_iterator siter_end = operations.get_statistics().end();
    for (statistics_type::const_iterator siter = operations.get_statistics().begin(); siter != siter_end; ++ siter) {
      os << siter->first
	 << ' ' << siter->second.count
	 << ' ' << siter->second.node
	 << ' ' << siter->second.edge;
      os << ' ';
      utils::encode_base64(siter->second.user_time, std::ostream_iterator<char>(os));
      os << ' ';
      utils::encode_base64(siter->second.cpu_time, std::ostream_iterator<char>(os));
      os << ' ';
      utils::encode_base64(siter->second.thread_time, std::ostream_iterator<char>(os));
      os << '\n';
    }
    os << '\n';
  }
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
      request_recv[rank] = MPI::COMM_WORLD.Irecv(0, 0, utils::mpi_traits<int>::data_type(), rank, notify_tag);
      request_send[rank] = MPI::COMM_WORLD.Isend(0, 0, utils::mpi_traits<int>::data_type(), rank, notify_tag);
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
    MPI::Request request_send = MPI::COMM_WORLD.Isend(0, 0, utils::mpi_traits<int>::data_type(), 0, notify_tag);
    MPI::Request request_recv = MPI::COMM_WORLD.Irecv(0, 0, utils::mpi_traits<int>::data_type(), 0, notify_tag);
    
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
  struct model_file_type
  {
    weight_set_type weights_;
    tree_rnn_type   theta_;
    path_type       path_weights_;
    path_type       path_theta_;

    model_file_type() : weights_(), theta_(), path_weights_(), path_theta_() {}
    model_file_type(const weight_set_type& weights,
		    const tree_rnn_type& theta,
		    const path_type& path_weights,
		    const path_type& path_theta)
      : weights_(weights),
	theta_(theta),
	path_weights_(path_weights),
	path_theta_(path_theta) {}
  };
  
  typedef utils::lockfree_list_queue<model_file_type, std::allocator<model_file_type> > queue_type;
  
  Dumper(queue_type& queue)
    : queue_(queue) {}
  
  void operator()()
  {
    model_file_type model_file;
    
    while (1) {
      queue_.pop(model_file);
      if (model_file.weights_.empty()) break;

      if (! model_file.path_weights_.empty()) {
	utils::compress_ostream os(model_file.path_weights_, 1024 * 1024);
	os.precision(20);
	os << model_file.weights_;
      }
      
      if (! model_file.path_theta_.empty())
	model_file.theta_.write(model_file.path_theta_);
    }
  }
  
  queue_type& queue_;
};

template <typename Learner>
struct Task
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef std::string update_encoded_type;
  typedef utils::lockfree_list_queue<update_encoded_type, std::allocator<update_encoded_type> > queue_type;
  
  typedef std::vector<size_t, std::allocator<size_t> > segment_set_type;

  typedef cicada::feature::FrontierBiTreeRNN::feature_name_set_type feature_name_set_type;

  typedef tree_rnn_type::tensor_type tensor_type;

  typedef typename Learner::candidate_type     candidate_type;
  typedef typename Learner::candidate_set_type candidate_set_type;
  typedef typename Learner::candidate_map_type candidate_map_type;
  
  typedef typename Learner::gradient_type gradient_type;
  
  struct Encoder
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    std::string operator()(const gradient_type& gradient)
    {
      buffer_.clear();

      boost::iostreams::filtering_ostream os;
      os.push(codec::lz4_compressor());
      os.push(boost::iostreams::back_insert_device<buffer_type>(buffer_));
      
      os << gradient.theta_;
      
      feature_set_type::const_iterator fiter_end = gradient.weights_.end();
      for (feature_set_type::const_iterator fiter = gradient.weights_.begin(); fiter != fiter_end; ++ fiter) {
	const size_type feature_size = fiter->first.size();
	
	os.write((char*) &feature_size, sizeof(size_type));
	os.write((char*) &(*fiter->first.begin()), feature_size);
	os.write((char*) &fiter->second, sizeof(feature_set_type::mapped_type));
      }

      os.reset();
      
      return std::string(buffer_.begin(), buffer_.end());
    }
    
    buffer_type buffer_;
  };
  
  struct Decoder
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;

    void operator()(const update_encoded_type& encoded, gradient_type& gradient)
    {
      gradient.clear();
      
      boost::iostreams::filtering_istream is;
      is.push(codec::lz4_decompressor());
      is.push(boost::iostreams::array_source(&(*encoded.begin()), encoded.size()));
      
      is >> gradient.theta_;
      
      size_type feature_size = 0;
      feature_set_type::mapped_type value;
      
      while (is) {
	if (! is.read((char*) &feature_size, sizeof(size_type))) break;
	buffer_.resize(feature_size);
	if (! is.read((char*) &(*buffer_.begin()), buffer_.size())) break;
	if (! is.read((char*) &value, sizeof(feature_set_type::data_type))) break;
	
	gradient.weights_[feature_type(buffer_.begin(), buffer_.end())] = value;
      }
    }
    
    buffer_type buffer_;
  };
  
  Task(const int rank,
       queue_type& queue_merge,
       queue_type& queue_bcast,
       const Learner& learner,
       operation_set_type& operations,
       const model_type& model,
       const event_set_type& events,
       const scorer_document_type& scorers,
       weight_set_type& weights,
       tree_rnn_type& theta)
    : rank_(rank),
      queue_merge_(queue_merge),
      queue_bcast_(queue_bcast),
      learner_(learner),
      operations_(operations),
      model_(model),
      events_(events),
      scorers_(scorers),
      weights_(weights),
      theta_(theta)
  {
    generator_.seed(utils::random_seed());
    
    // intialize segments_ 
    segments_.clear();
    for (size_t seg = 0; seg != events.size(); ++ seg)
      if (! events[seg].empty())
	segments_.push_back(seg);
    
    // check for the feature
    const cicada::feature::FrontierBiTreeRNN* tree_rnn_feature = 0;
    
    model_type::const_iterator miter_end = model.end();
    for (model_type::const_iterator miter = model.begin(); miter != miter_end; ++ miter)
      if (dynamic_cast<cicada::feature::FrontierBiTreeRNN*>(miter->get())) {
	tree_rnn_feature = dynamic_cast<cicada::feature::FrontierBiTreeRNN*>(miter->get());
	break;
      }
    
    if (! tree_rnn_feature)
      throw std::runtime_error("no tree-rnn feature?");
    
    if (&tree_rnn_feature->model() != &theta)
      throw std::runtime_error("different theta?");
    
    names_ = tree_rnn_feature->features();

    no_bos_eos_    = tree_rnn_feature->no_bos_eos();
    skip_sgml_tag_ = tree_rnn_feature->skip_sgml_tag();
  }
  
  const int rank_;

  queue_type& queue_merge_;
  queue_type& queue_bcast_;
  
  const Learner&  learner_;
  
  operation_set_type&           operations_;
  const model_type&             model_;
  const event_set_type&         events_;
  const scorer_document_type&   scorers_;
  weight_set_type&              weights_;
  tree_rnn_type&                theta_;
  
  KBestSentence   kbest_generator_;
  Oracle          oracle_generator_;
  
  Encoder         encoder_;
  Decoder         decoder_;
  
  boost::mt19937        generator_;
  segment_set_type      segments_;
  feature_name_set_type names_;
  
  bool no_bos_eos_;
  bool skip_sgml_tag_;
    
  score_ptr_type score_1best_;
  score_ptr_type score_oracle_;
  size_type      num_update_;
  
  struct history_type
  {
    history_type() {}
    
    segment_set_type   segments;
    candidate_map_type kbests;
    candidate_map_type oracles;
  };
  typedef std::deque<history_type, std::allocator<history_type> > history_set_type;
  
  void operator()()
  {
    candidate_set_type   kbests;
    
    segment_set_type     segments_batch;
    candidate_map_type   kbests_batch;
    candidate_map_type   oracles_batch;
    scorer_document_type scorers_batch(scorers_);

    history_set_type history;
    
    score_1best_.reset();
    score_oracle_.reset();
    num_update_ = 0;
    
    segment_set_type::const_iterator siter     = segments_.begin();
    segment_set_type::const_iterator siter_end = segments_.end();
    
    bool merge_finished = false;
    bool learn_finished = false;
    
    gradient_type       gradient(theta_.hidden_, theta_.embedding_);
    update_encoded_type encoded;

    tensor_type W;
    
    const_cast<Learner&>(learner_).initialize(names_,
					      no_bos_eos_,
					      skip_sgml_tag_,
					      weights_,
					      W,
					      theta_);
    
    int non_found_iter = 0;
    while (! merge_finished || ! learn_finished) {
      bool found = false;
      
      if (! merge_finished) {
	size_t num_updated = 0;
	
	while (queue_merge_.pop_swap(encoded, true)) {
	  if (encoded.empty()) {
	    merge_finished = true;
	    break;
	  }
	  
	  decoder_(encoded, gradient);
	  
	  const_cast<Learner&>(learner_).learn(weights_, W, theta_, gradient);

	  ++ num_updated;
	}
	
	if (num_updated && debug >= 2)
	  std::cerr << "rank: " << rank_ << " updated weights: " << num_updated << std::endl;

	found |= num_updated;
	num_update_ += num_updated;
      }
      
      if (! learn_finished) {
	segments_batch.clear();
	kbests_batch.clear();
	oracles_batch.clear();
	scorers_batch.clear();
	
	segment_set_type::const_iterator siter_last = std::min(siter + batch_size, siter_end);
	for (/**/; siter != siter_last; ++ siter) {
	  const size_t id = *siter;
	  
	  if (events_[id].empty() || ! scorers_[id]) continue;
	  
	  kbest_generator_(operations_, events_[id], weights_, kbests);
	  
	  if (kbests.empty()) continue;
	  
	  segments_batch.push_back(id);
	  kbests_batch.push_back(kbests);
	  scorers_batch.push_back(scorers_[id]);
	}
	
	// if we have segments!
	if (! segments_batch.empty()) {
	  // oracle computation
	  std::pair<score_ptr_type, score_ptr_type> scores = oracle_generator_(kbests_batch,
									       scorers_batch,
									       oracles_batch,
									       generator_);
	  
	  if (! score_1best_)
	    score_1best_ = scores.first;
	  else
	    *score_1best_ += *(scores.first);
	  
	  if (! score_oracle_)
	    score_oracle_ = scores.second;
	  else
	    *score_oracle_ += *(scores.second);
	    
	  if (debug >= 2)
	    std::cerr << "rank: " << rank_ << " batch 1best:  " << *scores.first << std::endl
		      << "rank: " << rank_ << " batch oracle: " << *scores.second << std::endl
		      << "rank: " << rank_ << " accumulated 1best:  " << *score_1best_ << std::endl
		      << "rank: " << rank_ << " accumulated oracle: " << *score_oracle_ << std::endl;
	  
	  gradient.clear();
	  
	  // encode into learner...
	  double objective = 0.0;
	  
	  if (! history.empty())
	    for (size_t j = 0; j != history.size(); ++ j)
	      for (size_t i = 0; i != history[j].segments.size(); ++ i)
		objective += const_cast<Learner&>(learner_).accumulate(history[j].segments[i],
								       history[j].kbests[i],
								       history[j].oracles[i],
								       weights_,
								       W,
								       theta_,
								       gradient);
	  
	  for (size_t i = 0; i != kbests_batch.size(); ++ i)
	    objective += const_cast<Learner&>(learner_).accumulate(segments_batch[i],
								   kbests_batch[i],
								   oracles_batch[i],
								   weights_,
								   W,
								   theta_,
								   gradient);
	  
	  // perform parameter updates...
	  const_cast<Learner&>(learner_).learn(weights_, W, theta_, gradient);
	  
	  if (debug >= 2)
	    std::cerr << "rank: " << rank_ << " objective: " << objective << " batch: " << kbests_batch.size() << std::endl;
	  
	  // here, we will bcast the updated amount to others...
	  if (gradient.count_) {
	    queue_bcast_.push(encoder_(gradient));
	    
	    ++ num_update_;
	  }

	  // update history
	  if (merge_history > 0) {
	    if (static_cast<int>(history.size()) >= merge_history)
	      history.erase(history.begin());
	    
	    history.push_back(history_type());
	    
	    history.back().segments.swap(segments_batch);
	    history.back().kbests.swap(kbests_batch);
	    history.back().oracles.swap(oracles_batch);
	  }
	}
	
	// signal finished!
	if (siter == siter_end) {
	  learn_finished = true;
	  queue_bcast_.push(update_encoded_type());
	}
	
	found = true;
      }
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
        
    // random shuffle the segments
    boost::random_number_generator<boost::mt19937> gen(generator_);
    
    std::random_shuffle(segments_.begin(), segments_.end(), gen);
  }
};

struct abs_sum
{
  double operator()(const double& x, const double& y) const
  {
    return x + std::fabs(y);
  }
};


template <typename Learner>
void cicada_learn(const Learner& learner,
		  operation_set_type& operations,
		  const model_type& model,
		  const event_set_type& events,
		  const scorer_document_type& scorers,
		  weight_set_type& weights,
		  tree_rnn_type& theta)
{
  typedef Dumper dumper_type;
  
  typedef Task<Learner> task_type;
  
  typedef typename task_type::update_encoded_type encoded_type;
  typedef typename dumper_type::model_file_type   model_file_type;
  
  typedef encoded_type                   buffer_type;
  typedef boost::shared_ptr<buffer_type> buffer_ptr_type;
  typedef std::deque<buffer_ptr_type, std::allocator<buffer_ptr_type> >  buffer_set_type;
  typedef std::vector<buffer_set_type, std::allocator<buffer_set_type> > buffer_map_type;
  
  typedef boost::shared_ptr<utils::mpi_ostream_simple> ostream_ptr_type;
  typedef boost::shared_ptr<utils::mpi_istream_simple> istream_ptr_type;
  
  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  buffer_type          buffer;
  buffer_map_type      buffers(mpi_size);
  ostream_ptr_set_type ostreams(mpi_size);
  istream_ptr_set_type istreams(mpi_size);
  
  typename task_type::queue_type queue_merge;
  typename task_type::queue_type queue_bcast;
  
  task_type task(mpi_rank,
		 queue_merge,
		 queue_bcast,
		 learner,
		 operations,
		 model,
		 events,
		 scorers,
		 weights,
		 theta);

  // prepare dumper for the root
  dumper_type::queue_type queue_dumper;
  std::auto_ptr<boost::thread> dumper(dump_weights_mode && mpi_rank == 0
				      ? new boost::thread(dumper_type(queue_dumper))
				      : 0);
  
  // first, bcast weights...
  bcast_weights(weights, theta);

  // start training
  for (int iter = 0; iter != iteration; ++ iter) {
    // perform learning
    
    if (debug && mpi_rank == 0)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    // prepare iostreams...
    for (int rank = 0; rank < mpi_size; ++ rank)
      if (rank != mpi_rank) {
	ostreams[rank].reset(new utils::mpi_ostream_simple(rank, update_tag, 4096));
	istreams[rank].reset(new utils::mpi_istream_simple(rank, update_tag, 4096));
      }
    
    // create thread!
    boost::thread worker(boost::ref(task));
    
    bool finished = false;
    
    int non_found_iter = 0;
    while (1) {
      bool found = false;
      
      // reduce samples...
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && istreams[rank] && istreams[rank]->test()) {
	  if (istreams[rank]->read(buffer)) {
	    buffer_type(buffer).swap(buffer);
	    queue_merge.push_swap(buffer);
	  } else
	    istreams[rank].reset();
	  
	  buffer.clear();
	  found = true;
	}
      
      // check termination...
      if (! finished && std::count(istreams.begin(), istreams.end(), istream_ptr_type()) == mpi_size) {
	queue_merge.push(buffer_type());
	finished = true;
      }
      
      // bcast samples...
      // first, get the queue...
      if (queue_bcast.pop_swap(buffer, true)) {
	buffer_ptr_type buffer_ptr;
	if (! buffer.empty()) {
	  buffer_ptr.reset(new buffer_type());
	  buffer_ptr->swap(buffer);
	  buffer.clear();
	}
	
	for (int rank = 0; rank != mpi_size; ++ rank) 
	  if (rank != mpi_rank)
	    buffers[rank].push_back(buffer_ptr);
	
	found = true;
      }
      
      // second, bcast...
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && ostreams[rank] && ostreams[rank]->test() && ! buffers[rank].empty()) {
	  if (! buffers[rank].front()) {
	    // termination!
	    if (! ostreams[rank]->terminated())
	      ostreams[rank]->terminate();
	    else {
	      ostreams[rank].reset();
	      buffers[rank].erase(buffers[rank].begin());
	    }
	  } else {
	    ostreams[rank]->write(*(buffers[rank].front()));
	    buffers[rank].erase(buffers[rank].begin());
	  }
	  
	  found = true;
	}
      
      // termination condition
      if (finished
	  && std::count(istreams.begin(), istreams.end(), istream_ptr_type()) == mpi_size
	  && std::count(ostreams.begin(), ostreams.end(), ostream_ptr_type()) == mpi_size) break;
      
      // a conventional loop...
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    // join worker!
    worker.join();
        
    if (debug) {
      score_ptr_type score_1best(task.score_1best_);
      score_ptr_type score_oracle(task.score_oracle_);
      
      if (debug >= 2 && mpi_rank == 0)
	std::cerr << "reducing evaluation scores" << std::endl;
      
      reduce_score_pair(score_1best, score_oracle);
      
      if (mpi_rank == 0)
	std::cerr << "total 1best:  " << (score_1best ? score_1best->description() : std::string("?")) << std::endl
		  << "total oracle: " << (score_oracle ? score_oracle->description() : std::string("?")) << std::endl;
    }
    
    if (mix_average_mode) {
      // +1 to avid zero...
      const size_t updated = task.num_update_ + 1;
      
      weights *= updated;
      theta   *= updated;
      
      reduce_weights(weights, theta);
      
      bcast_weights(weights, theta);
      
      size_t updated_total = 0;
      MPI::COMM_WORLD.Allreduce(&updated, &updated_total, 1, utils::mpi_traits<size_t>::data_type(), MPI::SUM);
      
      weights *= 1.0 / updated_total;
      theta   *= 1.0 / updated_total;
    } else if (mix_select_mode) {
      typedef std::vector<double, std::allocator<double> > buffer_type;
      
      const double l1 = std::accumulate(weights.begin(), weights.end(), double(0.0), abs_sum());
      
      buffer_type buffer_send(mpi_size, 0.0);
      buffer_type buffer_recv(mpi_size, 0.0);
      
      buffer_send[mpi_rank] = l1;
      buffer_recv[mpi_rank] = l1;
      
      MPI::COMM_WORLD.Reduce(&(*buffer_send.begin()), &(*buffer_recv.begin()), mpi_size, utils::mpi_traits<double>::data_type(), MPI::MAX, 0);
      
      int rank_min = (std::min_element(buffer_recv.begin(), buffer_recv.end()) - buffer_recv.begin());
      
      MPI::COMM_WORLD.Bcast(&rank_min, 1, utils::mpi_traits<int>::data_type(), 0);
      
      if (debug >= 2 && mpi_rank == 0)
	std::cerr << "minimum rank: " << rank_min << " L1: " << buffer_recv[rank_min] << std::endl;
      
      bcast_weights(rank_min, weights, theta);
    } 
    
    if (dump_weights_mode && mpi_rank == 0)
      queue_dumper.push(model_file_type(weights,
					theta,
					add_suffix(output_weights_file, "." + utils::lexical_cast<std::string>(iter + 1)),
					add_suffix(output_model_file, "." + utils::lexical_cast<std::string>(iter + 1))));
  }
  
  // finish dumper...
  if (dumper.get()) {
    queue_dumper.push(model_file_type());
    dumper->join();
  }
}

void reduce_score_pair(score_ptr_type& score_1best, score_ptr_type& score_oracle)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    std::string score_str_1best;
    std::string score_str_oracle;
    
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      is.push(utils::mpi_device_source(rank, score_tag, 256));
      
      is >> score_str_1best;
      is >> score_str_oracle;
      
      if (score_str_1best != "?") {
	if (! score_1best)
	  score_1best = scorer_type::score_type::decode(score_str_1best);
	else
	  *score_1best += *scorer_type::score_type::decode(score_str_1best);
      }
      
      if (score_str_oracle != "?") {
	if (! score_oracle)
	  score_oracle = scorer_type::score_type::decode(score_str_oracle);
	else
	  *score_oracle += *scorer_type::score_type::decode(score_str_oracle);      
      }
    }
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(utils::mpi_device_sink(0, score_tag, 256));
    os << (score_1best ? score_1best->encode() : "?");
    os << ' ';
    os << (score_oracle ? score_oracle->encode() : "?");
  }
}


void write_weights(std::ostream& os, const weight_set_type& weights)
{
  const size_t weights_size = weights.size();
  os.write((char*) &weights_size, sizeof(size_t));
  
  for (feature_type::id_type id = 0; id != weights_size; ++ id) {
    const feature_type feature(id);
    const size_t feature_size = feature.size();
    
    os.write((char*) &feature_size, sizeof(size_t));
    os.write((char*) &(*feature.begin()), feature_size);
    
    const weight_set_type::value_type value = weights[id];
    os.write((char*) &value, sizeof(weight_set_type::value_type));
  }
}

void read_weights(std::istream& is, weight_set_type& weights)
{
  std::vector<char, std::allocator<char> > buffer;
  weight_set_type::value_type              value;
  
  size_t weights_size = 0;
  is.read((char*) &weights_size, sizeof(size_t));
  
  for (size_t i = 0; i != weights_size; ++ i) {
    size_t feature_size = 0;
    is.read((char*) &feature_size, sizeof(size_t));
    
    buffer.resize(feature_size);
    is.read((char*) &(*buffer.begin()), feature_size);
    is.read((char*) &value, sizeof(weight_set_type::value_type));
    
    weights[feature_type(buffer.begin(), buffer.end())] = value;
  }
}

// reduce and perform "+"
void reduce_weights(std::istream& is, weight_set_type& weights)
{
  std::vector<char, std::allocator<char> > buffer;
  weight_set_type::value_type              value;
  
  size_t weights_size = 0;
  is.read((char*) &weights_size, sizeof(size_t));
  
  for (size_t i = 0; i != weights_size; ++ i) {
    size_t feature_size = 0;
    is.read((char*) &feature_size, sizeof(size_t));
    
    buffer.resize(feature_size);
    is.read((char*) &(*buffer.begin()), feature_size);
    is.read((char*) &value, sizeof(weight_set_type::value_type));
    
    weights[feature_type(buffer.begin(), buffer.end())] += value;
  }  
}

void recv_weights(const int rank, weight_set_type& weights, tree_rnn_type& theta)
{
  weights.clear();
  weights.allocate();
  
  boost::iostreams::filtering_istream is;
  is.push(codec::lz4_decompressor());
  is.push(utils::mpi_device_source(rank, weights_tag, 1024 * 1024));
  
  is >> theta;
  read_weights(is, weights);
}

void send_weights(const int rank, const weight_set_type& weights, const tree_rnn_type& theta)
{
  boost::iostreams::filtering_ostream os;
  os.push(codec::lz4_compressor());
  os.push(utils::mpi_device_sink(rank, weights_tag, 1024 * 1024));
  
  os << theta;
  write_weights(os, weights);
}

void reduce_weights(const int rank, weight_set_type& weights, tree_rnn_type& theta)
{
  tree_rnn_type theta_reduced;
  
  boost::iostreams::filtering_istream is;
  is.push(codec::lz4_decompressor());
  is.push(utils::mpi_device_source(rank, weights_tag, 1024 * 1024));
  
  is >> theta_reduced;
  reduce_weights(is, weights);
  
  theta += theta_reduced;
}

template <typename Iterator>
void reduce_weights(Iterator first, Iterator last, weight_set_type& weights, tree_rnn_type& theta)
{
  for (/**/; first != last; ++ first)
    reduce_weights(*first, weights, theta);
}

void reduce_weights(weight_set_type& weights, tree_rnn_type& theta)
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
	reduce_weights(ranks.front(), weights, theta);
      else
	reduce_weights(ranks.begin(), ranks.end(), weights, theta);
      
    } else
      send_weights(mpi_rank % reduce_size, weights, theta);
    
    merge_size = reduce_size;
  }
}

void bcast_weights(const int rank, weight_set_type& weights, tree_rnn_type& theta)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == rank) {
    boost::iostreams::filtering_ostream os;
    os.push(codec::lz4_compressor());
    os.push(utils::mpi_device_bcast_sink(rank, 1024 * 1024));

    os << theta;
    write_weights(os, weights);
  } else {
    weights.clear();
    weights.allocate();
    
    boost::iostreams::filtering_istream is;
    is.push(codec::lz4_decompressor());
    is.push(utils::mpi_device_bcast_source(rank, 1024 * 1024));
    
    is >> theta;
    read_weights(is, weights);
  }
}

void bcast_weights(weight_set_type& weights, tree_rnn_type& theta)
{
  bcast_weights(0, weights, theta);
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
    ("output-feature-function", po::value<path_type>(&output_feature),                                    "output feature function(s)")
    
    //operatins...
    ("operation",      po::value<op_set_type>(&ops)->composing(), "operations")
    ("operation-list", po::bool_switch(&op_list),                 "list of available operation(s)")
    
    ("scorer-list",    po::bool_switch(&scorer_list),    "list of available scorers")
    ("format-list",    po::bool_switch(&format_list),    "list of available formatters")
    ("signature-list", po::bool_switch(&signature_list), "list of available signatures")
    ("stemmer-list",   po::bool_switch(&stemmer_list),   "list of available stemmers")
    ("tokenizer-list", po::bool_switch(&tokenizer_list), "list of available tokenizers")
    ("matcher-list",   po::bool_switch(&matcher_list),   "list of available matchers")
    ;

  po::options_description opts_learn("learning options");
  opts_learn.add_options()
    ("refset",         po::value<path_type>(&refset_file),         "refset")
    ("output-weights", po::value<path_type>(&output_weights_file), "model (or weights) output")
    ("output-model",   po::value<path_type>(&output_model_file),   "rnn model output")
    ("weights",        po::value<path_type>(&weights_file),        "initial model (or weights)")
    
    ("scorer",           po::value<std::string>(&scorer_name)->default_value(scorer_name), "evaluation scorer")
    
    ("iteration",       po::value<int>(&iteration)->default_value(iteration),                "learning iterations")
    ("batch",           po::value<int>(&batch_size)->default_value(batch_size),              "batch (or batch, bin) size")
    ("kbest",           po::value<int>(&kbest_size)->default_value(kbest_size),              "kbest size")
    ("kbest-unique",    utils::true_false_switch(&kbest_unique_mode),                        "unique kbest")
    ("kbest-diversity", po::value<double>(&kbest_diversity)->default_value(kbest_diversity), "kbest diversity")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("lambda",           po::value<double>(&lambda)->default_value(lambda),      "regularization constant")
    ("eta0",             po::value<double>(&eta0)->default_value(eta0),          "\\eta_0 for decay")
    
    ("merge-history",       po::value<int>(&merge_history),         "merge history for decoded results")
    ("mix-none",            po::bool_switch(&mix_none_mode),        "no mixing")
    ("mix-average",         po::bool_switch(&mix_average_mode),     "mixing weights by averaging")
    ("mix-select",          po::bool_switch(&mix_select_mode),      "select weights by L1")
    ("dump-weights",        po::bool_switch(&dump_weights_mode),    "dump mode (or weights) during iterations")
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

