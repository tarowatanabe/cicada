//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// mini-batch kbest-based online asynchronous learning
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

#include "liblbfgs/lbfgs.hpp"
#include "cg_descent/cg.hpp"

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>
#include <boost/random.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file;
path_type oracle_file;

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
path_type output_file;  // weights output
path_type weights_file;

// scorers
std::string scorer_name = "bleu:order=4,exact=true";
bool yield_sentence   = false;
bool yield_alignment  = false;
bool yield_dependency = false;

// learning parameters
int iteration = 10;
int batch_size = 8;
int kbest_size = 1000;
bool kbest_diverse_mode = false;

// solver parameters
bool learn_xbleu = false;
bool learn_hinge = false;
bool learn_ohinge = false;
bool learn_pa = false;
bool learn_cw = false;
bool learn_arow = false;
bool learn_nherd = false;
bool learn_mira   = false;
bool learn_softmax    = false;
bool learn_osoftmax   = false;
bool learn_el     = false;
bool learn_oel     = false;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1e-3;
double oscar = 0.0;
double temperature = 0.0;
double scale = 1.0;
double alpha0 = 0.85;
double eta0 = 0.2;
int order = 4;

bool rate_simple = false;
bool rate_exponential = false;
bool rate_adagrad = false;

// additional misc parameters...
bool rda_mode = false;
bool loss_rank = false; // loss by rank
bool softmax_margin = false;
bool project_weight = false;
bool merge_oracle_mode = false;
int merge_history = 0;
bool mix_none_mode = false;
bool mix_average_mode = false;
bool mix_select_mode = false;
int  mix_kbest_features = 0;
bool dump_weights_mode   = false; // dump current weights... for debugging purpose etc.

int debug = 0;

#include "cicada_learn_asynchronous_kbest_impl.hpp"

// forward declarations...

void options(int argc, char** argv);
void merge_features();
void merge_statistics(const operation_set_type& operations, operation_set_type::statistics_type& statistics);

void cicada_learn(operation_set_type& operations,
		  const event_set_type& events,
		  const event_set_type& events_oracle,
		  const scorer_document_type& scorers,
		  weight_set_type& weights);
void synchronize();

void bcast_weights(const int rank, weight_set_type& weights);
void bcast_weights(weight_set_type& weights);
void reduce_weights(weight_set_type& weights);
void reduce_score_pair(score_ptr_type& score_1best, score_ptr_type& score_oracle);

void send_weights(const int rank, const weight_set_type& weights);
void recv_weights(const int rank, weight_set_type& weights);

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
    
    if (output_file.empty())
      throw std::runtime_error("no output?");
    
    if (int(yield_sentence) + yield_alignment + yield_dependency > 1)
      throw std::runtime_error("you can specify either --yield-{sentence,alignment,dependency}");
    if (int(yield_sentence) + yield_alignment + yield_dependency == 0)
      yield_sentence = true;
    
    if (int(learn_xbleu) + learn_mira + learn_softmax + learn_osoftmax + learn_el + learn_oel + learn_hinge + learn_ohinge + learn_pa + learn_cw + learn_arow + learn_nherd > 1)
      throw std::runtime_error("you can specify either --learn-{xbleu,mira,softmax,osoftmax,el,oel,hinge,ohinge,pa,cw,arow}");
    if (int(learn_xbleu) + learn_mira + learn_softmax + learn_osoftmax + learn_el + learn_oel + learn_hinge + learn_ohinge + learn_pa + learn_cw + learn_arow + learn_nherd== 0)
      learn_softmax = true;
    
    const bool regularize_oscar = (oscar > 0.0);
    
    if (int(regularize_l1) + regularize_l2 + regularize_oscar > 1)
      throw std::runtime_error("either L1 or L2, OSCAR regularization");
    if (int(regularize_l1) + regularize_l2 + regularize_oscar == 0)
      regularize_l2 = true;

    if (C <= 0.0)
      throw std::runtime_error("regularization constant must be positive: " + utils::lexical_cast<std::string>(C));
    if (scale <= 0.0)
      throw std::runtime_error("weight scale constant must be positive: " + utils::lexical_cast<std::string>(scale));

    if (batch_size <= 0)
      throw std::runtime_error("batch size must be possitive: " + utils::lexical_cast<std::string>(batch_size));
    if (kbest_size <= 0)
      throw std::runtime_error("kbest size must be possitive: " + utils::lexical_cast<std::string>(kbest_size));
    
    if (order <= 0)
      throw std::runtime_error("ngram order for xBLEU must be positive");

    if (int(mix_none_mode) + mix_average_mode + mix_select_mode + (mix_kbest_features > 0) > 1)
      throw std::runtime_error("you can specify only one of mix-{none,average,select,kbest-feautures}");
    
    if (int(mix_none_mode) + mix_average_mode + mix_select_mode + (mix_kbest_features > 0) == 0)
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

    event_set_type events_oracle;
    if (! oracle_file.empty()) {
      if (oracle_file != "-" && !boost::filesystem::exists(oracle_file))
	throw std::runtime_error("no oracle input? " + oracle_file.string());
      
      read_events(oracle_file, events_oracle, input_directory_mode, input_id_mode, mpi_rank, mpi_size);
    }
    
    // read reference data
    scorer_document_type scorers(scorer_name);
    read_refset(refset_file, scorers, mpi_rank, mpi_size);
    
    if (scorers.size() != events.size())
      throw std::runtime_error("training sample size and reference translation size does not match");
    
    if (! events_oracle.empty())
      if (scorers.size() != events_oracle.size())
	throw std::runtime_error("oracle size and reference translation size does not match");
    
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
    cicada_learn(operations, events, events_oracle, scorers, weights);
    
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_file, 1024 * 1024);
      os.precision(20);
      os << weights;
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

template <typename Learner, typename KBestGenerator, typename OracleGenerator, typename Generator>
struct Task
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef std::string update_encoded_type;
  typedef utils::lockfree_list_queue<update_encoded_type, std::allocator<update_encoded_type> > queue_type;
  
  typedef std::vector<size_t, std::allocator<size_t> > segment_set_type;
  
  //
  // revise this to use lz4 codec..
  //
  struct Encoder
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    template <typename Iterator>
    std::string operator()(Iterator first, Iterator last)
    {
      buffer.clear();

      boost::iostreams::filtering_ostream os;
      os.push(codec::lz4_compressor());
      os.push(boost::iostreams::back_insert_device<buffer_type>(buffer));
      
      for (/**/; first != last; ++ first) {
	const size_type feature_size = first->first.size();
	
	os.write((char*) &feature_size, sizeof(size_type));
	os.write((char*) &(*first->first.begin()), feature_size);
	os.write((char*) &first->second, sizeof(feature_set_type::data_type));
      }

      os.reset();
      
      return std::string(buffer.begin(), buffer.end());
    }
    
    buffer_type buffer;
  };
  
  struct Decoder
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;

    void operator()(const update_encoded_type& encoded, feature_set_type& updates)
    {
      updates.clear();
      
      boost::iostreams::filtering_istream is;
      is.push(codec::lz4_decompressor());
      is.push(boost::iostreams::array_source(&(*encoded.begin()), encoded.size()));
      
      size_type feature_size = 0;
      feature_set_type::data_type value;
      
      while (is) {
	if (! is.read((char*) &feature_size, sizeof(size_type))) break;
	buffer.resize(feature_size);
	if (! is.read((char*) &(*buffer.begin()), buffer.size())) break;
	if (! is.read((char*) &value, sizeof(feature_set_type::data_type))) break;
	
	updates[feature_type(buffer.begin(), buffer.end())] = value;
      }
    }
    
    buffer_type buffer;
  };
  
  Task(const int rank,
       queue_type& queue_merge,
       queue_type& queue_bcast,
       Learner& learner,
       KBestGenerator& kbest_generator,
       OracleGenerator& oracle_generator,
       operation_set_type& operations,
       const event_set_type& events,
       const event_set_type& events_oracle,
       const scorer_document_type& scorers,
       weight_set_type& weights,
       const segment_set_type& segments,
       Generator& generator)
    : rank_(rank),
      queue_merge_(queue_merge),
      queue_bcast_(queue_bcast),
      learner_(learner),
      kbest_generator_(kbest_generator),
      oracle_generator_(oracle_generator),
      operations_(operations),
      events_(events),
      events_oracle_(events_oracle),
      scorers_(scorers),
      weights_(weights),
      segments_(segments),
      generator_(generator)

  {}
  
  const int rank_;

  queue_type& queue_merge_;
  queue_type& queue_bcast_;

  Learner&         learner_;
  KBestGenerator&  kbest_generator_;
  OracleGenerator& oracle_generator_;
  
  operation_set_type&           operations_;
  const event_set_type&         events_;
  const event_set_type&         events_oracle_;
  const scorer_document_type&   scorers_;
  weight_set_type&              weights_;
  
  const segment_set_type&       segments_;
  Generator&      generator_;
  
  Encoder         encoder_;
  Decoder         decoder_;

  score_ptr_type score_1best_;
  score_ptr_type score_oracle_;
  size_type      num_update_;

  struct history_type
  {
    history_type() {}
    history_type(const std::string& parameter) { }
    
    segment_set_type    segments;
    hypothesis_map_type kbests;
    hypothesis_map_type oracles;
  };
  typedef std::deque<history_type, std::allocator<history_type> > history_set_type;

  
  void operator()()
  {
    hypothesis_set_type kbests;
    
    segment_set_type     segments_batch;
    hypothesis_map_type  kbests_batch;
    hypothesis_map_type  kbests_oracle_batch;
    hypothesis_map_type  oracles_batch;
    scorer_document_type scorers_batch(scorers_);

    history_set_type history;
    
    learner_.initialize(weights_);
    
    score_1best_.reset();
    score_oracle_.reset();
    num_update_ = 0;
    
    segment_set_type::const_iterator siter     = segments_.begin();
    segment_set_type::const_iterator siter_end = segments_.end();
    
    bool merge_finished = false;
    bool learn_finished = false;

    feature_set_type    updates;
    update_encoded_type encoded;
    
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
	  
	  decoder_(encoded, updates);
	  
	  learner_.update(weights_, updates);

	  ++ num_updated;
	}
	
	if (num_updated && debug >= 2)
	  std::cerr << "rank: " << rank_ << " updated weights: " << num_updated << std::endl;

	found |= num_updated;
	num_update_ += num_updated;
      }
      
      if (! learn_finished) {
	if (merge_history > 0) {
	  
	  if (static_cast<int>(history.size()) >= merge_history)
	    history.erase(history.begin());
	  
	  history.push_back(history_type(scorers_.parameter()));
	  
	  history.back().segments.swap(segments_batch);
	  history.back().kbests.swap(kbests_batch);
	  history.back().oracles.swap(oracles_batch);
	}
	
	segments_batch.clear();
	kbests_batch.clear();
	kbests_oracle_batch.clear();
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
	  
	  if (! events_oracle_.empty()) {
	    if (events_oracle_[id].empty())
	      throw std::runtime_error("no oracle? " + utils::lexical_cast<std::string>(id));
	      
	    kbest_generator_(operations_, events_oracle_[id], weights_, kbests);

	    if (kbests.empty())
	      throw std::runtime_error("no kbests for oracle? " + utils::lexical_cast<std::string>(id));
	      
	    if (merge_oracle_mode)
	      kbests_batch.back().insert(kbests_batch.back().end(), kbests.begin(), kbests.end());
	    else
	      kbests_oracle_batch.push_back(kbests);
	  }
	}
	  
	// if we have segments!
	if (! segments_batch.empty()) {
	  // oracle computation
	  std::pair<score_ptr_type, score_ptr_type> scores = (kbests_oracle_batch.empty()
							      ? oracle_generator_(kbests_batch, scorers_batch, oracles_batch, generator_)
							      : oracle_generator_(kbests_batch, kbests_oracle_batch, scorers_batch, oracles_batch, generator_));
	    
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
	    
	  // encode into learner...
	  if (! history.empty())
	    for (size_t j = 0; j != history.size(); ++ j)
	      for (size_t i = 0; i != history[j].segments.size(); ++ i)
		learner_.encode(history[j].segments[i], history[j].kbests[i], history[j].oracles[i]);

	  for (size_t i = 0; i != kbests_batch.size(); ++ i)
	    learner_.encode(segments_batch[i], kbests_batch[i], oracles_batch[i]);
	    
	  // perform learning...
	  const double objective = learner_.learn(weights_, updates);
	  
	  if (debug >= 2)
	    std::cerr << "rank: " << rank_ << " objective: " << objective << " batch: " << kbests_batch.size() << std::endl;
	    
	  // here, we will bcast the updated amount to others...
	  if (! updates.empty()) {
	    queue_bcast_.push(encoder_(updates.begin(), updates.end()));
	    
	    ++ num_update_;
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

    learner_.finalize(weights_);
    
    learner_.clear();
  }
};

struct abs_sum
{
  double operator()(const double& x, const double& y) const
  {
    return x + std::fabs(y);
  }
};

template <typename Learner, typename KBestGenerator, typename OracleGenerator>
void cicada_learn(Learner& learner,
		  KBestGenerator& kbest_generator,
		  OracleGenerator& oracle_generator,
		  operation_set_type& operations,
		  const event_set_type& events,
		  const event_set_type& events_oracle,
		  const scorer_document_type& scorers,
		  weight_set_type& weights);

template <typename Learner>
void cicada_learn_yield(Learner& learner,
			operation_set_type& operations,
			const event_set_type& events,
			const event_set_type& events_oracle,
			const scorer_document_type& scorers,
			weight_set_type& weights)
{
  if (yield_sentence) {
    KBestSentence kbest_generator;
    Oracle        oracle_generator;
    
    cicada_learn(learner, kbest_generator, oracle_generator, operations, events, events_oracle, scorers, weights);
  } else if (yield_alignment) {
    KBestAlignment kbest_generator;
    Oracle         oracle_generator;
    
    cicada_learn(learner, kbest_generator, oracle_generator, operations, events, events_oracle, scorers, weights);
  } else if (yield_dependency) {
    KBestDependency kbest_generator;
    Oracle          oracle_generator;
    
    cicada_learn(learner, kbest_generator, oracle_generator, operations, events, events_oracle, scorers, weights);    
  } else
    throw std::runtime_error("invalid yield");
}

void cicada_learn_learner(Regularize& regularizer,
			  Rate& rate,
			  operation_set_type& operations,
			  const event_set_type& events,
			  const event_set_type& events_oracle,
			  const scorer_document_type& scorers,
			  weight_set_type& weights)
{
  if (learn_pa) {
    LearnPA learner;
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_cw) {
    LearnCW learner;
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_arow) {
    LearnAROW learner;
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_nherd) {
    LearnNHERD learner;
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_mira) {
    LearnMIRA learner;
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_xbleu) {
    LearnXBLEU learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_softmax) {
    LearnSoftmax learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_osoftmax) {
    LearnOSoftmax learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_el) {
    LearnExpectedLoss learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_oel) {
    LearnOExpectedLoss learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_hinge) {
    LearnHinge learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else if (learn_ohinge) {
    LearnOHinge learner(regularizer, rate);
    
    cicada_learn_yield(learner, operations, events, events_oracle, scorers, weights);
  } else
    throw std::runtime_error("invalid learner");
  
}

void cicada_learn_regularizer(Rate& rate,
			      operation_set_type& operations,
			      const event_set_type& events,
			      const event_set_type& events_oracle,
			      const scorer_document_type& scorers,
			      weight_set_type& weights)
{
  const bool regularize_oscar = (oscar > 0.0);

  if (rda_mode) {
    if (regularize_l1) {
      RegularizeRDAL1 regularizer(C);
      
      cicada_learn_learner(regularizer, rate, operations, events, events_oracle, scorers, weights);
    } else if (regularize_l2) {
      RegularizeRDAL2 regularizer(C);
      
      cicada_learn_learner(regularizer, rate, operations, events, events_oracle, scorers, weights);
    } else if (regularize_oscar) {
      RegularizeRDAOSCAR regularizer(C, oscar);
      
      cicada_learn_learner(regularizer, rate, operations, events, events_oracle, scorers, weights);
    } else
      throw std::runtime_error("unsupported regularizer");
  } else {
    if (regularize_l1) {
      RegularizeL1 regularizer(C);
      
      cicada_learn_learner(regularizer, rate, operations, events, events_oracle, scorers, weights);
    } else if (regularize_l2) {
      RegularizeL2 regularizer(C);
      
      cicada_learn_learner(regularizer, rate, operations, events, events_oracle, scorers, weights);
    } else if (regularize_oscar) {
      RegularizeOSCAR regularizer(C, oscar);
      
      cicada_learn_learner(regularizer, rate, operations, events, events_oracle, scorers, weights);
    } else
      throw std::runtime_error("unsupported regularizer");
  }
}

void cicada_learn(operation_set_type& operations,
		  const event_set_type& events,
		  const event_set_type& events_oracle,
		  const scorer_document_type& scorers,
		  weight_set_type& weights)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
    
  size_t instances_rank = 0;
  for (size_t seg = 0; seg != events.size(); ++ seg)
    instances_rank += (! events[seg].empty());
  
  size_t instances = 0;
  MPI::COMM_WORLD.Allreduce(&instances_rank, &instances, 1, utils::mpi_traits<size_t>::data_type(), MPI::SUM);
  
  if (debug && mpi_rank == 0)
    std::cerr << "# of training instances: " << instances << std::endl;
  
  const size_t samples = (instances + batch_size - 1) / batch_size;

  if (rate_simple) {
    RateSimple rate(eta0);
    
    cicada_learn_regularizer(rate, operations, events, events_oracle, scorers, weights);
  } else if (rate_exponential) {
    RateExponential rate(alpha0, eta0, samples);
    
    cicada_learn_regularizer(rate, operations, events, events_oracle, scorers, weights);
  } else if (rate_adagrad) {
    RateAdaGrad rate(eta0);
    
    cicada_learn_regularizer(rate, operations, events, events_oracle, scorers, weights);
  } else
    throw std::runtime_error("unsupported learning rate");  
}

template <typename Learner, typename KBestGenerator, typename OracleGenerator>
void cicada_learn(Learner& learner,
		  KBestGenerator& kbest_generator,
		  OracleGenerator& oracle_generator,
		  operation_set_type& operations,
		  const event_set_type& events,
		  const event_set_type& events_oracle,
		  const scorer_document_type& scorers,
		  weight_set_type& weights)
{
  typedef Dumper dumper_type;

  typedef Task<Learner, KBestGenerator, OracleGenerator, boost::mt19937> task_type;
  
  typedef typename task_type::segment_set_type    segment_set_type;
  typedef typename task_type::update_encoded_type encoded_type;

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
  
  // random number generator...
  boost::mt19937 generator;
  generator.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> gen(generator);
  
  segment_set_type segments;
  for (size_t seg = 0; seg != events.size(); ++ seg)
    if (! events[seg].empty())
      segments.push_back(seg);
  
  segment_set_type(segments).swap(segments);
    
  typename task_type::queue_type queue_merge;
  typename task_type::queue_type queue_bcast;

  task_type task(mpi_rank,
		 queue_merge,
		 queue_bcast,
		 learner,
		 kbest_generator,
		 oracle_generator,
		 operations,
		 events,
		 events_oracle,
		 scorers,
		 weights,
		 segments,
		 generator);

  // prepare dumper for the root
  dumper_type::queue_type queue_dumper;
  std::auto_ptr<boost::thread> dumper(dump_weights_mode && mpi_rank == 0
				      ? new boost::thread(dumper_type(queue_dumper))
				      : 0);
  
  // first, bcast weights...
  bcast_weights(weights);

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
    
    // randomize..
    std::random_shuffle(segments.begin(), segments.end(), gen);
    
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
      
      reduce_weights(weights);
      
      bcast_weights(weights);
      
      size_t updated_total = 0;
      MPI::COMM_WORLD.Allreduce(&updated, &updated_total, 1, utils::mpi_traits<size_t>::data_type(), MPI::SUM);
      
      weights *= 1.0 / updated_total;
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
      
      bcast_weights(rank_min, weights);
    } else if (mix_kbest_features > 0) {
      // reduce column-L2 weights
      weight_set_type weights_l2(weights);
      
      std::transform(weights_l2.begin(), weights_l2.end(), weights_l2.begin(), weights_l2.begin(), std::multiplies<weight_set_type::value_type>());
      
      reduce_weights(weights_l2);
      
      // synchronize here...
      MPI::COMM_WORLD.Barrier();
      
      // reduced rank-averaged weights
      weights *= 1.0 / mpi_size;
      
      reduce_weights(weights);
      
      if (mpi_rank == 0) {
	typedef std::pair<double, feature_type::id_type> value_type;
	typedef std::vector<value_type, std::allocator<value_type> > heap_type;
	
	// compute k-best wrt column-L2
	
	heap_type heap;
	
	heap.reserve(weights_l2.size());
	
	for (feature_type::id_type id = 0; id != weights_l2.size(); ++ id)
	  if (! feature_type(id).empty() && weights_l2[id] != 0.0) {
	    heap.push_back(value_type(weights_l2[id], id));
	    std::push_heap(heap.begin(), heap.end(), std::less<value_type>());
	  }
	
	if (static_cast<int>(heap.size()) > mix_kbest_features) {
	  typedef std::vector<bool, std::allocator<bool> > survived_type;
	  
	  survived_type survived(utils::bithack::max(weights.size(), weights_l2.size()), false);
	  
	  heap_type::iterator iter_begin = heap.begin();
	  heap_type::iterator iter       = heap.end();
	  
	  // kbest features
	  size_t num_survived = 0;
	  for (int k = 0; k != mix_kbest_features && iter_begin != iter; ++ k, -- iter) {
	    survived[iter_begin->second] = true;
	    std::pop_heap(iter_begin, iter, std::less<value_type>());
	    ++ num_survived;
	  }
	  
	  // also keep the tied features...
	  if (iter != iter_begin && iter != heap.end()) {
	    const double threshold = iter->first;

	    for (/**/; iter_begin != iter && iter_begin->first == threshold; -- iter) {
	      survived[iter_begin->second] = true;
	      std::pop_heap(iter_begin, iter, std::less<value_type>());
	      ++ num_survived;
	    }
	  }

	  if (debug && mpi_rank == 0)
	    std::cerr << "survived: " << num_survived << " all: "<< heap.size() << std::endl;
	  
	  for (feature_type::id_type id = 0; id != weights.size(); ++ id)
	    if (! survived[id])
	      weights[id] = 0.0;
	}
      } 
      
      // bcast weights!
      bcast_weights(weights);
    }
    
    if (dump_weights_mode && mpi_rank == 0)
      queue_dumper.push(std::make_pair(add_suffix(output_file, "." + utils::lexical_cast<std::string>(iter + 1)), weights));
  }

  // finish dumper...
  if (dumper.get()) {
    queue_dumper.push(std::make_pair(path_type(), weight_set_type()));
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

void recv_weights(const int rank, weight_set_type& weights)
{
  weights.clear();
  weights.allocate();
  
  boost::iostreams::filtering_istream is;
  is.push(codec::lz4_decompressor());
  is.push(utils::mpi_device_source(rank, weights_tag, 1024 * 1024));
  
  std::string feature;
  std::string value;
  
  while ((is >> feature) && (is >> value))
    weights[feature] = utils::decode_base64<double>(value);
}

void send_weights(const int rank, const weight_set_type& weights)
{
  boost::iostreams::filtering_ostream os;
  os.push(codec::lz4_compressor());
  os.push(utils::mpi_device_sink(rank, weights_tag, 1024 * 1024));
  
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
  is.push(codec::lz4_decompressor());
  is.push(utils::mpi_device_source(rank, weights_tag, 1024 * 1024));
  
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
    device.push_back(device_ptr_type(new device_type(*first, weights_tag, 1024 * 1024)));
    stream.push_back(stream_ptr_type(new stream_type()));
    
    stream.back()->push(codec::lz4_decompressor());
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
    os.push(codec::lz4_compressor());
    os.push(utils::mpi_device_bcast_sink(rank, 1024 * 1024));
    
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
    is.push(codec::lz4_decompressor());
    is.push(utils::mpi_device_bcast_source(rank, 1024 * 1024));
    
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
    ("oracle", po::value<path_type>(&oracle_file)->default_value(oracle_file), "oracle file")
    
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
    ("refset",  po::value<path_type>(&refset_file),  "refset")
    ("output",  po::value<path_type>(&output_file),  "model (or weights) output")
    ("weights", po::value<path_type>(&weights_file), "initial model (or weights)")
    
    ("scorer",           po::value<std::string>(&scorer_name)->default_value(scorer_name), "evaluation scorer")
    ("yield-sentence",   po::bool_switch(&yield_sentence),                                 "sentence yield")
    ("yield-alignment",  po::bool_switch(&yield_alignment),                                "alignment yield")
    ("yield-dependency", po::bool_switch(&yield_dependency),                               "dependency yield")
    
    ("iteration",     po::value<int>(&iteration)->default_value(iteration),   "learning iterations")
    ("batch",         po::value<int>(&batch_size)->default_value(batch_size), "batch (or batch, bin) size")
    ("kbest",         po::value<int>(&kbest_size)->default_value(kbest_size), "kbest size")
    ("kbest-diverse", po::bool_switch(&kbest_diverse_mode),                   "non unique kbest")
    
    ("learn-xbleu",    po::bool_switch(&learn_xbleu),    "online SGD with xBLEU loss")
    ("learn-softmax",  po::bool_switch(&learn_softmax),  "online SGD with softmax loss")
    ("learn-osoftmax", po::bool_switch(&learn_osoftmax), "online optimized-SGD with softmax loss")
    ("learn-el",       po::bool_switch(&learn_el),       "online SGD with expected loss")
    ("learn-oel",      po::bool_switch(&learn_oel),      "online optimized-SGD with expected loss")    
    ("learn-hinge",    po::bool_switch(&learn_hinge),    "online SGD with hinge loss (Pegasos)")
    ("learn-ohinge",   po::bool_switch(&learn_ohinge),   "online optimized-SGD with hinge loss (optimized-Pegasos)")
    ("learn-mira",     po::bool_switch(&learn_mira),     "online MIRA algorithm")
    ("learn-pa",       po::bool_switch(&learn_pa),       "online PA algorithm")
    ("learn-cw",       po::bool_switch(&learn_cw),       "online CW algorithm")
    ("learn-arow",     po::bool_switch(&learn_arow),     "online AROW algorithm")
    ("learn-nherd",    po::bool_switch(&learn_nherd),    "online NHERD algorithm")
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-regularization")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-regularization")
    ("C",             po::value<double>(&C)->default_value(C),                     "regularization constant")
    ("oscar",         po::value<double>(&oscar)->default_value(oscar),             "OSCAR regularization constant")
    ("temperature",   po::value<double>(&temperature)->default_value(temperature), "temperature")
    ("scale",         po::value<double>(&scale)->default_value(scale),             "scaling for weight")
    ("alpha0",        po::value<double>(&alpha0)->default_value(alpha0),           "\\alpha_0 for decay")
    ("eta0",          po::value<double>(&eta0)->default_value(eta0),               "\\eta_0 for decay")
    ("order",         po::value<int>(&order)->default_value(order),                "ngram order for xBLEU")
    
    ("rate-exponential", po::bool_switch(&rate_exponential),  "exponential learning rate")
    ("rate-simple",      po::bool_switch(&rate_simple),       "simple learning rate")
    ("rate-adagrad",     po::bool_switch(&rate_adagrad),      "adaptive learning rate (AdaGrad)")
    
    ("rda", po::bool_switch(&rda_mode), "RDA method for optimization (regularized dual averaging method)")

    ("loss-rank",           po::bool_switch(&loss_rank),            "rank loss")
    ("softmax-margin",      po::bool_switch(&softmax_margin),       "softmax margin")
    ("project-weight",      po::bool_switch(&project_weight),       "project L2 weight")
    ("merge-oracle",        po::bool_switch(&merge_oracle_mode),    "merge oracle kbests")
    ("merge-history",       po::value<int>(&merge_history),         "merge history for decoded results")
    ("mix-none",            po::bool_switch(&mix_none_mode),        "no mixing")
    ("mix-average",         po::bool_switch(&mix_average_mode),     "mixing weights by averaging")
    ("mix-select",          po::bool_switch(&mix_select_mode),      "select weights by L1")
    ("mix-kbest-features",  po::value<int>(&mix_kbest_features),    "mix k-best features")
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

