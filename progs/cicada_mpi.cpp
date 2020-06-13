//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <cstdlib>

#include "cicada_impl.hpp"
#include "cicada_output_impl.hpp"

#include "cicada/eval/score.hpp"
#include "cicada/format.hpp"
#include "cicada/signature.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/tokenizer.hpp"
#include "cicada/matcher.hpp"

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/program_options.hpp"
#include "utils/filesystem.hpp"
#include "utils/base64.hpp"
#include "utils/space_separator.hpp"
#include "utils/random_seed.hpp"
#include "utils/getline.hpp"

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/thread.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file = "-";

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

int debug = 0;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

void cicada_file(operation_set_type& operations);
void cicada_directory(operation_set_type& operations);
void synchronize();
void merge_features();
void merge_statistics(const operation_set_type& operations, operation_set_type::statistics_type& statistics);

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
				  debug);

    if (mpi_rank == 0 && debug)
      std::cerr << "operations: " << operations.size() << std::endl;
    
    // make sure to synchronize here... otherwise, badthink may happen...
    if (mpi_rank == 0 && ! operations.get_output_data().directory.empty())
      prepare_directory(operations.get_output_data().directory);
    
    synchronize();
    
    ::sync();
    
    if (! operations.get_output_data().file.empty())
      cicada_file(operations);
    else
      cicada_directory(operations);
    
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
  sample_tag = 1000,
  result_tag,
  notify_tag,
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

template <typename Tp, typename Alloc=std::allocator<Tp> >
struct single_queue
{
  typedef utils::lockfree_list_queue<Tp, Alloc > queue_type;
  
  single_queue(size_t size=0) : queue(1), busy(0) {}
  
  bool empty() { return ! utils::atomicop::fetch_and_add(busy, int(0)) && queue.empty(); } 
  void ready() { busy = 0; }

  void wait_empty()
  {
    queue.wait_empty();
    
    for (;;) {
      for (int i = 0; i < 50; ++ i) {
	if (empty())
	  return;
	else
	  boost::thread::yield();
      }
      
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    }
  }

  void push(const Tp& x)
  {
    busy = 1;
    queue.push(x);
  }
  
  void push_swap(Tp& x)
  {
    busy = 1;
    queue.push_swap(x);
  }
  
  void pop(Tp& x)
  {
    queue.pop(x);
  }
  
  void pop_swap(Tp& x)
  {
    queue.pop_swap(x);
  }
  
private:
  queue_type queue;
  int busy;
};

struct MapFile
{
  typedef std::pair<std::string, bool> value_type;
  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;
  
  path_type   path;
  queue_type& queue;
  
  MapFile(const path_type& _path, queue_type& _queue)
    : path(_path), queue(_queue) {}
  
  void operator()()
  {
    if (input_directory_mode) {
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_input = path / file_name;
	
	if (! boost::filesystem::exists(path_input)) break;
	
	queue.push(std::make_pair(path_input.string(), false));
      }
    } else {
      utils::compress_istream is(path, 1024 * 1024);
      
      operation_set_type::operation_type::id_type id = 0;
      std::string line;
      
      while (utils::getline(is, line)) {
	if (input_id_mode) {
	  if (line.empty())
	    throw std::runtime_error("invalid empty input!");
	  
	  queue.push(std::make_pair(line, false));
	} else
	  queue.push(std::make_pair(utils::lexical_cast<std::string>(id) + " ||| " + line, false));
	
	
	++ id;
      }
    }
    
    queue.push(std::make_pair(std::string(), true));
  }
};


struct TaskFile
{
  typedef single_queue<std::string, std::allocator<std::string> > queue_single_type;
  //typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_single_type;
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  TaskFile(queue_single_type&   __queue_is,
	     queue_type&   __queue_os,
	     operation_set_type& __operations)
    : queue_is(__queue_is),
      queue_os(__queue_os),
      operations(__operations) {}

  void operator()()
  {
    if (input_directory_mode) {
      typedef boost::spirit::istream_iterator iter_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      std::string file;
      std::string line;
      
      while (1) {
	file.clear();
	queue_is.pop_swap(file);
	if (file.empty()) break;
	
	utils::compress_istream is(file, 1024 * 1024);
	
	if (utils::getline(is, line) && ! line.empty())
	  operations(line);
	else
	  throw std::runtime_error("invalid file? " + file);
	
	queue_is.ready();
	
	queue_os.push(utils::lexical_cast<std::string>(operations.get_data().id) + ' ' + operations.get_output_data().buffer);
      }
    } else {
      std::string line;
      
      while (1) {
	line.clear();
	queue_is.pop_swap(line);
	if (line.empty()) break;
	
	operations(line);
	
	queue_is.ready();
	
	queue_os.push(utils::lexical_cast<std::string>(operations.get_data().id) + ' ' + operations.get_output_data().buffer);
      }
    }
    
    operations.clear();
    const_cast<operation_set_type::data_type&>(operations.get_data()).clear();
    
    queue_os.push(std::string());
  }
  
  queue_single_type& queue_is;
  queue_type&        queue_os;
  operation_set_type& operations;
};

struct ReduceFile
{
  typedef TaskFile::queue_type queue_type;
  
  ReduceFile(queue_type& __queue, const path_type& __path)
    : queue(__queue), path(__path) {}
  
  void operator()()
  {
    typedef operation_set_type::operation_type::id_type id_type;
    typedef std::map<id_type, std::string, std::less<id_type>, std::allocator<std::pair<const id_type, std::string> > > buffer_map_type;
    
    buffer_map_type maps;
    std::string buffer;
    
    id_type     id = 0;
    
    const bool flush_output = (path == "-"
			       || (boost::filesystem::exists(path)
				   && ! boost::filesystem::is_regular_file(path)));

    utils::compress_ostream os(path, 1024 * 1024);
        
    for (;;) {
      buffer.clear();
      queue.pop_swap(buffer);
      
      if (buffer.empty()) break;

      bool dump = false;
      
      utils::piece buffer_piece(buffer);
      
      utils::piece::const_iterator iter = buffer_piece.begin();
      for (/**/; iter != buffer_piece.end() && ! std::isspace(*iter); ++ iter);
      
      // tokenize here...
      const id_type      buffer_id        = utils::lexical_cast<id_type>(buffer_piece.substr(0, iter - buffer_piece.begin()));
      const utils::piece buffer_tokenized = buffer_piece.substr(iter + 1 - buffer_piece.begin());
      
      if (buffer_id == id) {
	os << buffer_tokenized;
	dump = true;
	
	++ id;
      } else
	maps[buffer_id] = static_cast<std::string>(buffer_tokenized);
      
      for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; /**/) {
	os << iter->second;
	dump = true;
	maps.erase(iter ++);
	++ id;
      }
      
      if (dump && flush_output)
	os << std::flush;
    }
    
    for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; /**/) {
      os << iter->second;
      maps.erase(iter ++);
      ++ id;
    }
    
    // we will do twice, in case we have wrap-around for id...!
    if (! maps.empty())
      for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; /**/) {
	os << iter->second;
	maps.erase(iter ++);
	++ id;
      }
    
    if (flush_output)
      os << std::flush;
    
    if (! maps.empty())
      throw std::runtime_error("id mismatch! expecting: " + utils::lexical_cast<std::string>(id)
			       + " next: " + utils::lexical_cast<std::string>(maps.begin()->first)
			       + " renamining: " + utils::lexical_cast<std::string>(maps.size()));
  }
  
  queue_type& queue;
  path_type   path;
};

void cicada_file(operation_set_type& operations)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  typedef TaskFile   task_type;
  
  task_type::queue_single_type queue_is(1);
  task_type::queue_type        queue_os;
  
  boost::thread thread(task_type(queue_is, queue_os, operations));
  
  if (mpi_rank == 0) {
    typedef MapFile    map_type;
    typedef ReduceFile reduce_type;
    
    typedef utils::mpi_ostream        ostream_type;
    typedef utils::mpi_istream_simple istream_type;
    
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    
    typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;

    typedef map_type::queue_type queue_input_type;
    
    queue_input_type queue_input(mpi_size);
    
    boost::thread thread_reduce(reduce_type(queue_os, operations.get_output_data().file));
    boost::thread thread_map(map_type(input_file, queue_input));

    ostream_ptr_set_type ostream(mpi_size);
    istream_ptr_set_type istream(mpi_size);
    
    for (int rank = 1; rank < mpi_size; ++ rank) {
      ostream[rank].reset(new ostream_type(rank, sample_tag, 4096));
      istream[rank].reset(new istream_type(rank, result_tag, 4096));
    }

    std::string line;
    map_type::value_type line_input(std::string(), false);
    
    int non_found_iter = 0;
    
    while (! line_input.second) {
      bool found = false;
      
      for (int rank = 1; rank < mpi_size && ! line_input.second; ++ rank)
	if (ostream[rank]->test() && queue_input.pop(line_input, true) && ! line_input.second) {
	  ostream[rank]->write(line_input.first);
	  
	  found = true;
	}
      
      if (queue_is.empty() && queue_input.pop(line_input, true) && ! line_input.second) {
	queue_is.push(line_input.first);
	
	found = true;
      }
      
      // reduce from others...
      for (int rank = 1; rank < mpi_size; ++ rank)
	if (istream[rank] && istream[rank]->test()) {
	  if (istream[rank]->read(line))
	    queue_os.push_swap(line);
	  else
	    istream[rank].reset();
	  
	  found = true;
	}
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }

    for (;;) {
      bool found = false;
      
      // termination!
      for (int rank = 1; rank < mpi_size; ++ rank)
	if (ostream[rank] && ostream[rank]->test()) {
	  if (! ostream[rank]->terminated())
	    ostream[rank]->terminate();
	  else
	    ostream[rank].reset();
	  
	  found = true;
	}
      
      // reduce from others...
      for (int rank = 1; rank < mpi_size; ++ rank)
	if (istream[rank] && istream[rank]->test()) {
	  if (istream[rank]->read(line))
	    queue_os.push_swap(line);
	  else
	    istream[rank].reset();
	  
	  found = true;
	}

      // termination condition!
      if (std::count(istream.begin(), istream.end(), istream_ptr_type()) == mpi_size
	  && std::count(ostream.begin(), ostream.end(), ostream_ptr_type()) == mpi_size
	  && queue_is.empty()) {
	queue_is.push(std::string());
	
	break;
      }
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    thread_map.join();
    thread_reduce.join();
    
  } else {
    boost::shared_ptr<utils::mpi_istream>        is(new utils::mpi_istream(0, sample_tag, 4096));
    boost::shared_ptr<utils::mpi_ostream_simple> os(new utils::mpi_ostream_simple(0, result_tag, 4096));
    
    std::string line;
    bool terminated = false;
    
    int non_found_iter = 0;
    for (;;) {
      bool found = false;
      
      if (is && is->test() && queue_is.empty()) {
	if (is->read(line))
	  queue_is.push_swap(line);
	else {
	  queue_is.push(std::string());
	  is.reset();
	}
	found = true;
      }
      
      if (! terminated) {
	line.clear();
	if (os && os->test() && queue_os.pop_swap(line, true)) {
	  if (line.empty())
	    terminated = true;
	  else
	    os->write(line);
	  found = true;
	}
      } else {
	if (os && os->test()) {
	  if (! os->terminated())
	    os->terminate();
	  else
	    os.reset();
	  found = true;
	}
      }
      
      if (! is && ! os) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
  
  thread.join();
}

struct Task
{
  typedef single_queue<std::string, std::allocator<std::string> > queue_type;
  //typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  Task(queue_type&   __queue,
       operation_set_type& __operations)
    : queue(__queue),
      operations(__operations) {}

  void operator()()
  {
    if (input_directory_mode) {
      std::string file;
      std::string line;
      
      while (1) {
	file.clear();
	queue.pop_swap(file);
	if (file.empty()) break;

	utils::compress_istream is(file, 1024 * 1024);
	
	if (utils::getline(is, line) && ! line.empty())
	  operations(line);
	else
	  throw std::runtime_error("invalid file? " + file);
	
	queue.ready();
      }
    } else {
      std::string line;
      
      while (1) {
	line.clear();
	queue.pop_swap(line);
	if (line.empty()) break;

	operations(line);
	
	queue.ready();
      }
    }
    
    operations.clear();
    const_cast<operation_set_type::data_type&>(operations.get_data()).clear();
  }
  
  queue_type&   queue;
  operation_set_type& operations;
};

void cicada_directory(operation_set_type& operations)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  typedef Task  task_type;
  
  typedef task_type::queue_type queue_type;

  queue_type queue(1);
  
  boost::thread thread(task_type(queue, operations));
  
  if (mpi_rank == 0) {
    typedef utils::mpi_ostream ostream_type;
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    
    std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > stream(mpi_size);
    
    for (int rank = 1; rank < mpi_size; ++ rank)
      stream[rank].reset(new ostream_type(rank, sample_tag, 4096));
    
    if (input_directory_mode) {
      boost::filesystem::directory_iterator iter_end;
      boost::filesystem::directory_iterator iter(input_file);
      
      int non_found_iter = 0;
      while (iter != iter_end) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size && iter != iter_end; ++ rank)
	  if (stream[rank]->test()) {
	    const path_type path(*iter);
	    ++ iter;
	    
	    stream[rank]->write(path.string());
	    
	    found = true;
	  }
	
	if (queue.empty() && iter != iter_end) {
	  const path_type path(*iter);
	  ++ iter;
	  
	  queue.push(path.string());
	  
	  found = true;
	}
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }

    } else {
      utils::compress_istream is(input_file, 1024 * 1024);
      
      operation_set_type::operation_type::id_type id = 0;
      std::string line;
      
      int non_found_iter = 0;
      
      while (is) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size && is; ++ rank)
	  if (stream[rank]->test() && utils::getline(is, line)) {
	    if (input_id_mode) {
	      if (line.empty())
		throw std::runtime_error("invalid empty input!");
	      
	      stream[rank]->write(line);
	    } else
	      stream[rank]->write(utils::lexical_cast<std::string>(id) + " ||| " + line);
	    
	    ++ id;
	    
	    found = true;
	  }
	
	if (queue.empty() && utils::getline(is, line)) {
	  if (input_id_mode) {
	    if (line.empty())
	      throw std::runtime_error("invalid empty input!");
	    
	    queue.push(line);
	  } else
	    queue.push(utils::lexical_cast<std::string>(id) + " ||| " + line);
	  
	  ++ id;
	  
	  found = true;
	}
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
    }
    
    int non_found_iter = 0;
    bool terminated = false;
    while (1) {
      bool found = false;
      
      if (! terminated && queue.empty()) {
	queue.push(std::string());
	
	terminated = true;
	found = true;
      }
      
      for (int rank = 1; rank < mpi_size; ++ rank) 
	if (stream[rank] && stream[rank]->test()) {
	  if (! stream[rank]->terminated())
	    stream[rank]->terminate();
	  else
	    stream[rank].reset();
	  
	  found = true;
	}
      
      if (terminated && std::count(stream.begin(), stream.end(), ostream_ptr_type()) == mpi_size)
	break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
  } else {
    utils::mpi_istream is(0, sample_tag, 4096, true);
    
    std::string line;
    while (is.read(line)) {
      queue.push_swap(line);
      queue.wait_empty();
      is.ready();
    }
    queue.wait_empty();
    queue.push(std::string());
  }
  
  thread.join();
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
    ("feature-function",        po::value<feature_parameter_set_type >(&feature_parameters)->composing(), "feature function(s)")
    ("feature-function-list",   po::bool_switch(&feature_list),                                           "list of available feature function(s)")
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
