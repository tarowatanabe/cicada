//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unistd.h>

#include "cicada_impl.hpp"

#include "utils/mpi.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/program_options.hpp"

#include <boost/program_options.hpp>
#include <boost/thread.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file = "-";

bool input_id_mode = false;
bool input_bitext_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_span_mode = false;
bool input_directory_mode = false;

std::string symbol_goal         = vocab_type::S;
std::string symbol_non_terminal = vocab_type::X;
path_type   symbol_fallback_file;

grammar_file_set_type grammar_mutable_files;
grammar_file_set_type grammar_static_files;

bool grammar_glue_straight = false;
bool grammar_glue_inverted = false;
bool grammar_insertion = false;
bool grammar_deletion = false;

grammar_file_set_type tree_grammar_mutable_files;
grammar_file_set_type tree_grammar_static_files;

bool tree_grammar_fallback = false;

feature_parameter_set_type feature_parameters;
bool feature_list = false;

op_set_type ops;
bool op_list = false;

int debug = 0;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

void cicada_stdout(operation_set_type& operations);
void cicada_process(operation_set_type& operations);
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


    // read grammars...
    grammar_type grammar;
    const size_t grammar_static_size  = load_grammar<cicada::GrammarStatic>(grammar, grammar_static_files);
    const size_t grammar_mutable_size = load_grammar<cicada::GrammarMutable>(grammar, grammar_mutable_files);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "loaded static grammar: " << grammar_static_size << std::endl
		<< "loaded mutable grammar: " << grammar_mutable_size << std::endl;

    
    if (grammar_glue_straight || grammar_glue_inverted) {
      if (! symbol_fallback_file.empty()) {
	if (symbol_fallback_file != "-" && ! boost::filesystem::exists(symbol_fallback_file))
	  throw std::runtime_error("invalid fallback non-terminal file: " + symbol_fallback_file.file_string());
	
	utils::compress_istream is(symbol_fallback_file, 1024 * 1024);
	grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										    symbol_non_terminal,
										    std::istream_iterator<std::string>(is),
										    std::istream_iterator<std::string>(),
										    grammar_glue_straight,
										    grammar_glue_inverted)));
	
      } else
	grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										    symbol_non_terminal,
										    grammar_glue_straight,
										    grammar_glue_inverted)));
    }

    if (debug && mpi_rank == 0)
      std::cerr << "grammar: " << grammar.size() << std::endl;

    tree_grammar_type tree_grammar;
    const size_t tree_grammar_static_size  = load_grammar<cicada::TreeGrammarStatic>(tree_grammar, tree_grammar_static_files);
    const size_t tree_grammar_mutable_size = load_grammar<cicada::TreeGrammarMutable>(tree_grammar, tree_grammar_mutable_files);
    
    if (debug)
      std::cerr << "loaded static tree grammar: " << tree_grammar_static_size << std::endl
		<< "loaded mutable tree grammar: " << tree_grammar_mutable_size << std::endl;
    
    if (debug)
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
				  symbol_non_terminal,
				  grammar_insertion,
				  grammar_deletion,
				  tree_grammar_fallback,
				  true,
				  input_lattice_mode,
				  input_forest_mode,
				  input_span_mode,
				  input_bitext_mode,
				  true,
				  debug);
    
    // make sure to synchronize here... otherwise, badthink may happen...
    if (mpi_rank == 0 && ! operations.get_output_data().directory.empty()) {
      const path_type& directory = operations.get_output_data().directory;
      
      if (boost::filesystem::exists(directory) && ! boost::filesystem::is_directory(directory))
	boost::filesystem::remove_all(directory);
      
      boost::filesystem::create_directories(directory);
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(directory); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);

      ::sync();
    }

    MPI::COMM_WORLD.Barrier();
    
    if (! operations.get_output_data().file.empty())
      cicada_stdout(operations);
    else
      cicada_process(operations);
    
    operations.clear();
    
    ::sync();

    synchronize();
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
    MPI::COMM_WORLD.Send(0, 0, MPI::INT, 0, notify_tag);
    MPI::COMM_WORLD.Recv(0, 0, MPI::INT, 0, notify_tag);
  }
}

struct MapStdout
{
  typedef std::pair<std::string, bool> value_type;
  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;
  
  path_type   path;
  queue_type& queue;
  
  MapStdout(const path_type& _path, queue_type& _queue)
    : path(_path), queue(_queue) {}
  
  void operator()()
  {
    if (input_directory_mode) {
      std::string line;
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(path); iter != iter_end; ++ iter) {
	utils::compress_istream is(*iter, 1024 * 1024);
	
	if (std::getline(is, line))
	  queue.push(std::make_pair(line, false));
      }
      
    } else {
      size_t id = 0;
      utils::compress_istream is(path, 1024 * 1024);
      
      std::string line;
      while (std::getline(is, line)) {
	
	if (input_id_mode)
	  queue.push(std::make_pair(line, false));
	else
	  queue.push(std::make_pair(utils::lexical_cast<std::string>(id) + " ||| " + line, false));
	
	++ id;
      }
    }
    
    queue.push(std::make_pair(std::string(), true));
  }
};


struct TaskStdout
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  TaskStdout(queue_type&   __queue_is,
	     queue_type&   __queue_os,
	     operation_set_type& __operations)
    : queue_is(__queue_is),
      queue_os(__queue_os),
      operations(__operations) {}

  void operator()()
  {
    std::string line;
    while (1) {
      queue_is.pop_swap(line);
      if (line.empty()) break;
      
      operations(line);
      
      queue_os.push(utils::lexical_cast<std::string>(operations.get_data().id) + ' ' + operations.get_output_data().buffer);
    }

    queue_os.push(std::string());
  }
  
  
  queue_type&   queue_is;
  queue_type&   queue_os;
  operation_set_type& operations;
};

struct ReduceStdout
{
  typedef TaskStdout::queue_type queue_type;
  
  ReduceStdout(queue_type& __queue, const path_type& __path)
    : queue(__queue), path(__path) {}
  
  void operator()()
  {
    typedef size_t id_type;
    typedef std::map<id_type, std::string, std::less<id_type>, std::allocator<std::pair<const id_type, std::string> > > buffer_map_type;
    
    buffer_map_type maps;
    std::string buffer;
    
    id_type     id = 0;
    
    utils::compress_ostream os(path, 1024 * 1024);
    
    for (;;) {
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
      
      for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; ++ id) {
	os << iter->second;
	dump = true;
	maps.erase(iter ++);
      }
      
      if (dump)
	os << std::flush;
    }
    
    for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; ++ id) {
      os << iter->second;
      maps.erase(iter ++);
    }
    
    // we will do twice, in case we have wrap-around for id...!
    if (! maps.empty())
      for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; ++ id) {
	os << iter->second;
	maps.erase(iter ++);
      }
    
    os << std::flush;
    
    if (! maps.empty())
      throw std::runtime_error("id mismatch! expecting: " + utils::lexical_cast<std::string>(id)
			       + " next: " + utils::lexical_cast<std::string>(maps.begin()->first));
  }
  
  queue_type& queue;
  path_type   path;
};

void cicada_stdout(operation_set_type& operations)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  typedef TaskStdout   task_type;
  
  typedef task_type::queue_type queue_type;
  
  queue_type queue_is(1);
  queue_type queue_os;
  
  boost::thread thread(task_type(queue_is, queue_os, operations));
  
  if (mpi_rank == 0) {
    typedef MapStdout    map_type;
    typedef ReduceStdout reduce_type;
    
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
	  && queue_is.push(std::string(), true)) break;
      
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
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  Task(queue_type&   __queue,
       operation_set_type& __operations)
    : queue(__queue),
      operations(__operations) {}

  void operator()()
  {
    std::string       line;
    while (1) {
      queue.pop_swap(line);
      if (line.empty()) break;
      
      operations(line);
    }
  }
  
  queue_type&   queue;
  operation_set_type& operations;
};

void cicada_process(operation_set_type& operations)
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
      
      std::string line;
      
      int non_found_iter = 0;
      while (iter != iter_end) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size && iter != iter_end; ++ rank)
	  if (stream[rank]->test()) {
	    utils::compress_istream is(*iter, 1024 * 1024);
	    ++ iter;
	    
	    if (std::getline(is, line))
	      stream[rank]->write(line);
	    
	    found = true;
	  }
	
	if (queue.empty() && iter != iter_end) {
	  utils::compress_istream is(*iter, 1024 * 1024);
	  ++ iter;
	  
	  if (std::getline(is, line))
	    queue.push(line);
	  
	  found = true;
	}
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }

    } else {
      utils::compress_istream is(input_file, 1024 * 1024);
      
      size_t id = 0;
      std::string line;
      
      int non_found_iter = 0;
      while (is) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size && is; ++ rank)
	  if (stream[rank]->test() && std::getline(is, line)) {
	    if (input_id_mode)
	      stream[rank]->write(line);
	    else
	      stream[rank]->write(utils::lexical_cast<std::string>(id) + " ||| " + line);
	    
	    ++ id;
	    
	    found = true;
	  }
	
	if (queue.empty() && std::getline(is, line)) {
	  if (input_id_mode)
	    queue.push(line);
	  else
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
      
      if (! terminated && queue.push(std::string(), true)) {
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
    queue.push(std::string());
  }
  
  thread.join();
}



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
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-span",       po::bool_switch(&input_span_mode),       "span input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")
    
    // grammar
    ("goal",           po::value<std::string>(&symbol_goal)->default_value(symbol_goal),                 "goal symbol")
    ("non-terminal",   po::value<std::string>(&symbol_non_terminal)->default_value(symbol_non_terminal), "default non-terminal symbol")
    ("fallback",       po::value<path_type>(&symbol_fallback_file),                                      "fallback non-terminal list")
    ("grammar",        po::value<grammar_file_set_type >(&grammar_mutable_files)->composing(),           "grammar file(s)")
    ("grammar-static", po::value<grammar_file_set_type >(&grammar_static_files)->composing(),            "static binary grammar file(s)")
    
    // special handling
    ("grammar-glue-straight", po::bool_switch(&grammar_glue_straight), "add straight hiero glue rule")
    ("grammar-glue-inverted", po::bool_switch(&grammar_glue_inverted), "add inverted hiero glue rule")
    ("grammar-insertion",     po::bool_switch(&grammar_insertion),     "source-to-target transfer grammar")
    ("grammar-deletion",      po::bool_switch(&grammar_deletion),      "source-to-<epsilon> transfer grammar")
    
    // tree-grammar
    ("tree-grammar",          po::value<grammar_file_set_type >(&tree_grammar_mutable_files)->composing(), "tree grammar file(s)")
    ("tree-grammar-static",   po::value<grammar_file_set_type >(&tree_grammar_static_files)->composing(),  "static binary tree grammar file(s)")
    
    // special handling
    ("tree-grammar-fallback", po::bool_switch(&tree_grammar_fallback),                                     "source-to-target transfer tree grammar")


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
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  if (variables.count("config")) {
    const path_type path_config = variables["config"].as<path_type>();
    if (! boost::filesystem::exists(path_config))
      throw std::runtime_error("no config file: " + path_config.file_string());
	
    utils::compress_istream is(path_config);
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {

    if (mpi_rank == 0)
      std::cout << argv[0] << " [options]\n"
		<< desc_command << std::endl;

    MPI::Finalize();
    exit(0);
  }
}
