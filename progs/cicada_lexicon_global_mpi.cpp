// global lexicon learning...
//
// First, compile vocabulary
// Second, learn classifier for each target-word... (we will use "threads" for distributed computing)
//
// each thread will learn a classifier for each target-word and dump results
//    when dumping, we will threashold 
//    when dumping, we will dump in a directory, "target source parameter" format (similar to lexicon model of moses)
//
// supports: L1/L2 norm
//

// MPI version...!

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include "learn_global_lexicon_impl.hpp"

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"


path_type source_file = "-";
path_type target_file = "-";

path_type output_file;

int    max_iteration = 100;
bool   regularize_l1 = false;
bool   regularize_l2 = false;
double C = 1.0;

bool learn_maxent = false;
bool learn_sgd = false;
bool learn_mira = false;
bool learn_cw = false;
bool learn_arow = false;

path_type prog_name;

int debug = 0;

int getoptions(int argc, char** argv);

void learn(const path_type& path,
	   const bitext_set_type& bitexts,
	   const word_set_type& vocab);


int main(int argc, char** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;
    
    switch (int(learn_maxent) + learn_sgd + learn_mira + learn_cw + learn_arow) {
    case 0: learn_maxent = true; break;
    case 1: break;
    default: throw std::runtime_error("you can specily one of learn_{maxent,sgd,mira,cw,arow}");
    }
    
    // setup regularizer... default is L2
    if (int(regularize_l1) + regularize_l2 == 0)
      regularize_l2 = true;
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("you can specify either L1 or L2");
    
    if (C <= 0.0)
      throw std::runtime_error("C must be greater than zero");
    
    if (output_file.empty())
      throw std::runtime_error("no output path");
    
    srandom(getpid() * time(0));
    
    bitext_set_type bitexts;
    word_set_type   vocab;
    
    // read bitexts...
    read_bitexts(source_file, target_file, bitexts, vocab);
    
    learn(output_file, bitexts, vocab);
  }
  catch (std::exception& err) {
    if (mpi_rank == 0)
      std::cerr << "error: " << err.what() << std::endl;
    
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  return 0;
}


enum {
  word_tag = 1000,
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

void learn(const path_type& path,
	   const bitext_set_type& bitexts,
	   const word_set_type& vocab)
{
  typedef utils::repository repository_type;
  
  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::thread_type thread_type;
  typedef map_reduce_type::queue_type  queue_type;

  typedef boost::shared_ptr<thread_type> thread_ptr_type;

  const int mpi_size = MPI::COMM_WORLD.Get_size();
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  
  if (mpi_rank == 0) {
    repository_type rep(path, repository_type::write);
    
    rep["bitext-size"] = boost::lexical_cast<std::string>(bitexts.size());
    rep["vocab-size"]  = boost::lexical_cast<std::string>(vocab.size());
    rep["size"] = boost::lexical_cast<std::string>(mpi_size);
  }

  if (mpi_rank == 0 && debug)
    std::cerr << "bitexts: " << bitexts.size() << std::endl
	      << "vocabulary: " << vocab.size() << std::endl;

  MPI::COMM_WORLD.Barrier();
  
  while (! boost::filesystem::exists(path))
    boost::thread::yield();
  
  repository_type rep(path, repository_type::read);
  
  queue_type      queue(1);
  const path_type path_lexicon = rep.path(std::string("lexicon.") + boost::lexical_cast<std::string>(mpi_rank) + ".gz");
  
  thread_ptr_type mapper;
  if (learn_maxent)
    mapper.reset(new thread_type(Mapper<OptimizeLBFGS>(bitexts,
						       queue,
						       path_lexicon)));
  if (learn_sgd)
    mapper.reset(new thread_type(Mapper<OptimizeSGD>(bitexts,
						     queue,
						     path_lexicon)));
  else if (learn_mira)
    mapper.reset(new thread_type(Mapper<OptimizeMIRA>(bitexts,
						      queue,
						      path_lexicon)));
  else if (learn_cw)
    mapper.reset(new thread_type(Mapper<OptimizeCW>(bitexts,
						    queue,
						    path_lexicon)));
  else if (learn_arow)
    mapper.reset(new thread_type(Mapper<OptimizeAROW>(bitexts,
						      queue,
						      path_lexicon)));

  if (mpi_rank == 0) {
    typedef utils::mpi_ostream ostream_type;
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
    
    ostream_ptr_set_type stream(mpi_size);
    for (int rank = 1; rank < mpi_size; ++ rank)
      stream[rank].reset(new ostream_type(rank, word_tag, 4096));
    
    word_set_type::const_iterator viter = vocab.begin();
    word_set_type::const_iterator viter_end = vocab.end();
    
    int non_found_iter = 0;
    size_t num_processed = 0;
    
    while (viter != viter_end) {
      bool found = false;
      
      for (int rank = 1; rank < mpi_size && viter != viter_end; ++ rank)
	if (stream[rank]->test()) {
	  stream[rank]->write(*viter);
	  ++ viter;
	  
	  if (debug && num_processed > 0) {
	    const int prev = int(100.0 * num_processed / vocab.size());
	    const int curr = int(100.0 * (num_processed + 1) / vocab.size());
	    
	    if (prev != curr)
	      std::cerr << '.';
	  }
	  ++ num_processed;
	  
	  found = true;
	}
      
      if (queue.empty() && viter != viter_end) {
	queue.push(*viter);
	++ viter;

	if (debug && num_processed > 0) {
	  const int prev = int(100.0 * num_processed / vocab.size());
	  const int curr = int(100.0 * (num_processed + 1) / vocab.size());
	  
	  if (prev != curr)
	    std::cerr << '.';
	}
	++ num_processed;
	
	found = true;
      }
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    if (debug)
      std::cerr << '\n';
    
    //terminating...
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
      
      if (std::count(stream.begin(), stream.end(), ostream_ptr_type()) == mpi_size) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
  } else {
    utils::mpi_istream is(0, word_tag, 4096, true);
    
    std::string word;
    while (is.read(word)) {
      queue.push(word);
      is.ready();
    }
  }
  
  queue.push(word_type());
  mapper->join();
}


int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source language corpus")
    ("target", po::value<path_type>(&target_file), "target language corpus")
    
    ("output", po::value<path_type>(&output_file), "output model")

    ("learn-maxent", po::bool_switch(&learn_maxent),  "maximum entropy (default)")
    ("learn-sgd",    po::bool_switch(&learn_sgd),     "SGD (Pegasos for L2, FOBOS-cummulative for L1)")
    ("learn-mira",   po::bool_switch(&learn_mira),    "MIRA")
    ("learn-cw",     po::bool_switch(&learn_cw),      "Confidence-Weighted")
    ("learn-arow",   po::bool_switch(&learn_arow),    "Adaptive-Regularization")
    
    ("max-iteration", po::value<int>(&max_iteration),  "maximum iteration")
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-norm (only for logistic-loss, maxent/SGD)")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-norm (default)")
    ("C",             po::value<double>(&C),           "C for regularization")
    
    ("prog", po::value<path_type>(&prog_name), "this prog's path")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}



