//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//
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

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>

#include "cicada_lexicon_global_impl.hpp"

#include "utils/lexical_cast.hpp"

path_type source_file = "-";
path_type target_file = "-";

path_type output_file;

int    vocab_size = 0;
int    max_iteration = 100;
bool   regularize_l1 = false;
bool   regularize_l2 = false;
double C = 1.0;

bool learn_maxent = false;
bool learn_sgd = false;
bool learn_mira = false;
bool learn_cw = false;
bool learn_arow = false;

int threads = 1;
int debug = 0;

int getoptions(int argc, char** argv);

void learn(const path_type& path,
	   const bitext_set_type& bitexts,
	   const word_set_type& vocab);


int main(int argc, char** argv)
{
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
    
    threads = (threads <= 0 ? 1 : threads);
    
    bitext_set_type bitexts;
    word_set_type   vocab;
    
    // read bitexts...
    read_bitexts(source_file, target_file, bitexts, vocab, vocab_size);
    
    // learn... do we dump parameters at the same time? (in single presision float?)
    learn(output_file, bitexts, vocab);
    
    // perform indexing...?
    
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
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
  typedef std::vector<thread_ptr_type, std::allocator<thread_ptr_type> > thread_ptr_set_type;

  {
    repository_type rep(path, repository_type::write);
    
    rep["bitext-size"] = utils::lexical_cast<std::string>(bitexts.size());
    rep["vocab-size"]  = utils::lexical_cast<std::string>(vocab.size());
    rep["size"] = utils::lexical_cast<std::string>(threads);
  }
  
  if (debug)
    std::cerr << "bitexts: " << bitexts.size() << std::endl
	      << "vocabulary: " << vocab.size() << std::endl;
  
  repository_type rep(path, repository_type::read);
  
  queue_type queue(threads);
  thread_ptr_set_type mapper(threads);
  
  for (int i = 0; i < mapper.size(); ++ i) {
    const path_type path_lexicon = rep.path(std::string("lexicon.") + utils::lexical_cast<std::string>(i) + ".gz");
    
    if (learn_maxent)
      mapper[i].reset(new thread_type(Mapper<OptimizeLBFGS>(bitexts,
							    queue,
							    path_lexicon)));
    if (learn_sgd)
      mapper[i].reset(new thread_type(Mapper<OptimizeSGD>(bitexts,
							  queue,
							  path_lexicon)));
    else if (learn_mira)
      mapper[i].reset(new thread_type(Mapper<OptimizeMIRA>(bitexts,
							   queue,
							   path_lexicon)));
    else if (learn_cw)
      mapper[i].reset(new thread_type(Mapper<OptimizeCW>(bitexts,
							 queue,
							 path_lexicon)));
    else if (learn_arow)
      mapper[i].reset(new thread_type(Mapper<OptimizeAROW>(bitexts,
							   queue,
							   path_lexicon)));
  }
  
  size_t num_processed = 0;
  for (word_set_type::const_iterator viter = vocab.begin(); viter != vocab.end(); ++ viter, ++ num_processed) {
    queue.push(*viter);
    
    if (debug && num_processed > 0) {
      const int prev = int(100.0 * num_processed / vocab.size());
      const int curr = int(100.0 * (num_processed + 1) / vocab.size());
      
      if (prev != curr)
	std::cerr << '.';
    }
  }
  if (debug)
    std::cerr << '\n';

  
  for (int i = 0; i < mapper.size(); ++ i)
    queue.push(word_type());
  
  for (int i = 0; i < mapper.size(); ++ i)
    mapper[i]->join();
}


int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source language corpus")
    ("target", po::value<path_type>(&target_file), "target language corpus")
    
    ("output", po::value<path_type>(&output_file), "output model")
    
    ("vocab", po::value<int>(&vocab_size), "vocabulary size")

    ("learn-maxent", po::bool_switch(&learn_maxent),  "maximum entropy (default)")
    ("learn-sgd",    po::bool_switch(&learn_sgd),     "SGD (Pegasos for L2, FOBOS-cummulative for L1)")
    ("learn-mira",   po::bool_switch(&learn_mira),    "MIRA")
    ("learn-cw",     po::bool_switch(&learn_cw),      "Confidence-Weighted")
    ("learn-arow",   po::bool_switch(&learn_arow),    "Adaptive-Regularization")
    
    ("max-iteration", po::value<int>(&max_iteration),  "maximum iteration")
    ("regularize-l1", po::bool_switch(&regularize_l1), "L1-norm (only for logistic-loss, maxent/SGD")
    ("regularize-l2", po::bool_switch(&regularize_l2), "L2-norm (default)")
    ("C",             po::value<double>(&C),           "C for regularization")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
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



