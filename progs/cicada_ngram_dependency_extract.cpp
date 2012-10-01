//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// dependency ngram language model learning... actually, this will collect counts only
//

#include "cicada_ngram_dependency_extract_impl.hpp"

#include <iostream>
#include <sstream>

#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/functional/hash.hpp>
#include <boost/thread.hpp>

#include <utils/compress_stream.hpp>

typedef DependencyCounts::path_type     path_type;
typedef DependencyCounts::path_set_type path_set_type;

typedef DependencyCounts::count_type    count_type;
typedef DependencyCounts::word_type     word_type;
typedef DependencyCounts::vocab_type    vocab_type;

typedef DependencyCounts::count_path_set_type count_path_set_type;

path_type corpus_file;
path_type corpus_list_file;
path_type output_file;

bool map_line = false;
int threads = 2;
double max_malloc = 1.0; // 1G bytes

int debug = 0;

void accumulate(const path_set_type& corpus_files,
		const path_type& output_file,
		count_path_set_type& count_paths);

int getoptions(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;
    
    if (output_file.empty())
      throw std::runtime_error("no output file");
    
    threads = utils::bithack::max(1, threads);
    
    path_set_type corpus_files;
    
    if (corpus_file == "-" || boost::filesystem::exists(corpus_file))
      corpus_files.push_back(corpus_file);
    
    if (corpus_list_file == "-" || boost::filesystem::exists(corpus_list_file)) {
      utils::compress_istream is(corpus_list_file);
      std::string line;
      while (std::getline(is, line)) {
	boost::algorithm::trim(line);
	if (! line.empty()) {
	  if (boost::filesystem::exists(line))
	    corpus_files.push_back(line);
	  else if (boost::filesystem::exists(corpus_list_file.parent_path() / line))
	    corpus_files.push_back(corpus_list_file.parent_path() / line);
	  else
	    throw std::runtime_error(std::string("no file? ") + line);
	}
      }
    }
    
    if (corpus_files.empty()) 
      throw std::runtime_error("no corpus files");
    
    DependencyCounts::preprocess(output_file);
    
    count_path_set_type count_paths;
    
    accumulate(corpus_files, output_file, count_paths);
    
    DependencyCounts::postprocess(output_file, count_paths);
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void accumulate(const path_set_type& corpus_files,
		const path_type& output_file,
		count_path_set_type& count_paths)
{
  if (map_line) {
    typedef DependencyCounts::TaskLine data_type;
    
    data_type::line_set_type lines;
    data_type::queue_type    queue(threads * 8);
    
    boost::thread_group workers;
    std::vector<data_type, std::allocator<data_type> > data(threads, data_type(queue, output_file, max_malloc));
    
    for (int shard = 0; shard != threads; ++ shard)
      workers.add_thread(new boost::thread(data[shard]));
    
    path_set_type::const_iterator fiter_end = corpus_files.end();
    for (path_set_type::const_iterator fiter = corpus_files.begin(); fiter != fiter_end; ++ fiter) {
      if (*fiter != "-" && ! boost::filesystem::exists(*fiter))
	throw std::runtime_error(std::string("no file? ") + fiter->string());
      
      if (debug)
	std::cerr << "file: " << fiter->string() << std::endl;
      
      utils::compress_istream is(*fiter, 1024 * 1024);
      std::string line;
      
      while (std::getline(is, line)) {
	if (line.empty()) continue;
	
	lines.push_back(line);
	
	if (lines.size() >= 1024 * 8) {
	  equeue.push_swap(lines);
	  lines.clear();
	}
      }
    }
    
    if (! lines.empty())
      queue.push_swap(lines);
    
    for (int shard = 0; shard != threads; ++ shard)
      queue.push(data_type::line_set_type());
    
    workers.join_all();
  } else {
    typedef DependencyCounts::TaskFile data_type;
    
    data_type::line_set_type lines;
    data_type::queue_type    queue(threads);
    
    boost::thread_group workers;
    std::vector<data_type, std::allocator<data_type> > data(threads, data_type(queue, output_file, max_malloc));
    
    for (int shard = 0; shard != threads; ++ shard)
      workers.add_thread(new boost::thread(data[shard]));
    
    path_set_type::const_iterator fiter_end = corpus_files.end();
    for (path_set_type::const_iterator fiter = corpus_files.begin(); fiter != fiter_end; ++ fiter) {
      if (*fiter != "-" && ! boost::filesystem::exists(*fiter))
	throw std::runtime_error(std::string("no file? ") + fiter->string());
      
      if (debug)
	std::cerr << "file: " << fiter->string() << std::endl;
      
      queue.push(*fiter);
    }
    
    for (int shard = 0; shard != threads; ++ shard)
      queue.push(path_type());
    
    workers.join_all();
    
    for (int shard = 0; shard != threads; ++ shard)
      count_paths += data[shard].paths;
  }
}

int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("corpus",       po::value<path_type>(&corpus_file),      "corpus file")
    ("corpus-list",  po::value<path_type>(&corpus_list_file), "corpus list file")
    
    ("output",       po::value<path_type>(&output_file),      "output directory")
    ("map-line",   po::bool_switch(&map_line),     "map by lines, not by files")
    ("threads",    po::value<int>(&threads),       "# of threads")
    ("max-malloc", po::value<double>(&max_malloc), "maximum malloc in GB")
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

