// -*- encoding: utf-8 -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cstdio>

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/algorithm/string.hpp>

#include <utils/compress_stream.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type> path_set_type;

path_set_type input_files;
int threads = 0;
int debug = 0;

int getoptions(int argc, char** argv);

void run_command(const std::string& command)
{
  static const size_t buffer_size = 1024;
  
  char buffer[buffer_size];
  
  FILE* fp = ::popen(command.c_str(), "r");
  
  for (;;) {
    const size_t size = ::fread(buffer, 1, buffer_size, fp);
    if (size == 0) {
      if (::feof(fp))
	break;
      else if (::ferror(fp))
	throw std::runtime_error("error?");
    } else
      std::cout.write(buffer, size);
  }
  ::pclose(fp);
}

struct Task
{
  typedef boost::thread thread_type;
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;
  
  queue_type& queue;
  int id;
  
  Task(queue_type& __queue, const int __id)
    : queue(__queue), id(__id) {}

  void operator()()
  {
    std::string command;
    
    while (1) {
      queue.pop_swap(command);
      if (command.empty()) break;

      if (debug)
	std::cerr << "id: " << id << " " << command << std::endl;
      
      run_command(command);
    }
  }
};


int main(int argc, char** argv)
{
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;

    typedef Task task_type;
    typedef task_type::queue_type  queue_type;

    threads = utils::bithack::max(threads, 1);
    
    queue_type queue(threads);
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(task_type(queue, i)));
    
    if (input_files.empty())
      input_files.push_back("-");
    
    for (path_set_type::const_iterator piter = input_files.begin(); piter != input_files.end(); ++ piter) {
      if (debug)
	std::cerr << "file: " << piter->string() << std::endl;
      
      utils::compress_istream is(*piter);
      
      std::string command;
      while (std::getline(is, command)) {
	boost::algorithm::trim(command);
	
	if (command.empty()) continue;
	
	queue.push_swap(command);
      }
    }
    
    for (int i = 0; i < threads; ++ i)
      queue.push(std::string());
    
    workers.join_all();
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  catch (...) {
    std::cerr << "error... " << std::endl;
    return -1;
  }
  return 0;
}

int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("threads", po::value<int>(&threads), "# of threads")
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::options_description hidden;
  hidden.add_options()
    ("input-file", po::value<path_set_type>(&input_files), "input file");

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);
  
  po::positional_options_description pos;
  pos.add("input-file", -1); // all the files

  po::command_line_parser parser(argc, argv);
  parser.style(po::command_line_style::unix_style & (~po::command_line_style::allow_guessing));
  parser.options(cmdline_options);
  parser.positional(pos);
  
  po::variables_map vm;
  po::store(parser.run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options] [file(s) listing shell commands]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}
