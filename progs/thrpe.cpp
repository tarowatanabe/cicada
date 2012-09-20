// -*- encoding: utf-8 -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include <utils/subprocess.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/async_device.hpp>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type> path_set_type;

path_type input_file = "-";
path_type output_file = "-";
std::string command;
int threads = 1;
int debug = 0;

struct MapReduce
{
  typedef size_t id_type;
  typedef std::pair<id_type, std::string> value_type;
  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;
  typedef utils::lockfree_list_queue<id_type, std::allocator<id_type> > queue_id_type;

  typedef utils::subprocess subprocess_type;
};

namespace std
{
  inline
  void swap(MapReduce::value_type& x, MapReduce::value_type& y)
  {
    std::swap(x.first,  y.first);
    std::swap(x.second, y.second);
  }
};

struct Mapper
{
  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::id_type         id_type;
  typedef map_reduce_type::value_type      value_type;
  typedef map_reduce_type::queue_type      queue_type;
  typedef map_reduce_type::queue_id_type   queue_id_type;
  typedef map_reduce_type::subprocess_type subprocess_type;

  Mapper(queue_type& __queue,
	 queue_id_type& __queue_id,
	 subprocess_type& __subprocess)
    : queue(__queue),
      queue_id(__queue_id),
      subprocess(__subprocess)
  {}

  void operator()()
  {
    boost::iostreams::filtering_ostream os;
    os.push(utils::async_sink(subprocess.desc_write(), true));
    subprocess.desc_write() = -1;
    
    value_type value;
    
    while (1) {
      queue.pop_swap(value);
      if (value.first == id_type(-1)) break;
      
      queue_id.push(value.first);
      
      os << value.second << '\n';
    }
    
    queue_id.push(id_type(-1));
  }

  queue_type&      queue;
  queue_id_type&   queue_id;
  subprocess_type& subprocess;
};

struct Reducer
{
  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::id_type         id_type;
  typedef map_reduce_type::value_type      value_type;
  typedef map_reduce_type::queue_type      queue_type;
  typedef map_reduce_type::queue_id_type   queue_id_type;
  typedef map_reduce_type::subprocess_type subprocess_type;
  
  Reducer(queue_type& __queue,
	  queue_id_type& __queue_id,
	  subprocess_type& __subprocess)
    : queue(__queue),
      queue_id(__queue_id),
      subprocess(__subprocess)
  {}

  void operator()()
  {
    boost::iostreams::filtering_istream is;
    is.push(utils::async_source(subprocess.desc_read(), true));
    subprocess.desc_read() = -1;
    
    value_type value;
    
    while (1) {
      queue_id.pop(value.first);
      if (value.first == id_type(-1)) break;
      
      if (! std::getline(is, value.second))
	throw std::runtime_error("invalid lines?");
      
      queue.push_swap(value);
    }
    
    value.first  = id_type(-1);
    value.second = std::string();
    
    queue.push_swap(value);
  }
  
  queue_type&      queue;
  queue_id_type&   queue_id;
  subprocess_type& subprocess;
};

struct Merger
{
  typedef MapReduce map_reduce_type;
 
  typedef map_reduce_type::id_type         id_type;
  typedef map_reduce_type::value_type      value_type;
  typedef map_reduce_type::queue_type      queue_type;
  
  Merger(queue_type& __queue,
	 std::ostream& __os,
	 const int __workers)
    : queue(__queue),
      os(__os),
      workers(__workers) {}

  void operator()()
  {
    typedef std::map<id_type, std::string, std::less<id_type>, std::allocator<std::pair<const id_type, std::string> > > value_set_type;

    value_type value;
    value_set_type values;
    
    id_type curr = 0;

    int consumed_size = 0;
    
    while (1) {
      queue.pop_swap(value);
      if (value.first == id_type(-1)) {
	++ consumed_size;
	if (consumed_size == workers)
	  break;
	else
	  continue;
      }
      
      if (curr == value.first) {
	os << value.second << '\n';
	++ curr;
      } else 
	values.insert(value);
      
      while (! values.empty() && values.begin()->first == curr) {
	os << values.begin()->second << '\n';
	values.erase(values.begin());
	++ curr;
      }
    }
    
    if (! values.empty() && curr != values.begin()->first) {
      throw std::runtime_error("invalid id: expected: " + utils::lexical_cast<std::string>(curr) + " current: " + utils::lexical_cast<std::string>(values.begin()->first));
    }
    
    value_set_type::const_iterator viter_end = values.end();
    for (value_set_type::const_iterator viter = values.begin(); viter != viter_end; ++ viter)
      os << viter->second << '\n';
  }
  
  queue_type&   queue;
  std::ostream& os;
  int workers;
};

int getoptions(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;
    
    threads = utils::bithack::max(threads, 1);

    if (command.empty())
      throw std::runtime_error("no command?");

    typedef MapReduce map_reduce_type;
    typedef Mapper    mapper_type;
    typedef Reducer   reducer_type;
    typedef Merger    merger_type;
    

    // spawn processes first...
    std::vector<boost::shared_ptr<map_reduce_type::subprocess_type>, std::allocator<boost::shared_ptr<map_reduce_type::subprocess_type> > > subprocess(threads);
    for (size_t shard = 0; shard != subprocess.size(); ++ shard)
      subprocess[shard].reset(new map_reduce_type::subprocess_type(command));
    
    map_reduce_type::queue_type queue_mapper(threads);
    map_reduce_type::queue_type queue_reducer(threads);
    
    boost::thread_group mappers;
    boost::thread_group reducers;
    std::vector<map_reduce_type::queue_id_type, std::allocator<map_reduce_type::queue_id_type> > queues(threads);

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    std::auto_ptr<boost::thread> merger(new boost::thread(merger_type(queue_reducer, os, threads)));
    
    for (size_t shard = 0; shard != queues.size(); ++ shard) {
      reducers.add_thread(new boost::thread(reducer_type(queue_reducer, queues[shard], *subprocess[shard])));
      
      mappers.add_thread(new boost::thread(mapper_type(queue_mapper, queues[shard], *subprocess[shard])));
    }
    
    map_reduce_type::id_type id = 0;
    std::string line;
    while (std::getline(is, line)) {
      queue_mapper.push(std::make_pair(id, line));
      ++ id;
    }
    
    for (size_t shard = 0; shard != queues.size(); ++ shard)
      queue_mapper.push(std::make_pair(map_reduce_type::id_type(-1), std::string()));
    
    merger->join();
    reducers.join_all();
    mappers.join_all();
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << argv[0] << " "<< err.what() << std::endl;
    return 1;
  }
  return 0;
}

int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",   po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",  po::value<path_type>(&output_file)->default_value(output_file), "output file")
    ("command", po::value<std::string>(&command),                               "command")
    
    ("threads", po::value<int>(&threads),                                       "# of threads")
    
    ("debug",   po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  //po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}

