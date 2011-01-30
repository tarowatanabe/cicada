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

#include <utils/mpi.hpp>
#include <utils/mpi_stream.hpp>
#include <utils/mpi_stream_simple.hpp>
#include <utils/async_device.hpp>
#include <utils/compress_stream.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/subprocess.hpp>


struct MapReduce
{
  typedef std::pair<int, std::string> value_type;
  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_id_type;

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
      if (value.first < 0) break;
      
      queue_id.push(value.first);
      
      os << value.second << '\n';
    }
    
    queue_id.push(-1);
  }

  queue_type&      queue;
  queue_id_type&   queue_id;
  subprocess_type& subprocess;
};

struct Reducer
{
  typedef MapReduce map_reduce_type;
  
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
      if (value.first < 0) break;
      
      if (! std::getline(is, value.second))
	throw std::runtime_error("invalid lines?");
      
      queue.push_swap(value);
    }
    
    value.first = -1;
    value.second = std::string();
    
    queue.push_swap(value);
  }
  
  queue_type&      queue;
  queue_id_type&   queue_id;
  subprocess_type& subprocess;
};

struct Consumer
{
  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::value_type      value_type;
  typedef map_reduce_type::queue_type      queue_type;

  Consumer(queue_type& __queue,
	   std::istream& __is)
    : queue(__queue),
      is(__is) {}
  
  void operator()()
  {
    value_type value(0, std::string());
    
    while (std::getline(is, value.second)) {
      queue.push(value);
      ++ value.first;
    }
    
    value.first = -1;
    value.second = std::string();
    
    queue.push_swap(value);
  }
  
  queue_type&   queue;
  std::istream& is;
};

struct Dumper
{
  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::value_type      value_type;
  typedef map_reduce_type::queue_type      queue_type;
  
  Dumper(queue_type& __queue,
	 std::ostream& __os,
	 const int __mpi_size)
    : queue(__queue),
      os(__os),
      mpi_size(__mpi_size) {}

  void operator()()
  {
    typedef std::map<int, std::string, std::less<int>, std::allocator<std::pair<const int, std::string> > > value_set_type;

    value_type value;
    value_set_type values;
    
    int curr = 0;

    int consumed_size = 0;
    
    while (1) {
      queue.pop_swap(value);
      if (value.first < 0) {
	++ consumed_size;
	if (consumed_size == mpi_size)
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
      throw std::runtime_error("invalid id: expected: " + boost::lexical_cast<std::string>(curr) + " current: " + boost::lexical_cast<std::string>(values.begin()->first));
    }
    
    value_set_type::const_iterator viter_end = values.end();
    for (value_set_type::const_iterator viter = values.begin(); viter != viter_end; ++ viter)
      os << viter->second << '\n';
  }
  
  queue_type&   queue;
  std::ostream& os;
  int mpi_size;
};

enum {
  line_tag = 1000,
};

inline
void tokenize(const std::string& buffer, MapReduce::value_type& value)
{
#if 0
  boost::iostreams::filtering_istream is_buffer;
  is_buffer.push(boost::iostreams::array_source(buffer.c_str(), buffer.size()));
  
  is_buffer >> value.first;
  std::getline(is_buffer, value.second);
#endif
  
  std::string::const_iterator iter = buffer.begin();

  for (/**/; iter != buffer.end() && ! std::isspace(*iter); ++ iter);
  
  value.first = boost::lexical_cast<int>(buffer.substr(0, iter - buffer.begin()));
  value.second = buffer.substr(iter + 1 - buffer.begin());
}


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

typedef boost::filesystem::path path_type;

path_type input_file = "-";
path_type output_file = "-";
std::string command;

bool even = false;
int debug = 0;

int getoptions(int argc, char** argv);


int main(int argc, char** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;

    if (command.empty())
      throw std::runtime_error("no command?");

    typedef MapReduce map_reduce_type;
      
    typedef map_reduce_type::value_type      value_type;
    typedef map_reduce_type::queue_type      queue_type;
    typedef map_reduce_type::queue_id_type   queue_id_type;
    typedef map_reduce_type::subprocess_type subprocess_type;

    typedef Mapper  mapper_type;
    typedef Reducer reducer_type;
    
    typedef Consumer consumer_type;
    typedef Dumper   dumper_type;
    
    if (mpi_rank == 0) {
      subprocess_type subprocess(command);
      
      queue_type    queue_is(mpi_size);
      queue_type    queue_send(1);
      queue_id_type queue_id;
      queue_type    queue_recv;

      utils::compress_istream is(input_file, 1024 * 1024);
      utils::compress_ostream os(output_file, output_file == "-" ? 4096 : 1024 * 1024);

      boost::thread consumer(consumer_type(queue_is, is));
      boost::thread dumper(dumper_type(queue_recv, os, mpi_size));
      
      boost::thread mapper(mapper_type(queue_send, queue_id, subprocess));
      boost::thread reducer(reducer_type(queue_recv, queue_id, subprocess));
      
      
      typedef utils::mpi_ostream        ostream_type;
      typedef utils::mpi_istream_simple istream_type;

      typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
      typedef boost::shared_ptr<istream_type> istream_ptr_type;
      
      typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
      typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
      
      ostream_ptr_set_type ostream(mpi_size);
      istream_ptr_set_type istream(mpi_size);
      
      for (int rank = 1; rank < mpi_size; ++ rank) {
	ostream[rank].reset(new ostream_type(rank, line_tag, 4096));
	istream[rank].reset(new istream_type(rank, line_tag, 4096));
      }
      
      std::string line;
      value_type value(0, std::string());
      value_type value_recv(0, std::string());
      
      int non_found_iter = 0;
      
      if (even) {
	
	int id = 0;

	while (value.first >= 0) {
	  bool found = false;

	  if ((id % mpi_size == 0 ? queue_send.empty() : ostream[id % mpi_size]->test()) && queue_is.pop(value, true) && value.first >= 0) {
	    const int rank = value.first % mpi_size;
	    
	    if (rank == 0)
	      queue_send.push(value);
	    else
	      ostream[rank]->write(boost::lexical_cast<std::string>(value.first) + ' ' + value.second);
	    
	    ++ id;
	    
	    found = true;
	  }
	  
	  // reduce...
	  for (int rank = 1; rank < mpi_size; ++ rank)
	    if (istream[rank] && istream[rank]->test()) {
	      if (istream[rank]->read(line)) {
		tokenize(line, value_recv);
		
		queue_recv.push_swap(value_recv);
	      } else
		istream[rank].reset();
	      
	      found = true;
	    }
	  
	  non_found_iter = loop_sleep(found, non_found_iter);
	}
	
      } else {
	while (value.first >= 0) {
	  bool found = false;
	  
	  for (int rank = 1; rank < mpi_size && value.first >= 0; ++ rank)
	    if (ostream[rank]->test() && queue_is.pop(value, true) && value.first >= 0) {
	      ostream[rank]->write(boost::lexical_cast<std::string>(value.first) + ' ' + value.second);
	    
	      found = true;
	    }
	
	  if (queue_send.empty() && queue_is.pop(value, true) && value.first >= 0) {
	    queue_send.push(value);
	  
	    found = true;
	  }
	
	  // reduce...
	  for (int rank = 1; rank < mpi_size; ++ rank)
	    if (istream[rank] && istream[rank]->test()) {
	      if (istream[rank]->read(line)) {
		tokenize(line, value_recv);
	      
		queue_recv.push_swap(value_recv);
	      } else {
		queue_recv.push(std::make_pair(-1, std::string()));
		istream[rank].reset();
	      }
	    
	      found = true;
	    }
	  
	  non_found_iter = loop_sleep(found, non_found_iter);
	}
      }
      
      bool terminated = false;

      for (;;) {
	bool found = false;

	if (! terminated && queue_send.push(std::make_pair(-1, std::string()), true)) {
	  terminated = true;
	  found = true;
	}
	
	// termination...
	for (int rank = 1; rank < mpi_size; ++ rank)
	  if (ostream[rank] && ostream[rank]->test()) {
	    if (! ostream[rank]->terminated())
	      ostream[rank]->terminate();
	    else
	      ostream[rank].reset();
	    
	    found = true;
	  }
	
	
	// reduce...
	for (int rank = 1; rank < mpi_size; ++ rank)
	  if (istream[rank] && istream[rank]->test()) {
	    if (istream[rank]->read(line)) {
	      tokenize(line, value_recv);
	      
	      queue_recv.push_swap(value_recv);
	    } else {
	      queue_recv.push(std::make_pair(-1, std::string()));
	      istream[rank].reset();
	    }
		
	    found = true;
	  }
	
	// termination condition!
	if (std::count(istream.begin(), istream.end(), istream_ptr_type()) == mpi_size
	    && std::count(ostream.begin(), ostream.end(), ostream_ptr_type()) == mpi_size
	    && terminated) break;
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
      
      
      mapper.join();
      reducer.join();
      consumer.join();
      dumper.join();
    } else {
      subprocess_type subprocess(command);
      
      queue_type    queue_send(1);
      queue_id_type queue_id;
      queue_type    queue_recv;
      
      boost::thread mapper(mapper_type(queue_send, queue_id, subprocess));
      boost::thread reducer(reducer_type(queue_recv, queue_id, subprocess));

      typedef utils::mpi_istream        istream_type;
      typedef utils::mpi_ostream_simple ostream_type;
      
      boost::shared_ptr<istream_type> is(new istream_type(0, line_tag, 4096));
      boost::shared_ptr<ostream_type> os(new ostream_type(0, line_tag, 4096));
      
      std::string line;
      value_type value;

      bool terminated = false;
      
      int non_found_iter = 0;
      for (;;) {
	bool found = false;
	
	if (is && is->test() && queue_send.empty()) {
	  if (is->read(line))
	    tokenize(line, value);
	  else {
	    value.first = -1;
	    value.second = std::string();
	    
	    is.reset();
	  }
	  
	  queue_send.push_swap(value);
	  
	  found = true;
	}
	
	if (! terminated) {
	  if (os && os->test() && queue_recv.pop_swap(value, true)) {
	    if (value.first < 0)
	      terminated = true;
	    else
	      os->write(boost::lexical_cast<std::string>(value.first) + ' ' + value.second);
	    
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
      
      mapper.join();
      reducer.join();
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

int getoptions(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    ("command", po::value<std::string>(&command),                              "command")
    ("even",   po::bool_switch(&even),                                          "evenly split data")
    ("debug",  po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  //po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    
    if (mpi_rank == 0)
      std::cout << argv[0] << " [options] [file(s) listing shell commands]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}
