// -*- encoding: utf-8 -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include <utils/mpi_device.hpp>
#include <utils/mpi_stream.hpp>
#include <utils/mpi_stream_simple.hpp>
#include <utils/compress_stream.hpp>
#include <utils/subprocess.hpp>


enum {
  command_tag = 1000,
  line_tag,
  notify_tag,
  sync_tag,
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

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
typedef std::vector<std::string, std::allocator<std::string> > command_set_type;

path_type input_file = "-";
path_set_type command_files;

path_type prog_name;

bool even = false;
int debug = 0;

int getoptions(int argc, char** argv);


int main(int argc, char** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    std::vector<const char*, std::allocator<const char*> > args;
    args.reserve(argc);
    for (int i = 1; i < argc; ++ i)
      args.push_back(argv[i]);
    args.push_back(0);
    
    if (getoptions(argc, argv) != 0)
      return 1;

    if (! prog_name.empty() && ! boost::filesystem::exists(prog_name))
      throw std::runtime_error(std::string("no binary? ") + prog_name.string());

    if (command_files.empty())
      command_files.push_back("-");
    
    if (MPI::Comm::Get_parent() != MPI::COMM_NULL) {
      utils::mpi_intercomm comm_parent(MPI::Comm::Get_parent());
      
      std::string command;
      {
	utils::mpi_istream_simple is(comm_parent.comm, 0, command_tag, 4096);
	is.read(command);
      }
      
      FILE* fp = ::popen(command.c_str(), "w");
      
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::zlib_decompressor());
      is.push(utils::mpi_device_source(comm_parent.comm, 0, line_tag, 1024 * 1024));
      
      std::string line;
      while (std::getline(is, line)) {
	::fwrite(line.c_str(), 1, line.size(), fp);
	::fwrite("\n", 1, 1, fp);
      }
      ::pclose(fp);
      
      comm_parent.comm.Send(0, 0, MPI::INT, 0, notify_tag);

      MPI::COMM_WORLD.Barrier();
      
    } else {
      command_set_type commands;
      
      if (mpi_rank == 0)
	for (path_set_type::const_iterator piter = command_files.begin(); piter != command_files.end(); ++ piter) {
	  if (debug)
	    std::cerr << "file: " << piter->string() << std::endl;
	  
	  utils::compress_istream is(*piter);
	  
	  std::string command;
	  while (std::getline(is, command)) {
	    boost::algorithm::trim(command);
	    if (! command.empty())
	      commands.push_back(command);
	  }
	}

      int commands_size = commands.size();
      MPI::COMM_WORLD.Bcast(&commands_size, 1, MPI::INT, 0);
      
      if (commands_size == 0)
	throw std::runtime_error("no commands to run?");
      
      const int mpi_child_size = commands_size;
      
      const std::string name = (boost::filesystem::exists(prog_name) ? prog_name.string() : std::string(argv[0]));
      utils::mpi_intercomm comm_child(MPI::COMM_WORLD.Spawn(name.c_str(), &(*args.begin()), mpi_child_size, MPI::INFO_NULL, 0));
      
      if (mpi_rank == 0) {
	for (int rank = 0; rank != mpi_child_size; ++ rank) {
	  utils::mpi_ostream_simple stream(comm_child.comm, rank, command_tag, 4096);
	  stream.write(commands[rank]);
	}

	typedef boost::iostreams::filtering_ostream ostream_type;
	typedef utils::mpi_device_sink              odevice_type;
	
	typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
	typedef boost::shared_ptr<odevice_type> odevice_ptr_type;

	typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
	typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;
	
	ostream_ptr_set_type stream(mpi_child_size);
	odevice_ptr_set_type device(mpi_child_size);
	
	for (int rank = 0; rank != mpi_child_size; ++ rank) {
	  device[rank].reset(new odevice_type(comm_child.comm, rank, line_tag, 1024 * 1024, false, true));
	  
	  stream[rank].reset(new ostream_type());
	  stream[rank]->push(boost::iostreams::zlib_compressor());
	  stream[rank]->push(*device[rank]);
	}

	utils::compress_istream is(input_file, 1024 * 1024);
	
	std::string line;
	int non_found_iter = 0;
	
	if (even) {
	  for (;;) {
	    bool found = false;

	    if (is)
	      for (int rank = 0; rank != mpi_child_size && std::getline(is, line); ++ rank)
		if (stream[rank]) {
		  *stream[rank] << line << '\n';
		  
		  found = true;
		}
	    
	    if (! is) {
	      for (int rank = 0; rank != mpi_child_size; ++ rank)
		if (stream[rank] && device[rank] && device[rank]->test()) {
		  stream[rank].reset();
		  
		  found = true;
		}
	      
	      found |= utils::mpi_terminate_devices(stream, device);
	      
	      if (std::count(device.begin(), device.end(), odevice_ptr_type()) == device.size()) break;
	    }
	  
	    non_found_iter = loop_sleep(found, non_found_iter);
	  }
	  
	} else {
	  for (;;) {
	    bool found = false;
	  
	    for (int rank = 0; rank != mpi_child_size && is; ++ rank)
	      if (stream[rank] && device[rank] && device[rank]->test() && device[rank]->flush(true) == 0 && std::getline(is, line)) {
		*stream[rank] << line << '\n';
	      
		found = true;
	      }
	  
	    if (! is) {
	      for (int rank = 0; rank != mpi_child_size; ++ rank)
		if (stream[rank] && device[rank] && device[rank]->test()) {
		  stream[rank].reset();
		
		  found = true;
		}
	    
	      found |= utils::mpi_terminate_devices(stream, device);
	    
	      if (std::count(device.begin(), device.end(), odevice_ptr_type()) == device.size()) break;
	    }
	  
	    non_found_iter = loop_sleep(found, non_found_iter);
	  }
	}
	
	// termination...
	
	std::vector<MPI::Request, std::allocator<MPI::Request> > request(mpi_child_size);
	std::vector<bool, std::allocator<bool> > terminated(mpi_child_size, false);
	
	for (int rank = 0; rank != mpi_child_size; ++ rank)
	  request[rank] = comm_child.comm.Irecv(0, 0, MPI::INT, rank, notify_tag);
	
	for (;;) {
	  bool found = false;
	  
	  for (int rank = 0; rank != mpi_child_size; ++ rank)
	    if (! terminated[rank] && request[rank].Test()) {
	      terminated[rank] = true;
	      found = true;
	    }
	  
	  if (std::count(terminated.begin(), terminated.end(), true) == mpi_child_size) break;
	  
	  non_found_iter = loop_sleep(found, non_found_iter);
	}
      }
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
    ("input", po::value<path_type>(&input_file)->default_value(input_file), "input file")
    ("prog",  po::value<path_type>(&prog_name),  "this binary")
    ("even",  po::bool_switch(&even),            "evenly split data")
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::options_description hidden;
  hidden.add_options()
    ("command-file", po::value<path_set_type>(&command_files), "command files");

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);
  
  po::positional_options_description pos;
  pos.add("command-file", -1); // all the files
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    
    if (mpi_rank == 0)
      std::cout << argv[0] << " [options] [file(s) listing shell commands]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}
