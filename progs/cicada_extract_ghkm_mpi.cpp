//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada_extract_ghkm_impl.hpp"

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <utils/resource.hpp>

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"

typedef cicada::Sentence  sentence_type;
typedef cicada::Alignment alignment_type;

typedef Bitext bitext_type;

typedef Task task_type;
typedef task_type::queue_type queue_type;
typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

typedef boost::filesystem::path path_type;

path_type source_file;
path_type target_file;
path_type alignment_file;

path_type output_file;

std::string non_terminal = "[x]";
int max_nodes = 15;
int max_height = 3;

bool exhaustive = false;
bool inverse = false;
bool swap_source_target = false;

double max_malloc = 8; // 8 GB

int debug = 0;


void options(int argc, char** argv);

enum {
  bitext_tag = 1000,
  file_tag
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

int main(int argc, char** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    options(argc, argv);
    
    if (source_file.empty() || (! boost::filesystem::exists(source_file) && source_file != "-"))
      throw std::runtime_error("no source file? " + source_file.file_string());
    if (target_file.empty() || (! boost::filesystem::exists(target_file) && target_file != "-"))
      throw std::runtime_error("no target file? " + target_file.file_string());
    if (alignment_file.empty() || (! boost::filesystem::exists(alignment_file) && alignment_file != "-"))
      throw std::runtime_error("no alignment file? " + alignment_file.file_string());
    if (output_file.empty())
      throw std::runtime_error("no output directory?");
    

    if (mpi_rank == 0) {
      // create directories for output
      if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
	boost::filesystem::remove_all(output_file);
      
      boost::filesystem::create_directories(output_file);
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
    }

    static const size_t queue_size = 64;
    
    queue_type queue(queue_size);
    task_type task(queue, output_file, non_terminal, max_nodes, max_height, exhaustive, inverse, swap_source_target, max_malloc);
    boost::thread worker(boost::ref(task));

    if (mpi_rank == 0) {
      typedef boost::iostreams::filtering_ostream ostream_type;
      typedef utils::mpi_device_sink              odevice_type;
      
      typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
      typedef boost::shared_ptr<odevice_type> odevice_ptr_type;
      
      typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
      typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;
      
      ostream_ptr_set_type stream(mpi_size);
      odevice_ptr_set_type device(mpi_size);
      
      for (int rank = 1; rank < mpi_size; ++ rank) {
	stream[rank].reset(new ostream_type());
	device[rank].reset(new odevice_type(rank, bitext_tag, 4096, false, true));
	
	stream[rank]->push(boost::iostreams::gzip_compressor());
	stream[rank]->push(*device[rank]);
      }

      utils::resource start_extract;
      
      utils::compress_istream is_src(source_file, 1024 * 1024);
      utils::compress_istream is_trg(target_file, 1024 * 1024);
      utils::compress_istream is_alg(alignment_file, 1024 * 1024);
      
      bitext_type bitext;

      int non_found_iter = 0;
      size_t num_samples = 0;
      while (is_src && is_trg && is_alg) {
	bool found = false;
	
	for (int rank = 1; rank != mpi_size && is_src && is_trg && is_alg; ++ rank) 
	  if (device[rank]->test()) {
	    
	    while (is_src && is_trg && is_alg) {
	      is_src >> bitext.source;
	      is_trg >> bitext.target;
	      is_alg >> bitext.alignment;
	      
	      if (bitext.source.is_valid() && ! bitext.target.empty()) break;
	    }
	    
	    if (! is_src || ! is_trg || ! is_alg) break;
	    
	    *stream[rank] << bitext << '\n';
	    ++ num_samples;
	    if (debug) {
	      if (num_samples % 10000 == 0)
		std::cerr << '.';
	      if (num_samples % 1000000 == 0)
		std::cerr << std::endl;
	    }
	    found = true;
	  }
	
	if (is_src && is_trg && is_alg && queue.size() < queue_size) {
	  while (is_src && is_trg && is_alg) {
	    is_src >> bitext.source;
	    is_trg >> bitext.target;
	    is_alg >> bitext.alignment;
	    
	    if (bitext.source.is_valid() && ! bitext.target.empty()) break;
	  }
	  
	  if (! is_src || ! is_trg || ! is_alg) break;
	  
	  queue.push_swap(bitext);
	  ++ num_samples;
	  if (debug) {
	    if (num_samples % 10000 == 0)
	      std::cerr << '.';
	    if (num_samples % 1000000 == 0)
	      std::cerr << std::endl;
	  }
	  found = true;
	}
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
      
      if (is_src || is_trg || is_alg)
	throw std::runtime_error("# of lines do not match");

      if (debug)
	std::cerr << std::endl;
      if (debug)
	std::cerr << "# of samples: " << num_samples << std::endl;
      
      for (int rank = 1; rank != mpi_size; ++ rank)
	stream[rank].reset();
      
      for (;;) {
	const bool found = utils::mpi_terminate_devices(stream, device);
	
	if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }

      utils::resource end_extract;

      if (debug)
	std::cerr << "extract counts cpu time:  " << end_extract.cpu_time() - start_extract.cpu_time() << std::endl
		  << "extract counts user time: " << end_extract.user_time() - start_extract.user_time() << std::endl;
      
    } else {
      boost::iostreams::filtering_istream stream;
      stream.push(boost::iostreams::gzip_decompressor());
      stream.push(utils::mpi_device_source(0, bitext_tag, 4096));
      
      bitext_type bitext;
      
      while (stream >> bitext)
	queue.push_swap(bitext);
    }
    
    // termination...
    queue.push(bitext_type());
    worker.join();
    
    if (mpi_rank == 0) {
      utils::compress_ostream os(output_file / "files");
      
      task_type::path_set_type::const_iterator piter_end = task.paths.end();
      for (task_type::path_set_type::const_iterator piter = task.paths.begin(); piter != piter_end; ++ piter) {
	utils::tempfile::erase(*piter);
	os << piter->filename() << '\n';
      }
      
      for (int rank = 1; rank != mpi_size; ++ rank) {
	boost::iostreams::filtering_istream is;
	is.push(utils::mpi_device_source(rank, file_tag, 4096));
	
	std::string line;
	while (std::getline(is, line)) {
	  if (! boost::filesystem::exists(line))
	    throw std::runtime_error("no file?");
	  
	  os << path_type(line).filename() << '\n';
	}
      }
      
    } else {
      boost::iostreams::filtering_ostream os;
      os.push(utils::mpi_device_sink(0, file_tag, 4096));
      
      task_type::path_set_type::const_iterator piter_end = task.paths.end();
      for (task_type::path_set_type::const_iterator piter = task.paths.begin(); piter != piter_end; ++ piter) {
	utils::tempfile::erase(*piter);
	os << piter->file_string() << '\n';
      }
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("source",    po::value<path_type>(&source_file),    "source file")
    ("target",    po::value<path_type>(&target_file),    "target file")
    ("alignment", po::value<path_type>(&alignment_file), "alignment file")
    ("output",    po::value<path_type>(&output_file),    "output directory")
    
    ("non-terminal", po::value<std::string>(&non_terminal)->default_value(non_terminal), "non-terminal for target side")
    
    ("max-nodes",  po::value<int>(&max_nodes)->default_value(max_nodes),   "maximum # of nodes in a rule")
    ("max-height", po::value<int>(&max_height)->default_value(max_height), "maximum height of a rule")
    
    ("exhaustive", po::bool_switch(&exhaustive),                           "exhausive extraction")
    ("inverse",    po::bool_switch(&inverse),                              "inversed word alignment")
    ("swap",       po::bool_switch(&swap_source_target),                   "swap source/target")
    
    ("max-malloc", po::value<double>(&max_malloc), "maximum malloc in GB")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    if (mpi_rank == 0)
      std::cout << argv[0] << " [options]\n"
		<< desc_command << std::endl;
    MPI::Finalize();
    exit(0);
  }
}
