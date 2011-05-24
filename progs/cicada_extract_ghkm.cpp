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
bool constrained = false;
bool inverse = false;
bool swap_source_target = false;

double max_malloc = 8; // 8 GB
int threads = 1;

int debug = 0;


void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (source_file.empty() || (! boost::filesystem::exists(source_file) && source_file != "-"))
      throw std::runtime_error("no source file? " + source_file.string());
    if (target_file.empty() || (! boost::filesystem::exists(target_file) && target_file != "-"))
      throw std::runtime_error("no target file? " + target_file.string());
    if (alignment_file.empty() || (! boost::filesystem::exists(alignment_file) && alignment_file != "-"))
      throw std::runtime_error("no alignment file? " + alignment_file.string());
    if (output_file.empty())
      throw std::runtime_error("no output directory?");

    threads = utils::bithack::max(threads, 1);
    
    // create directories for output
    if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
      boost::filesystem::remove_all(output_file);
    
    boost::filesystem::create_directories(output_file);
    
    boost::filesystem::directory_iterator iter_end;
    for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
      boost::filesystem::remove_all(*iter);

    utils::resource start_extract;
    
    queue_type queue(1024 * threads);
    task_set_type tasks(threads, task_type(queue, output_file, non_terminal, max_nodes, max_height, exhaustive, constrained, inverse, swap_source_target, max_malloc));
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    utils::compress_istream is_alg(alignment_file, 1024 * 1024);
    
    bitext_type bitext;
    
    size_t num_samples = 0;
    for (;;) {
      is_src >> bitext.source;
      is_trg >> bitext.target;
      is_alg >> bitext.alignment;
      
      if (! is_src || ! is_trg || ! is_alg) break;
      
      if (bitext.source.is_valid() && ! bitext.target.empty()) {
	queue.push_swap(bitext);
	++ num_samples;
	
	if (debug) {
	  if (num_samples % 10000 == 0)
	    std::cerr << '.';
	  if (num_samples % 1000000 == 0)
	    std::cerr << std::endl;
	}
      }
    }
    if (debug)
      std::cerr << std::endl;
    if (debug)
      std::cerr << "# of samples: " << num_samples << std::endl;
    
    if (is_src || is_trg || is_alg)
      throw std::runtime_error("# of lines do not match");
    
    for (int i = 0; i != threads; ++ i) {
      bitext.clear();
      queue.push_swap(bitext);
    }
    
    workers.join_all();

    utils::resource end_extract;
    
    if (debug)
      std::cerr << "extract counts cpu time:  " << end_extract.cpu_time() - start_extract.cpu_time() << std::endl
		<< "extract counts user time: " << end_extract.user_time() - start_extract.user_time() << std::endl;
    
    utils::compress_ostream os(output_file / "files");
    for (int i = 0; i != threads; ++ i) {
      task_type::path_set_type::const_iterator piter_end = tasks[i].paths.end();
      for (task_type::path_set_type::const_iterator piter = tasks[i].paths.begin(); piter != piter_end; ++ piter) {
	utils::tempfile::erase(*piter);
	os << path_type(piter->filename()).string() << '\n';
      }
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("source",    po::value<path_type>(&source_file),    "source hypergraph file")
    ("target",    po::value<path_type>(&target_file),    "target file")
    ("alignment", po::value<path_type>(&alignment_file), "alignment file")
    ("output",    po::value<path_type>(&output_file),    "output directory")
    
    ("non-terminal", po::value<std::string>(&non_terminal)->default_value(non_terminal), "non-terminal for target side")
    
    ("max-nodes",  po::value<int>(&max_nodes)->default_value(max_nodes),   "maximum # of nodes in a rule")
    ("max-height", po::value<int>(&max_height)->default_value(max_height), "maximum height of a rule")
    
    ("exhaustive",  po::bool_switch(&exhaustive),                           "exhausive extraction")
    ("constrained", po::bool_switch(&constrained),                          "constrained minimum extraction")
    ("inverse",     po::bool_switch(&inverse),                              "inversed word alignment")
    ("swap",        po::bool_switch(&swap_source_target),                   "swap source/target")
    
    ("max-malloc", po::value<double>(&max_malloc), "maximum malloc in GB")
    ("threads",    po::value<int>(&threads),       "# of threads")
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
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
