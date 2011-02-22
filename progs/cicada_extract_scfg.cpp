//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada_extract_scfg_impl.hpp"

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
path_type spans_source_file;
path_type spans_target_file;

path_type output_file;

int max_length = 7;
int max_fertility = 10;
int max_span = 15;
int min_hole = 1;
bool ternary = false;
bool sentential = false;
bool inverse = false;

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
    task_set_type tasks(threads, task_type(queue, output_file, max_length, max_fertility, max_span, min_hole, ternary, sentential, inverse, max_malloc));
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    utils::compress_istream is_alg(alignment_file, 1024 * 1024);
    
    std::auto_ptr<std::istream> is_span_src(boost::filesystem::exists(spans_source_file)
					    ? new utils::compress_istream(spans_source_file, 1024 * 1024) : 0);
    std::auto_ptr<std::istream> is_span_trg(boost::filesystem::exists(spans_target_file)
					    ? new utils::compress_istream(spans_target_file, 1024 * 1024) : 0);
    
    bitext_type bitext;
    
    size_t num_samples = 0;
    for (;;) {
      is_src >> bitext.source;
      is_trg >> bitext.target;
      is_alg >> bitext.alignment;

      if (is_span_src.get())
	*is_span_src >> bitext.spans_source;

      if (is_span_trg.get())
	*is_span_trg >> bitext.spans_target;
      
      if (! is_src || ! is_trg || ! is_alg || (is_span_src.get() && ! *is_span_src) || (is_span_trg.get() && ! *is_span_trg)) break;
      
      if (! bitext.source.empty() && ! bitext.target.empty()) {
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
    
    if (is_src || is_trg || is_alg || (is_span_src.get() && *is_span_src) || (is_span_trg.get() && *is_span_trg))
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
	os << piter->filename() << '\n';
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
    ("source",       po::value<path_type>(&source_file),       "source file")
    ("target",       po::value<path_type>(&target_file),       "target file")
    ("alignment",    po::value<path_type>(&alignment_file),    "alignment file")
    ("spans-source", po::value<path_type>(&spans_source_file), "source span file")
    ("spans-target", po::value<path_type>(&spans_target_file), "target span file")
    
    ("output",    po::value<path_type>(&output_file),    "output directory")
    
    ("max-length",    po::value<int>(&max_length)->default_value(max_length),       "maximum terminal length")
    ("max-fertility", po::value<int>(&max_fertility)->default_value(max_fertility), "maximum terminal fertility ratio")
    ("max-span",      po::value<int>(&max_span)->default_value(max_span),           "maximum span for rule")
    ("min-hole",      po::value<int>(&min_hole)->default_value(min_hole),           "minimum hole for antecedent non-terminals")
    ("ternary",       po::bool_switch(&ternary),                                    "extract ternary rules")
    ("sentential",    po::bool_switch(&sentential),                                 "extract sentential rules")
    ("inverse",       po::bool_switch(&inverse),                                    "inversed word alignment")
    
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
