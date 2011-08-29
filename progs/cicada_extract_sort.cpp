//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada_extract_score_impl.hpp"

#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <utils/resource.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>

typedef PhrasePair       rule_pair_type;
typedef PhrasePairParser rule_pair_parser_type;

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type input_files;
path_type output_file;

double max_malloc = 8;
int threads = 1;

int debug = 0;

struct Task
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  Task(queue_type& __queue,
       const path_type& __output,
       const size_t __malloc_threshold)
    : queue(__queue), output(__output), malloc_threshold(__malloc_threshold) {}

  queue_type&   queue;
  path_type     output;
  path_set_type paths;
  const size_t  malloc_threshold;
  
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
				  std::allocator<rule_pair_type> > rule_pair_set_type;
#else
  typedef sgi::hash_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
			std::allocator<rule_pair_type> > rule_pair_set_type;
#endif

  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };

  void dump(const rule_pair_set_type& rule_pairs)
  {
    typedef std::vector<const rule_pair_type*, std::allocator<const rule_pair_type*> > sorted_type;
    
    if (rule_pairs.empty()) return;
    
    sorted_type sorted(rule_pairs.size());
    {
      sorted_type::iterator siter = sorted.begin();
      rule_pair_set_type::const_iterator citer_end = rule_pairs.end();
      for (rule_pair_set_type::const_iterator citer = rule_pairs.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    
    std::sort(sorted.begin(), sorted.end(), less_ptr<rule_pair_type>());
    
    const path_type path_tmp = utils::tempfile::file_name(output / "counts-XXXXXX");
    utils::tempfile::insert(path_tmp);
    const path_type path = path_tmp.string() + ".gz";
    utils::tempfile::insert(path);
    
    utils::compress_ostream os(path, 1024 * 1024);
    os.precision(20);
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << *(*siter) << '\n';
    
    paths.push_back(path);
  }
  
  void operator()()
  {
    rule_pair_parser_type parser;

    rule_pair_type        rule_pair;
    rule_pair_set_type    rule_pairs;
    
    std::string line;
    
    const size_t iter_mask = (1 << 10) - 1;
    
    for (size_t iter = 0; /**/; ++ iter) {
      queue.pop_swap(line);
      if (line.empty()) break;

      if (! parser(line, rule_pair)) continue;
      
      std::pair<rule_pair_set_type::iterator, bool> result = rule_pairs.insert(rule_pair);
      if (! result.second)
	const_cast<rule_pair_type&>(*result.first).increment(rule_pair.counts.begin(), rule_pair.counts.end());
      
      if ((iter & iter_mask) == iter_mask && utils::malloc_stats::used() > malloc_threshold) {
	dump(rule_pairs);
	rule_pairs.clear();
      }
    }
    
    dump(rule_pairs);
  }
};

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (input_files.empty())
      input_files.push_back("-");
    
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

    typedef Task task_type;
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
    
    task_type::queue_type queue(16 * 1024 * threads);
    task_set_type tasks(threads, task_type(queue, output_file, max_malloc * 1024 * 1024 * 1024));
    
    utils::resource start;
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    std::string line;
    
    path_set_type::const_iterator fiter_end = input_files.end();
    for (path_set_type::const_iterator fiter = input_files.begin(); fiter != fiter_end; ++ fiter) {
      
      if (boost::filesystem::is_directory(*fiter)) {
	if (! boost::filesystem::exists(*fiter / "files"))
	  throw std::runtime_error("no files? " +  (*fiter / "files").string());
	
	std::string list;
	utils::compress_istream is_list(*fiter / "files");
	while (std::getline(is_list, list)) 
	  if (! list.empty()) {
	    const path_type path = *fiter / list;
	    
	    if (! boost::filesystem::exists(path))
	      throw std::runtime_error("no count files? " + path.string());
	    
	    utils::compress_istream is(path, 1024 * 1024);
	    
	    while (std::getline(is, line)) 
	      if (! line.empty())
		queue.push_swap(line);
	  }
      } else {
	utils::compress_istream is(*fiter, 1024 * 1024);
	
	while (std::getline(is, line)) 
	  if (! line.empty())
	    queue.push_swap(line);
      }
    }

    for (int i = 0; i != threads; ++ i)
      queue.push(std::string());
    
    workers.join_all();
    
    utils::resource end;
    
    if (debug)
      std::cerr << "sort counts cpu time:  " << end.cpu_time() - start.cpu_time() << std::endl
		<< "sort counts user time: " << end.user_time() - start.user_time() << std::endl;
    
    utils::compress_ostream os(output_file / "files");
    for (int i = 0; i != threads; ++ i) {
      path_set_type::const_iterator piter_end = tasks[i].paths.end();
      for (path_set_type::const_iterator piter = tasks[i].paths.begin(); piter != piter_end; ++ piter) {
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
  
  po::options_description opts_config("options");
  
  opts_config.add_options()
    ("input",      po::value<path_set_type>(&input_files)->multitoken(), "input file(s)")
    ("output",     po::value<path_type>(&output_file),                   "output file")
    ("max-malloc", po::value<double>(&max_malloc),                       "maximum malloc in GB")
    ("threads",    po::value<int>(&threads),                             "# of threads")

    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc;
  desc.add(opts_config);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}


