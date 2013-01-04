//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada_extract_score_impl.hpp"
#include "cicada_output_impl.hpp"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>

#include <utils/filesystem.hpp>
#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/lexical_cast.hpp>

typedef boost::filesystem::path                                    path_type;
typedef std::vector<path_type, std::allocator<path_type> >         path_set_type;
typedef std::vector<path_set_type, std::allocator<path_set_type> > path_map_type;

typedef RootCount root_count_type;
typedef utils::unordered_set<root_count_type, boost::hash<root_count_type>, std::equal_to<root_count_type>,
			     std::allocator<root_count_type> >::type root_count_set_type;
typedef std::vector<root_count_set_type, std::allocator<root_count_set_type> > root_count_map_type;

struct less_file_size
{
  bool operator()(const path_type& x, const path_type& y) const
  {
    return boost::filesystem::file_size(x) < boost::filesystem::file_size(y);
  }
};

struct greater_file_size
{
  bool operator()(const path_type& x, const path_type& y) const
  {
    return boost::filesystem::file_size(x) > boost::filesystem::file_size(y);
  }
};

path_set_type input_files;
path_type output_file = "";

bool score_phrase = false;
bool score_scfg   = false;
bool score_ghkm   = false;

double max_malloc = 8; // 8 GB
int    threads = 1;

int debug = 0;

void merge_counts(path_set_type& counts_files);

template <typename Extractor>
void source_counts(const path_set_type& counts_files,
		   path_set_type& source_files,
		   root_count_set_type& root_joint,
		   root_count_set_type& root_sources);

void reverse_counts(const path_set_type& counts_files,
		   path_map_type& reversed_files);
template <typename Extractor>
void target_counts(const path_map_type& reversed_files,
		   path_map_type& target_files,
		   root_count_set_type& root_targets);
void score_counts(const path_type& output_file,
		  const path_set_type& counts_files,
		  const path_set_type& source_files,
		  const path_map_type& target_files);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (output_file.empty())
      throw std::runtime_error("no output file?");

    if (int(score_phrase) + score_scfg + score_ghkm != 1)
      throw std::runtime_error("specify either one of --score-phrase|scfg|ghkm");
    
    threads = utils::bithack::max(1, threads);

    path_set_type counts_files;
    
    for (path_set_type::const_iterator fiter = input_files.begin(); fiter != input_files.end(); ++ fiter) {
      if (! boost::filesystem::exists(*fiter))
	throw std::runtime_error("no file? " + fiter->string());
      
      if (boost::filesystem::is_directory(*fiter)) {
	const path_type path = *fiter / "files";
	
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error("no files? " + path.string());
	
	utils::compress_istream is(path, 1024 * 1024);
	std::string line;
	while (std::getline(is, line))
	  if (! line.empty()) {
	    const path_type path(*fiter / line);
	    
	    if (! boost::filesystem::exists(path))
	      throw std::runtime_error("no file? " + path.string());
	    
	    counts_files.push_back(path);
	  }
      } else
	counts_files.push_back(*fiter);
    }
    
    std::sort(counts_files.begin(), counts_files.end(), greater_file_size());
    
    if (counts_files.size() > 128) {
      if (debug)
	std::cerr << "merge counts: " << counts_files.size() << std::endl;
      
      utils::resource start_merge;
      
      merge_counts(counts_files);
      
      utils::resource end_merge;
      
      if (debug)
	std::cerr << "merge counts: " << counts_files.size()
		  << " cpu time:  " << end_merge.cpu_time() - start_merge.cpu_time()
		  << " user time: " << end_merge.user_time() - start_merge.user_time()
		  << std::endl;
    }

    
    // reverse counts...
    path_set_type source_files;
    path_map_type reversed_files;
    path_map_type target_files;
    root_count_set_type root_joint;
    root_count_set_type root_sources;
    root_count_set_type root_targets;
    
    utils::resource start_source;
    if (score_phrase)
      source_counts<ExtractRootPhrase>(counts_files, source_files, root_joint, root_sources);
    else if (score_scfg)
      source_counts<ExtractRootSCFG>(counts_files, source_files, root_joint, root_sources);
    else
      source_counts<ExtractRootGHKM>(counts_files, source_files, root_joint, root_sources);
    utils::resource end_source;
    
    if (debug)
      std::cerr << "source counts"
		<< " cpu time:  " << end_source.cpu_time() - start_source.cpu_time()
		<< " user time: " << end_source.user_time() - start_source.user_time()
		<< std::endl;
    
    utils::resource start_reverse;
    reverse_counts(counts_files, reversed_files);
    utils::resource end_reverse;
    
    if (debug)
      std::cerr << "reverse counts"
		<< " cpu time:  " << end_reverse.cpu_time() - start_reverse.cpu_time()
		<< " user time: " << end_reverse.user_time() - start_reverse.user_time()
		<< std::endl;
    
    utils::resource start_target;
    if (score_phrase)
      target_counts<ExtractRootPhrase>(reversed_files, target_files, root_targets);
    else if (score_scfg)
      target_counts<ExtractRootSCFG>(reversed_files, target_files, root_targets);
    else
      target_counts<ExtractRootGHKM>(reversed_files, target_files, root_targets);
    utils::resource end_target;
    
    if (debug)
      std::cerr << "target counts"
		<< " cpu time:  " << end_target.cpu_time() - start_target.cpu_time()
		<< " user time: " << end_target.user_time() - start_target.user_time()
		<< std::endl;
    
    // scoring...
    utils::resource start_score;
    score_counts(output_file, counts_files, source_files, target_files);
    utils::resource end_score;
    if (debug)
      std::cerr << "score counts"
		<< " cpu time:  " << end_score.cpu_time() - start_score.cpu_time()
		<< " user time: " << end_score.user_time() - start_score.user_time()
		<< std::endl;
    
    // finally, dump files, root-sources and root-targets...
    {
      utils::compress_ostream os_file(output_file / "files");
      utils::compress_ostream os_joint(output_file / "root-joint.gz");
      utils::compress_ostream os_src(output_file / "root-source.gz");
      utils::compress_ostream os_trg(output_file / "root-target.gz");
      
      os_joint.precision(20);
      os_src.precision(20);
      os_trg.precision(20);
      
      for (int shard = 0; shard != threads; ++ shard)
	os_file << (utils::lexical_cast<std::string>(shard) + ".gz") << '\n';

      root_count_set_type::const_iterator jiter_end = root_joint.end();
      for (root_count_set_type::const_iterator jiter = root_joint.begin(); jiter != jiter_end; ++ jiter)
	os_joint << *jiter << '\n';
      
      root_count_set_type::const_iterator siter_end = root_sources.end();
      for (root_count_set_type::const_iterator siter = root_sources.begin(); siter != siter_end; ++ siter)
	os_src << *siter << '\n';
      
      root_count_set_type::const_iterator titer_end = root_targets.end();
      for (root_count_set_type::const_iterator titer = root_targets.begin(); titer != titer_end; ++ titer)
	os_trg << *titer << '\n';
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct TaskMerge
{
  TaskMerge(path_set_type& __files,
	    const path_type& __prefix,
	    const size_t& __size) : files(__files), prefix(__prefix), size(__size) {}
  
  void operator()()
  {
    typedef PhrasePair       rule_pair_type;
    typedef PhrasePairParser rule_pair_parser_type;

    typedef utils::unordered_set<path_type, boost::hash<path_type>, std::equal_to<path_type>,
				 std::allocator<path_type> >::type path_temporary_type;

    typedef std::pair<size_t, path_type> size_path_type;
    typedef std::vector<size_path_type, std::allocator<size_path_type> > size_path_set_type;
    
    size_path_set_type size_files;
    
    path_set_type::const_iterator fiter_end = files.end();
    for (path_set_type::const_iterator fiter = files.begin(); fiter != fiter_end; ++ fiter)
      size_files.push_back(size_path_type(boost::filesystem::file_size(*fiter), *fiter));
    
    rule_pair_parser_type parser;
    
    path_temporary_type temp;
    
    while (size_files.size() > size && size_files.size() >= 2) {
      std::sort(size_files.begin(), size_files.end(), std::greater<size_path_type>());
      
      const path_type file1 = size_files.back().second;
      size_files.pop_back();
      
      const path_type file2 = size_files.back().second;
      size_files.pop_back();
      
      const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.merged.XXXXXX");
      utils::tempfile::insert(counts_file_tmp);
      const path_type counts_file = counts_file_tmp.string() + ".gz";
      utils::tempfile::insert(counts_file);
      
      temp.insert(counts_file);

      {
	utils::compress_istream is1(file1, 1024 * 1024);
	utils::compress_istream is2(file2, 1024 * 1024);
      
	utils::compress_ostream os(counts_file, 1024 * 1024);
	os.exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
      
	std::string line1;
	std::string line2;

	rule_pair_type rule1;
	rule_pair_type rule2;
      
	bool parsed1 = std::getline(is1, line1) && parser(line1, rule1);
	bool parsed2 = std::getline(is2, line2) && parser(line2, rule2);
      
	while (parsed1 && parsed2) {
	  if (rule1 < rule2) {
	    os << rule1 << '\n';
	    parsed1 = std::getline(is1, line1) && parser(line1, rule1);
	  } else if (rule2 < rule1) {
	    os << rule2 << '\n';
	    parsed2 = std::getline(is2, line2) && parser(line2, rule2);
	  } else {
	    rule1.increment(rule2.counts.begin(), rule2.counts.end());
	    os << rule1 << '\n';
	    parsed1 = std::getline(is1, line1) && parser(line1, rule1);
	    parsed2 = std::getline(is2, line2) && parser(line2, rule2);
	  }
	}
      
	// dump remaining...
	while (parsed1) {
	  os << rule1 << '\n';
	  parsed1 = std::getline(is1, line1) && parser(line1, rule1);
	}
      
	while (parsed2) {
	  os << rule2 << '\n';
	  parsed2 = std::getline(is2, line2) && parser(line2, rule2);
	}
      }
      
      if (temp.find(file1) != temp.end()) {
	boost::filesystem::remove(file1);
	utils::tempfile::erase(file1);
      }
      
      if (temp.find(file2) != temp.end()) {
	boost::filesystem::remove(file2);
	utils::tempfile::erase(file2);
      }
      
      size_files.push_back(size_path_type(boost::filesystem::file_size(counts_file), counts_file));
    }
    
    files.clear();
    
    size_path_set_type::const_iterator siter_end = size_files.end();
    for (size_path_set_type::const_iterator siter = size_files.begin(); siter != siter_end; ++ siter)
      files.push_back(siter->second);
  }
  
  path_set_type& files;
  path_type      prefix;
  size_t         size;
};

void merge_counts(path_set_type& counts_files)
{
  typedef TaskMerge task_type;

  typedef std::vector<path_set_type, std::allocator<path_set_type> > path_map_type;
  
  std::sort(counts_files.begin(), counts_files.end(), greater_file_size());
  
  path_map_type mapped_files(threads);
  for (size_t i = 0; i != counts_files.size(); ++ i)
    mapped_files[i % threads].push_back(counts_files[i]);
  
  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(task_type(mapped_files[i], utils::tempfile::tmp_dir(), 128 / threads)));
  
  workers.join_all();
  
  counts_files.clear();
  for (int i = 0; i != threads; ++ i)
    counts_files.insert(counts_files.end(), mapped_files[i].begin(), mapped_files[i].end());

  std::sort(counts_files.begin(), counts_files.end(), greater_file_size());
}


void score_counts(const path_type& output_file,
		  const path_set_type& counts_files,
		  const path_set_type& source_files,
		  const path_map_type& target_files)
{
  typedef PhrasePairScore map_reduce_type;
  
  typedef PhrasePairScoreMapper  mapper_type;
  typedef PhrasePairScoreReducer reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef std::vector<queue_ptr_set_type, std::allocator<queue_ptr_set_type> > queue_ptr_map_type;

  typedef utils::compress_ostream ostream_type;
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;

  if (static_cast<int>(source_files.size()) != threads)
    throw std::runtime_error("# of threads differ");  
  if (static_cast<int>(target_files.size()) != threads)
    throw std::runtime_error("# of threads differ");
  
  prepare_directory(output_file);
  
  path_map_type mapped_files(threads);
  for (size_t i = 0; i != counts_files.size(); ++ i)
    mapped_files[i % threads].push_back(counts_files[i]);
  
  queue_ptr_map_type   queues_mapper(threads, queue_ptr_set_type(threads));
  queue_ptr_map_type   queues_reducer(threads, queue_ptr_set_type(threads));
  ostream_ptr_set_type ostreams(threads);
  
  // construct queue matrix...
  for (int i = 0; i != threads; ++ i)
    for (int j = 0; j != threads; ++ j) {
      queues_mapper[i][j].reset(new queue_type(1024));
      queues_reducer[j][i] = queues_mapper[i][j];
    }
  
  boost::thread_group reducers;
  for (int shard = 0; shard != threads; ++ shard) {
    const path_type path = output_file / (utils::lexical_cast<std::string>(shard) + ".gz");
    
    ostreams[shard].reset(new utils::compress_ostream(path, 1024 * 1024));
    
    reducers.add_thread(new boost::thread(reducer_type(source_files[shard],
						       target_files[shard],
						       queues_reducer[shard],
						       *ostreams[shard],
						       debug)));
  }
  
  boost::thread_group mappers;
  for (int shard = 0; shard != threads; ++ shard)
    mappers.add_thread(new boost::thread(mapper_type(mapped_files[shard],
						     queues_mapper[shard],
						     max_malloc,
						     debug)));
  
  mappers.join_all();
  reducers.join_all();
}


template <typename Extractor>
void target_counts(const path_map_type& reversed_files,
		   path_map_type& target_files,
		   root_count_set_type& root_targets)
{
  typedef PhrasePairTarget map_reduce_type;
  
  typedef PhrasePairTargetMapper<Extractor>  mapper_type;
  typedef PhrasePairTargetReducer            reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  if (static_cast<int>(reversed_files.size()) != threads)
    throw std::runtime_error("# of threads differ");

  target_files.clear();
  target_files.reserve(threads);
  target_files.resize(threads);
  
  queue_ptr_set_type  queues(threads);
  root_count_map_type root_counts(threads);
  
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    queues[shard].reset(new queue_type(1024 * threads));
  
  boost::thread_group reducers;
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    reducers.add_thread(new boost::thread(reducer_type(*queues[shard],
						       utils::tempfile::tmp_dir(),
						       target_files[shard],
						       threads,
						       max_malloc,
						       debug)));
  
  boost::thread_group mappers;
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    mappers.add_thread(new boost::thread(mapper_type(reversed_files[shard],
						     queues,
						     root_counts[shard],
						     max_malloc,
						     debug)));
  
  reducers.join_all();
  mappers.join_all();
  
  // merge root counts...
  for (size_t shard = 0; shard != root_counts.size(); ++ shard) {
    root_count_set_type::const_iterator citer_end = root_counts[shard].end();
    for (root_count_set_type::const_iterator citer = root_counts[shard].begin(); citer != citer_end; ++ citer) {
      std::pair<root_count_set_type::iterator, bool> result = root_targets.insert(*citer);
      if (! result.second) {
	const_cast<root_count_type&>(*result.first).increment(citer->counts.begin(), citer->counts.end());
	const_cast<root_count_type&>(*result.first).observed += citer->observed;
      }
    }
  }
}

template <typename Extractor>
void source_counts(const path_set_type& counts_files,
		   path_set_type& source_files,
		   root_count_set_type& root_joint,
		   root_count_set_type& root_sources)
{
  typedef PhrasePairSource map_reduce_type;
  
  typedef PhrasePairSourceMapper             mapper_type;
  typedef PhrasePairSourceReducer<Extractor> reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef std::vector<queue_ptr_set_type, std::allocator<queue_ptr_set_type> > queue_ptr_map_type;
  
  path_map_type mapped_files(threads);
  for (size_t i = 0; i != counts_files.size(); ++ i) {
    if (! boost::filesystem::exists(counts_files[i]))
      throw std::runtime_error("no file? " + counts_files[i].string());
    
    mapped_files[i % threads].push_back(counts_files[i]);
  }

  source_files.clear();
  source_files.reserve(threads);
  source_files.resize(threads);
  
  root_count_map_type joint_counts(threads);
  root_count_map_type source_counts(threads);
  
  queue_ptr_map_type   queues_mapper(threads, queue_ptr_set_type(threads));
  queue_ptr_map_type   queues_reducer(threads, queue_ptr_set_type(threads));

  // construct queue matrix...
  for (int i = 0; i != threads; ++ i)
    for (int j = 0; j != threads; ++ j) {
      queues_mapper[i][j].reset(new queue_type(1024));
      queues_reducer[j][i] = queues_mapper[i][j];
    }

  boost::thread_group reducers;
  for (int shard = 0; shard != threads; ++ shard)
    reducers.add_thread(new boost::thread(reducer_type(queues_reducer[shard],
						       utils::tempfile::tmp_dir(),
						       source_files[shard],
						       joint_counts[shard],
						       source_counts[shard],
						       max_malloc,
						       debug)));  
  
  boost::thread_group mappers;
  for (int shard = 0; shard != threads; ++ shard)
    mappers.add_thread(new boost::thread(mapper_type(mapped_files[shard],
						     queues_mapper[shard],
						     max_malloc,
						     debug)));

  mappers.join_all();
  reducers.join_all();
  
  // merge root_joint and root_sources...
  for (size_t shard = 0; shard != source_counts.size(); ++ shard) {
    root_count_set_type::const_iterator jiter_end = joint_counts[shard].end();
    for (root_count_set_type::const_iterator jiter = joint_counts[shard].begin(); jiter != jiter_end; ++ jiter) {
      std::pair<root_count_set_type::iterator, bool> result = root_joint.insert(*jiter);
      if (! result.second) {
	const_cast<root_count_type&>(*result.first).increment(jiter->counts.begin(), jiter->counts.end());
	const_cast<root_count_type&>(*result.first).observed += jiter->observed;
      }
    }
    
    root_count_set_type::const_iterator citer_end = source_counts[shard].end();
    for (root_count_set_type::const_iterator citer = source_counts[shard].begin(); citer != citer_end; ++ citer) {
      std::pair<root_count_set_type::iterator, bool> result = root_sources.insert(*citer);
      if (! result.second) {
	const_cast<root_count_type&>(*result.first).increment(citer->counts.begin(), citer->counts.end());
	const_cast<root_count_type&>(*result.first).observed += citer->observed;
      }
    }
  }
}


void reverse_counts(const path_set_type& counts_files,
		    path_map_type& reversed_files)
{
  typedef PhrasePairReverse map_reduce_type;
  
  typedef PhrasePairReverseMapper  mapper_type;
  typedef PhrasePairReverseReducer reducer_type;

  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  path_map_type mapped_files(threads);
  for (size_t i = 0; i != counts_files.size(); ++ i)
    mapped_files[i % threads].push_back(counts_files[i]);
  
  reversed_files.clear();
  reversed_files.reserve(threads);
  reversed_files.resize(threads);

  queue_ptr_set_type  queues(threads);
  
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    queues[shard].reset(new queue_type(1024 * threads));
  
  boost::thread_group reducers;
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    reducers.add_thread(new boost::thread(reducer_type(*queues[shard],
						       utils::tempfile::tmp_dir(),
						       reversed_files[shard],
						       threads,
						       max_malloc,
						       debug)));

  boost::thread_group mappers;
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    mappers.add_thread(new boost::thread(mapper_type(mapped_files[shard],
						     queues,
						     max_malloc,
						     debug)));
  
  reducers.join_all();
  mappers.join_all();
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",                  po::value<path_set_type>(&input_files)->multitoken(), "input files")
    ("output",                 po::value<path_type>(&output_file),                   "output directory")
    
    ("score-phrase", po::bool_switch(&score_phrase), "score phrase pair counts")
    ("score-scfg",   po::bool_switch(&score_scfg),   "score synchronous-CFG counts")
    ("score-ghkm",   po::bool_switch(&score_ghkm),   "score ghkm fragment counts")
    
    ("max-malloc", po::value<double>(&max_malloc), "maximum malloc in GB")
    ("threads", po::value<int>(&threads), "# of threads")
    
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

