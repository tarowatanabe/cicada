//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada_extract_score_impl.hpp"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>

#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/lexical_cast.hpp>

typedef boost::filesystem::path                                    path_type;
typedef std::vector<path_type, std::allocator<path_type> >         path_set_type;
typedef std::vector<path_set_type, std::allocator<path_set_type> > path_map_type;

typedef RootCount root_count_type;
typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
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

path_set_type counts_files;
path_type     list_file;

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

path_type output_file = "";

bool score_phrase = false;
bool score_scfg   = false;
bool score_ghkm   = false;

double max_malloc = 8; // 8 GB
int    threads = 1;

int debug = 0;


void modify_counts(const path_set_type& counts_files,
		   path_map_type& modified_files);
template <typename Extractor>
void reverse_counts(const path_map_type& modified_files,
		    path_map_type& reversed_files,
		    root_count_set_type& root_targets);
template <typename Extractor, typename Lexicon>
void score_counts(const path_type& output_file,
		  const path_set_type& counts_files,
		  const path_map_tyep& reversed_files,
		  root_count_set_type& root_sources,
		  const Lexicon& lexicon);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (counts_files.empty() && (list_file.empty() || (! boost::filesystem::exists(list_file) && list_file != "-")))
      throw std::runtime_error("no count files?");
    if (output_file.empty())
      throw std::runtime_error("no output file?");
    if (lexicon_source_target_file.empty() || ! boost::filesystem::exists(lexicon_source_target_file))
      throw std::runtime_error("no lexicon model for lex(target | source): " + lexicon_source_target_file.string());
    if (lexicon_target_source_file.empty() || ! boost::filesystem::exists(lexicon_target_source_file))
      throw std::runtime_error("no lexicon model for lex(source | target): " + lexicon_target_source_file.string());

    if (int(score_phrase) + score_scfg + score_ghkm != 1)
      throw std::runtime_error("specify either one of --score-phrase|scfg|ghkm");
    
    threads = utils::bithack::max(1, threads);

    // sort input files by its size...
    if (! list_file.empty()) {
      const path_type path = (boost::filesystem::is_directory(list_file) ? list_file / "files" : list_file);
      
      if (! boost::filesystem::exists(path) && path != "-")
	throw std::runtime_error("no file? " + path.string());

      const path_type dirname = (path == "-" ? path_type() : path.parent_path());
      
      utils::compress_istream is(path, 1024 * 1024);
      std::string line;
      while (std::getline(is, line)) 
	if (! line.empty()) {
	  const path_type path(line);
	  
	  if (boost::filesystem::exists(path))
	    counts_files.push_back(path);
	  else if (! dirname.empty() && boost::filesystem::exists(dirname / path))
	    counts_files.push_back(dirname / path);
	}

      if (counts_files.empty())
	throw std::runtime_error("no count files?");
    }
    
    std::sort(counts_files.begin(), counts_files.end(), greater_file_size());
    
    // modify counts...
    path_map_type modified_files;
    root_count_set_type root_sources;
    root_count_set_type root_targets;
    
    utils::resource start_modify;
    modify_counts(counts_files, modified_files);
    utils::resource end_modify;
    
    if (debug)
      std::cerr << "modify counts cpu time:  " << end_modify.cpu_time() - start_modify.cpu_time() << std::endl
		<< "modify counts user time: " << end_modify.user_time() - start_modify.user_time() << std::endl;
    
    // indexing...
    modified_counts_set_type modified_counts(threads);
    
    utils::resource start_index;
    if (score_phrase)
      index_counts<ExtractRootPhrase>(modified_files, modified_counts, root_targets);
    else if (score_scfg)
      index_counts<ExtractRootSCFG>(modified_files, modified_counts, root_targets);
    else
      index_counts<ExtractRootGHKM>(modified_files, modified_counts, root_targets);
    utils::resource end_index;
 
    if (debug)
      std::cerr << "index counts cpu time:  " << end_index.cpu_time() - start_index.cpu_time() << std::endl
		<< "index counts user time: " << end_index.user_time() - start_index.user_time() << std::endl;
   
    // scoring...
    const LexiconModel lexicon_source_target(lexicon_source_target_file);
    const LexiconModel lexicon_target_source(lexicon_target_source_file);
    
    utils::resource start_score;
    if (score_phrase)
      score_counts<ExtractRootPhrase, LexiconPhrase>(output_file,
						     counts_files,
						     modified_counts,
						     root_sources,
						     LexiconPhrase(lexicon_source_target, lexicon_target_source));
      
    else if (score_scfg)
      score_counts<ExtractRootSCFG, LexiconSCFG>(output_file,
						 counts_files,
						 modified_counts,
						 root_sources,
						 LexiconSCFG(lexicon_source_target, lexicon_target_source));
    else
      score_counts<ExtractRootGHKM, LexiconGHKM>(output_file,
						 counts_files,
						 modified_counts,
						 root_sources,
						 LexiconGHKM(lexicon_source_target, lexicon_target_source));
    utils::resource end_score;
    if (debug)
      std::cerr << "score counts cpu time:  " << end_score.cpu_time() - start_score.cpu_time() << std::endl
		<< "score counts user time: " << end_score.user_time() - start_score.user_time() << std::endl;
    
    // finally, dump files, root-sources and root-targets...
    {
      utils::compress_ostream os_file(output_file / "files");
      utils::compress_ostream os_src(output_file / "root-source.gz");
      utils::compress_ostream os_trg(output_file / "root-target.gz");
      
      os_src.precision(20);
      os_trg.precision(20);

      for (int shard = 0; shard != threads; ++ shard)
	os_file << (utils::lexical_cast<std::string>(shard) + ".gz") << '\n';
      
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


template <typename Extractor, typename Lexicon>
void score_counts(const path_type& output_file,
		  const path_set_type& counts_files,
		  const modified_counts_set_type& modified,
		  root_count_set_type& root_sources,
		  const Lexicon& lexicon)
{
  typedef PhrasePairScore map_reduce_type;
  
  typedef PhrasePairScoreMapper             mapper_type;
  typedef PhrasePairScoreReducer<Extractor, Lexicon> reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef std::vector<queue_ptr_set_type, std::allocator<queue_ptr_set_type> > queue_ptr_map_type;

  typedef utils::compress_ostream ostream_type;
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  
  if (static_cast<int>(modified.size()) != threads)
    throw std::runtime_error("# of threads differ");
  
  // create directories for output
  if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
    boost::filesystem::remove_all(output_file);
  
  boost::filesystem::create_directories(output_file);
  
  boost::filesystem::directory_iterator iter_end;
  for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
    boost::filesystem::remove_all(*iter);

  path_map_type mapped_files(threads);
  for (size_t i = 0; i != counts_files.size(); ++ i)
    mapped_files[i % threads].push_back(counts_files[i]);
  
  root_count_map_type root_counts(threads);
  
  boost::thread_group mappers;
  boost::thread_group reducers;
  
  queue_ptr_map_type   queues_mapper(threads, queue_ptr_set_type(threads));
  queue_ptr_map_type   queues_reducer(threads, queue_ptr_set_type(threads));
  ostream_ptr_set_type ostreams(threads);
  
  // construct queue matrix...
  for (int i = 0; i != threads; ++ i)
    for (int j = 0; j != threads; ++ j) {
      queues_mapper[i][j].reset(new queue_type(1024));
      queues_reducer[j][i] = queues_mapper[i][j];
    }
  
  for (int shard = 0; shard != threads; ++ shard)
    mappers.add_thread(new boost::thread(mapper_type(mapped_files[shard],
						     queues_mapper[shard],
						     max_malloc,
						     debug)));
  for (int shard = 0; shard != threads; ++ shard) {
    const path_type path = output_file / (utils::lexical_cast<std::string>(shard) + ".gz");
    
    // TODO: ZERO BUFFER!
    ostreams[shard].reset(new utils::compress_ostream(path, 0));
    
    reducers.add_thread(new boost::thread(reducer_type(modified,
						       root_counts[shard],
						       Extractor(),
						       lexicon,
						       queues_reducer[shard],
						       *ostreams[shard],
						       debug)));
  }
  
  mappers.join_all();
  reducers.join_all();
  
  
  // merge root_sources...
  for (size_t shard = 0; shard != root_counts.size(); ++ shard) {
    root_count_set_type::const_iterator citer_end = root_counts[shard].end();
    for (root_count_set_type::const_iterator citer = root_counts[shard].begin(); citer != citer_end; ++ citer) {
      std::pair<root_count_set_type::iterator, bool> result = root_sources.insert(*citer);
      if (! result.second) {
	const_cast<root_count_type&>(*result.first).increment(citer->counts.begin(), citer->counts.end());
	const_cast<root_count_type&>(*result.first).observed_joint += citer->observed_joint;
	const_cast<root_count_type&>(*result.first).observed       += citer->observed;
      }
    }
  }
}


template <typename Extractor>
struct IndexTask
{
  const path_set_type&  paths;
  modified_counts_type& counts;

  IndexTask(const path_set_type& __paths,
	    modified_counts_type& __counts)
    : paths(__paths), counts(__counts) {}
  
  void operator()()
  {
    counts.open(paths, Extractor());
  }
};

template <typename Extractor>
void index_counts(const path_map_type& modified_files,
		  modified_counts_set_type& modified_counts,
		  root_count_set_type& root_targets)
{
  if (static_cast<int>(modified_files.size()) != threads)
    throw std::runtime_error("# of threads differ");
  if (static_cast<int>(modified_counts.size()) != threads)
    throw std::runtime_error("# of threads differ");
  
  boost::thread_group workers;
  
  for (size_t shard = 0; shard != modified_counts.size(); ++ shard)
    workers.add_thread(new boost::thread(IndexTask<Extractor>(modified_files[shard], modified_counts[shard])));
  
  workers.join_all();
  
  // merge root counts...
  for (size_t shard = 0; shard != modified_counts.size(); ++ shard) {
    root_count_set_type::const_iterator citer_end = modified_counts[shard].root_counts.end();
    for (root_count_set_type::const_iterator citer = modified_counts[shard].root_counts.begin(); citer != citer_end; ++ citer) {
      std::pair<root_count_set_type::iterator, bool> result = root_targets.insert(*citer);
      if (! result.second) {
	const_cast<root_count_type&>(*result.first).increment(citer->counts.begin(), citer->counts.end());
	const_cast<root_count_type&>(*result.first).observed_joint += citer->observed_joint;
	const_cast<root_count_type&>(*result.first).observed       += citer->observed;
      }
    }
  }
}

void modify_counts(const path_set_type& counts_files,
		   path_map_type& modified_files)
{
  typedef PhrasePairModify map_reduce_type;
  
  typedef PhrasePairModifyMapper  mapper_type;
  typedef PhrasePairModifyReducer reducer_type;

  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  path_map_type mapped_files(threads);
  for (size_t i = 0; i != counts_files.size(); ++ i)
    mapped_files[i % threads].push_back(counts_files[i]);
  
  modified_files.clear();
  modified_files.reserve(threads);
  modified_files.resize(threads);

  boost::thread_group mappers;
  boost::thread_group reducers;
  queue_ptr_set_type  queues(threads);
  
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    queues[shard].reset(new queue_type(128));
  
  for (size_t shard = 0; shard != queues.size(); ++ shard)
    reducers.add_thread(new boost::thread(reducer_type(*queues[shard],
						       modified_files[shard],
						       threads,
						       max_malloc,
						       debug)));
  
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
    ("counts",                 po::value<path_set_type>(&counts_files)->multitoken(), "counts files")
    ("list",                   po::value<path_type>(&list_file),                      "count file list")
    ("lexicon-source-target",  po::value<path_type>(&lexicon_source_target_file),     "lexicon model for lex(target | source)")
    ("lexicon-target-source",  po::value<path_type>(&lexicon_target_source_file),     "lexicon model for lex(source | target)")
    ("output",                 po::value<path_type>(&output_file),                    "output directory")
    
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

