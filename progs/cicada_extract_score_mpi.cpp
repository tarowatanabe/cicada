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

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"

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


typedef PhraseCounts modified_counts_type;
typedef std::vector<modified_counts_type, std::allocator<modified_counts_type> > modified_counts_set_type;

path_set_type counts_files;
path_type     list_file;

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

path_type output_file = "";

bool score_phrase = false;
bool score_scfg   = false;
bool score_ghkm   = false;

double max_malloc = 8; // 8 GB
path_type prog_name;

int debug = 0;

void score_counts_mapper(utils::mpi_intercomm& reducer,
			 const path_set_type& counts_files);
template <typename Extractor, typename Lexicon>
void score_counts_reducer(utils::mpi_intercomm& mapper,
			  const path_type& output_file,
			  const modified_counts_set_type& modified,
			  root_count_set_type& root_sources,
			  const Lexicon& lexicon);
template <typename Extractor>
void index_counts(const path_set_type& modified_files,
		  modified_counts_set_type& modified_counts,
		  root_count_set_type& root_targets);
void modify_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files);
void modify_counts_reducer(utils::mpi_intercomm& mapper,
			   path_set_type& modified_files);
void synchronize_mapper(utils::mpi_intercomm& reducer);
void synchronize_reducer(utils::mpi_intercomm& mapper);

void options(int argc, char** argv);

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

    if (! prog_name.empty() && ! boost::filesystem::exists(prog_name))
      throw std::runtime_error(std::string("no binary? ") + prog_name.string());
    
    if (MPI::Comm::Get_parent() != MPI::COMM_NULL) {
      utils::mpi_intercomm comm_parent(MPI::Comm::Get_parent());
      
      path_set_type modified_files;
      root_count_set_type root_sources;
      root_count_set_type root_targets;
      
      utils::resource start_modify;
      modify_counts_reducer(comm_parent, modified_files);
      utils::resource end_modify;
      
      if (debug && mpi_rank == 0)
	std::cerr << "modify counts reducer cpu time:  " << end_modify.cpu_time() - start_modify.cpu_time() << std::endl
		  << "modify counts reducer  user time: " << end_modify.user_time() - start_modify.user_time() << std::endl;
            
      modified_counts_set_type modified_counts(mpi_size);
      
      utils::resource start_index;
      if (score_phrase)
	index_counts<ExtractRootPhrase>(modified_files, modified_counts, root_targets);
      else if (score_scfg)
	index_counts<ExtractRootSCFG>(modified_files, modified_counts, root_targets);
      else
	index_counts<ExtractRootGHKM>(modified_files, modified_counts, root_targets);
      
      utils::resource end_index;
      
      if (debug && mpi_rank == 0)
	std::cerr << "index counts cpu time:  " << end_index.cpu_time() - start_index.cpu_time() << std::endl
		  << "index counts user time: " << end_index.user_time() - start_index.user_time() << std::endl;
      
      // synchronize here...
      synchronize_reducer(comm_parent);
      
      // scoring...
      const LexiconModel lexicon_source_target(lexicon_source_target_file);
      const LexiconModel lexicon_target_source(lexicon_target_source_file);
      
      utils::resource start_score;
      if (score_phrase)
	score_counts_reducer<ExtractRootPhrase, LexiconPhrase>(comm_parent,
							       output_file,
							       modified_counts,
							       root_sources,
							       LexiconPhrase(lexicon_source_target, lexicon_target_source));
      else if (score_scfg)
	score_counts_reducer<ExtractRootSCFG, LexiconSCFG>(comm_parent,
							   output_file,
							   modified_counts,
							   root_sources,
							   LexiconSCFG(lexicon_source_target, lexicon_target_source));
      else
	score_counts_reducer<ExtractRootGHKM, LexiconGHKM>(comm_parent,
							   output_file,
							   modified_counts,
							   root_sources,
							   LexiconGHKM(lexicon_source_target, lexicon_target_source));
      utils::resource end_score;
      if (debug && mpi_rank == 0)
	std::cerr << "score counts reducer cpu time:  " << end_score.cpu_time() - start_score.cpu_time() << std::endl
		  << "score counts reducer user time: " << end_score.user_time() - start_score.user_time() << std::endl;
      
      // finally, dump root-sources and root-targets...
      if (mpi_rank == 0) {
	utils::compress_ostream os_file(output_file / "files");
	utils::compress_ostream os_src(output_file / "root-source.gz");
	utils::compress_ostream os_trg(output_file / "root-target.gz");
      
	os_src.precision(20);
	os_trg.precision(20);
      
	for (int shard = 0; shard != mpi_size; ++ shard)
	  os_file << (utils::lexical_cast<std::string>(shard) + ".gz") << '\n';
	
	root_count_set_type::const_iterator siter_end = root_sources.end();
	for (root_count_set_type::const_iterator siter = root_sources.begin(); siter != siter_end; ++ siter)
	  os_src << *siter << '\n';
	
	root_count_set_type::const_iterator titer_end = root_targets.end();
	for (root_count_set_type::const_iterator titer = root_targets.begin(); titer != titer_end; ++ titer)
	  os_trg << *titer << '\n';
      }
      
      // synchronize here...
      synchronize_reducer(comm_parent);
    } else {
      const std::string name = (boost::filesystem::exists(prog_name) ? prog_name.string() : std::string(argv[0]));
      utils::mpi_intercomm comm_child(MPI::COMM_WORLD.Spawn(name.c_str(), &(*args.begin()), mpi_size, MPI::INFO_NULL, 0));
      
      // sort input files by its size...
      if (mpi_rank == 0) {
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
	
	boost::iostreams::filtering_ostream os;
	os.push(boost::iostreams::gzip_compressor());
	os.push(utils::mpi_device_bcast_sink(0, 4096));
	
	path_set_type::const_iterator citer_end = counts_files.end();
	for (path_set_type::const_iterator citer = counts_files.begin(); citer != citer_end; ++ citer)
	  os << citer->string() << '\n';
      } else {
	counts_files.clear();
	
	boost::iostreams::filtering_istream is;
	is.push(boost::iostreams::gzip_decompressor());
	is.push(utils::mpi_device_bcast_source(0, 4096));
	
	std::string line;
	while (std::getline(is, line))
	  counts_files.push_back(line);
      }

      if (debug && mpi_rank == 0)
	std::cerr << "count files: " << counts_files.size() << std::endl;
      
      utils::resource start_modify;
      modify_counts_mapper(comm_child, counts_files);
      utils::resource end_modify;
      
      if (debug && mpi_rank == 0)
	std::cerr << "modify counts mapper cpu time:  " << end_modify.cpu_time() - start_modify.cpu_time() << std::endl
		  << "modify counts mapper user time: " << end_modify.user_time() - start_modify.user_time() << std::endl;
      
      // synchronize here...
      synchronize_mapper(comm_child);
      
      utils::resource start_score;
      score_counts_mapper(comm_child, counts_files);
      utils::resource end_score;
      
      if (debug && mpi_rank == 0)
	std::cerr << "score counts mapper cpu time:  " << end_score.cpu_time() - start_score.cpu_time() << std::endl
		  << "score counts mapper user time: " << end_score.user_time() - start_score.user_time() << std::endl;

      // synchronize here...
      synchronize_mapper(comm_child);
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  return 0;
}

enum {
  root_count_tag = 2000,
  modified_tag,
  phrase_pair_tag,
  notify_tag,
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

void synchronize_mapper(utils::mpi_intercomm& reducer)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  std::vector<MPI::Request, std::allocator<MPI::Request> > request(mpi_size);
  std::vector<bool, std::allocator<bool> > terminated(mpi_size, false);
  
  for (int rank = 0; rank != mpi_size; ++ rank)
    request[rank] = reducer.comm.Irecv(0, 0, MPI::INT, rank, notify_tag);
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    for (int rank = 0; rank != mpi_size; ++ rank)
      if (! terminated[rank] && request[rank].Test()) {
	terminated[rank] = true;
	found = true;
      }
    
    if (std::count(terminated.begin(), terminated.end(), true) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
}

void synchronize_reducer(utils::mpi_intercomm& mapper)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  std::vector<MPI::Request, std::allocator<MPI::Request> > request(mpi_size);
  std::vector<bool, std::allocator<bool> > terminated(mpi_size, false);
  
  for (int rank = 0; rank != mpi_size; ++ rank)
    request[rank] = mapper.comm.Isend(0, 0, MPI::INT, rank, notify_tag);
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    for (int rank = 0; rank != mpi_size; ++ rank)
      if (! terminated[rank] && request[rank].Test()) {
	terminated[rank] = true;
	found = true;
      }
    
    if (std::count(terminated.begin(), terminated.end(), true) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
}


void score_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files)
{
  typedef boost::iostreams::filtering_ostream ostream_type;
  typedef utils::mpi_device_sink              odevice_type;
  
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef boost::shared_ptr<odevice_type> odevice_ptr_type;
  
  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;
  
  typedef PhrasePairScore map_reduce_type;
  
  typedef PhrasePairScoreMapper  mapper_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  
  typedef PhrasePairGenerator generator_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  path_set_type mapped_files;
  for (size_t i = 0; i != counts_files.size(); ++ i)
    if (static_cast<int>(i % mpi_size) == mpi_rank)
      mapped_files.push_back(counts_files[i]);
  
  ostream_ptr_set_type stream(mpi_size);
  odevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    stream[rank].reset(new ostream_type());
    device[rank].reset(new odevice_type(reducer.comm, rank, phrase_pair_tag, 1024 * 1024, false, true));
    
    stream[rank]->push(boost::iostreams::gzip_compressor());
    stream[rank]->push(*device[rank]);
    
    queues[rank].reset(new queue_type(1024));
    
    stream[rank]->precision(20);
  }
  
  phrase_pair_type phrase_pair;
  generator_type   generator;
  
  boost::thread_group mapper;
  mapper.add_thread(new boost::thread(mapper_type(mapped_files, queues, max_malloc, debug)));
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    for (int rank = 0; rank != mpi_size; ++ rank)
      if (stream[rank] && device[rank]) {

	bool non_sleep = false;
	if (! device[rank]->test() || device[rank]->flush(true) != 0)
	  boost::thread::yield();
	else
	  non_sleep = true;
	
	if (queues[rank]->pop_swap(phrase_pair, true)) {
	  if (! phrase_pair.source.empty())
	    generator(*stream[rank], phrase_pair) << '\n';
	  else
	    stream[rank].reset();
	  
	  found |= non_sleep;
	}
      }
    
    found |= utils::mpi_terminate_devices(stream, device);
    
    if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  mapper.join_all();
}

template <typename Extractor, typename Lexicon>
void score_counts_reducer(utils::mpi_intercomm& mapper,
			  const path_type& output_file,
			  const modified_counts_set_type& modified,
			  root_count_set_type& root_sources,
			  const Lexicon& lexicon)
{
  typedef boost::iostreams::filtering_istream istream_type;
  typedef utils::mpi_device_source            idevice_type;

  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  typedef boost::shared_ptr<idevice_type> idevice_ptr_type;
  
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  typedef std::vector<idevice_ptr_type, std::allocator<idevice_ptr_type> > idevice_ptr_set_type;

  typedef PhrasePairScore map_reduce_type;
  
  typedef PhrasePairScoreReducer<Extractor, Lexicon> reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  
  typedef PhrasePairParser parser_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  // create directories for output
  if (mpi_rank == 0) {
    if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
      boost::filesystem::remove_all(output_file);
    
    boost::filesystem::create_directories(output_file);
    
    boost::filesystem::directory_iterator iter_end;
    for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
      boost::filesystem::remove_all(*iter);
  }

  MPI::COMM_WORLD.Barrier();
  
  istream_ptr_set_type stream(mpi_size);
  idevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  const size_t queue_size = 1024 * 128;
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    stream[rank].reset(new istream_type());
    device[rank].reset(new idevice_type(mapper.comm, rank, phrase_pair_tag, 1024 * 1024));

    queues[rank].reset(new queue_type(queue_size));
    
    stream[rank]->push(boost::iostreams::gzip_decompressor());
    stream[rank]->push(*device[rank]);
  }
  
  utils::compress_ostream os(output_file / (utils::lexical_cast<std::string>(mpi_rank) + ".gz"), 1024 * 1024);
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(modified,
						    root_sources,
						    Extractor(),
						    lexicon,
						    queues,
						    os,
						    debug)));
  
  phrase_pair_type phrase_pair;
  parser_type      parser;
  std::string line;
   
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    for (int rank = 0; rank < mpi_size; ++ rank)
      while (stream[rank] && device[rank] && device[rank]->test() && queues[rank]->size() < queue_size) {
	if (std::getline(*stream[rank], line)) {
	  
	  if (parser(line, phrase_pair))
	    queues[rank]->push_swap(phrase_pair);
	  else
	    std::cerr << "failed phrase-pair parsing: " << line << std::endl;
	} else {
	  phrase_pair.clear();
	  queues[rank]->push_swap(phrase_pair);
	  
	  stream[rank].reset();
	  device[rank].reset();
	}
	
	found = true;
      }
    
    if (std::count(device.begin(), device.end(), idevice_ptr_type()) == mpi_size)
      break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  reducer.join_all();
  
  // merge root-sources...
  if (mpi_rank == 0) {
    RootCountParser parser;
    root_count_type root_count;
    std::string line;
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::gzip_decompressor());
      is.push(utils::mpi_device_source(rank, root_count_tag, 1024 * 1024));
      
      while (std::getline(is, line)) {
	if (! parser(line, root_count)) {
	  std::cerr << "warning: root-count parsing failed: " << line << std::endl;
	  continue;
	}
	
	std::pair<root_count_set_type::iterator, bool> result = root_sources.insert(root_count);
	if (! result.second) {
	  const_cast<root_count_type&>(*result.first).increment(root_count.counts.begin(), root_count.counts.end());
	  const_cast<root_count_type&>(*result.first).observed_joint += root_count.observed_joint;
	  const_cast<root_count_type&>(*result.first).observed       += root_count.observed;
	}
      }
    }
    
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::gzip_compressor());
    os.push(utils::mpi_device_sink(0, root_count_tag, 1024 * 1024));
    
    RootCountGenerator generator;
    root_count_set_type::const_iterator siter_end = root_sources.end();
    for (root_count_set_type::const_iterator siter = root_sources.begin(); siter != siter_end; ++ siter)
      generator(os, *siter) << '\n';
  }
}

template <typename Extractor>
void index_counts(const path_set_type& modified_files,
		  modified_counts_set_type& modified_counts,
		  root_count_set_type& root_targets)
{
  typedef utils::repository repository_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  modified_counts[mpi_rank].open(modified_files, Extractor());
  
  path_type path_index;
  if (mpi_rank == 0) {
    path_index = utils::tempfile::directory_name(std::string("cicada.extract.index.XXXXXX"));
    repository_type rep(path_index, repository_type::write);
    
    boost::iostreams::filtering_ostream os;
    os.push(utils::mpi_device_bcast_sink(0, 4096));
    
    os << path_index.string() << '\n';
  } else {
    boost::iostreams::filtering_istream is;
    is.push(utils::mpi_device_bcast_source(0, 4096));
    
    std::string line;
    if (std::getline(is, line))
      path_index = line;
  }
  
  utils::tempfile::insert(path_index);
  
  if (path_index.empty())
    throw std::runtime_error("no index path? " + path_index.string());

  while (! repository_type::exists(path_index))
    boost::thread::yield();
  
  repository_type rep(path_index, repository_type::read);
  
  modified_counts[mpi_rank].write(rep.path(utils::lexical_cast<std::string>(mpi_rank)));
  
  ::sync();
  
  // merge root-targets.
  if (mpi_rank == 0) {
    root_targets = modified_counts[0].root_counts;
    
    RootCountParser parser;
    root_count_type root_count;
    std::string line;
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::gzip_decompressor());
      is.push(utils::mpi_device_source(rank, root_count_tag, 1024 * 1024));
      
      while (std::getline(is, line)) {
	if (! parser(line, root_count)) {
	  std::cerr << "warning: root-count parsing failed: " << line << std::endl;
	  continue;
	}
	
	std::pair<root_count_set_type::iterator, bool> result = root_targets.insert(root_count);
	if (! result.second) {
	  const_cast<root_count_type&>(*result.first).increment(root_count.counts.begin(), root_count.counts.end());
	  const_cast<root_count_type&>(*result.first).observed_joint += root_count.observed_joint;
	  const_cast<root_count_type&>(*result.first).observed       += root_count.observed;
	}
      }
    }
  } else {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::gzip_compressor());
    os.push(utils::mpi_device_sink(0, root_count_tag, 1024 * 1024));
    
    RootCountGenerator generator;
    root_count_set_type::const_iterator titer_end = modified_counts[mpi_rank].root_counts.end();
    for (root_count_set_type::const_iterator titer = modified_counts[mpi_rank].root_counts.begin(); titer != titer_end; ++ titer)
      generator(os, *titer) << '\n';
  }
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    MPI::COMM_WORLD.Barrier();
    
    const path_type path = rep.path(utils::lexical_cast<std::string>(rank));
    
    while (! modified_counts[rank].exists(path))
      boost::thread::yield();
    
    modified_counts[rank].open(path);
  }
}


void modify_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files)
{
  typedef boost::iostreams::filtering_ostream ostream_type;
  typedef utils::mpi_device_sink              odevice_type;
  
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef boost::shared_ptr<odevice_type> odevice_ptr_type;

  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;

  typedef PhrasePairModify map_reduce_type;
  
  typedef PhrasePairModifyMapper  mapper_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::modified_type     modified_type;
  typedef map_reduce_type::modified_set_type modified_set_type;

  typedef PhrasePairModifiedGenerator modified_generator_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  path_set_type mapped_files;
  for (size_t i = 0; i != counts_files.size(); ++ i)
    if (static_cast<int>(i % mpi_size) == mpi_rank)
      mapped_files.push_back(counts_files[i]);
  
  ostream_ptr_set_type stream(mpi_size);
  odevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new odevice_type(reducer.comm, rank, modified_tag, 1024 * 1024, false, true));
    
    stream[rank].reset(new ostream_type());
    stream[rank]->push(boost::iostreams::gzip_compressor());
    stream[rank]->push(*device[rank]);
    stream[rank]->precision(20);
    
    queues[rank].reset(new queue_type(128));
  }

  if (debug >= 2)
    std::cerr << "modify counts: rank: " << mpi_rank << " files: " << mapped_files.size() << std::endl;

  boost::thread_group mapper;
  mapper.add_thread(new boost::thread(mapper_type(mapped_files, queues, max_malloc, debug)));
  
  modified_set_type       modified;
  modified_generator_type generator;

  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    for (int rank = 0; rank != mpi_size; ++ rank)
      if (stream[rank] && device[rank]) {

	bool non_sleep = false;
	
	if (! device[rank]->test() || device[rank]->flush(true) != 0)
	  boost::thread::yield();
	else
	  non_sleep = true;
	
	if (queues[rank]->pop_swap(modified, true)) {
	  if (! modified.empty()) {
	    if (debug >= 5)
	      std::cerr << "modify counts mapper: " << modified.size() << std::endl;
	    
	    modified_set_type::const_iterator citer_end = modified.end();
	    for (modified_set_type::const_iterator citer = modified.begin(); citer != citer_end; ++ citer) {
	      generator(*stream[rank], *citer) << '\n';
	      
	      if (! device[rank]->test() || device[rank]->flush(true) != 0)
		boost::thread::yield();
	      else
		non_sleep = true;
	    }
	    
	    modified.clear();
	    modified_set_type(modified).swap(modified);
	  } else
	    stream[rank].reset();
	  
	  found |= non_sleep;
	}
      }
    
    found |= utils::mpi_terminate_devices(stream, device);
    
    if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  mapper.join_all();
}

void modify_counts_reducer(utils::mpi_intercomm& mapper,
			   path_set_type& modified_files)
{
  typedef boost::iostreams::filtering_istream istream_type;
  typedef utils::mpi_device_source            idevice_type;
  
  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  typedef boost::shared_ptr<idevice_type> idevice_ptr_type;
  
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  typedef std::vector<idevice_ptr_type, std::allocator<idevice_ptr_type> > idevice_ptr_set_type;
  
  typedef PhrasePairModify map_reduce_type;
  
  typedef PhrasePairModifyReducer reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::modified_type     modified_type;
  typedef map_reduce_type::modified_set_type modified_set_type;
  
  typedef PhrasePairModifiedParser modified_parser_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  istream_ptr_set_type stream(mpi_size);
  idevice_ptr_set_type device(mpi_size);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new idevice_type(mapper.comm, rank, modified_tag, 1024 * 1024));
    
    stream[rank].reset(new istream_type());
    stream[rank]->push(boost::iostreams::gzip_decompressor());
    stream[rank]->push(*device[rank]);
  }

  const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
  
  const size_t queue_size = mpi_size * 128;
  queue_type queue(queue_size);
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(queue, modified_files, 1, max_malloc, debug)));
  
  modified_set_type modified;
  modified_type     parsed;
  
  modified_parser_type parser;
  std::string line;
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    if (queue.size() < queue_size) {
      for (int rank = 0; rank != mpi_size; ++ rank) {
	for (int iter = 0; iter != 64 && stream[rank] && device[rank] && device[rank]->test(); ++ iter) {
	  if (std::getline(*stream[rank], line)) {
	    if (parser(line, parsed))
	      modified.push_back(parsed);
	    else
	      std::cerr << "failed modified phrase parsing: " << line << std::endl;
	  } else {
	    stream[rank].reset();
	    device[rank].reset();
	  }
	  
	  found = true;
	}
	
	if (found && utils::malloc_stats::used() > malloc_threshold) {
	  boost::thread::yield();
	  found = false;
	}
      }
    }
    
    const size_t modified_size = modified.size();
    
    if (! modified.empty() && queue.push_swap(modified, true)) {
      
      if (debug >= 5)
	std::cerr << "modify counts reducer: " << modified_size << std::endl;
      
      modified.clear();
      modified_set_type(modified).swap(modified);
      
      found = true;
    }
    
    if (modified.empty() && std::count(device.begin(), device.end(), idevice_ptr_type()) == mpi_size)
      break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  
  modified.clear();
  queue.push_swap(modified);
  
  reducer.join_all();
}


void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

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
    
    ("max-malloc", po::value<double>(&max_malloc),   "maximum malloc in GB")
    ("prog",       po::value<path_type>(&prog_name), "this binary");
  
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

