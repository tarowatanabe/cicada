//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada_extract_score_impl.hpp"
#include "cicada_output_impl.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/random.hpp>

#include <utils/filesystem.hpp>
#include <utils/resource.hpp>
#include <utils/bithack.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/random_seed.hpp>

#include <codec/lz4.hpp>

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"

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
path_type temporary_dir = "";

bool score_phrase = false;
bool score_scfg   = false;
bool score_ghkm   = false;

double max_malloc = 8; // 8 GB
path_type prog_name;

int debug = 0;

void score_counts_mapper(utils::mpi_intercomm& reducer,
			 const path_set_type& counts_files);
void score_counts_reducer(utils::mpi_intercomm& mapper,
			  const path_type& output_file,
			  const path_type& source_file,
			  const path_set_type& target_files);

void source_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files);
template <typename Extractor>
void source_counts_reducer(utils::mpi_intercomm& mapper,
			   path_type& source_file,
			   root_count_set_type& root_joint,
			   root_count_set_type& root_source);

template <typename Extractor>
void target_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files,
			  root_count_set_type& root_counts);
void target_counts_reducer(utils::mpi_intercomm& mapper,
			   path_set_type& target_files);

void reverse_counts_mapper(utils::mpi_intercomm& reducer,
			   const path_set_type& counts_files,
			   path_set_type& reversed_files);
void reverse_counts_reducer(utils::mpi_intercomm& mapper,
			    const path_type& output_file);

void synchronize_mapper(utils::mpi_intercomm& reducer);
void synchronize_reducer(utils::mpi_intercomm& mapper);

void reduce_root_counts(root_count_set_type& root_counts);

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
    
    if (! temporary_dir.empty())
      ::setenv("TMPDIR_SPEC", temporary_dir.string().data(), 1);
    
    if (output_file.empty())
      throw std::runtime_error("no output file?");
    
    if (int(score_phrase) + score_scfg + score_ghkm != 1)
      throw std::runtime_error("specify either one of --score-phrase|scfg|ghkm");

    if (! prog_name.empty() && ! boost::filesystem::exists(prog_name))
      throw std::runtime_error(std::string("no binary? ") + prog_name.string());
    
    if (MPI::Comm::Get_parent() != MPI::COMM_NULL) {
      utils::mpi_intercomm comm_parent(MPI::Comm::Get_parent());
      
      path_type     source_file;
      path_set_type target_files;
      root_count_set_type root_joint;
      root_count_set_type root_sources;
      
      utils::resource start_source;
      if (score_phrase)
	source_counts_reducer<ExtractRootPhrase>(comm_parent, source_file, root_joint, root_sources);
      else if (score_scfg)
	source_counts_reducer<ExtractRootSCFG>(comm_parent, source_file, root_joint, root_sources);
      else
	source_counts_reducer<ExtractRootGHKM>(comm_parent, source_file, root_joint, root_sources);
      utils::resource end_source;
      
      if (debug && mpi_rank == 0)
	std::cerr << "source counts reducer"
		  << " cpu time:  " << end_source.cpu_time() - start_source.cpu_time() 
		  << " user time: " << end_source.user_time() - start_source.user_time()
		  << std::endl;
      
      if (mpi_rank == 0)
	prepare_directory(output_file);
      
      MPI::COMM_WORLD.Barrier();
      
      utils::resource start_reverse;
      reverse_counts_reducer(comm_parent, output_file);
      utils::resource end_reverse;
      
      if (debug && mpi_rank == 0)
	std::cerr << "reverse counts reducer"
		  << " cpu time:  " << end_reverse.cpu_time() - start_reverse.cpu_time()
		  << " user time: " << end_reverse.user_time() - start_reverse.user_time()
		  << std::endl;
      
      utils::resource start_target;
      target_counts_reducer(comm_parent, target_files);
      utils::resource end_target;
      
      if (debug && mpi_rank == 0)
	std::cerr << "target counts reducer"
		  << " cpu time:  " << end_target.cpu_time() - start_target.cpu_time()
		  << " user time: " << end_target.user_time() - start_target.user_time()
		  << std::endl;
      
      // scoring...
      utils::resource start_score;
      score_counts_reducer(comm_parent, output_file, source_file, target_files);
      utils::resource end_score;
      
      if (debug && mpi_rank == 0)
	std::cerr << "score counts reducer"
		  << " cpu time:  " << end_score.cpu_time() - start_score.cpu_time()
		  << " user time: " << end_score.user_time() - start_score.user_time()
		  << std::endl;
      
      reduce_root_counts(root_joint);
      
      MPI::COMM_WORLD.Barrier();
      
      reduce_root_counts(root_sources);
      
      // finally, dump root-sources and root-targets...
      if (mpi_rank == 0) {
	utils::compress_ostream os_file(output_file / "files");
	utils::compress_ostream os_joint(output_file / "root-joint.gz");
	utils::compress_ostream os_src(output_file / "root-source.gz");
	
	os_joint.precision(20);
	os_src.precision(20);
	
	for (int shard = 0; shard != mpi_size; ++ shard)
	  os_file << (utils::lexical_cast<std::string>(shard) + ".gz") << '\n';

	root_count_set_type::const_iterator jiter_end = root_joint.end();
	for (root_count_set_type::const_iterator jiter = root_joint.begin(); jiter != jiter_end; ++ jiter)
	  os_joint << *jiter << '\n';
	
	root_count_set_type::const_iterator siter_end = root_sources.end();
	for (root_count_set_type::const_iterator siter = root_sources.begin(); siter != siter_end; ++ siter)
	  os_src << *siter << '\n';
      }
      
      // synchronize here...
      synchronize_reducer(comm_parent);
    } else {
      std::vector<int, std::allocator<int> > error_codes(mpi_size, MPI_SUCCESS);
      const std::string name = (boost::filesystem::exists(prog_name) ? prog_name.string() : std::string(argv[0]));
      utils::mpi_intercomm comm_child(MPI::COMM_WORLD.Spawn(name.c_str(), &(*args.begin()), mpi_size, MPI::INFO_NULL, 0, &(*error_codes.begin())));
      
      for (size_t i = 0; i != error_codes.size(); ++ i)
	if (error_codes[i] != MPI_SUCCESS)
	  throw std::runtime_error("one of children failed to launch!");
      
      path_set_type counts_files;
      
      // sort input files by its size...
      if (mpi_rank == 0) {
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
	
	boost::iostreams::filtering_ostream os;
	//os.push(boost::iostreams::zlib_compressor());
	os.push(codec::lz4_compressor());
	os.push(utils::mpi_device_bcast_sink(0, 4096));
	
	path_set_type::const_iterator citer_end = counts_files.end();
	for (path_set_type::const_iterator citer = counts_files.begin(); citer != citer_end; ++ citer)
	  os << citer->string() << '\n';
      } else {
	boost::iostreams::filtering_istream is;
	//is.push(boost::iostreams::zlib_decompressor());
	is.push(codec::lz4_decompressor());
	is.push(utils::mpi_device_bcast_source(0, 4096));
	
	std::string line;
	while (std::getline(is, line))
	  counts_files.push_back(line);
      }
      
      if (debug && mpi_rank == 0)
	std::cerr << "count files: " << counts_files.size() << std::endl;
      
      utils::resource start_source;
      source_counts_mapper(comm_child, counts_files);
      utils::resource end_source;
      
      if (debug && mpi_rank == 0)
	std::cerr << "source counts mapper"
		  << " cpu time:  " << end_source.cpu_time() - start_source.cpu_time()
		  << " user time: " << end_source.user_time() - start_source.user_time()
		  << std::endl;
      
      path_set_type reversed_files;
      root_count_set_type root_targets;
      
      utils::resource start_reverse;
      reverse_counts_mapper(comm_child, counts_files, reversed_files);
      utils::resource end_reverse;
      
      if (debug && mpi_rank == 0)
	std::cerr << "reverse counts mapper"
		  << " cpu time:  " << end_reverse.cpu_time() - start_reverse.cpu_time()
		  << " user time: " << end_reverse.user_time() - start_reverse.user_time()
		  << std::endl;
      
      utils::resource start_target;
      if (score_phrase)
	target_counts_mapper<ExtractRootPhrase>(comm_child, reversed_files, root_targets);
      else if (score_scfg)
	target_counts_mapper<ExtractRootSCFG>(comm_child, reversed_files, root_targets);
      else
	target_counts_mapper<ExtractRootGHKM>(comm_child, reversed_files, root_targets);
      utils::resource end_target;
      
      if (debug && mpi_rank == 0)
	std::cerr << "target counts mapper"
		  << " cpu time:  " << end_target.cpu_time() - start_target.cpu_time()
		  << " user time: " << end_target.user_time() - start_target.user_time()
		  << std::endl;
      
      // remove all the reversed files
      for (path_set_type::const_iterator miter = reversed_files.begin(); miter != reversed_files.end(); ++ miter) {
	boost::filesystem::remove(*miter);
	utils::tempfile::erase(*miter);

	const path_type tmp_file = miter->parent_path() / miter->stem();
	if (boost::filesystem::exists(tmp_file)) {
	  boost::filesystem::remove(tmp_file);
	  utils::tempfile::erase(tmp_file);
	}
      }
      
      utils::resource start_score;
      score_counts_mapper(comm_child, counts_files);
      utils::resource end_score;
      
      if (debug && mpi_rank == 0)
	std::cerr << "score counts mapper"
		  << " cpu time:  " << end_score.cpu_time() - start_score.cpu_time()
		  << " user time: " << end_score.user_time() - start_score.user_time()
		  << std::endl;

      reduce_root_counts(root_targets);
      
      if (mpi_rank == 0) {
	utils::compress_ostream os_trg(output_file / "root-target.gz");
	os_trg.precision(20);
	
	root_count_set_type::const_iterator titer_end = root_targets.end();
	for (root_count_set_type::const_iterator titer = root_targets.begin(); titer != titer_end; ++ titer)
	  os_trg << *titer << '\n';
      }
      
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
  reverse_tag,
  source_tag,
  target_tag,
  phrase_pair_tag,
  file_tag,
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

void reduce_root_counts(root_count_set_type& root_counts)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    RootCountParser parser;
    root_count_type root_count;
    
    std::string line;
    for (int rank = 1; rank != mpi_size; ++ rank) {
      boost::iostreams::filtering_istream is;
      //is.push(boost::iostreams::zlib_decompressor());
      is.push(codec::lz4_decompressor());
      is.push(utils::mpi_device_source(rank, root_count_tag, 4096));
      
      while (std::getline(is, line)) {
	if (! parser(line, root_count)) {
	  std::cerr << "warning: root-count parsing failed: " << line << std::endl;
	  continue;
	}
	
	std::pair<root_count_set_type::iterator, bool> result = root_counts.insert(root_count);
	if (! result.second) {
	  const_cast<root_count_type&>(*result.first).increment(root_count.counts.begin(), root_count.counts.end());
	  const_cast<root_count_type&>(*result.first).observed += root_count.observed;
	}
      }
    }
  } else {
    boost::iostreams::filtering_ostream os;
    //os.push(boost::iostreams::zlib_compressor());
    os.push(codec::lz4_compressor());
    os.push(utils::mpi_device_sink(0, root_count_tag, 4096));
    
    RootCountGenerator generator;
    root_count_set_type::const_iterator riter_end = root_counts.end();
    for (root_count_set_type::const_iterator riter = root_counts.begin(); riter != riter_end; ++ riter)
      generator(os, *riter) << '\n';
  }
}

template <typename Mapper>
struct progress_mapper : public Mapper
{
  progress_mapper(const Mapper& mapper) : Mapper(mapper) {}

  void operator()()
  {
    Mapper::operator()(progress_type());
  }
  
  struct progress_type
  {
    progress_type() : i(0) {}
    
    void operator()() const
    {
      ++ const_cast<size_t&>(i);
      
      if (i % 100000 == 0)
	std::cerr << '.';
      if (i % 10000000 == 0)
	std::cerr << std::endl;
    }
    
    void final() const
    {
      if ((i % 100000) % 100)
	std::cerr << std::endl;
    }
  
    size_t i;
  };
};

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

  typedef std::vector<int, std::allocator<int> > rank_set_type;

  static const size_t buffer_size     = 1024 * 1024 * 4;
  static const size_t buffer_size_max = buffer_size << 2;
  static const size_t queue_size      = 1024;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  path_set_type mapped_files;
  for (size_t i = 0; i != counts_files.size(); ++ i)
    if (static_cast<int>(i % mpi_size) == mpi_rank)
      mapped_files.push_back(counts_files[i]);
  
  ostream_ptr_set_type stream(mpi_size);
  odevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  rank_set_type ranks(mpi_size);
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    stream[rank].reset(new ostream_type());
    device[rank].reset(new odevice_type(reducer.comm, rank, phrase_pair_tag, buffer_size, false, true));
    
    //stream[rank]->push(boost::iostreams::zlib_compressor());
    stream[rank]->push(codec::lz4_compressor());
    stream[rank]->push(*device[rank], buffer_size);
    stream[rank]->precision(20);
    
    queues[rank].reset(new queue_type(queue_size));
    
    ranks[rank] = rank;
  }
  
  phrase_pair_type phrase_pair;
  generator_type   generator;
  
  boost::thread_group mapper;
  if (mpi_rank == 0 && debug)
    mapper.add_thread(new boost::thread(progress_mapper<mapper_type>(mapper_type(mapped_files, queues, max_malloc, debug))));
  else
    mapper.add_thread(new boost::thread(mapper_type(mapped_files, queues, max_malloc, debug)));

  const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    rank_set_type::const_iterator riter_end = ranks.end();
    for (rank_set_type::const_iterator riter = ranks.begin(); riter != riter_end; ++ riter) {
      const int rank = *riter;
      
      if (stream[rank] && device[rank]) {
	
	if (device[rank]->test() && device[rank]->flush(true))
	  found = true;

	const size_t committed = device[rank]->committed();
	
	if (committed < buffer_size) {
	  while (static_cast<size_t>(device[rank]->committed()) < buffer_size && queues[rank]->pop_swap(phrase_pair, true)) {
	    found = true;
	    
	    if (! phrase_pair.source.empty())
	      generator(*stream[rank], phrase_pair) << '\n';
	    else {
	      stream[rank].reset();
	      break;
	    }
	  }
	} else if (committed < buffer_size_max || utils::malloc_stats::used() < malloc_threshold) {
	  if (queues[rank]->pop_swap(phrase_pair, true)) {
	    found |= (committed < buffer_size_max);
	    
	    if (! phrase_pair.source.empty())
	      generator(*stream[rank], phrase_pair) << '\n';
	    else
	      stream[rank].reset();
	  }
	}
      }
    }
    
    if (found)
      std::random_shuffle(ranks.begin(), ranks.end(), rgen);
    
    found |= utils::mpi_terminate_devices(stream, device);
    
    if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  mapper.join_all();
}

void score_counts_reducer(utils::mpi_intercomm& mapper,
			  const path_type& output_file,
			  const path_type& source_file,
			  const path_set_type& target_files)
{
  typedef boost::iostreams::filtering_istream istream_type;
  typedef utils::mpi_device_source            idevice_type;

  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  typedef boost::shared_ptr<idevice_type> idevice_ptr_type;
  
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  typedef std::vector<idevice_ptr_type, std::allocator<idevice_ptr_type> > idevice_ptr_set_type;

  typedef PhrasePairScore map_reduce_type;
  
  typedef PhrasePairScoreReducer reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  
  typedef PhrasePairParser parser_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  const size_t queue_size  = 1024 * 16;
  const size_t buffer_size = 1024 * 1024 * 4;
  
  istream_ptr_set_type stream(mpi_size);
  idevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    stream[rank].reset(new istream_type());
    device[rank].reset(new idevice_type(mapper.comm, rank, phrase_pair_tag, buffer_size));
    
    queues[rank].reset(new queue_type(queue_size));
    
    //stream[rank]->push(boost::iostreams::zlib_decompressor());
    stream[rank]->push(codec::lz4_decompressor());
    stream[rank]->push(*device[rank]);
  }
  
  utils::compress_ostream os(output_file / (utils::lexical_cast<std::string>(mpi_rank) + ".gz"), 1024 * 1024);
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(source_file,
						    target_files,
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
}

void source_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files)
{
  typedef boost::iostreams::filtering_ostream ostream_type;
  typedef utils::mpi_device_sink              odevice_type;
  
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef boost::shared_ptr<odevice_type> odevice_ptr_type;

  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;
  
  typedef PhrasePairSource       map_reduce_type;
  typedef PhrasePairSourceMapper mapper_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef map_reduce_type::simple_type simple_type;
  
  typedef PhrasePairSimpleGenerator simple_generator_type;
  
  typedef std::vector<int, std::allocator<int> > rank_set_type;
  
  static const size_t buffer_size     = 1024 * 1024 * 4;
  static const size_t buffer_size_max = buffer_size << 2;
  static const size_t queue_size      = 1024;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  path_set_type mapped_files;
  for (size_t i = 0; i != counts_files.size(); ++ i)
    if (static_cast<int>(i % mpi_size) == mpi_rank)
      mapped_files.push_back(counts_files[i]);
  
  ostream_ptr_set_type stream(mpi_size);
  odevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  rank_set_type ranks(mpi_size);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new odevice_type(reducer.comm, rank, source_tag, buffer_size, false, true));
    
    stream[rank].reset(new ostream_type());
    //stream[rank]->push(boost::iostreams::zlib_compressor());
    stream[rank]->push(codec::lz4_compressor());
    stream[rank]->push(*device[rank], buffer_size);
    stream[rank]->precision(20);
    
    queues[rank].reset(new queue_type(queue_size));
    
    ranks[rank] = rank;
  }
  
  if (debug >= 2)
    std::cerr << "source counts: rank: " << mpi_rank << " files: " << counts_files.size() << std::endl;

  boost::thread_group mapper;
  if (mpi_rank == 0 && debug)
    mapper.add_thread(new boost::thread(progress_mapper<mapper_type>(mapper_type(mapped_files, queues, max_malloc, debug))));
  else
    mapper.add_thread(new boost::thread(mapper_type(mapped_files, queues, max_malloc, debug)));
  
  simple_type simple;
  simple_generator_type generator;
  
  const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    rank_set_type::const_iterator riter_end = ranks.end();
    for (rank_set_type::const_iterator riter = ranks.begin(); riter != riter_end; ++ riter) {
      const int rank = *riter;
      
      if (stream[rank] && device[rank]) {
	
	if (device[rank]->test() && device[rank]->flush(true))
	  found = true;

	const size_t committed = device[rank]->committed();
	
	if (committed < buffer_size) {
	  while (static_cast<size_t>(device[rank]->committed()) < buffer_size && queues[rank]->pop_swap(simple, true)) {
	    found = true;
	    
	    if (! simple.source.empty())
	      generator(*stream[rank], simple) << '\n';
	    else {
	      stream[rank].reset();
	      break;
	    }
	  }
	} else if (committed < buffer_size_max || utils::malloc_stats::used() < malloc_threshold) {
	  if (queues[rank]->pop_swap(simple, true)) {
	    found |= (committed < buffer_size_max);
	    
	    if (! simple.source.empty())
	      generator(*stream[rank], simple) << '\n';
	    else 
	      stream[rank].reset();
	  }
	}
      }
    }
    
    if (found)
      std::random_shuffle(ranks.begin(), ranks.end(), rgen);
    
    found |= utils::mpi_terminate_devices(stream, device);
    
    if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  mapper.join_all();
}

template <typename Extractor>
void source_counts_reducer(utils::mpi_intercomm& mapper,
			   path_type& source_file,
			   root_count_set_type& root_joint,
			   root_count_set_type& root_source)
{
  typedef boost::iostreams::filtering_istream istream_type;
  typedef utils::mpi_device_source            idevice_type;
  
  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  typedef boost::shared_ptr<idevice_type> idevice_ptr_type;
  
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  typedef std::vector<idevice_ptr_type, std::allocator<idevice_ptr_type> > idevice_ptr_set_type;
  
  typedef PhrasePairSource                   map_reduce_type;
  typedef PhrasePairSourceReducer<Extractor> reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef map_reduce_type::simple_type     simple_type;

  typedef PhrasePairSimpleParser simple_parser_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  const size_t queue_size  = 1024 * 16;
  const size_t buffer_size = 1024 * 1024 * 4;
  
  istream_ptr_set_type stream(mpi_size);
  idevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new idevice_type(mapper.comm, rank, source_tag, buffer_size));
    stream[rank].reset(new istream_type());
    
    queues[rank].reset(new queue_type(queue_size));
    
    //stream[rank]->push(boost::iostreams::zlib_decompressor());
    stream[rank]->push(codec::lz4_decompressor());
    stream[rank]->push(*device[rank]);
  }
  
  queue_type queue(queue_size);
  boost::thread reducer(reducer_type(queues,
				     utils::tempfile::tmp_dir(),
				     source_file,
				     root_joint,
				     root_source,
				     max_malloc,
				     debug));
  
  simple_type        source;
  simple_parser_type parser;
  std::string line;
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    for (int rank = 0; rank < mpi_size; ++ rank)
      while (stream[rank] && device[rank] && device[rank]->test() && queues[rank]->size() < queue_size) {
	if (std::getline(*stream[rank], line)) {
	  if (parser(line, source))
	    queues[rank]->push_swap(source);
	  else
	    std::cerr << "failed phrase-pair-simple parsing: " << line << std::endl;
	} else {
	  source.clear();
	  queues[rank]->push_swap(source);
	  
	  stream[rank].reset();
	  device[rank].reset();
	}
	
	found = true;
      }
    
    if (std::count(device.begin(), device.end(), idevice_ptr_type()) == mpi_size)
      break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  reducer.join();
}


template <typename Extractor>
void target_counts_mapper(utils::mpi_intercomm& reducer,
			  const path_set_type& counts_files,
			  root_count_set_type& root_counts)
{
  typedef boost::iostreams::filtering_ostream ostream_type;
  typedef utils::mpi_device_sink              odevice_type;
  
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef boost::shared_ptr<odevice_type> odevice_ptr_type;

  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;

  typedef PhrasePairTarget map_reduce_type;
  
  typedef PhrasePairTargetMapper<Extractor>  mapper_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef map_reduce_type::simple_type     simple_type;
  
  typedef PhrasePairSimpleGenerator simple_generator_type;
  
  typedef std::vector<int, std::allocator<int> > rank_set_type;
  
  static const size_t buffer_size     = 1024 * 1024 * 4;
  static const size_t buffer_size_max = buffer_size << 2;
  static const size_t queue_size      = 1024;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  ostream_ptr_set_type stream(mpi_size);
  odevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  rank_set_type ranks(mpi_size);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new odevice_type(reducer.comm, rank, target_tag, buffer_size, false, true));
    
    stream[rank].reset(new ostream_type());
    //stream[rank]->push(boost::iostreams::zlib_compressor());
    stream[rank]->push(codec::lz4_compressor());
    stream[rank]->push(*device[rank], buffer_size);
    stream[rank]->precision(20);
    
    queues[rank].reset(new queue_type(queue_size));
    
    ranks[rank] = rank;
  }
  
  if (debug >= 2)
    std::cerr << "target counts: rank: " << mpi_rank << " files: " << counts_files.size() << std::endl;
  
  boost::thread_group mapper;
  if (mpi_rank == 0 && debug)
    mapper.add_thread(new boost::thread(progress_mapper<mapper_type>(mapper_type(counts_files, queues, root_counts, max_malloc, debug))));
  else
    mapper.add_thread(new boost::thread(mapper_type(counts_files, queues, root_counts, max_malloc, debug)));
  
  simple_type target;
  simple_generator_type generator;
  
  const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    rank_set_type::const_iterator riter_end = ranks.end();
    for (rank_set_type::const_iterator riter = ranks.begin(); riter != riter_end; ++ riter) {
      const int rank = *riter;
      
      if (stream[rank] && device[rank]) {
	
	if (device[rank]->test() && device[rank]->flush(true))
	  found = true;

	const size_t committed = static_cast<size_t>(device[rank]->committed());
	
	if (committed < buffer_size) {
	  for (int i = 0; i != 128 && static_cast<size_t>(device[rank]->committed()) < buffer_size && queues[rank]->pop_swap(target, true); ++ i) {
	    found = true;
	    
	    if (! target.source.empty())
	      generator(*stream[rank], target) << '\n';
	    else {
	      stream[rank].reset();
	      break;
	    }
	  }
	} else if (committed < buffer_size_max || utils::malloc_stats::used() < malloc_threshold) {
	  if (queues[rank]->pop_swap(target, true)) {
	    found |= (committed < buffer_size_max);
	    
	    if (! target.source.empty())
	      generator(*stream[rank], target) << '\n';
	    else 
	      stream[rank].reset();
	  }
	} 
      }
    }

    if (found)
      std::random_shuffle(ranks.begin(), ranks.end(), rgen);
    
    found |= utils::mpi_terminate_devices(stream, device);
    
    if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  mapper.join_all();
}

void target_counts_reducer(utils::mpi_intercomm& mapper,
			   path_set_type& target_files)
{
  typedef boost::iostreams::filtering_istream istream_type;
  typedef utils::mpi_device_source            idevice_type;
  
  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  typedef boost::shared_ptr<idevice_type> idevice_ptr_type;
  
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  typedef std::vector<idevice_ptr_type, std::allocator<idevice_ptr_type> > idevice_ptr_set_type;
  
  typedef PhrasePairTarget map_reduce_type;
  
  typedef PhrasePairTargetReducer reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef map_reduce_type::simple_type     simple_type;
  typedef map_reduce_type::simple_set_type simple_set_type;
  
  typedef PhrasePairSimpleParser simple_parser_type;

  typedef std::vector<int, std::allocator<int> > rank_set_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const size_t queue_size  = 1024 * 16;
  const size_t buffer_size = 1024 * 1024 * 4;
  
  istream_ptr_set_type stream(mpi_size);
  idevice_ptr_set_type device(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  rank_set_type ranks(mpi_size);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new idevice_type(mapper.comm, rank, target_tag, buffer_size));
    
    stream[rank].reset(new istream_type());
    //stream[rank]->push(boost::iostreams::zlib_decompressor());
    stream[rank]->push(codec::lz4_decompressor());
    stream[rank]->push(*device[rank]);
    
    ranks[rank] = rank;
  }
  
  queue_type queue(queue_size);
  boost::thread reducer(reducer_type(queue, utils::tempfile::tmp_dir(), target_files, 1, max_malloc, debug));
  
  simple_type     target;
  
  simple_parser_type parser;
  std::string line;
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    rank_set_type::const_iterator riter_end = ranks.end();
    for (rank_set_type::const_iterator riter = ranks.begin(); riter != riter_end; ++ riter) {
      const int rank = *riter;
      
      for (int i = 0; i != 128 && stream[rank] && device[rank] && device[rank]->test() && queue.size() < queue_size; ++ i) {
	if (std::getline(*stream[rank], line)) {
	  if (parser(line, target))
	    queue.push_swap(target);
	  else
	    std::cerr << "failed simple phrase parsing: " << line << std::endl;
	} else {
	  stream[rank].reset();
	  device[rank].reset();
	}
	
	found = true;
      }
    }
    
    if (found)
      std::random_shuffle(ranks.begin(), ranks.end(), rgen);
        
    if (std::count(device.begin(), device.end(), idevice_ptr_type()) == mpi_size)
      break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  target.clear();
  queue.push_swap(target);
  
  reducer.join();
}

void reverse_counts_mapper(utils::mpi_intercomm& reducer,
			   const path_set_type& counts_files,
			   path_set_type& reversed_files)
{
  typedef boost::iostreams::filtering_ostream ostream_type;
  typedef utils::mpi_device_sink              odevice_type;
  
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef boost::shared_ptr<odevice_type> odevice_ptr_type;

  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;

  typedef PhrasePairReverse map_reduce_type;
  
  typedef PhrasePairReverseMapper  mapper_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::simple_type     simple_type;
  typedef map_reduce_type::simple_set_type simple_set_type;

  typedef PhrasePairSimpleGenerator simple_generator_type;
  
  typedef std::vector<int, std::allocator<int> > rank_set_type;
  
  static const size_t buffer_size     = 1024 * 1024 * 4;
  static const size_t buffer_size_max = buffer_size << 2;
  static const size_t queue_size      = 1024;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
    
  path_set_type mapped_files;
  for (size_t i = 0; i != counts_files.size(); ++ i)
    if (static_cast<int>(i % mpi_size) == mpi_rank)
      mapped_files.push_back(counts_files[i]);
  
  ostream_ptr_set_type stream(mpi_size);
  odevice_ptr_set_type device(mpi_size);
  queue_ptr_set_type   queues(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  rank_set_type ranks(mpi_size);

  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new odevice_type(reducer.comm, rank, reverse_tag, buffer_size, false, true));
    
    stream[rank].reset(new ostream_type());
    //stream[rank]->push(boost::iostreams::zlib_compressor());
    stream[rank]->push(codec::lz4_compressor());
    stream[rank]->push(*device[rank], buffer_size);
    stream[rank]->precision(20);
    
    queues[rank].reset(new queue_type(queue_size));
    
    ranks[rank] = rank;
  }

  if (debug >= 2)
    std::cerr << "reverse counts: rank: " << mpi_rank << " files: " << mapped_files.size() << std::endl;

  boost::thread_group mapper;
  if (mpi_rank == 0 && debug)
    mapper.add_thread(new boost::thread(progress_mapper<mapper_type>(mapper_type(mapped_files, queues, max_malloc, debug))));
  else
    mapper.add_thread(new boost::thread(mapper_type(mapped_files, queues, max_malloc, debug)));
  
  simple_type           reversed;
  simple_generator_type generator;

  const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);

  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    rank_set_type::const_iterator riter_end = ranks.end();
    for (rank_set_type::const_iterator riter = ranks.begin(); riter != riter_end; ++ riter) {
      const int rank = *riter;
      
      if (stream[rank] && device[rank]) {
	
	if (device[rank]->test() && device[rank]->flush(true))
	  found = true;

	const size_t committed = device[rank]->committed();
	
	if (committed < buffer_size) {
	  for (int i = 0; i != 128 && static_cast<size_t>(device[rank]->committed()) < buffer_size && queues[rank]->pop_swap(reversed, true); ++ i) {
	    found = true;
	    
	    if (! reversed.source.empty())
	      generator(*stream[rank], reversed) << '\n';
	    else {
	      stream[rank].reset();
	      break;
	    }
	  }
	} else if (committed < buffer_size_max || utils::malloc_stats::used() < malloc_threshold) {
	  if (queues[rank]->pop_swap(reversed, true)) {
	    found |= (committed < buffer_size_max);
	    
	    if (! reversed.source.empty())
	      generator(*stream[rank], reversed) << '\n';
	    else
	      stream[rank].reset();
	  }
	}
      }
    }
    
    if (found)
      std::random_shuffle(ranks.begin(), ranks.end(), rgen);
    
    found |= utils::mpi_terminate_devices(stream, device);
    
    if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  mapper.join_all();
  
  // receive reversed files from the reducer sharing the same rank
  boost::iostreams::filtering_istream is;
  is.push(utils::mpi_device_source(reducer.comm, mpi_rank, file_tag, 4096));
  
  std::string line;
  while (std::getline(is, line)) {
    const path_type path(line);
    
    if (! boost::filesystem::exists(path))
      throw std::runtime_error("no reversed counts? " + line);
    
    reversed_files.push_back(path);
    utils::tempfile::insert(path);
    
    const path_type tmp_file = path.parent_path() / path.stem();
    if (boost::filesystem::exists(tmp_file))
      utils::tempfile::insert(tmp_file);
  }
}

void reverse_counts_reducer(utils::mpi_intercomm& mapper,
			    const path_type& output_file)
{
  typedef utils::repository repository_type;

  typedef boost::iostreams::filtering_istream istream_type;
  typedef utils::mpi_device_source            idevice_type;
  
  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  typedef boost::shared_ptr<idevice_type> idevice_ptr_type;
  
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  typedef std::vector<idevice_ptr_type, std::allocator<idevice_ptr_type> > idevice_ptr_set_type;
  
  typedef PhrasePairReverse map_reduce_type;
  
  typedef PhrasePairReverseReducer reducer_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::simple_type     simple_type;
  typedef map_reduce_type::simple_set_type simple_set_type;
  
  typedef PhrasePairSimpleParser simple_parser_type;

  typedef std::vector<int, std::allocator<int> > rank_set_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  const size_t queue_size  = 1024 * 16;
  const size_t buffer_size = 1024 * 1024 * 4;
  
  path_set_type reversed_files;
  
  istream_ptr_set_type stream(mpi_size);
  idevice_ptr_set_type device(mpi_size);
  
  boost::mt19937 gen;
  gen.seed(utils::random_seed());
  boost::random_number_generator<boost::mt19937> rgen(gen);
  
  rank_set_type ranks(mpi_size);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    device[rank].reset(new idevice_type(mapper.comm, rank, reverse_tag, buffer_size));
    
    stream[rank].reset(new istream_type());
    //stream[rank]->push(boost::iostreams::zlib_decompressor());
    stream[rank]->push(codec::lz4_decompressor());
    stream[rank]->push(*device[rank]);
    
    ranks[rank] = rank;
  }
  
  queue_type queue(queue_size);
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(queue, output_file, reversed_files, 1, max_malloc, debug)));
  
  simple_type     reversed;
  
  simple_parser_type parser;
  std::string line;
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    rank_set_type::const_iterator riter_end = ranks.end();
    for (rank_set_type::const_iterator riter = ranks.begin(); riter != riter_end; ++ riter) {
      const int rank = *riter;
      
      for (int i = 0; i != 128 && stream[rank] && device[rank] && device[rank]->test() && queue.size() < queue_size; ++ i) {
	if (std::getline(*stream[rank], line)) {
	  if (parser(line, reversed))
	    queue.push_swap(reversed);
	  else
	    std::cerr << "failed reversed phrase parsing: " << line << std::endl;
	} else {
	  stream[rank].reset();
	  device[rank].reset();
	}

	found = true;
      }
    }
    
    if (found)
      std::random_shuffle(ranks.begin(), ranks.end(), rgen);
    
    if (std::count(device.begin(), device.end(), idevice_ptr_type()) == mpi_size)
      break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  reversed.clear();
  queue.push_swap(reversed);
  
  reducer.join_all();
  
  {
    // send reversed files to mapper sharing the same rank
    boost::iostreams::filtering_ostream os;
    os.push(utils::mpi_device_sink(mapper.comm, mpi_rank, file_tag, 4096));
    
    for (path_set_type::const_iterator fiter = reversed_files.begin(); fiter != reversed_files.end(); ++ fiter)
      os << fiter->string() << '\n';
  }
  
  for (path_set_type::const_iterator fiter = reversed_files.begin(); fiter != reversed_files.end(); ++ fiter) {
    utils::tempfile::erase(*fiter);
    
    const path_type tmp_file = fiter->parent_path() / fiter->stem();
    if (boost::filesystem::exists(tmp_file))
      utils::tempfile::erase(tmp_file);
  }
}


void options(int argc, char** argv)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",                  po::value<path_set_type>(&input_files)->multitoken(), "input files")
    ("output",                 po::value<path_type>(&output_file),                   "output directory")
    ("temporary",              po::value<path_type>(&temporary_dir),                 "temporary directory")
    
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

