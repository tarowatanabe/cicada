//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// an implementation for neural network ngram language model
//
// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by Pegasos
// SGD with L2 regularizer inspired by AdaGrad (default)
// SGD with L2/L2 regularizer from RDA (TODO)
//
// an implementation of NCE estimate for ngram LM
//

#include <cstdlib>
#include <cmath>
#include <climits>

#include <deque>
#include <memory>

#include "cicada_nn_ngram_impl.hpp"

#include "utils/piece.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/unordered_map.hpp"
#include "utils/program_options.hpp"
#include "utils/random_seed.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"
#include "utils/mpi_traits.hpp"

#include "codec/lz4.hpp"

#include <boost/algorithm/string/trim.hpp>
#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/progress.hpp>

typedef boost::filesystem::path path_type;

typedef Data data_type;

typedef Model    model_type;
typedef Unigram  unigram_type;

typedef uint64_t count_type;
typedef utils::unordered_map<unigram_type::word_type, count_type,
			     boost::hash<unigram_type::word_type>, std::equal_to<unigram_type::word_type>,
			     std::allocator<std::pair<const unigram_type::word_type, count_type> > >::type word_set_type;

path_type input_file;
path_type list_file;
path_type embedding_file;
path_type output_model_file;

int dimension_embedding = 32;
int dimension_hidden = 256;
int order = 5;

bool optimize_sgd = false;
bool optimize_adagrad = false;

bool mix_simple = false;
bool mix_average = false;

int iteration = 10;
int batch_size = 64;
int samples = 100;
int cutoff = 3;
double lambda = 0;
double eta0 = 0.1;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const data_type& data,
		  const unigram_type& unigram,
		  model_type& theta);
void bcast_model(model_type& theta);
void read_data(const path_type& input_file,
	       const path_type& list_file,
	       data_type& data,
	       word_set_type& words);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  try {
    options(argc, argv);

    if (dimension_embedding <= 0)
      throw std::runtime_error("dimension must be positive");
    if (dimension_hidden <= 0)
      throw std::runtime_error("dimension must be positive");
    if (order <= 1)
      throw std::runtime_error("order size should be positive");

    if (samples <= 0)
      throw std::runtime_error("invalid sample size");
    if (batch_size <= 0)
      throw std::runtime_error("invalid batch size");
        
    if (int(optimize_sgd) + optimize_adagrad > 1)
      throw std::runtime_error("either one of optimize-{sgd,adagrad}");
    
    if (int(optimize_sgd) + optimize_adagrad == 0)
      optimize_sgd = true;
    
    if (int(mix_simple) + mix_average > 1)
      throw std::runtime_error("either one of mix-{simple,average}");
    
    // srand is used in Eigen
    std::srand(utils::random_seed());
  
    // this is optional, but safe to set this
    ::srandom(utils::random_seed());
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());

    if (input_file.empty() && list_file.empty())
      throw std::runtime_error("no data?");
    
    if (! input_file.empty())
      if (input_file != "-" && ! boost::filesystem::exists(input_file))
	throw std::runtime_error("no input file? " + input_file.string());
    
    if (! list_file.empty())
      if (list_file != "-" && ! boost::filesystem::exists(list_file))
	throw std::runtime_error("no list file? " + list_file.string());
    
    data_type     data(order);
    word_set_type words;
    
    read_data(input_file, list_file, data, words);

    if (debug && mpi_rank == 0)
      std::cerr << "# of ngrams: " << data.size() << std::endl
		<< "vocabulary: " << (words.size() - 1) << std::endl;
    
    unigram_type unigram(words.begin(), words.end());
    
    model_type theta(dimension_embedding, dimension_hidden, order, unigram, generator);

    bcast_model(theta);
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, order, lambda, eta0), data, unigram, theta);
      else
	learn_online(LearnSGD(lambda, eta0), data, unigram, theta);
    }
    
    if (mpi_rank == 0 && ! output_model_file.empty())
      theta.write(output_model_file);
    
  } catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  
  return 0;
}

enum {
  sentence_tag = 1000,
  ngram_tag,
  word_count_tag,
  model_tag,
  log_likelihood_tag,
  file_tag,
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

template <typename Learner>
struct TaskAccumulate
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef NGram ngram_type;

  typedef ngram_type::log_likelihood_type log_likelihood_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef std::string encoded_type;
  
  typedef utils::lockfree_list_queue<encoded_type, std::allocator<encoded_type> > queue_type;
  
  TaskAccumulate(const Learner& learner,
		 const data_type& data,
		 const unigram_type& unigram,
		 const size_type& samples,
		 model_type& theta,
		 queue_type& mapper,
		 queue_type& reducer,
		 size_type batch_size)
    : learner_(learner),
      data_(data),
      theta_(theta),
      mapper_(mapper),
      reducer_(reducer),
      ngram_(unigram, samples),
      gradient_(theta.dimension_embedding_, theta.dimension_hidden_, theta.order_),
      log_likelihood_(),
      batch_size_(batch_size)
  {
    generator_.seed(utils::random_seed());
  }
  
  void operator()()
  {
    clear();
    
    size_type batch = 0;

    encoded_type buffer;
    
    bool merge_finished = false;
    bool learn_finished = batch != data_.size();
    
    int non_found_iter = 0;
    
    while (! merge_finished || ! learn_finished) {
      bool found = false;
      
      if (! merge_finished)
	while (reducer_.pop(buffer, true)) {

	  if (buffer.empty())
	    merge_finished = true;
	  else {
	    gradient_.decode(buffer);
	    
	    learner_(theta_, gradient_);
	  }
	  
	  found = true;
	}
      
      if (! learn_finished) {
	found = true;
	
	gradient_.clear();
	
	const size_type last = utils::bithack::min(batch + batch_size_, data_.size());
	for (/**/; batch != last; ++ batch)
	  log_likelihood_ += ngram_.learn(data_.begin(batch), data_.end(batch), theta_, gradient_, generator_);
	
	learn_finished = (batch == data_.size());
	
	learner_(theta_, gradient_);
	
	gradient_.encode(buffer);
	
	mapper_.push_swap(buffer);
      }
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    // rescale current model
    theta_.rescale();
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
  
  void clear()
  {
    log_likelihood_ = log_likelihood_type();
  }
  
  Learner          learner_;
  const data_type& data_;
  model_type&      theta_;
  queue_type&      mapper_;
  queue_type&      reducer_;
  
  ngram_type ngram_;
  
  gradient_type       gradient_;
  log_likelihood_type log_likelihood_;
  
  size_type batch_size_;

  boost::mt19937 generator_;
};

inline
path_type add_suffix(const path_type& path, const std::string& suffix)
{
  bool has_suffix_gz  = false;
  bool has_suffix_bz2 = false;
  
  path_type path_added = path;
  
  if (path.extension() == ".gz") {
    path_added = path.parent_path() / path.stem();
    has_suffix_gz = true;
  } else if (path.extension() == ".bz2") {
    path_added = path.parent_path() / path.stem();
    has_suffix_bz2 = true;
  }
  
  path_added = path_added.string() + suffix;
  
  if (has_suffix_gz)
    path_added = path_added.string() + ".gz";
  else if (has_suffix_bz2)
    path_added = path_added.string() + ".bz2";
  
  return path_added;
}

template <typename Learner>
void learn_online(const Learner& learner,
		  const data_type& data,
		  const unigram_type& unigram,
		  model_type& theta)
{
  typedef TaskAccumulate<Learner> task_type;
  
  typedef typename task_type::size_type           size_type;
  typedef typename task_type::log_likelihood_type log_likelihood_type;
  typedef typename task_type::encoded_type        buffer_type;
  
  typedef typename task_type::queue_type queue_type;
  
  typedef boost::shared_ptr<buffer_type> buffer_ptr_type;
  typedef std::deque<buffer_ptr_type, std::allocator<buffer_ptr_type> >  buffer_set_type;
  typedef std::vector<buffer_set_type, std::allocator<buffer_set_type> > buffer_map_type;

  typedef boost::shared_ptr<utils::mpi_ostream_simple> ostream_ptr_type;
  typedef boost::shared_ptr<utils::mpi_istream_simple> istream_ptr_type;
  
  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  buffer_type          buffer;
  buffer_map_type      buffers(mpi_size);
  ostream_ptr_set_type ostreams(mpi_size);
  istream_ptr_set_type istreams(mpi_size);
  
  queue_type mapper;
  queue_type reducer;
  
  task_type task(learner,
		 data,
		 unigram,
		 samples,
		 theta,
		 mapper,
		 reducer,
		 batch_size);
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug && mpi_rank == 0)
      std::cerr << "iteration: " << (t + 1) << std::endl;
    
    utils::resource start;
    
    // create thread!
    boost::thread worker(boost::ref(task));
    
    bool finished = false;
    
    int non_found_iter = 0;
    for (;;) {
      bool found = false;
      
      // reduce samples...
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && istreams[rank] && istreams[rank]->test()) {
	  if (istreams[rank]->read(buffer))
	    reducer.push_swap(buffer);
	  else
	    istreams[rank].reset();
	  
	  buffer.clear();
	  found = true;
	}

      // check termination...
      if (! finished && std::count(istreams.begin(), istreams.end(), istream_ptr_type()) == mpi_size) {
	reducer.push(buffer_type());
	finished = true;
      }
      
      // bcast...
      // first, get the encoded buffer from mapper
      if (mapper.pop_swap(buffer, true)) {
	buffer_ptr_type buffer_ptr;
	
	if (! buffer.empty()) {
	  buffer_ptr.reset(new buffer_type());
	  buffer_ptr->swap(buffer);
	  buffer.clear();
	}
	
	for (int rank = 0; rank != mpi_size; ++ rank) 
	  if (rank != mpi_rank)
	    buffers[rank].push_back(buffer_ptr);
	
	found = true;
      }
      
      // second, bcast...
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && ostreams[rank] && ostreams[rank]->test() && ! buffers[rank].empty()) {
	  if (! buffers[rank].front()) {
	    // termination!
	    if (! ostreams[rank]->terminated())
	      ostreams[rank]->terminate();
	    else {
	      ostreams[rank].reset();
	      buffers[rank].erase(buffers[rank].begin());
	    }
	  } else {
	    ostreams[rank]->write(*(buffers[rank].front()));
	    buffers[rank].erase(buffers[rank].begin());
	  }
	  
	  found = true;
	}
      
      // termination condition
      if (finished
	  && std::count(istreams.begin(), istreams.end(), istream_ptr_type()) == mpi_size
	  && std::count(ostreams.begin(), ostreams.end(), ostream_ptr_type()) == mpi_size) break;
      
      // a conventional loop...
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    worker.join();
    
    utils::resource end;
    
    // merge log-likelihood...
    log_likelihood_type log_likelihood = task.log_likelihood_;
    
    if (mpi_rank == 0) {
      for (int rank = 1; rank != mpi_size; ++ rank) {
	log_likelihood_type ll;
	
	boost::iostreams::filtering_istream is;
	is.push(utils::mpi_device_source(rank, log_likelihood_tag, 4096));
	is.read((char*) &ll, sizeof(log_likelihood_type));
	
	log_likelihood += ll;
      }
    } else {
      boost::iostreams::filtering_ostream os;
      os.push(utils::mpi_device_sink(0, log_likelihood_tag, 4096));
      os.write((char*) &task.log_likelihood_, sizeof(log_likelihood_type));
    }
    
    if (debug && mpi_rank == 0)
      std::cerr << "log-likelihood: " << static_cast<double>(log_likelihood) << std::endl
		<< "perplexity: " << std::exp(- static_cast<double>(log_likelihood)) << std::endl;
    
    if (debug && mpi_rank == 0)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
  }
  
  // finalize model...
  theta.finalize();
}

struct Reader
{
  typedef cicada::Sentence sentence_type;

  Reader(data_type& data, word_set_type& words) 
    : data_(data), words_(words) {}

  void operator()(const std::string& line)
  {
    sentence_.assign(line);
    
    operator()(sentence_);
  }
  
  void operator()(const sentence_type& sentence)
  {
    typedef cicada::Vocab vocab_type;
    
    if (sentence.empty()) return;
    
    data_.insert(sentence);
    
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter)
      ++ words_[*siter];
    
    ++ words_[vocab_type::EOS];
  }

  sentence_type sentence_;
  
  data_type&     data_;
  word_set_type& words_;
};

struct ReaderFile : public Reader
{
  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
  
  ReaderFile(queue_type& queue, data_type& data, word_set_type& words)
    : Reader(data, words), queue_(queue) {}
  
  void operator()()
  {
    sentence_type sentence;
    path_type path;
    
    for (;;) {
      queue_.pop_swap(path);
      
      if (path.empty()) break;
      
      utils::compress_istream is(path, 1024 * 1024);
      
      while (is >> sentence)
	Reader::operator()(sentence);
    }
  }
  
  queue_type& queue_;
};

void bcast_model(model_type& theta)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    boost::iostreams::filtering_ostream os;
    os.push(codec::lz4_compressor());
    os.push(utils::mpi_device_bcast_sink(0, 1024 * 1024));
    
    os << theta;
  } else {
    boost::iostreams::filtering_istream is;
    is.push(codec::lz4_decompressor());
    is.push(utils::mpi_device_bcast_source(0, 1024 * 1024));

    is >> theta;
  }
}

void read_data(const path_type& input_file,
	       const path_type& list_file,
	       data_type& data,
	       word_set_type& words)
{
  typedef cicada::Vocab vocab_type;
  typedef size_t size_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  data.clear();
  words.clear();
  
  // read data from a file
  if (! input_file.empty()) {
    
    if (input_file != "-" && ! boost::filesystem::exists(input_file))
      throw std::runtime_error("no input file? " + input_file.string());
    
    Reader reader(data, words);
    
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
	device[rank].reset(new odevice_type(rank, sentence_tag, 1024 * 1024, false, true));
	
	stream[rank]->push(codec::lz4_compressor());
	stream[rank]->push(*device[rank]);
      }
      
      std::string line;
      utils::compress_istream is(input_file, 1024 * 1024);
      
      while (is) {
	for (int rank = 1; rank != mpi_size && is; ++ rank) {
	  if (! std::getline(is, line)) break;
	  
	  *stream[rank] << line << '\n';
	}
	
	if (std::getline(is, line))
	  reader(line);
      }
      
      // asynchronous termination...
      int non_found_iter = 0;
      for (;;) {
	bool found = false;
	
	// flush streams...
	for (int rank = 1; rank != mpi_size; ++ rank)
	  if (stream[rank] && device[rank]->test() && device[rank]->flush(true) == 0) {
	    stream[rank].reset();
	    found = true;
	  }
	
	found |= utils::mpi_terminate_devices(stream, device);
	
	if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
    } else {
      boost::iostreams::filtering_istream is;
      is.push(codec::lz4_decompressor());
      is.push(utils::mpi_device_source(0, sentence_tag, 1024 * 1024));
      
      std::string line;
      
      while (std::getline(is, line))
	reader(line);
    }
  }
  
  // read data from list
  if (! list_file.empty()) {
    if (list_file != "-" && ! boost::filesystem::exists(list_file))
	throw std::runtime_error("no list file? " + list_file.string());
    
    ReaderFile::queue_type queue(1);
    boost::thread worker(ReaderFile(queue, data, words));
    
    if (mpi_rank == 0) {
      typedef utils::mpi_ostream ostream_type;
      
      typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
      
      typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
      
      path_set_type paths;

      // first, we will read all the file names
      
      std::string line;
      utils::compress_istream is(list_file, 1024 * 1024);
      
      while (std::getline(is, line)) {
	boost::algorithm::trim(line);
	
	if (line.empty()) continue;
	
	if (boost::filesystem::exists(line))
	  paths.push_back(line);
	else if (boost::filesystem::exists(list_file.parent_path() / line))
	  paths.push_back(list_file.parent_path() / line);
	else
	  throw std::runtime_error(std::string("no file? ") + line);
      }
      
      std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > stream(mpi_size);
      for (int rank = 1; rank < mpi_size; ++ rank)
	stream[rank].reset(new ostream_type(rank, file_tag, 4096));
      
      int non_found_iter = 0;
      path_set_type::const_iterator piter_end = paths.end();
      path_set_type::const_iterator piter     = paths.begin();

      while (piter != piter_end) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size && piter != piter_end; ++ rank) 
	  if (stream[rank]->test()) {
	    stream[rank]->write(piter->string());
	    ++ piter;

	    found = true;
	  }
	
	if (piter != piter_end && queue.empty()) {
	  queue.push(*piter);
	  ++ piter;
	  
	  found = true;
	}

	non_found_iter = loop_sleep(found, non_found_iter);
      }
      
      // termination
      while (1) {
	bool found = false;
	
	for (int rank = 1; rank < mpi_size; ++ rank) 
	  if (stream[rank] && stream[rank]->test()) {
	    if (! stream[rank]->terminated())
	      stream[rank]->terminate();
	    else
	      stream[rank].reset();
	    
	    found = true;
	  }
	
	if (std::count(stream.begin(), stream.end(), ostream_ptr_type()) == mpi_size)
	  break;
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
      
    } else {
      utils::mpi_istream is(0, file_tag, 4096, true);
      
      std::string file;
      
      while (is.read(file)) {
	queue.push(file);
	queue.wait_empty();
	is.ready();
      }
    }
    
    queue.push(path_type());
    worker.join();
  }
  
  // balancing...  it is very stupid, but probably easier to implement...
  size_type data_size = data.size();
  size_type data_size_min = data.size();
  MPI::COMM_WORLD.Allreduce(&data_size, &data_size_min, 1, utils::mpi_traits<size_type>::data_type(), MPI::MIN);
  
  for (int rank = 0; rank != mpi_size; ++ rank) {
    if (rank == mpi_rank) {
      typedef boost::iostreams::filtering_ostream ostream_type;
      typedef utils::mpi_device_sink              odevice_type;
      
      typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
      typedef boost::shared_ptr<odevice_type> odevice_ptr_type;
      
      typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
      typedef std::vector<odevice_ptr_type, std::allocator<odevice_ptr_type> > odevice_ptr_set_type;
      
      const size_type diff = data_size - data_size_min;
      const size_type send_size = diff / mpi_size;
      
      ostream_ptr_set_type stream(mpi_size);
      odevice_ptr_set_type device(mpi_size);
      
      size_type pos = data.size() - send_size * (mpi_size - 1);
      
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank) {
	  stream[rank].reset(new ostream_type());
	  device[rank].reset(new odevice_type(rank, ngram_tag, 1024 * 1024, false, true));
	  
	  stream[rank]->push(codec::lz4_compressor());
	  stream[rank]->push(*device[rank]);
	  
	  std::copy(data.begin(pos), data.end(pos + send_size),
		    std::ostream_iterator<data_type::word_type>(*stream[rank], " "));
	  
	  pos += send_size;
	}
      
      // asynchronous termination...
      int non_found_iter = 0;
      for (;;) {
	bool found = false;
	
	// flush streams...
	for (int rank = 0; rank != mpi_size; ++ rank)
	  if (stream[rank] && device[rank]->test() && device[rank]->flush(true) == 0) {
	    stream[rank].reset();
	    found = true;
	  }
	
	found |= utils::mpi_terminate_devices(stream, device);
	
	if (std::count(device.begin(), device.end(), odevice_ptr_type()) == mpi_size) break;
	
	non_found_iter = loop_sleep(found, non_found_iter);
      }
      
      data.erase(data.begin(data.size() - send_size * (mpi_size - 1)), data.end());
    } else {
      // receive ngrams from rank
      boost::iostreams::filtering_istream is;
      is.push(codec::lz4_decompressor());
      is.push(utils::mpi_device_source(rank, ngram_tag, 1024 * 1024));
      
      data_type::word_type word;
      while (is >> word)
	data.data_.push_back(word);
    }
  }
  
  // bcast word + counts to everybody...
  {
    const word_set_type words_rank(words);
    
    for (int rank = 0; rank != mpi_size; ++ rank) {
      if (rank == mpi_rank) {
	namespace karma = boost::spirit::karma;
	namespace standard = boost::spirit::standard;

	boost::iostreams::filtering_ostream os;
	os.push(codec::lz4_compressor());
	os.push(utils::mpi_device_bcast_sink(rank, 1024 * 1024));

	std::ostream_iterator<char> iter(os);

	karma::uint_generator<count_type> count;
	
	word_set_type::const_iterator witer_end = words_rank.end();
	for (word_set_type::const_iterator witer = words_rank.begin(); witer != witer_end; ++ witer)
	  karma::generate(iter,
			  standard::string << karma::lit(' ') << count << karma::lit('\n'),
			  witer->first, witer->second);
      } else {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	typedef boost::spirit::istream_iterator iter_type;
	typedef standard::blank_type blank_type;
	
	qi::rule<iter_type, std::string(), blank_type> parse_word;
	qi::uint_parser<count_type>                    parse_count;
	
	parse_word %= qi::lexeme[+(standard::char_ - standard::space)];
	
	boost::iostreams::filtering_istream is;
	is.push(codec::lz4_decompressor());
	is.push(utils::mpi_device_bcast_source(rank, 1024 * 1024));
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;

	std::string word;
	count_type  count;
	
	while (iter != iter_end) {
	  word.clear();
	  
	  if (! qi::phrase_parse(iter, iter_end,
				 parse_word >> parse_count >> (qi::eol | qi::eoi),
				 standard::blank, word, count))				 
	    if (iter != iter_end)
	      throw std::runtime_error("parsing failed");
	  
	  words[word] += count;
	}
      }
    }
  }
  
  if (cutoff > 1) {
    word_set_type words_new;
    count_type count_unk = 0;
    
    word_set_type::const_iterator witer_end = words.end();
    for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer)
      if (witer->second >= cutoff)
	words_new.insert(*witer);
      else
	count_unk += witer->second;
    
    words_new[vocab_type::UNK] = count_unk;
    words_new.swap(words);
    words_new.clear();
    
    // enumerate data and replace by UNK
    data_type::iterator diter_end = data.end();
    for (data_type::iterator diter = data.begin(); diter != diter_end; ++ diter)
      if (*diter != vocab_type::BOS
	  && *diter != vocab_type::EOS
	  && *diter != vocab_type::EPSILON
	  && words.find(*diter) == words.end())
	*diter = vocab_type::UNK;
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("input",    po::value<path_type>(&input_file),  "input file")
    ("list",    po::value<path_type>(&list_file),    "list file")
    
    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("order",     po::value<int>(&order)->default_value(order),         "context order size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD fixed rate optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")

    ("mix-simple",  po::bool_switch(&mix_simple),  "mixing by selection")
    ("mix-average", po::bool_switch(&mix_average), "mixing by averaging")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("samples",           po::value<int>(&samples)->default_value(samples),       "# of NCE samples")
    ("cutoff",            po::value<int>(&cutoff)->default_value(cutoff),         "cutoff count for vocabulary (<= 1 to keep all)")
    ("lambda",            po::value<double>(&lambda)->default_value(lambda),      "regularization constant")
    ("eta0",              po::value<double>(&eta0)->default_value(eta0),          "\\eta_0 for decay")

    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  desc_command.add(opts_command);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options] [operations]\n"
	      << opts_command << std::endl;
    exit(0);
  }
}
