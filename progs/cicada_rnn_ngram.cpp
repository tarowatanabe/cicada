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

#include "cicada_rnn_ngram_impl.hpp"

#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/unordered_map.hpp"
#include "utils/program_options.hpp"
#include "utils/random_seed.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"

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

int dimension_embedding = 64;
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

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const data_type& data,
		  const unigram_type& unigram,
		  model_type& theta);
void read_data(const path_type& input_file,
	       const path_type& list_file,
	       data_type& data,
	       word_set_type& words);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (dimension_embedding <= 0)
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

    if (int(mix_simple) + mix_average == 0)
      mix_simple = true;
    
    threads = utils::bithack::max(threads, 1);
    
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

    if (debug)
      std::cerr << "# of ngrams: " << data.size() << std::endl
		<< "vocabulary: " << (words.size() - 1) << std::endl;
    
    unigram_type unigram(words.begin(), words.end());
    
    model_type theta(dimension_embedding, order, unigram, generator);
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, order, lambda, eta0), data, unigram, theta);
      else
	learn_online(LearnSGD(lambda, eta0), data, unigram, theta);
    }
    
    if (! output_model_file.empty())
      theta.write(output_model_file);
    
  } catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  
  return 0;
}

// We perform parallelization inspired by
//
// @InProceedings{zhao-huang:2013:NAACL-HLT,
//   author    = {Zhao, Kai  and  Huang, Liang},
//   title     = {Minibatch and Parallelization for Online Large Margin Structured Learning},
//   booktitle = {Proceedings of the 2013 Conference of the North American Chapter of the Association for Computational Linguistics: Human Language Technologies},
//   month     = {June},
//   year      = {2013},
//   address   = {Atlanta, Georgia},
//   publisher = {Association for Computational Linguistics},
//   pages     = {370--379},
//   url       = {http://www.aclweb.org/anthology/N13-1038}
// }
//
// which is a strategy very similar to those used in pialign.
//
// Basically, we split data into mini batch, and compute gradient only over the minibatch
//

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
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_mapper_type;
  typedef utils::lockfree_list_queue<gradient_type*, std::allocator<gradient_type*> > queue_merger_type;
  typedef std::vector<queue_merger_type, std::allocator<queue_merger_type> > queue_merger_set_type;
  
  typedef std::deque<gradient_type, std::allocator<gradient_type> > gradient_set_type;
  
  TaskAccumulate(const Learner& learner,
		 const data_type& data,
		 const unigram_type& unigram,
		 const size_type& samples,
		 const model_type& theta,
		 queue_mapper_type& mapper,
		 queue_merger_set_type& mergers,
		 size_type batch_size)
    : learner_(learner),
      data_(data),
      theta_(theta),
      mapper_(mapper),
      mergers_(mergers),
      ngram_(unigram, samples),
      log_likelihood_(),
      shard_(0),
      batch_size_(batch_size)
  {
    generator_.seed(utils::random_seed());
  }
  
  void operator()()
  {
    clear();
    
    const size_type shard_size = mergers_.size();
    
    size_type batch = 0;
    gradient_type* grad = 0;
    
    size_type merge_finished = 0;
    bool learn_finished = false;
    
    int non_found_iter = 0;
    
    while (merge_finished != shard_size || ! learn_finished) {
      bool found = false;
      
      if (merge_finished != shard_size)
	while (mergers_[shard_].pop(grad, true)) {
	  if (! grad)
	    ++ merge_finished;
	  else {
	    learner_(theta_, *grad);
	    grad->increment();
	  }
	  
	  found = true;
	}
      
      if (! learn_finished && mapper_.pop(batch, true)) {
	
	found = true;
	
	if (batch == size_type(-1)) {
	  // send termination!
	  for (size_type i = 0; i != shard_size; ++ i)
	    mergers_[i].push(0);
	  
	  learn_finished = true;
	} else {
	  gradient_type* grad = 0;
	  
	  for (size_type j = 0; j != gradients_.size(); ++ j)
	    if (gradients_[j].shared() >= shard_size || gradients_[j].shared() == 0) {
	      if (! grad)
		grad = &gradients_[j];
	      else
		gradients_[j].clear();
	    }
	  
	  if (! grad) {
	    gradients_.push_back(gradient_type(theta_.dimension_, theta_.order_));
	    grad = &gradients_.back();
	  }
	  
	  grad->clear();
	  
	  // use of batch from batch * batch_size and (batch + 1) * batch_size
	  
	  const size_type first = batch * batch_size_;
	  const size_type last  = utils::bithack::min(first + batch_size_, data_.size());
	  
	  for (size_type id = first; id != last; ++ id)
	    log_likelihood_ += ngram_.learn(data_.begin(id), data_.end(id), theta_, *grad, generator_);
	  
	  learner_(theta_, *grad);
	  grad->increment();
	  
	  for (size_type i = 0; i != shard_size; ++ i)
	    if (i != shard_)
	      mergers_[i].push(grad);
	}
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
    //gradients_.clear();
    log_likelihood_ = log_likelihood_type();
  }
  
  Learner                learner_;
  const data_type&       data_;
  model_type             theta_;
  queue_mapper_type&     mapper_;
  queue_merger_set_type& mergers_;
  
  ngram_type ngram_;
  
  gradient_set_type   gradients_;
  log_likelihood_type log_likelihood_;

  int shard_;
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
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef typename task_type::size_type           size_type;
  typedef typename task_type::log_likelihood_type log_likelihood_type;

  typedef typename task_type::queue_mapper_type     queue_mapper_type;
  typedef typename task_type::queue_merger_set_type queue_merger_set_type;
  
  typedef std::vector<size_type, std::allocator<size_type> > batch_set_type;
  
  const size_type batches_size = (data.size() + batch_size - 1) / batch_size;
  
  batch_set_type batches(batches_size);
  for (size_type batch = 0; batch != batches_size; ++ batch)
    batches[batch] = batch;
  
  queue_mapper_type     mapper(threads);
  queue_merger_set_type mergers(threads);
  
  task_set_type tasks(threads, task_type(learner,
					 data,
					 unigram,
					 samples,
					 theta,
					 mapper,
					 mergers,
					 batch_size));
  
  // assign shard id
  for (size_type shard = 0; shard != tasks.size(); ++ shard)
    tasks[shard].shard_ = shard;
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;
    
    std::auto_ptr<boost::progress_display> progress(debug
						    ? new boost::progress_display(batches_size, std::cerr, "", "", "")
						    : 0);
    
    utils::resource start;
    
    boost::thread_group workers;
    
    for (size_type i = 0; i != tasks.size(); ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));

    typename batch_set_type::const_iterator biter_end = batches.end();
    for (typename batch_set_type::const_iterator biter = batches.begin(); biter != biter_end; ++ biter) {
      mapper.push(*biter);
      
      if (debug)
	++ (*progress);
    }

    // termination
    for (size_type i = 0; i != tasks.size(); ++ i)
      mapper.push(size_type(-1));
    
    workers.join_all();
    
    utils::resource end;
    
    log_likelihood_type log_likelihood;
    for (size_type i = 0; i != tasks.size(); ++ i)
      log_likelihood += tasks[i].log_likelihood_;
    
    if (debug)
      std::cerr << "log-likelihood: " << static_cast<double>(log_likelihood) << std::endl
		<< "perplexity: " << std::exp(- static_cast<double>(log_likelihood)) << std::endl;
    
    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle ngrams!
    std::random_shuffle(batches.begin(), batches.end());
    
    if (mix_average) {
      for (size_type i = 1; i != tasks.size(); ++ i)
	tasks.front().theta_ += tasks[i].theta_;
      
      tasks.front().theta_ *= (1.0 / tasks.size());
      
      for (size_type i = 1; i != tasks.size(); ++ i)
	tasks[i].theta_ = tasks.front().theta_;
    } else if (mix_simple)
      for (size_type i = 1; i != tasks.size(); ++ i)
	tasks[i].theta_ = tasks.front().theta_;
  }
  
  // copy model!
  theta = tasks.front().theta_;

  // finalize model...
  theta.finalize();
}

struct Reader
{
  typedef cicada::Sentence sentence_type;

  Reader(const int order)  : data_(order) {}

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
  
  data_type     data_;
  word_set_type words_;
};

struct ReaderFile : public Reader
{
  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
  
  ReaderFile(queue_type& queue, const int order) : Reader(order), queue_(queue) {}
  
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
      
struct ReaderLines : public Reader
{
  typedef std::vector<std::string, std::allocator<std::string> > line_set_type;
  typedef utils::lockfree_list_queue<line_set_type, std::allocator<line_set_type> > queue_type;
  
  ReaderLines(queue_type& queue, const int order) : Reader(order), queue_(queue) {}
  
  void operator()() 
  {
    line_set_type lines;
    sentence_type sentence;
    path_type path;
    
    for (;;) {
      queue_.pop_swap(lines);
      
      if (lines.empty()) break;
      
      line_set_type::const_iterator liter_end = lines.end();
      for (line_set_type::const_iterator liter = lines.begin(); liter != liter_end; ++ liter) {
	sentence.assign(*liter);
	
	Reader::operator()(sentence);
      }
    }
  }
  
  queue_type& queue_;
};

void read_data(const path_type& input_file,
	       const path_type& list_file,
	       data_type& data,
	       word_set_type& words)
{
  typedef cicada::Vocab vocab_type;

  data.clear();
  words.clear();
  
  if (! input_file.empty()) {
    if (input_file != "-" && ! boost::filesystem::exists(input_file))
      throw std::runtime_error("no input file? " + input_file.string());

    ReaderLines::queue_type queue;

    std::vector<ReaderLines> tasks(threads, ReaderLines(queue, order));
    
    boost::thread_group workers;
    for (size_t i = 0; i != tasks.size(); ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    std::string                line;
    ReaderLines::line_set_type lines;
    
    utils::compress_istream is(input_file, 1024 * 1024);
    
    while (std::getline(is, line)) {
      if (line.empty()) continue;
      
      lines.push_back(line);

      if (lines.size() == 1024) {
	queue.push_swap(lines);
	lines.clear();
      }
    }
    
    if (! lines.empty())
      queue.push_swap(lines);
    lines.clear();
    
    // termination
    for (size_t i = 0; i != tasks.size(); ++ i)
      queue.push(ReaderLines::line_set_type());
    
    workers.join_all();

    size_t data_size = data.size();
    for (size_t i = 0; i != tasks.size(); ++ i)
      data_size += tasks[i].data_.size();
    data.reserve(data_size);
    
    // join data...
    for (size_t i = 0; i != tasks.size(); ++ i) {
      data += tasks[i].data_;
      
      word_set_type::const_iterator witer_end = tasks[i].words_.end();
      for (word_set_type::const_iterator witer = tasks[i].words_.begin(); witer != witer_end; ++ witer)
	words[witer->first] += witer->second;
    }
  }
  
  if (! list_file.empty()) {
    if (list_file != "-" && ! boost::filesystem::exists(list_file))
	throw std::runtime_error("no list file? " + list_file.string());
    
    ReaderFile::queue_type queue;
    
    std::vector<ReaderFile> tasks(threads, ReaderFile(queue, order));
    
    boost::thread_group workers;
    for (size_t i = 0; i != tasks.size(); ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    std::string line;
    
    utils::compress_istream is(list_file, 1024 * 1024);
    
    while (std::getline(is, line)) {
      boost::algorithm::trim(line);
      
      if (line.empty()) continue;

      if (boost::filesystem::exists(line))
	queue.push(line);
      else if (boost::filesystem::exists(list_file.parent_path() / line))
	queue.push(list_file.parent_path() / line);
      else
	throw std::runtime_error(std::string("no file? ") + line);
    }
    
     // termination
    for (size_t i = 0; i != tasks.size(); ++ i)
      queue.push(path_type());
    
    workers.join_all();
    
    size_t data_size = data.size();
    for (size_t i = 0; i != tasks.size(); ++ i)
      data_size += tasks[i].data_.size();
    data.reserve(data_size);
    
    // join data...
    for (size_t i = 0; i != tasks.size(); ++ i) {
      data += tasks[i].data_;
      
      word_set_type::const_iterator witer_end = tasks[i].words_.end();
      for (word_set_type::const_iterator witer = tasks[i].words_.begin(); witer != witer_end; ++ witer)
	words[witer->first] += witer->second;
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
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension")
    ("order",               po::value<int>(&order)->default_value(order),                             "context order size")
    
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

    ("threads", po::value<int>(&threads), "# of threads")
    
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
