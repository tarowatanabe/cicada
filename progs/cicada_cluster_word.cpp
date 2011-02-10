//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// word clustering...
//

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <deque>
#include <map>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/functional/hash.hpp>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>

#include <utils/alloc_vector.hpp>
#include <utils/compress_stream.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/simple_vector.hpp>
#include <utils/space_separator.hpp>
#include <utils/mathop.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>

#include <google/dense_hash_map>

typedef int64_t          count_type;
typedef int32_t          cluster_id_type;

typedef cicada::Symbol   word_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

typedef boost::filesystem::path path_type;

struct WordClassCount
{
  typedef std::pair<word_type, count_type> word_count_pair_type;
  typedef utils::simple_vector<word_count_pair_type, std::allocator<word_count_pair_type> > word_count_type;
  
  count_type      count;
  word_count_type words;
  word_type       word;
  cluster_id_type cluster;
  
  WordClassCount() : count(0), words(), word(), cluster(-1) {}
};
typedef WordClassCount word_class_count_type;

inline
bool operator<(const word_class_count_type& x, const word_class_count_type& y)
{
  return x.word < y.word;
}

template <typename Tp>
struct lessp
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return *x < *y;
  }
};

typedef std::vector<word_class_count_type, std::allocator<word_class_count_type > > word_class_count_set_type;

struct Cluster
{
  typedef google::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type> > word_count_type;
  
  word_count_type words;
  count_type count;
  count_type size;
  
  Cluster() : words(), count(0), size(0)
  {
    words.set_empty_key(word_type());
    words.set_deleted_key(word_type(word_type::id_type(-1)));
  }

  void clear()
  {
    words.clear();
    count = 0;
    size = 0;
  }
};
typedef Cluster cluster_type;

typedef std::vector<cluster_type, std::allocator<cluster_type> > cluster_set_type;

path_type input_file = "-";
path_type output_file = "-";

int num_cluster = 50;
int max_iteration = 20;

int num_thread = 2;

int debug = 0;

int getoptions(int argc, char** argv);

double likelihood_cluster(const cluster_set_type& clusters);
void read_corpus(const path_type& file,
		 word_class_count_set_type& words);
void initial_cluster(const word_class_count_set_type& words,
		     cluster_set_type& clusters);
template <typename Generator>
double cluster_words(const word_class_count_set_type& words,
		     cluster_set_type& clusters,
		     Generator& generator,
		     const bool randomize=false);
void dump_clusters(const path_type& file,
		   const word_class_count_set_type& words,
		   const cluster_set_type& clusters);

int main(int argc, char** argv)
{
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;
    
    boost::mt19937 gen;
    gen.seed(time(0) * getpid());
    boost::random_number_generator<boost::mt19937> generator(gen);
    
    word_class_count_set_type words;
    
    read_corpus(input_file, words);
    
    cluster_set_type clusters(num_cluster);
    
    initial_cluster(words, clusters);
    
    if (debug)
      std::cerr << "objective: " << likelihood_cluster(clusters) << std::endl;

    for (int i = 0; i < max_iteration; ++ i) {
      const double delta = cluster_words(words, clusters, generator, i != 0);
      
      if (debug)
	std::cerr << "iteration: " << i 
		  << " objective: " << likelihood_cluster(clusters) 
		  << " delta: " << delta
		  << std::endl;
      
      if (delta == 0.0) break;
    }
    
    
    dump_clusters(output_file, words, clusters);
  }  
  catch (std::exception& e) {
    std::cerr << "error: " << e.what() << std::endl;
    return -1;
  }
  return 0;
}

double likelihood_cluster(const cluster_set_type& clusters)
{
  double likelihood = 0.0;
  
  for (size_t i = 0; i < clusters.size(); ++ i) {
    const cluster_type& cluster = clusters[i];
    
    likelihood -= utils::mathop::log(double(cluster.count)) * cluster.count;
    
    cluster_type::word_count_type::const_iterator witer_end = cluster.words.end();
    for (cluster_type::word_count_type::const_iterator witer = cluster.words.begin(); witer != witer_end; ++ witer)
      likelihood += utils::mathop::log(double(witer->second)) * witer->second;
  }
  return likelihood;
}

double remove_word(const word_class_count_type& word_class, const cluster_type& cluster, cluster_type& cluster_local, const bool tentative=true)
{
  double delta = utils::mathop::log(double(cluster.count + cluster_local.count)) * (cluster.count + cluster_local.count);
  const count_type count_local_new = cluster_local.count - word_class.count;
  delta -= utils::mathop::log(double(cluster.count + count_local_new)) * (cluster.count + count_local_new);
  if (! tentative)
    cluster_local.count = count_local_new;
  
  // bigram
  word_class_count_type::word_count_type::const_iterator witer_end = word_class.words.end();
  for (word_class_count_type::word_count_type::const_iterator witer = word_class.words.begin(); witer != witer_end; ++ witer) {
    
    cluster_type::word_count_type::const_iterator citer = cluster.words.find(witer->first);
    const count_type word_count = (citer == cluster.words.end() ? count_type(0) : citer->second);
    
    if (tentative) {
      cluster_type::word_count_type::const_iterator citer_local = cluster_local.words.find(witer->first);
      count_type word_count_local(citer_local == cluster_local.words.end() ? count_type(0) : citer_local->second);
      
      if (word_count + word_count_local > 0)
	delta -= utils::mathop::log(double(word_count + word_count_local)) * (word_count + word_count_local);
      word_count_local -= witer->second;
      delta += utils::mathop::log(double(word_count + word_count_local)) * (word_count + word_count_local);
    } else {
      count_type& word_count_local = cluster_local.words[witer->first];
      
      if (word_count + word_count_local > 0)
	delta -= utils::mathop::log(double(word_count + word_count_local)) * (word_count + word_count_local);
      word_count_local -= witer->second;
      delta += utils::mathop::log(double(word_count + word_count_local)) * (word_count + word_count_local);
      
      if (word_count_local == 0)
	cluster_local.words.erase(witer->first);
    }
  }
  
  return delta;
}

double remove_word(const word_class_count_type& word_class, cluster_type& cluster, const bool tentative=true)
{
  // class->word mapping...
  double delta = utils::mathop::log(double(cluster.count)) * cluster.count;
  const count_type count_new = cluster.count - word_class.count;
  delta -= utils::mathop::log(double(count_new)) * count_new;
  if (! tentative)
    cluster.count = count_new;
  
  // bigram
  word_class_count_type::word_count_type::const_iterator witer_end = word_class.words.end();
  for (word_class_count_type::word_count_type::const_iterator witer = word_class.words.begin(); witer != witer_end; ++ witer) {

    if (tentative) {
      cluster_type::word_count_type::const_iterator citer = cluster.words.find(witer->first);
      count_type word_count = (citer == cluster.words.end() ? count_type(0) : citer->second);
      if (word_count > 0)
	delta -= utils::mathop::log(double(word_count)) * word_count;
      word_count -= witer->second;
      delta += utils::mathop::log(double(word_count)) * word_count;
    } else {
      count_type& word_count = cluster.words[witer->first];
      if (word_count > 0)
	delta -= utils::mathop::log(double(word_count)) * word_count;
      word_count -= witer->second;
      delta += utils::mathop::log(double(word_count)) * word_count;
      if (word_count == 0)
	cluster.words.erase(witer->first);
    }
  }

  if (! tentative)
    -- cluster.size;
  
  return delta;
}

double move_word(const word_class_count_type& word_class, const cluster_type& cluster, cluster_type& cluster_local, const bool tentative=true)
{
  double delta = (cluster.count + cluster_local.count) * utils::mathop::log(double(cluster.count + cluster_local.count));
  const count_type count_local_new = cluster_local.count + word_class.count;
  delta -= (count_local_new + cluster.count) * std::log(count_local_new + cluster.count);
  if (! tentative)
    cluster_local.count = count_local_new;
  
  // bigram
  word_class_count_type::word_count_type::const_iterator witer_end = word_class.words.end();
  for (word_class_count_type::word_count_type::const_iterator witer = word_class.words.begin(); witer != witer_end; ++ witer) {
    
    cluster_type::word_count_type::const_iterator citer = cluster.words.find(witer->first);
    const count_type word_count = (citer == cluster.words.end() ? count_type(0) : citer->second);
    
    if (tentative) {
      cluster_type::word_count_type::const_iterator citer_local = cluster_local.words.find(witer->first);
      count_type word_count_local = (citer_local == cluster_local.words.end() ? count_type(0) : citer_local->second);
      
      if (word_count + word_count_local > 0)
	delta -= (word_count + word_count_local) * utils::mathop::log(double(word_count + word_count_local));
      word_count_local += witer->second;
      delta += (word_count + word_count_local) * utils::mathop::log(double(word_count + word_count_local));
    } else {
      count_type& word_count_local = cluster_local.words[witer->first];
      
      if (word_count + word_count_local > 0)
	delta -= (word_count + word_count_local) * utils::mathop::log(double(word_count + word_count_local));
      word_count_local += witer->second;
      delta += (word_count + word_count_local) * utils::mathop::log(double(word_count + word_count_local));
    }
  }
  
  return delta;
}

double move_word(const word_class_count_type& word_class, const int id, cluster_type& cluster, const bool tentative=true)
{
  // class->word mapping...
  double delta = cluster.count * utils::mathop::log(double(cluster.count));
  const count_type count_new = cluster.count + word_class.count;
  delta -= count_new * std::log(count_new);
  if (! tentative)
    cluster.count = count_new;
  
  // bigram
  word_class_count_type::word_count_type::const_iterator witer_end = word_class.words.end();
  for (word_class_count_type::word_count_type::const_iterator witer = word_class.words.begin(); witer != witer_end; ++ witer) {
    
    if (tentative) {
      cluster_type::word_count_type::const_iterator citer = cluster.words.find(witer->first);
      count_type word_count = (citer == cluster.words.end() ? count_type(0) : citer->second);
      if (word_count > 0)
	delta -= word_count * utils::mathop::log(double(word_count));
      word_count += witer->second;
      delta += word_count * utils::mathop::log(double(word_count));
    } else {
      count_type& word_count = cluster.words[witer->first];
      if (word_count > 0)
	delta -= word_count * utils::mathop::log(double(word_count));
      word_count += witer->second;
      delta += word_count * utils::mathop::log(double(word_count));
    }
  }
  if (! tentative) {
    ++ cluster.size;
    const_cast<word_class_count_type&>(word_class).cluster = id;
  }
  
  return delta;
}

struct Task
{
  typedef utils::lockfree_list_queue<size_t, std::allocator<size_t> > queue_type;
  typedef std::vector<int, std::allocator<int> > move_set_type;
  
  queue_type& queue;
  const word_class_count_set_type& words;
  cluster_set_type& clusters;
  move_set_type& moves;
  
  Task(queue_type& _queue,
       const word_class_count_set_type& _words,
       cluster_set_type& _clusters,
       move_set_type& _moves)
    : queue(_queue), words(_words), clusters(_clusters), moves(_moves) {}
  
  void operator()() throw()
  {
    
    cluster_set_type clusters_local(clusters.size());
    
    while (1) {
      size_t pos(size_t(-1));
      queue.pop(pos);
      if (pos == size_t(-1)) break;
      
#if 0
      for (cluster_set_type::iterator citer = clusters_local.begin(); citer != clusters_local.end(); ++ citer)
	citer->clear();
#endif
      
      const word_class_count_type& word_class = words[pos];
      
      const int cluster_from = word_class.cluster;
      int cluster_max = cluster_from;
      double delta_max = 0;
      
      const double delta_remove = remove_word(word_class, clusters[cluster_from], clusters_local[cluster_from]);
      
      for (int cluster_to = 0; cluster_to != static_cast<int>(clusters.size()); ++ cluster_to)
	if (cluster_from != cluster_to) {
	  const double delta_move = move_word(word_class, clusters[cluster_to], clusters_local[cluster_to]);
	  
	  if (delta_remove + delta_move > delta_max) {
	    cluster_max = cluster_to;
	    delta_max = delta_remove + delta_move;
	  }
	}
      
      if (cluster_max != cluster_from) {
	remove_word(word_class, clusters[cluster_from], clusters_local[cluster_from], false);
	move_word(word_class, clusters[cluster_max], clusters_local[cluster_max], false);
	moves[pos] = cluster_max;
	
	if (debug >= 3)
	  std::cerr << "move word: " << word_class.word 
		    << " from: " << cluster_from
		    << " to: " << cluster_max
		    << " delta: " << delta_max
		    << std::endl;
      }
    }
  }
};

template <typename Generator>
double cluster_words(const word_class_count_set_type& words,
		     cluster_set_type& clusters,
		     Generator& generator,
		     const bool randomize)
{
  typedef uint32_t pos_type;
  typedef std::vector<pos_type, std::allocator<pos_type> > pos_set_type;
  
  pos_set_type positions(words.size());
  for (size_t pos = 0; pos < positions.size(); ++ pos)
    positions[pos] = pos;
  if (randomize)
    std::random_shuffle(positions.begin(), positions.end(), generator);
  
  if (num_thread > 1) {
    // evenly split data...
    typedef Task task_type;
    typedef task_type::queue_type queue_type;
    typedef task_type::move_set_type move_set_type;
    typedef boost::thread thread_type;
    
    queue_type queue(1024 * 1024);
    move_set_type moves(words.size(), -1);
    
    std::vector<thread_type*, std::allocator<thread_type*> > threads(num_thread);
    for (size_t i = 0; i != threads.size(); ++ i)
      threads[i] = new thread_type(task_type(queue, words, clusters, moves));
    
    for (size_t pos = 0; pos != words.size(); ++ pos)
      queue.push(positions[pos]);
    
    for (int i = 0; i < num_thread; ++ i)
      queue.push(size_t(-1));
    
    for (size_t i = 0; i != threads.size(); ++ i) {
      threads[i]->join();
      delete threads[i];
    }
    threads.clear();
    
    // move according to delta queue...
    double delta_sum = 0.0;
    std::vector<size_t, std::allocator<size_t> > moved(clusters.size());
    
    pos_set_type::const_iterator piter_end = positions.end();
    for (pos_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter)
      if (moves[*piter] >= 0) {
	const word_class_count_type& word_class = words[*piter];
	
	const int cluster_from = word_class.cluster;
	const int cluster_to = moves[*piter];
	
	const double delta = remove_word(word_class, clusters[word_class.cluster]) + move_word(word_class, cluster_to, clusters[cluster_to]);
	if (delta > 0.0) {
	  delta_sum += remove_word(word_class, clusters[word_class.cluster], false);
	  delta_sum += move_word(word_class, cluster_to, clusters[cluster_to], false);
	  ++ moved[cluster_to];
	}
      }
    
    if (debug >= 2)
      for (size_t cluster = 0; cluster != moved.size(); ++ cluster)
	std::cerr << "cluster: " << cluster
		  << " size: " << clusters[cluster].size
		  << " moved: " << moved[cluster] << std::endl;
    
    return delta_sum;
  } else {
    
    double delta_sum = 0.0;
    std::vector<size_t, std::allocator<size_t> > moved(clusters.size());
    
    pos_set_type::const_iterator piter_end = positions.end();
    for (pos_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
      const word_class_count_type& word_class = words[*piter];
      
      const int cluster_from = word_class.cluster;
      
      int cluster_max = cluster_from;
      double delta_max = 0;
      
      const double delta_remove = remove_word(word_class, clusters[cluster_from]);
      
      for (int cluster_to = 0; cluster_to < static_cast<int>(clusters.size()); ++ cluster_to)
	if (cluster_from != cluster_to) {
	  const double delta_move = move_word(word_class, cluster_to, clusters[cluster_to]);
	  
	  if (delta_remove + delta_move > delta_max) {
	    cluster_max = cluster_to;
	    delta_max = delta_remove + delta_move;
	  }
	}
      
      // move word...
      if (cluster_max != cluster_from) {
	remove_word(word_class, clusters[cluster_from], false);
	move_word(word_class, cluster_max, clusters[cluster_max], false);
	  
	delta_sum += delta_max;
	
	++ moved[cluster_max];
	
	if (debug >= 3)
	  std::cerr << "move word: " << word_class.word 
		    << " from: " << cluster_from
		    << " to: " << cluster_max
		    << " delta: " << delta_max
		    << std::endl;
      }
    }
    
    if (debug >= 2)
      for (size_t cluster = 0; cluster != moved.size(); ++ cluster)
	std::cerr << "cluster: " << cluster
		  << " size: " << clusters[cluster].size
		  << " moved: " << moved[cluster] << std::endl;
    
    
    return delta_sum;
  }
}

void dump_clusters(const path_type& file,
		   const word_class_count_set_type& words,
		   const cluster_set_type& clusters)
{
  typedef std::vector<std::string, std::allocator<std::string> > class_set_type;

  utils::compress_ostream os(file, 1024 * 1024);
  
  // sri-mode classes...
  class_set_type classes(clusters.size());
  for (size_t i = 0; i != clusters.size(); ++ i)
    classes[i] = "<class-" + boost::lexical_cast<std::string>(i+1) + ">";
  
  word_class_count_set_type::const_iterator witer_end = words.end();
  for (word_class_count_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer)
    os << classes[witer->cluster] << '\t' << witer->word << '\n';
}

template <typename Tp>
struct greater_countp
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->count > y->count;
  }
};

void initial_cluster(const word_class_count_set_type& words,
		     cluster_set_type& clusters)
{
  
  size_t num_word = 0;
  word_class_count_set_type::const_iterator witer_end = words.end();
  for (word_class_count_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer, ++ num_word) {
    //const int cluster_id = (num_word < clusters.size() ? num_word % clusters.size() : int(double(clusters.size()) * random() / (RAND_MAX + 1.0)));
    const int cluster_id = num_word % clusters.size();
    
    cluster_type& cluster = clusters[cluster_id];
    cluster.count += witer->count;
    ++ cluster.size;
    
    word_class_count_type::word_count_type::const_iterator citer_end = witer->words.end();
    for (word_class_count_type::word_count_type::const_iterator citer = witer->words.begin(); citer != citer_end; ++ citer)
      cluster.words[citer->first] += citer->second;
    
    const_cast<word_class_count_type&>(*witer).cluster = cluster_id;
  }
}

struct WordCount
{
  typedef google::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type> > word_count_type;
  
  word_count_type words;
  count_type count;
  
  WordCount() : words(), count(0) { words.set_empty_key(word_type()); }
};

void read_corpus(const path_type& file,
		 word_class_count_set_type& words)
{
  // we don't care EOS! since it is one-sided word-class...
  
  typedef WordCount word_count_type;
  typedef utils::alloc_vector<word_count_type, std::allocator<word_count_type> > word_count_set_type;
  typedef std::multimap<const word_count_type*, word_type, greater_countp<word_count_type>, std::allocator<std::pair<const word_count_type* const, word_type> > > sorted_type;

  
  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
  
  utils::compress_istream is(file, 1024 * 1024);
  
  std::string line;
  sentence_type sent;
  
  word_count_set_type word_counts;
    
  while (std::getline(is, line)) {
    utils::piece line_piece(line);
    tokenizer_type tokenizer(line_piece);
    
    sent.clear();
    sent.push_back(vocab_type::BOS);
    sent.insert(sent.end(), tokenizer.begin(), tokenizer.end());
    
    sentence_type::const_iterator siter_end = sent.end();
    for (sentence_type::const_iterator siter = sent.begin() + 1; siter != siter_end; ++ siter) {
      // we consider *(siter - 1) and *siter
      
      word_count_type& word_class = word_counts[siter->id()];
      
      ++ word_class.count;
      ++ word_class.words[(siter - 1)->id()];
    }
  }
  
  sorted_type sorted;
  word_count_set_type::const_iterator citer_begin = word_counts.begin();
  word_count_set_type::const_iterator citer_end = word_counts.end();
  for (word_count_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) 
    if (*citer)
      sorted.insert(std::make_pair(&(*(*citer)), word_type(word_type::id_type(citer - citer_begin))));

  words.clear();
  words.reserve(sorted.size());
  words.resize(sorted.size());
  
  word_class_count_set_type::iterator witer = words.begin();
  sorted_type::const_iterator siter_end = sorted.end();
  for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter, ++ witer) {
    witer->word = siter->second;
    witer->count = siter->first->count;
    witer->words = word_class_count_type::word_count_type(siter->first->words.begin(), siter->first->words.end());
    const_cast<word_count_type*>(siter->first)->words.clear();
  }
}

int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",  po::value<path_type>(&input_file),  "input file")
    ("output", po::value<path_type>(&output_file), "output file")
    
    ("cluster",   po::value<int>(&num_cluster),   "# of clusters")
    ("iteration", po::value<int>(&max_iteration), "# of iterations")
    
    ("threads", po::value<int>(&num_thread), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}
