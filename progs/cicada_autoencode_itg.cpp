//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// 
// ITG-based autoencoder
//
// X -> f/e : [vec_f, vec_e]
// X -> [X1, X2] : W_s [vec_x1_f, vec_x2_f, vec_x1_e, vec_x2_e] + B_s
// X -> <X1, X2> : W_i [vec_x1_f, vec_x2_f, vec_x2_e, vec_x1_e] + B_i
//
// - Use ITG beam search of Saers et al., (2009)
// - Learning is performed by autoencoding: recovering representation of children vectors
//
// word vector format: word dim1 dim2 dim3 ...
// tensor: [[dim1, dim2, ...], [dim1, dim2, ...], ... ] (make it compatible with neuron package...?)
// 

// currently, we learn and output the last derivation in a format similar to pialign:
//
// < [ ((( word |||  word ))) ] <  ((( word ||| word ))) > >
//
// [ ] indicates straight and < > indicates inversion.

// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by Pegasos (default)
// SGD with L2 regularizer inspired by AdaGrad
// SGD with L2/L2 regularizer from RDA
//
// + batch algorithm using LBFGS
//

#include <cstdlib>
#include <cmath>
#include <climits>

#include <set>
#include <deque>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"
#include "cicada/alignment.hpp"
#include "cicada/bitext.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/mathop.hpp"
#include "utils/program_options.hpp"
#include "utils/random_seed.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"

#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/progress.hpp>

#include "cicada_autoencode_itg_impl.hpp"

typedef boost::filesystem::path path_type;

typedef cicada::Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;

typedef Model model_type;
typedef Dictionary dictionary_type;

path_type source_file;
path_type target_file;

path_type embedding_source_file;
path_type embedding_target_file;

path_type derivation_file;
path_type alignment_source_target_file;
path_type alignment_target_source_file;
path_type output_model_file;

double alpha = 0.99;
double beta = 0.01;
int dimension_embedding = 32;
int dimension_hidden = 128;
int span = 8;
int window = 0;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int baby_steps = 1;
int batch_size = 4;
int beam = 10;
double lambda = 0;
double eta0 = 0.1;
int cutoff = 3;

bool moses_mode = false;
bool giza_mode = false;

bool dump_mode = false;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta);
void derivation(const bitext_set_type& bitexts,
		const dictionary_type& dict_source_target,
		const dictionary_type& dict_target_source,
		const model_type& theta);
void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       dictionary_type& dict_source_target,
	       dictionary_type& dict_target_source);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (dimension_embedding <= 0)
      throw std::runtime_error("dimension must be positive");
    if (dimension_hidden <= 0)
      throw std::runtime_error("dimension must be positive");
    if (span < 0)
      throw std::runtime_error("span context size should be positive or zero");
    if (window < 0)
      throw std::runtime_error("window size should be positive");
    
    if (alpha < 0.0)
      throw std::runtime_error("alpha should be >= 0.0");
    if (beta < 0.0)
      throw std::runtime_error("beta should be >= 0.0");

    if (beam <= 0)
      throw std::runtime_error("beam width should be positive");
    
    if (int(giza_mode) + moses_mode > 1)
      throw std::runtime_error("either giza style output or moses style output");

    if (int(giza_mode) + moses_mode == 0)
      moses_mode = true;

    if (int(optimize_sgd) + optimize_adagrad > 1)
      throw std::runtime_error("either one of optimize-{sgd,adagrad}");
    
    if (int(optimize_sgd) + optimize_adagrad == 0)
      optimize_sgd = true;
    
    threads = utils::bithack::max(threads, 1);
    
    // srand is used in Eigen
    std::srand(utils::random_seed());
  
    // this is optional, but safe to set this
    ::srandom(utils::random_seed());

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
        
    if (source_file.empty())
      throw std::runtime_error("no source data?");
    if (target_file.empty())
      throw std::runtime_error("no target data?");
    
    bitext_set_type bitexts;
    
    dictionary_type dict_source_target;
    dictionary_type dict_target_source;

    read_data(source_file, target_file, bitexts, dict_source_target, dict_target_source);
    
    const dictionary_type::dict_type::word_set_type& sources = dict_target_source[cicada::Vocab::EPSILON].words_;
    const dictionary_type::dict_type::word_set_type& targets = dict_source_target[cicada::Vocab::EPSILON].words_;

    if (debug)
      std::cerr << "# of unique source words: " << sources.size() << std::endl
		<< "# of unique target words: " << targets.size() << std::endl
		<< "# of sentences: " << bitexts.size() << std::endl;
    
    model_type theta(dimension_embedding, dimension_hidden, span, window, sources, targets, generator);
    
    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());
      
      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta.read_embedding(embedding_source_file, embedding_target_file);
    }
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, span, window, lambda, eta0),
		     bitexts,
		     dict_source_target,
		     dict_target_source,
		     theta);
      else
	learn_online(LearnSGD(lambda, eta0),
		     bitexts,
		     dict_source_target,
		     dict_target_source,
		     theta);
    }
    
    if (! derivation_file.empty() || ! alignment_source_target_file.empty() || ! alignment_target_source_file.empty())
      derivation(bitexts, dict_source_target, dict_target_source, theta);
    
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

struct OutputMapReduce
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef boost::filesystem::path path_type;
  
  typedef cicada::Bitext bitext_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;

  typedef ITG itg_type;

  typedef itg_type::vocab_type vocab_type;
  
  typedef itg_type::span_type       span_type;
  typedef itg_type::span_pair_type  span_pair_type;
  typedef itg_type::hyperedge_type  hyperedge_type;
  typedef itg_type::derivation_type derivation_type;

  struct bitext_derivation_type
  {
    size_type       id_;
    derivation_type derivation_;
    
    bitext_derivation_type() : id_(size_type(-1)), derivation_() {}
    bitext_derivation_type(const size_type& id,
			   const derivation_type& derivation)
      : id_(id), derivation_(derivation) {}
    
    void swap(bitext_derivation_type& x)
    {
      std::swap(id_, x.id_);
      derivation_.swap(x.derivation_);
    }

    void clear()
    {
      id_ = size_type(-1);
      derivation_.clear();
    }
  };
  
  typedef bitext_derivation_type value_type;

  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;

  struct compare_value
  {
    bool operator()(const value_type& x, const value_type& y) const
    {
      return x.id_ < y.id_;
    }
  };
  typedef std::set<value_type, compare_value, std::allocator<value_type> > bitext_reduced_type;

};

namespace std
{
  inline
  void swap(OutputMapReduce::value_type& x,
	    OutputMapReduce::value_type& y)
  {
    x.swap(y);
  }
};

struct OutputDerivation : OutputMapReduce
{
  typedef std::vector<std::string, std::allocator<std::string> > stack_type;
	  
  OutputDerivation(const bitext_set_type& bitexts,
		   const path_type& path,
		   queue_type& queue)
    : bitexts_(bitexts), path_(path), queue_(queue) {}
  
  void operator()()
  {
    if (path_.empty()) {
      bitext_derivation_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_reduced_type bitexts;
      bitext_derivation_type bitext;
      size_type id = 0;
      
      utils::compress_ostream os(path_, 1024 * 1024);
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
	
	if (bitext.id_ == id) {
	  write(os, bitext);
	  ++ id;
	} else
	  bitexts.insert(bitext);
	
	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  write(os, *bitexts.begin());
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	write(os, *bitexts.begin());
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while writing derivation output?");
    }
  }
  
  void write(std::ostream& os, const value_type& bitext)
  {
    stack_.clear();
    
    derivation_type::const_iterator diter_end = bitext.derivation_.end();
    for (derivation_type::const_iterator diter = bitext.derivation_.begin(); diter != diter_end; ++ diter) {
      if (diter->terminal()) {
	const word_type& source = (! diter->span_.source_.empty()
				   ? bitexts_[bitext.id_].source_[diter->span_.source_.first_]
				   : vocab_type::EPSILON);
	const word_type& target = (! diter->span_.target_.empty()
				   ? bitexts_[bitext.id_].target_[diter->span_.target_.first_]
				   : vocab_type::EPSILON);
	
	os << "((( " << source << " ||| " << target << " )))";
	
	while (! stack_.empty() && stack_.back() != " ") {
	  os << stack_.back();
	  stack_.pop_back();
	}
	
	if (! stack_.empty() && stack_.back() == " ") {
	  os << stack_.back();
	  stack_.pop_back();
	}
      } else if (diter->straight()) {
	os << "[ ";
	stack_.push_back(" ]");
	stack_.push_back(" ");
      } else {
	os << "< ";
	stack_.push_back(" >");
	stack_.push_back(" ");
      }
    }
    
    os << '\n';
  }
  
  const bitext_set_type& bitexts_;
  path_type              path_;
  queue_type&            queue_;
  
  stack_type stack_;
};

struct OutputAlignment : OutputMapReduce
{
  typedef cicada::Alignment alignment_type;

  OutputAlignment(const bitext_set_type& bitexts,
		  const path_type& path_source_target,
		  const path_type& path_target_source,
		  queue_type& queue)
    : bitexts_(bitexts),
      path_source_target_(path_source_target),
      path_target_source_(path_target_source),
      queue_(queue) {}
  
  void operator()()
  {
    if (path_source_target_.empty() && path_target_source_.empty()) {
      bitext_derivation_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_reduced_type bitexts;
      bitext_derivation_type bitext;
      size_type id = 0;
      
      std::auto_ptr<std::ostream> os_source_target(! path_source_target_.empty()
						   ? new utils::compress_ostream(path_source_target_, 1024 * 1024)
						   : 0);
      std::auto_ptr<std::ostream> os_target_source(! path_target_source_.empty()
						   ? new utils::compress_ostream(path_target_source_, 1024 * 1024)
						   : 0);
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
	
	if (bitext.id_ == id) {
	  write(os_source_target.get(), os_target_source.get(), bitext);
	  ++ id;
	} else
	  bitexts.insert(bitext);
	
	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  write(os_source_target.get(), os_target_source.get(), *bitexts.begin());
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	write(os_source_target.get(), os_target_source.get(), *bitexts.begin());
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while writing derivation output?");
    }
  }

  void write(std::ostream* os_source_target, std::ostream* os_target_source, const value_type& bitext)
  {
    alignment_.clear();

    derivation_type::const_iterator diter_end = bitext.derivation_.end();
    for (derivation_type::const_iterator diter = bitext.derivation_.begin(); diter != diter_end; ++ diter)
      if (diter->terminal() && diter->aligned())
	alignment_.push_back(std::make_pair(diter->span_.source_.first_, diter->span_.target_.first_));
    
    if (os_source_target) {
      std::sort(alignment_.begin(), alignment_.end());

      if (moses_mode)
	*os_source_target << alignment_ << '\n';
      else
	output(*os_source_target, bitext.id_, bitexts_[bitext.id_].source_, bitexts_[bitext.id_].target_, alignment_);
    }
    
    if (os_target_source) {
      alignment_.inverse();
      
      std::sort(alignment_.begin(), alignment_.end());

      if (moses_mode)
	*os_target_source << alignment_ << '\n';
      else
	output(*os_target_source, bitext.id_, bitexts_[bitext.id_].target_, bitexts_[bitext.id_].source_, alignment_);
    }
  }

  typedef int index_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > align_set_type;
  typedef std::set<index_type, std::less<index_type>, std::allocator<index_type> > align_none_type;
  
  align_set_type  aligns_;
  align_none_type aligns_none_;

  void output(std::ostream& os,
	      const size_type& id,
	      const sentence_type& source,
	      const sentence_type& target,
	      const alignment_type& alignment)
  {
    os << "# Sentence pair (" << (id + 1) << ')'
       << " source length " << source.size()
       << " target length " << target.size()
       << " alignment score : " << 0 << '\n';
    os << target << '\n';
    
    if (source.empty() || target.empty()) {
      os << "NULL ({ })";
      sentence_type::const_iterator siter_end = source.end();
      for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
	os << ' ' << *siter << " ({ })";
      os << '\n';
    } else {
      aligns_.clear();
      aligns_.resize(source.size());
      
      aligns_none_.clear();
      for (size_type trg = 0; trg != target.size(); ++ trg)
	aligns_none_.insert(trg + 1);
      
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
	aligns_[aiter->source].push_back(aiter->target + 1);
	aligns_none_.erase(aiter->target + 1);
      }
      
      os << "NULL";
      os << " ({ ";
      std::copy(aligns_none_.begin(), aligns_none_.end(), std::ostream_iterator<index_type>(os, " "));
      os << "})";
      
      for (size_type src = 0; src != source.size(); ++ src) {
	os << ' ' << source[src];
	os << " ({ ";
	std::copy(aligns_[src].begin(), aligns_[src].end(), std::ostream_iterator<index_type>(os, " "));
	os << "})";
      }
      os << '\n';
    }
  }

  const bitext_set_type& bitexts_;
  path_type              path_source_target_;
  path_type              path_target_source_;
  queue_type&            queue_;
  
  alignment_type alignment_;
};


template <typename Learner>
struct TaskAccumulate
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef ITG itg_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  
  typedef itg_type::vocab_type vocab_type;
  
  typedef Average loss_type;
  
  typedef OutputMapReduce output_map_reduce_type;
  
  typedef output_map_reduce_type::bitext_derivation_type bitext_derivation_type;
  typedef output_map_reduce_type::queue_type             queue_derivation_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_mapper_type;
  typedef utils::lockfree_list_queue<gradient_type*, std::allocator<gradient_type*> > queue_merger_type;
  typedef std::vector<queue_merger_type, std::allocator<queue_merger_type> > queue_merger_set_type;
  
  typedef std::deque<gradient_type, std::allocator<gradient_type> > gradient_set_type;
  
  TaskAccumulate(const Learner& learner,
		 const bitext_set_type& bitexts,
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 const model_type& theta,
		 const int& beam,
		 const size_type batch_size,
		 queue_mapper_type& mapper,
		 queue_merger_set_type& mergers,
		 queue_derivation_type& queue_derivation,
		 queue_derivation_type& queue_alignment)
    : learner_(learner),
      bitexts_(bitexts),
      theta_(theta),
      mapper_(mapper),
      mergers_(mergers),
      queue_derivation_(queue_derivation),
      queue_alignment_(queue_alignment),
      itg_(dict_source_target, dict_target_source, beam),
      parsed_(0),
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
    
    bitext_derivation_type bitext_derivation;

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
	    if (gradients_[j].shared() == shard_size) {
	      grad = &gradients_[j];
	      break;
	    }
	  
	  if (! grad) {
	    gradients_.push_back(gradient_type(theta_.embedding_, theta_.hidden_, theta_.span_, theta_.window_));
	    grad = &gradients_.back();
	  }
	  
	  grad->clear();
	  
	  const size_type first = batch * batch_size_;
	  const size_type last  = utils::bithack::min(first + batch_size_, bitexts_.size());
	  
	  for (size_type id = first; id != last; ++ id) {
	    const sentence_type& source = bitexts_[id].source_;
	    const sentence_type& target = bitexts_[id].target_;
	    
	    bitext_derivation.id_ = id;
	    bitext_derivation.derivation_.clear();
	    
	    if (! source.empty() && ! target.empty()) {
#if 0
	      std::cerr << "source: " << source << std::endl
			<< "target: " << target << std::endl;
#endif
	      
	      const double error = itg_.forward(source, target, theta_);
	      
	      const bool parsed = (error != std::numeric_limits<double>::infinity());
	      
	      //std::cerr << "error: " << error << std::endl;
	      
	      if (parsed) {
		itg_.backward(source, target, theta_, *grad, generator_);
		
		loss_ += error;
		++ parsed_;
		
		itg_.derivation(source, target, bitext_derivation.derivation_);
	      } else
		std::cerr << "failed parsing: " << std::endl
			  << "source: " << source << std::endl
			  << "target: " << target << std::endl;
	    }
	    
	    queue_derivation_.push(bitext_derivation);
	    queue_alignment_.push(bitext_derivation);
	  }
	  
	  learner_(theta_, *grad);
	  grad->increment();
	  
	  for (size_type i = 0; i != shard_size; ++ i)
	    if (i != shard_)
	      mergers_[i].push(grad);
	}
      }
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    theta_.finalize();
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
    loss_ = loss_type();
    parsed_ = 0;
  }

  Learner                learner_;
  const bitext_set_type& bitexts_;
  model_type             theta_;

  queue_mapper_type&     mapper_;
  queue_merger_set_type& mergers_;  
  queue_derivation_type& queue_derivation_;
  queue_derivation_type& queue_alignment_;
  
  itg_type itg_;

  gradient_set_type gradients_;
  loss_type         loss_;
  size_type         parsed_;
  
  int            shard_;
  size_type      batch_size_;
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

template <typename Lengths>
struct less_lengths
{
  typedef typename Lengths::size_type size_type;
  
  less_lengths(const Lengths& lengths)
    : lengths_(lengths) {}

  bool operator()(const size_type& x, const size_type& y) const
  {
    return lengths_[x] < lengths_[y];
  }

  const Lengths& lengths_;
};

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta)
{
  typedef TaskAccumulate<Learner> task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef typename task_type::size_type size_type;

  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputDerivation output_derivation_type;
  typedef OutputAlignment  output_alignment_type;

  typedef typename task_type::queue_mapper_type     queue_mapper_type;
  typedef typename task_type::queue_merger_set_type queue_merger_set_type;

  typedef typename task_type::loss_type loss_type;

  typedef std::vector<size_type, std::allocator<size_type> > batch_set_type;

  const size_type batches_size = (bitexts.size() + batch_size - 1) / batch_size;
  
  batch_set_type batches(batches_size);
  batch_set_type lengths(batches_size);
  for (size_type batch = 0; batch != batches_size; ++ batch) {
    batches[batch] = batch;
    
    const size_type first = batch * batch_size;
    const size_type last  = utils::bithack::min(first + batch_size, bitexts.size());
    
    lengths[batch] = 0;
    for (size_type pos = first; pos != last; ++ pos)
      lengths[batch] += bitexts[pos].source_.size() + bitexts[pos].target_.size();
  }
  
  queue_mapper_type     mapper(threads);
  queue_merger_set_type mergers(threads);

  typename output_map_reduce_type::queue_type queue_derivation;
  typename output_map_reduce_type::queue_type queue_alignment;
  
  task_set_type tasks(threads, task_type(learner,
					 bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta,
					 beam,
					 batch_size,
					 mapper,
					 mergers,
					 queue_derivation,
					 queue_alignment));

  // assign shard id
  for (size_type shard = 0; shard != tasks.size(); ++ shard)
    tasks[shard].shard_ = shard;

  // iterations for baby-steps
  int baby_iter = 0;
  const int baby_last = utils::bithack::branch(baby_steps > 0, baby_steps, 0);
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;
    
    // baby-steps...
    bool baby_finished = true;
    if (baby_iter != baby_last) {
      ++ baby_iter;
      baby_finished = false;
    }

    if (! baby_finished) {
      // sort bitexts...
      typename batch_set_type::iterator biter     = batches.begin();
      typename batch_set_type::iterator biter_end = batches.end();
      
      while (biter < biter_end) {
	typename batch_set_type::iterator iter_end = std::min(biter + utils::bithack::max(4096 / batch_size, 1), biter_end);
	
	std::sort(biter, iter_end, less_lengths<batch_set_type>(lengths));
	biter = iter_end;
      }
    }

    std::auto_ptr<boost::progress_display> progress(debug
						    ? new boost::progress_display(batches_size, std::cerr, "", "", "")
						    : 0);
    
    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);

    boost::thread output_derivation(output_derivation_type(bitexts,
							   ! derivation_file.empty() && dump_mode
							   ? add_suffix(derivation_file, iter_tag)
							   : path_type(),
							   queue_derivation));
    boost::thread output_alignment(output_alignment_type(bitexts,
							 ! alignment_source_target_file.empty() && dump_mode
							 ? add_suffix(alignment_source_target_file, iter_tag)
							 : path_type(),
							 ! alignment_target_source_file.empty() && dump_mode
							 ? add_suffix(alignment_target_source_file, iter_tag)
							 : path_type(),
							 queue_alignment));

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
    
    queue_derivation.push(typename output_map_reduce_type::value_type());
    queue_alignment.push(typename output_map_reduce_type::value_type());
    
    utils::resource end;
    
    loss_type loss;
    size_type parsed = 0;
    
    for (size_type i = 0; i != tasks.size(); ++ i) {
      loss   += tasks[i].loss_;
      parsed += tasks[i].parsed_;
    }

    if (debug)
      std::cerr << "reconstruction error: " << static_cast<double>(loss) << std::endl
		<< "parsed: " << parsed << std::endl;
    
    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;

    // shuffle bitexts!
    {
      typename batch_set_type::iterator biter     = batches.begin();
      typename batch_set_type::iterator biter_end = batches.end();
      
      while (biter < biter_end) {
	typename batch_set_type::iterator iter_end = std::min(biter + utils::bithack::max(4096 / batch_size, 1), biter_end);
	
	std::random_shuffle(biter, iter_end);
	biter = iter_end;
      }
    }

    // mixing
    for (size_type i = 1; i != tasks.size(); ++ i)
      tasks[i].theta_ = tasks.front().theta_;
    
    output_derivation.join();
    output_alignment.join();
  }
  
  // copy model!
  theta = tasks.front().theta_;
}

struct TaskDerivation
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef ITG itg_type;
  
  typedef itg_type::vocab_type vocab_type;
  typedef itg_type::bitext_type bitext_type;

  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  
  typedef OutputMapReduce output_map_reduce_type;

  typedef output_map_reduce_type::bitext_derivation_type bitext_derivation_type;
  typedef output_map_reduce_type::queue_type queue_derivation_type;

  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_type;
  
  TaskDerivation(const bitext_set_type& bitexts,
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 const model_type& theta,
		 const int& beam,
		 queue_type& queue,
		 queue_derivation_type& queue_derivation,
		 queue_derivation_type& queue_alignment)
    : bitexts_(bitexts),
      theta_(theta),
      queue_(queue),
      queue_derivation_(queue_derivation),
      queue_alignment_(queue_alignment),
      itg_(dict_source_target, dict_target_source, beam) {}

  
  void operator()()
  {
    bitext_derivation_type bitext_derivation;
    
    size_type bitext_id;
    for (;;) {
      queue_.pop(bitext_id);
      
      if (bitext_id == size_type(-1)) break;
      
      const sentence_type& source = bitexts_[bitext_id].source_;
      const sentence_type& target = bitexts_[bitext_id].target_;
      
      bitext_derivation.id_ = bitext_id;
      bitext_derivation.derivation_.clear();
      
      if (! source.empty() && ! target.empty()) {

#if 0
	std::cerr << "source: " << source << std::endl
		  << "target: " << target << std::endl;
#endif
	
	const double error = itg_.forward(source, target, theta_);
	
	const bool parsed = (error != std::numeric_limits<double>::infinity());
	
	if (parsed)
	  itg_.derivation(source, target, bitext_derivation.derivation_);
	else
	  std::cerr << "failed parsing: " << std::endl
		    << "source: " << source << std::endl
		    << "target: " << target << std::endl;
      }
      
      queue_derivation_.push(bitext_derivation);
      queue_alignment_.push(bitext_derivation);
    }
  }

  const bitext_set_type& bitexts_;
  const model_type& theta_;
  
  queue_type&            queue_;
  queue_derivation_type& queue_derivation_;
  queue_derivation_type& queue_alignment_;
  
  itg_type itg_;
};

void derivation(const bitext_set_type& bitexts,
		const dictionary_type& dict_source_target,
		const dictionary_type& dict_target_source,
		const model_type& theta)
{
  typedef TaskDerivation task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  
  typedef task_type::size_type size_type;
  
  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputDerivation output_derivation_type;
  typedef OutputAlignment  output_alignment_type;

  task_type::queue_type   mapper(8 * threads);
  output_map_reduce_type::queue_type queue_derivation;
  output_map_reduce_type::queue_type queue_alignment;
  
  task_set_type tasks(threads, task_type(bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta,
					 beam,
					 mapper,
					 queue_derivation,
					 queue_alignment));

  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  boost::thread_group workers_dump;
  workers_dump.add_thread(new boost::thread(output_derivation_type(bitexts,
								   ! derivation_file.empty()
								   ? derivation_file
								   : path_type(),
								   queue_derivation)));
  workers_dump.add_thread(new boost::thread(output_alignment_type(bitexts,
								  ! alignment_source_target_file.empty()
								  ? alignment_source_target_file
								  : path_type(),
								  ! alignment_target_source_file.empty()
								  ? alignment_target_source_file
								  : path_type(),
								  queue_alignment)));

  if (debug)
    std::cerr << "max derivation" << std::endl;

  std::auto_ptr<boost::progress_display> progress(debug
						  ? new boost::progress_display(bitexts.size(), std::cerr, "", "", "")
						  : 0);

  utils::resource start;
  
  for (size_type i = 0; i != bitexts.size(); ++ i) {
    mapper.push(i);
    
    if (debug)
      ++ (*progress);
  }
  
  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    mapper.push(size_type(-1));
  
  workers.join_all();

  utils::resource end;

  if (debug)
    std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
	      << "user time:   " << end.user_time() - start.user_time() << std::endl;

  queue_derivation.push(output_map_reduce_type::value_type());
  queue_alignment.push(output_map_reduce_type::value_type());
  
  workers_dump.join_all();
}

void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       dictionary_type& dict_source_target,
	       dictionary_type& dict_target_source)
{
  typedef cicada::Vocab vocab_type;
  typedef cicada::Symbol word_type;
  typedef bitext_type::sentence_type sentence_type;

  bitexts.clear();
  dict_source_target.clear();
  dict_target_source.clear();

  utils::compress_istream src(source_file, 1024 * 1024);
  utils::compress_istream trg(target_file, 1024 * 1024);
  
  sentence_type source;
  sentence_type target;
  
  while (src && trg) {
    src >> source;
    trg >> target;
    
    if (! src || ! trg) break;
    
    bitexts.push_back(bitext_type(source, target));

    if (source.empty() || target.empty()) continue;
    
    sentence_type::const_iterator siter_begin = source.begin();
    sentence_type::const_iterator siter_end   = source.end();
    sentence_type::const_iterator titer_begin = target.begin();
    sentence_type::const_iterator titer_end   = target.end();
    
    {
      dictionary_type::dict_type& dict = dict_source_target[vocab_type::EPSILON];
      
      for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	++ dict[*titer];
      
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	dictionary_type::dict_type& dict = dict_source_target[*siter];
	
	for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	  ++ dict[*titer];
      }
    }

    {
      dictionary_type::dict_type& dict = dict_target_source[vocab_type::EPSILON];
      
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	++ dict[*siter];
      
      for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	dictionary_type::dict_type& dict = dict_target_source[*titer];
	
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  ++ dict[*siter];
      }
    }
  }
  
  if (src || trg)
    throw std::runtime_error("# of sentnces do not match");

  if (cutoff > 1) {
    typedef dictionary_type::dict_type::count_set_type word_set_type;
    
    word_set_type words_source;
    word_set_type words_target;
    
    const word_set_type& counts_source = dict_target_source[vocab_type::EPSILON].counts_;
    const word_set_type& counts_target = dict_source_target[vocab_type::EPSILON].counts_;
    
    word_set_type::const_iterator siter_end = counts_source.end();
    for (word_set_type::const_iterator siter = counts_source.begin(); siter != siter_end; ++ siter)
      if (siter->second >= cutoff)
	words_source.insert(*siter);
    
    word_set_type::const_iterator titer_end = counts_target.end();
    for (word_set_type::const_iterator titer = counts_target.begin(); titer != titer_end; ++ titer)
      if (titer->second >= cutoff)
	words_target.insert(*titer);
    
    dictionary_type dict_source_target_new;
    dictionary_type dict_target_source_new;
    
    for (word_type::id_type i = 0; i != dict_source_target.dicts_.size(); ++ i)
      if (dict_source_target.dicts_.exists(i)) {
	word_type source(i);
	
	if (source != vocab_type::EPSILON && words_source.find(source) == words_source.end())
	  source = vocab_type::UNK;
	
	dictionary_type::dict_type& dict = dict_source_target_new[source];
	
	word_set_type::const_iterator titer_end = dict_source_target[i].counts_.end();
	for (word_set_type::const_iterator titer = dict_source_target[i].counts_.begin(); titer != titer_end; ++ titer)
	  if (words_target.find(titer->first) == words_target.end())
	    dict[vocab_type::UNK] += titer->second;
	  else
	    dict[titer->first] += titer->second;
      }
    
    for (word_type::id_type i = 0; i != dict_target_source.dicts_.size(); ++ i)
      if (dict_target_source.dicts_.exists(i)) {
	word_type target(i);
	
	if (target != vocab_type::EPSILON && words_target.find(target) == words_target.end())
	  target = vocab_type::UNK;
	
	dictionary_type::dict_type& dict = dict_target_source_new[target];
	
	word_set_type::const_iterator siter_end = dict_target_source[i].counts_.end();
	for (word_set_type::const_iterator siter = dict_target_source[i].counts_.begin(); siter != siter_end; ++ siter)
	  if (words_source.find(siter->first) == words_source.end())
	    dict[vocab_type::UNK] += siter->second;
	  else
	    dict[siter->first] += siter->second;
      }

    dict_source_target.swap(dict_source_target_new);
    dict_target_source.swap(dict_target_source_new);
    
    bitext_set_type::iterator biter_end = bitexts.end();
    for (bitext_set_type::iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {

      sentence_type::iterator siter_end = biter->source_.end();
      for (sentence_type::iterator siter = biter->source_.begin(); siter != siter_end; ++ siter)
	if (words_source.find(*siter) == words_source.end())
	  *siter = vocab_type::UNK;

      sentence_type::iterator titer_end = biter->target_.end();
      for (sentence_type::iterator titer = biter->target_.begin(); titer != titer_end; ++ titer)
	if (words_target.find(*titer) == words_target.end())
	  *titer = vocab_type::UNK;	
    }
    
  }

  dict_source_target[vocab_type::BOS][vocab_type::BOS] = 1;
  dict_source_target[vocab_type::EOS][vocab_type::EOS] = 1;

  dict_target_source[vocab_type::BOS][vocab_type::BOS] = 1;
  dict_target_source[vocab_type::EOS][vocab_type::EOS] = 1;
  
  dict_source_target.initialize();
  dict_target_source.initialize();
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("source",    po::value<path_type>(&source_file),    "source file")
    ("target",    po::value<path_type>(&target_file),    "target file")
    
    ("embedding-source", po::value<path_type>(&embedding_source_file), "initial source embedding")
    ("embedding-target", po::value<path_type>(&embedding_target_file), "initial target embedding")

    ("derivation",               po::value<path_type>(&derivation_file),               "output derivation")
    ("alignment-source-target",  po::value<path_type>(&alignment_source_target_file),  "output alignemnt for P(target | source)")
    ("alignment-target-source",  po::value<path_type>(&alignment_target_source_file),  "output alignemnt for P(source | target)")

    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("alpha",     po::value<double>(&alpha)->default_value(alpha),      "parameter for reconstruction error")
    ("beta",      po::value<double>(&beta)->default_value(beta),        "parameter for classificaiton error")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("span",                po::value<int>(&span)->default_value(span),                               "span context size")
    ("window",              po::value<int>(&window)->default_value(window),                           "context window size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("baby-steps",        po::value<int>(&baby_steps)->default_value(baby_steps), "# of baby steps")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("cutoff",            po::value<int>(&cutoff)->default_value(cutoff),         "cutoff count for vocabulary (<= 1 to keep all)")
    ("beam",              po::value<int>(&beam)->default_value(beam),             "beam width for parsing")
    ("lambda",            po::value<double>(&lambda)->default_value(lambda),      "regularization constant")
    ("eta0",              po::value<double>(&eta0)->default_value(eta0),          "\\eta_0 for decay")

    ("moses", po::bool_switch(&moses_mode), "dump alignment in Moses format")
    ("giza",  po::bool_switch(&giza_mode),  "dump alignment in Giza format")
    ("dump",  po::bool_switch(&dump_mode),  "dump intermediate derivations and alignments")
    
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
