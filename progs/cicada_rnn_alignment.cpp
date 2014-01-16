//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// an implementation for neural network alignment model with NCE estimate...
// 
// we will try SGD with L2 regularizer inspired by AdaGrad (default)
//

#include <cstdlib>
#include <cmath>
#include <climits>

#include <set>

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

#include "cicada_rnn_alignment_impl.hpp"

typedef boost::filesystem::path path_type;

typedef cicada::Sentence sentence_type;
typedef cicada::Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;

typedef Model      model_type;
typedef Dictionary dictionary_type;

path_type source_file;
path_type target_file;

path_type embedding_source_file;
path_type embedding_target_file;

path_type output_source_target_file;
path_type output_target_source_file;
path_type alignment_source_target_file;
path_type alignment_target_source_file;

int dimension_embedding = 32;
int dimension_hidden = 128;
int window = 0;
int alignment = 8;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int baby_steps = 1;
int batch_size = 4;
int sample_size = 10;
int beam_size = 50;
int cutoff = 3;
double lambda = 0;
double lambda2 = 0;
double eta0 = 0.1;

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
		  model_type& theta_source_target,
		  model_type& theta_target_source);
void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       dictionary_type& dict_source_target,
	       dictionary_type& dict_target_source);
void viterbi(const bitext_set_type& bitexts,
	     const dictionary_type& dict_source_target,
	     const dictionary_type& dict_target_source,
	     const model_type& theta_source_target,
	     const model_type& theta_target_source);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (dimension_embedding <= 0)
      throw std::runtime_error("dimension must be positive");
    if (dimension_hidden <= 0)
      throw std::runtime_error("dimension must be positive");
    if (alignment <= 1)
      throw std::runtime_error("order size should be positive");

    if (sample_size <= 0)
      throw std::runtime_error("invalid sample size");
    if (beam_size <= 0)
      throw std::runtime_error("invalid beam size");
    if (batch_size <= 0)
      throw std::runtime_error("invalid batch size");

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

    if (source_file.empty() || target_file.empty())
      throw std::runtime_error("no data?");

    if (source_file != "-" && ! boost::filesystem::exists(source_file))
      throw std::runtime_error("no source file? " + source_file.string());
    
    if (target_file != "-" && ! boost::filesystem::exists(target_file))
      throw std::runtime_error("no target file? " + target_file.string());
    

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

    model_type theta_source_target(dimension_embedding, dimension_hidden, alignment, sources, targets, generator);
    model_type theta_target_source(dimension_embedding, dimension_hidden, alignment, targets, sources, generator);

    const size_t cols = utils::bithack::min(utils::bithack::min(theta_source_target.source_.cols(),
								theta_source_target.target_.cols()),
					    utils::bithack::min(theta_target_source.source_.cols(),
								theta_target_source.target_.cols()));
    
    theta_source_target.source_.block(0, 0, dimension_embedding, cols)
      = theta_target_source.target_.block(0, 0, dimension_embedding, cols);
    theta_source_target.target_.block(0, 0, dimension_embedding, cols)
      = theta_target_source.source_.block(0, 0, dimension_embedding, cols);

    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());
      
      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta_source_target.read_embedding(embedding_source_file, embedding_target_file);
      theta_target_source.read_embedding(embedding_target_file, embedding_source_file);
    }
        
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, alignment, lambda, lambda2, eta0),
		     bitexts,
		     dict_source_target,
		     dict_target_source,
		     theta_source_target,
		     theta_target_source);
      else
	learn_online(LearnSGD(lambda, lambda2, eta0),
		     bitexts,
		     dict_source_target,
		     dict_target_source,
		     theta_source_target,
		     theta_target_source);
    }

    if (! alignment_source_target_file.empty() || ! alignment_target_source_file.empty())
      viterbi(bitexts,
	      dict_source_target,
	      dict_target_source,
	      theta_source_target,
	      theta_target_source);
    
    if (! output_source_target_file.empty())
      theta_source_target.write(output_source_target_file);
    
    if (! output_target_source_file.empty())
      theta_target_source.write(output_target_source_file);
    
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
  
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Alignment alignment_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;

  struct bitext_alignment_type
  {
    size_type       id_;
    bitext_type     bitext_;
    alignment_type  alignment_;
    
    bitext_alignment_type() : id_(size_type(-1)), bitext_(), alignment_() {}
    bitext_alignment_type(const size_type& id,
			   const bitext_type& bitext,
			   const alignment_type& alignment)
      : id_(id), bitext_(bitext), alignment_(alignment) {}
    
    void swap(bitext_alignment_type& x)
    {
      std::swap(id_, x.id_);
      bitext_.swap(x.bitext_);
      alignment_.swap(x.alignment_);
    }

    void clear()
    {
      id_ = size_type(-1);
      bitext_.clear();
      alignment_.clear();
    }
  };
  
  typedef bitext_alignment_type value_type;

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

struct OutputAlignment : OutputMapReduce
{
  OutputAlignment(const path_type& path,
		  queue_type& queue)
    : path_(path),
      queue_(queue) {}
  
  void operator()()
  {
    if (path_.empty()) {
      bitext_alignment_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_reduced_type bitexts;
      bitext_alignment_type bitext;
      size_type id = 0;
      
      utils::compress_ostream os(path_, 1024 * 1024);
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
	
	// sort
	std::sort(bitext.alignment_.begin(), bitext.alignment_.end());

	if (bitext.id_ == id) {
	  if (moses_mode)
	    os << bitext.alignment_ << '\n';
	  else
	    output(os, bitext);
	  ++ id;
	} else
	  bitexts.insert(bitext);

	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  if (moses_mode)
	    os << bitexts.begin()->alignment_ << '\n';
	  else
	    output(os, *bitexts.begin());
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }

      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	if (moses_mode)
	  os << bitexts.begin()->alignment_ << '\n';
	else
	  output(os, *bitexts.begin());
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while writing alignment output?");
    }
  }

  typedef int index_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > align_set_type;
  typedef std::set<index_type, std::less<index_type>, std::allocator<index_type> > align_none_type;
  
  align_set_type  aligns_;
  align_none_type aligns_none_;

  void output(std::ostream& os, const bitext_alignment_type& bitext)
  {
    os << "# Sentence pair (" << (bitext.id_ + 1) << ')'
       << " source length " << bitext.bitext_.source_.size()
       << " target length " << bitext.bitext_.target_.size()
       << " alignment score : " << 0 << '\n';
    os << bitext.bitext_.target_ << '\n';
    
    if (bitext.bitext_.source_.empty() || bitext.bitext_.target_.empty()) {
      os << "NULL ({ })";
      sentence_type::const_iterator siter_end = bitext.bitext_.source_.end();
      for (sentence_type::const_iterator siter = bitext.bitext_.source_.begin(); siter != siter_end; ++ siter)
	os << ' ' << *siter << " ({ })";
      os << '\n';
    } else {
      aligns_.clear();
      aligns_.resize(bitext.bitext_.source_.size());
      
      aligns_none_.clear();
      for (size_type trg = 0; trg != bitext.bitext_.target_.size(); ++ trg)
	aligns_none_.insert(trg + 1);
      
      alignment_type::const_iterator aiter_end = bitext.alignment_.end();
      for (alignment_type::const_iterator aiter = bitext.alignment_.begin(); aiter != aiter_end; ++ aiter) {
	aligns_[aiter->source].push_back(aiter->target + 1);
	aligns_none_.erase(aiter->target + 1);
      }
      
      os << "NULL";
      os << " ({ ";
      std::copy(aligns_none_.begin(), aligns_none_.end(), std::ostream_iterator<index_type>(os, " "));
      os << "})";
      
      for (size_type src = 0; src != bitext.bitext_.source_.size(); ++ src) {
	os << ' ' << bitext.bitext_.source_[src];
	os << " ({ ";
	std::copy(aligns_[src].begin(), aligns_[src].end(), std::ostream_iterator<index_type>(os, " "));
	os << "})";
      }
      os << '\n';
    }
  }

  path_type   path_;
  queue_type& queue_;
};

template <typename Learner>
struct TaskAccumulate
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef HMM hmm_type;

  typedef hmm_type::log_likelihood_type log_likelihood_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef OutputMapReduce output_map_reduce_type;
  
  typedef output_map_reduce_type::queue_type queue_alignment_type;
  typedef output_map_reduce_type::value_type bitext_alignment_type;

  typedef std::pair<gradient_type*, gradient_type*> gradient_pair_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_mapper_type;
  typedef utils::lockfree_list_queue<gradient_pair_type, std::allocator<gradient_pair_type> > queue_merger_type;
  typedef std::vector<queue_merger_type, std::allocator<queue_merger_type> > queue_merger_set_type;
  
  typedef std::deque<gradient_type, std::allocator<gradient_type> > gradient_set_type;

  TaskAccumulate(const Learner& learner,
		 const bitext_set_type& bitexts,
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 const model_type& theta_source_target,
		 const model_type& theta_target_source,
		 const size_type batch_size,
		 queue_mapper_type& mapper,
		 queue_merger_set_type& mergers,
		 queue_alignment_type& queue_source_target,
		 queue_alignment_type& queue_target_source)
    : learner_source_target_(learner),
      learner_target_source_(learner),
      embedding_source_target_(theta_source_target.embedding_),
      embedding_target_source_(theta_target_source.embedding_),
      bitexts_(bitexts),
      theta_source_target_(theta_source_target),
      theta_target_source_(theta_target_source),
      mapper_(mapper),
      mergers_(mergers),
      queue_source_target_(queue_source_target),
      queue_target_source_(queue_target_source),
      hmm_source_target_(dict_source_target, sample_size, beam_size),
      hmm_target_source_(dict_target_source, sample_size, beam_size),
      batch_size_(batch_size)
  {
    generator_.seed(utils::random_seed());
  }

  void operator()()
  {
    clear();
    
    const size_type shard_size = mergers_.size();
    const size_type embedding_size = theta_source_target_.embedding_;
    const size_type hidden_size    = theta_source_target_.hidden_;
    const size_type alignment_size = theta_source_target_.alignment_;
    
    size_type batch = 0;
    gradient_pair_type grads;
    
    size_type merge_finished = 0;
    bool learn_finished = false;
    
    int non_found_iter = 0;
    
    bitext_alignment_type bitext_source_target;
    bitext_alignment_type bitext_target_source;
    
    while (merge_finished != shard_size || ! learn_finished) {
      bool found = false;
      
      if (merge_finished != shard_size)
	while (mergers_[shard_].pop(grads, true)) {
	  if (! grads.first)
	    ++ merge_finished;
	  else {
	    embedding_source_target_.assign(*grads.first,  theta_source_target_);
	    embedding_target_source_.assign(*grads.second, theta_target_source_);
	    
	    learner_source_target_(theta_source_target_, *grads.first,  embedding_target_source_);
	    learner_target_source_(theta_target_source_, *grads.second, embedding_source_target_);
	    
	    grads.first->increment();
	    grads.second->increment();
	  }
	  
	  found = true;
	}
      
      if (! learn_finished && mapper_.pop(batch, true)) {
	found = true;
	
	if (batch == size_type(-1)) {
	  // send termination!
	  for (size_type i = 0; i != shard_size; ++ i)
	    mergers_[i].push(std::make_pair(static_cast<gradient_type*>(0), static_cast<gradient_type*>(0)));
	  
	  learn_finished = true;
	} else {
	  gradient_type* grad_source_target = 0;
	  gradient_type* grad_target_source = 0;
	  
	  for (size_type j = 0; j != gradients_.size(); ++ j)
	    if (gradients_[j].shared() == shard_size) {
	      if (! grad_source_target)
		grad_source_target = &gradients_[j];
	      else if (! grad_target_source)
		grad_target_source = &gradients_[j];
	      
	      if (grad_source_target && grad_target_source) break;
	    }
	  
	  if (! grad_source_target) {
	    gradients_.push_back(gradient_type(embedding_size, hidden_size, alignment_size));
	    grad_source_target = &gradients_.back();
	  }
	  
	  if (! grad_target_source) {
	    gradients_.push_back(gradient_type(embedding_size, hidden_size, alignment_size));
	    grad_target_source = &gradients_.back();
	  }
	  
	  grad_source_target->clear();
	  grad_target_source->clear();
	  
	  const size_type first = batch * batch_size_;
	  const size_type last  = utils::bithack::min(first + batch_size_, bitexts_.size());
	  
	  for (size_type id = first; id != last; ++ id) {
	    const bitext_type& bitext = bitexts_[id];
	    
	    bitext_source_target.id_ = id;
	    bitext_source_target.bitext_.source_ = bitext.source_;
	    bitext_source_target.bitext_.target_ = bitext.target_;
	    bitext_source_target.alignment_.clear();
	    
	    bitext_target_source.id_ = id;
	    bitext_target_source.bitext_.source_ = bitext.target_;
	    bitext_target_source.bitext_.target_ = bitext.source_;
	    bitext_target_source.alignment_.clear();
	    
	    if (! bitext.source_.empty() && ! bitext.target_.empty()) {
	      hmm_source_target_.forward(bitext.source_, bitext.target_, theta_source_target_, bitext_source_target.alignment_);
	      hmm_target_source_.forward(bitext.target_, bitext.source_, theta_target_source_, bitext_target_source.alignment_);
	      
	      log_likelihood_source_target_
		+= hmm_source_target_.backward(bitext.source_,
					       bitext.target_,
					       theta_source_target_,
					       *grad_source_target,
					       generator_);
	      
	      log_likelihood_target_source_
		+= hmm_target_source_.backward(bitext.target_,
					       bitext.source_,
					       theta_target_source_,
					       *grad_target_source,
					       generator_);
	    }
	    
	    // reduce alignment
	    queue_source_target_.push_swap(bitext_source_target);
	    queue_target_source_.push_swap(bitext_target_source);
	  }
	  
	  embedding_source_target_.assign(*grad_source_target, theta_source_target_);
	  embedding_target_source_.assign(*grad_target_source, theta_target_source_);
	  
	  learner_source_target_(theta_source_target_, *grad_source_target, embedding_target_source_);
	  learner_target_source_(theta_target_source_, *grad_target_source, embedding_source_target_);
	  
	  grad_source_target->increment();
	  grad_target_source->increment();
	  
	  for (size_type i = 0; i != shard_size; ++ i)
	    if (i != shard_)
	      mergers_[i].push(std::make_pair(grad_source_target, grad_target_source));
	}
      }
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    theta_source_target_.finalize();
    theta_target_source_.finalize();
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
    log_likelihood_source_target_ = log_likelihood_type();
    log_likelihood_target_source_ = log_likelihood_type();
  }
  
  Learner   learner_source_target_;
  Learner   learner_target_source_;
  Embedding embedding_source_target_;
  Embedding embedding_target_source_;
  
  const bitext_set_type& bitexts_;
  model_type             theta_source_target_;
  model_type             theta_target_source_;

  queue_mapper_type&     mapper_;
  queue_merger_set_type& mergers_;
  queue_alignment_type& queue_source_target_;
  queue_alignment_type& queue_target_source_;
  
  hmm_type hmm_source_target_;
  hmm_type hmm_target_source_;
  
  gradient_set_type   gradients_;
  log_likelihood_type log_likelihood_source_target_;
  log_likelihood_type log_likelihood_target_source_;

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
		  model_type& theta_source_target,
		  model_type& theta_target_source)
{
  typedef TaskAccumulate<Learner> task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef typename task_type::size_type size_type;

  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  typedef typename task_type::queue_mapper_type     queue_mapper_type;
  typedef typename task_type::queue_merger_set_type queue_merger_set_type;
  
  typedef typename task_type::log_likelihood_type log_likelihood_type;

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
  
  typename output_map_reduce_type::queue_type queue_source_target;
  typename output_map_reduce_type::queue_type queue_target_source;
  
  task_set_type tasks(threads, task_type(learner,
					 bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta_source_target,
					 theta_target_source,
					 batch_size, 
					 mapper,
					 mergers,
					 queue_source_target,
					 queue_target_source));
  
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
    
    boost::thread output_source_target(output_alignment_type(! alignment_source_target_file.empty() && dump_mode
							     ? add_suffix(alignment_source_target_file, iter_tag)
							     : path_type(),
							     queue_source_target));
    boost::thread output_target_source(output_alignment_type(! alignment_target_source_file.empty() && dump_mode
							     ? add_suffix(alignment_target_source_file, iter_tag)
							     : path_type(),
							     queue_target_source));
    
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
    
    queue_source_target.push(typename output_map_reduce_type::value_type());
    queue_target_source.push(typename output_map_reduce_type::value_type());
    
    utils::resource end;
    
    log_likelihood_type log_likelihood_source_target;
    log_likelihood_type log_likelihood_target_source;
    
    for (size_type i = 0; i != tasks.size(); ++ i) {
      log_likelihood_source_target += tasks[i].log_likelihood_source_target_;
      log_likelihood_target_source += tasks[i].log_likelihood_target_source_;
    }
    
    if (debug)
      std::cerr << "log-likelihood P(target | source): " << static_cast<double>(log_likelihood_source_target) << std::endl
		<< "entropy        P(target | source): " << std::exp(- static_cast<double>(log_likelihood_source_target)) << std::endl
		<< "perplexity     P(target | source): " << (- static_cast<double>(log_likelihood_source_target) / std::log(2.0)) << std::endl
		<< "log-likelihood P(source | target): " << static_cast<double>(log_likelihood_target_source) << std::endl
		<< "entropy        P(source | target): " << std::exp(- static_cast<double>(log_likelihood_target_source)) << std::endl
		<< "perplexity     P(source | target): " << (- static_cast<double>(log_likelihood_target_source) / std::log(2.0)) << std::endl;

    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle bitexts!
    {
      boost::random_number_generator<boost::mt19937> gen(tasks.front().generator_);
      
      typename batch_set_type::iterator biter     = batches.begin();
      typename batch_set_type::iterator biter_end = batches.end();
      
      while (biter < biter_end) {
	typename batch_set_type::iterator iter_end = std::min(biter + utils::bithack::max(4096 / batch_size, 1), biter_end);
	
	std::random_shuffle(biter, iter_end, gen);
	biter = iter_end;
      }
    }
    
    //mixing
    for (size_type i = 1; i != tasks.size(); ++ i) {
      tasks[i].theta_source_target_ = tasks.front().theta_source_target_;
      tasks[i].theta_target_source_ = tasks.front().theta_target_source_;
    }
    
    output_source_target.join();
    output_target_source.join();
  }
  
  // copy models
  theta_source_target = tasks.front().theta_source_target_;
  theta_target_source = tasks.front().theta_target_source_;
}

struct TaskViterbi
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;

  typedef HMM hmm_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef OutputMapReduce output_map_reduce_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_type;

  typedef output_map_reduce_type::queue_type queue_alignment_type;
  typedef output_map_reduce_type::value_type bitext_alignment_type;

  TaskViterbi(const bitext_set_type& bitexts,
	      const dictionary_type& dict_source_target,
	      const dictionary_type& dict_target_source,
	      const model_type& theta_source_target,
	      const model_type& theta_target_source,
	      queue_type& queue,
	      queue_alignment_type& queue_source_target,
	      queue_alignment_type& queue_target_source)
    : bitexts_(bitexts),
      theta_source_target_(theta_source_target),
      theta_target_source_(theta_target_source),
      queue_(queue),
      queue_source_target_(queue_source_target),
      queue_target_source_(queue_target_source),
      hmm_source_target_(dict_source_target, sample_size, beam_size),
      hmm_target_source_(dict_target_source, sample_size, beam_size) {}

  void operator()()
  {
    bitext_alignment_type bitext_source_target;
    bitext_alignment_type bitext_target_source;
    
    size_type sentence_id;
    for (;;) {
      queue_.pop(sentence_id);
      
      if (sentence_id == size_type(-1)) break;
      
      const bitext_type& bitext = bitexts_[sentence_id];
      
      bitext_source_target.id_ = sentence_id;
      bitext_source_target.bitext_.source_ = bitext.source_;
      bitext_source_target.bitext_.target_ = bitext.target_;
      bitext_source_target.alignment_.clear();
      
      bitext_target_source.id_ = sentence_id;
      bitext_target_source.bitext_.source_ = bitext.target_;
      bitext_target_source.bitext_.target_ = bitext.source_;
      bitext_target_source.alignment_.clear();
      
      if (! bitext.source_.empty() && ! bitext.target_.empty()) {
	hmm_source_target_.viterbi(bitext.source_, bitext.target_, theta_source_target_, bitext_source_target.alignment_);
	
	hmm_target_source_.viterbi(bitext.target_, bitext.source_, theta_target_source_, bitext_target_source.alignment_);
      }
      
      // reduce alignment
      queue_source_target_.push_swap(bitext_source_target);
      queue_target_source_.push_swap(bitext_target_source);
    }
  }

  void clear()
  {
  }

  const bitext_set_type& bitexts_;
  const model_type& theta_source_target_;
  const model_type& theta_target_source_;
  
  queue_type&           queue_;
  queue_alignment_type& queue_source_target_;
  queue_alignment_type& queue_target_source_;
  
  hmm_type hmm_source_target_;
  hmm_type hmm_target_source_;
};

void viterbi(const bitext_set_type& bitexts,
	     const dictionary_type& dict_source_target,
	     const dictionary_type& dict_target_source,
	     const model_type& theta_source_target,
	     const model_type& theta_target_source)
{
  typedef TaskViterbi task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  task_type::queue_type   mapper(64 * threads);
  
  output_map_reduce_type::queue_type queue_source_target;
  output_map_reduce_type::queue_type queue_target_source;
  
  task_set_type tasks(threads, task_type(bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta_source_target,
					 theta_target_source,
					 mapper,
					 queue_source_target,
					 queue_target_source));

  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  
  boost::thread output_source_target(output_alignment_type(! alignment_source_target_file.empty()
							   ? alignment_source_target_file
							   : path_type(),
							   queue_source_target));
  boost::thread output_target_source(output_alignment_type(! alignment_target_source_file.empty()
							   ? alignment_target_source_file
							   : path_type(),
							   queue_target_source));

  if (debug)
    std::cerr << "Viterbi alignment" << std::endl;
  
  std::auto_ptr<boost::progress_display> progress(debug
						  ? new boost::progress_display(bitexts.size(), std::cerr, "", "", "")
						  : 0);
  
  utils::resource start;
  
  // actually run...
  for (size_type id = 0; id != bitexts.size(); ++ id) {
    mapper.push(id);
    
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
  
  queue_source_target.push(output_map_reduce_type::value_type());
  queue_target_source.push(output_map_reduce_type::value_type());

  output_source_target.join();
  output_target_source.join();
}

void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       dictionary_type& dict_source_target,
	       dictionary_type& dict_target_source)
{
  typedef cicada::Symbol word_type;
  typedef cicada::Vocab  vocab_type;
  
  bitexts.clear();
  dict_source_target.clear();
  dict_target_source.clear();

  utils::compress_istream src(source_file, 1024 * 1024);
  utils::compress_istream trg(target_file, 1024 * 1024);
  
  sentence_type source;
  sentence_type target;
  
  for (;;) {
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
    throw std::runtime_error("# of sentences does not match");
  
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
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    
    ("embedding-source", po::value<path_type>(&embedding_source_file), "initial source embedding")
    ("embedding-target", po::value<path_type>(&embedding_target_file), "initial target embedding")
    
    ("output-source-target", po::value<path_type>(&output_source_target_file), "output model parameter for P(target | source)")
    ("output-target-source", po::value<path_type>(&output_target_source_file), "output model parameter for P(source | target)")

    ("alignment-source-target", po::value<path_type>(&alignment_source_target_file), "output alignment for P(target | source)")
    ("alignment-target-source", po::value<path_type>(&alignment_target_source_file), "output alignment for P(source | target)")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("window",              po::value<int>(&window)->default_value(window),                           "context window size")
    ("alignment",           po::value<int>(&alignment)->default_value(alignment),                     "alignment model size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),     "max # of iterations")
    ("baby-steps",        po::value<int>(&baby_steps)->default_value(baby_steps),   "# of baby steps")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size),   "mini-batch size")
    ("sample",            po::value<int>(&sample_size)->default_value(sample_size), "sampling size")
    ("beam",              po::value<int>(&beam_size)->default_value(beam_size),     "histogram beam size")
    ("cutoff",            po::value<int>(&cutoff)->default_value(cutoff),           "cutoff count for vocabulary (<= 1 to keep all)")
    ("lambda",            po::value<double>(&lambda)->default_value(lambda),        "regularization constant")
    ("lambda2",           po::value<double>(&lambda2)->default_value(lambda2),      "regularization constant for bilingual agreement")
    ("eta0",              po::value<double>(&eta0)->default_value(eta0),            "\\eta_0 for decay")

    ("moses",      po::bool_switch(&moses_mode),       "dump alignment in Moses format")
    ("giza",       po::bool_switch(&giza_mode),        "dump alignment in Giza format")
    ("dump",       po::bool_switch(&dump_mode),        "dump intermediate alignments")
    
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
