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

#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"
#include "utils/mpi_traits.hpp"

#include "codec/lz4.hpp"

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

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta_source_target,
		  model_type& theta_target_source);
void bcast_model(model_type& theta);
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
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
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

    if (debug && mpi_rank == 0)
      std::cerr << "# of unique source words: " << sources.size() << std::endl
		<< "# of unique target words: " << targets.size() << std::endl
		<< "# of sentences: " << bitexts.size() << std::endl;

    model_type theta_source_target(dimension_embedding, dimension_hidden, alignment, sources, targets, generator);
    model_type theta_target_source(dimension_embedding, dimension_hidden, alignment, targets, sources, generator);
    
    if (mpi_rank == 0) {
      const size_t cols = utils::bithack::min(utils::bithack::min(theta_source_target.source_.cols(),
								  theta_source_target.target_.cols()),
					      utils::bithack::min(theta_target_source.source_.cols(),
								  theta_target_source.target_.cols()));
      
      theta_source_target.source_.block(0, 0, dimension_embedding, cols)
	= theta_target_source.target_.block(0, 0, dimension_embedding, cols);
      theta_source_target.target_.block(0, 0, dimension_embedding, cols)
	= theta_target_source.source_.block(0, 0, dimension_embedding, cols);
    }
    
    if (mpi_rank == 0)
      if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
	if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	  throw std::runtime_error("no embedding: " + embedding_source_file.string());
	
	if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	  throw std::runtime_error("no embedding: " + embedding_target_file.string());
	
	theta_source_target.read_embedding(embedding_source_file, embedding_target_file);
	theta_target_source.read_embedding(embedding_target_file, embedding_source_file);
      }
    
    bcast_model(theta_source_target);
    bcast_model(theta_target_source);
    
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
    
    if (mpi_rank == 0 && ! output_source_target_file.empty())
      theta_source_target.write(output_source_target_file);
    
    if (mpi_rank == 0 && ! output_target_source_file.empty())
      theta_target_source.write(output_target_source_file);
    
  } catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  
  return 0;
}

enum {
  bitext_tag = 1000,
  model_tag,
  gradient_tag,
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

struct MapReduce
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
    alignment_type  alignment_source_target_;
    alignment_type  alignment_target_source_;
    
    bitext_alignment_type()
      : id_(size_type(-1)),
	bitext_(),
	alignment_source_target_(),
	alignment_target_source_() {}
    bitext_alignment_type(const size_type& id,
			  const bitext_type& bitext,
			  const alignment_type& alignment_source_target,
			  const alignment_type& alignment_target_source)
      : id_(id),
	bitext_(bitext),
	alignment_source_target_(alignment_source_target),
	alignment_target_source_(alignment_target_source) {}
    
    void swap(bitext_alignment_type& x)
    {
      std::swap(id_, x.id_);
      bitext_.swap(x.bitext_);
      alignment_source_target_.swap(x.alignment_source_target_);
      alignment_target_source_.swap(x.alignment_target_source_);
    }

    void clear()
    {
      id_ = size_type(-1);
      bitext_.clear();
      alignment_source_target_.clear();
      alignment_target_source_.clear();
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
  
  
  struct codec_type
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    buffer_type buffer_;

    void encode(const bitext_alignment_type& bitext, std::string& encoded) const
    {
      buffer_type& buffer = const_cast<buffer_type&>(buffer_);
      buffer.clear();
      
      {
	boost::iostreams::filtering_ostream os;
	os.push(boost::iostreams::back_insert_device<buffer_type>(buffer));

	// id
	os.write((char*) &bitext.id_, sizeof(size_type));
	
	// source
	const size_type source_size = bitext.bitext_.source_.size();
	os.write((char*) &source_size, sizeof(size_type));

	sentence_type::const_iterator siter_end = bitext.bitext_.source_.end();
	for (sentence_type::const_iterator siter = bitext.bitext_.source_.begin(); siter != siter_end; ++ siter) {
	  const size_type word_size = siter->size();
	  os.write((char*) &word_size, sizeof(size_type));

	  os.write((char*) &(*siter->begin()), word_size);
	}
	
	// target
	const size_type target_size = bitext.bitext_.target_.size();
	os.write((char*) &target_size, sizeof(size_type));
	
	sentence_type::const_iterator titer_end = bitext.bitext_.target_.end();
	for (sentence_type::const_iterator titer = bitext.bitext_.target_.begin(); titer != titer_end; ++ titer) {
	  const size_type word_size = titer->size();
	  os.write((char*) &word_size, sizeof(size_type));
	  
	  os.write((char*) &(*titer->begin()), word_size);
	}
	
	// alignment
	write(os, bitext.alignment_source_target_);
	write(os, bitext.alignment_target_source_);
      }
      
      encoded = std::string(buffer.begin(), buffer.end());
    }

    void decode(bitext_alignment_type& bitext, const std::string& encoded) const
    {
      buffer_type& buffer = const_cast<buffer_type&>(buffer_);
      
      bitext.clear();
      
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::array_source(&(*encoded.begin()), encoded.size()));
      
      // id
      is.read((char*) &bitext.id_, sizeof(size_type));
      
      // source
      size_type source_size = 0;
      is.read((char*) &source_size, sizeof(size_type));
      
      for (size_type src = 0; src != source_size; ++ src) {
	size_type word_size = 0;
	is.read((char*) &word_size, sizeof(size_type));
	
	buffer.resize(word_size);
	is.read((char*) &(*buffer.begin()), word_size);
	
	bitext.bitext_.source_.push_back(word_type(buffer.begin(), buffer.end()));
      }
      
      // target
      size_type target_size = 0;
      is.read((char*) &target_size, sizeof(size_type));
      
      for (size_type trg = 0; trg != target_size; ++ trg) {
	size_type word_size = 0;
	is.read((char*) &word_size, sizeof(size_type));
	
	buffer.resize(word_size);
	is.read((char*) &(*buffer.begin()), word_size);
	
	bitext.bitext_.target_.push_back(word_type(buffer.begin(), buffer.end()));
      }
      
      // alignment
      read(is, bitext.alignment_source_target_);
      read(is, bitext.alignment_target_source_);
    }

    void write(std::ostream& os, const alignment_type& alignment) const
    {
      const size_type alignment_size = alignment.size();
      
      os.write((char*) &alignment_size, sizeof(size_type));
      os.write((char*) &(*alignment.begin()), sizeof(alignment_type::value_type) * alignment_size);
    }

    void read(std::istream& is, alignment_type& alignment) const
    {
      size_type alignment_size = 0;
      is.read((char*) &alignment_size, sizeof(size_type));
      
      alignment.resize(alignment_size);
      is.read((char*) &(*alignment.begin()), sizeof(alignment_type::value_type) * alignment_size);
    }
  };
};

namespace std
{
  inline
  void swap(MapReduce::value_type& x,
	    MapReduce::value_type& y)
  {
    x.swap(y);
  }
};

struct OutputAlignment : MapReduce
{
  OutputAlignment(const path_type& path_source_target,
		  const path_type& path_target_source,
		  queue_type& queue)
    : path_source_target_(path_source_target),
      path_target_source_(path_target_source),
      queue_(queue) {}
  
  void operator()()
  {
    if (path_source_target_.empty() && path_target_source_.empty()) {
      bitext_alignment_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_reduced_type bitexts;
      bitext_alignment_type bitext;
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
	
	// sort
	std::sort(bitext.alignment_source_target_.begin(), bitext.alignment_source_target_.end());
	std::sort(bitext.alignment_target_source_.begin(), bitext.alignment_target_source_.end());

	if (bitext.id_ == id) {
	  
	  if (os_source_target.get()) {
	    if (moses_mode)
	      *os_source_target << bitext.alignment_source_target_ << '\n';
	    else
	      output(*os_source_target,
		     bitext.id_, bitext.bitext_.source_, bitext.bitext_.target_, bitext.alignment_source_target_);
	  }
	  
	  if (os_target_source.get()) {
	    if (moses_mode)
	      *os_target_source << bitext.alignment_target_source_ << '\n';
	    else
	      output(*os_target_source,
		     bitext.id_, bitext.bitext_.target_, bitext.bitext_.source_, bitext.alignment_target_source_);
	  }

	  ++ id;
	} else
	  bitexts.insert(bitext);

	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  const bitext_alignment_type& bitext = *(bitexts.begin());
	  
	  if (os_source_target.get()) {
	    if (moses_mode)
	      *os_source_target << bitext.alignment_source_target_ << '\n';
	    else
	      output(*os_source_target,
		     bitext.id_, bitext.bitext_.source_, bitext.bitext_.target_, bitext.alignment_source_target_);
	  }
	  
	  if (os_target_source.get()) {
	    if (moses_mode)
	      *os_target_source << bitext.alignment_target_source_ << '\n';
	    else
	      output(*os_target_source,
		     bitext.id_, bitext.bitext_.target_, bitext.bitext_.source_, bitext.alignment_target_source_);
	  }
	  
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	const bitext_alignment_type& bitext = *(bitexts.begin());
	
	if (os_source_target.get()) {
	  if (moses_mode)
	    *os_source_target << bitext.alignment_source_target_ << '\n';
	  else
	    output(*os_source_target,
		   bitext.id_, bitext.bitext_.source_, bitext.bitext_.target_, bitext.alignment_source_target_);
	}
	
	if (os_target_source.get()) {
	  if (moses_mode)
	    *os_target_source << bitext.alignment_target_source_ << '\n';
	  else
	    output(*os_target_source,
		   bitext.id_, bitext.bitext_.target_, bitext.bitext_.source_, bitext.alignment_target_source_);
	}
	
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
      os << "NULL ({ ";
      for (size_type trg = 0; trg != target.size(); ++ trg)
	os << (trg + 1) << ' ';
      os << "})";
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

  path_type   path_source_target_;
  path_type   path_target_source_;
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

  typedef MapReduce map_reduce_type;
  
  typedef map_reduce_type::queue_type queue_bitext_type;
  typedef map_reduce_type::value_type bitext_alignment_type;
  
  typedef std::string encoded_type;
  
  typedef utils::lockfree_list_queue<encoded_type, std::allocator<encoded_type> > queue_gradient_type;

  typedef std::vector<char, std::allocator<char> > buffer_type;

  TaskAccumulate(const Learner& learner,
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 model_type& theta_source_target,
		 model_type& theta_target_source,
		 const size_type batch_size,
		 queue_bitext_type& bitext_mapper,
		 queue_bitext_type& bitext_reducer,
		 queue_gradient_type& gradient_mapper,
		 queue_gradient_type& gradient_reducer,
		 const int rank)
    : learner_source_target_(learner),
      learner_target_source_(learner),
      embedding_source_target_(theta_source_target.embedding_),
      embedding_target_source_(theta_target_source.embedding_),
      theta_source_target_(theta_source_target),
      theta_target_source_(theta_target_source),
      bitext_mapper_(bitext_mapper),
      bitext_reducer_(bitext_reducer),
      gradient_mapper_(gradient_mapper),
      gradient_reducer_(gradient_reducer),
      hmm_source_target_(dict_source_target, sample_size, beam_size),
      hmm_target_source_(dict_target_source, sample_size, beam_size),
      gradient_source_target_(theta_source_target.embedding_, theta_source_target.hidden_, theta_source_target.alignment_),
      gradient_target_source_(theta_target_source.embedding_, theta_target_source.hidden_, theta_target_source.alignment_),
      gradient_source_target_batch_(theta_source_target.embedding_, theta_source_target.hidden_, theta_source_target.alignment_),
      gradient_target_source_batch_(theta_target_source.embedding_, theta_target_source.hidden_, theta_target_source.alignment_),
      batch_size_(batch_size),
      rank_(rank)
  {
    generator_.seed(utils::random_seed());
  }

  void operator()()
  {
    clear();
    
    bitext_alignment_type bitext;
    size_type batch_learn = 0;
    
    encoded_type encoded;
    
    bool merge_finished = false;
    bool learn_finished = false;
    
    int non_found_iter = 0;
    
    while (! merge_finished || ! learn_finished) {
      bool found = false;
      
      if (! merge_finished)
	while (gradient_reducer_.pop(encoded, true)) {
	  if (encoded.empty())
	    merge_finished = true;
	  else {
	    boost::iostreams::filtering_istream is;
	    is.push(codec::lz4_decompressor());
	    is.push(boost::iostreams::array_source(&(*encoded.begin()), encoded.size()));
	    
	    is >> gradient_source_target_;
	    is >> gradient_target_source_;
	    
	    embedding_source_target_.assign(gradient_source_target_, theta_source_target_);
	    embedding_target_source_.assign(gradient_target_source_, theta_target_source_);
	    
	    learner_source_target_(theta_source_target_, gradient_source_target_, embedding_target_source_);
	    learner_target_source_(theta_target_source_, gradient_target_source_, embedding_source_target_);
	  }
	  
	  found = true;
	}
      
      if (! learn_finished && bitext_mapper_.pop(bitext, true)) {
	found = true;
	
	if (bitext.id_ != size_type(-1)) {
	  const sentence_type& source = bitext.bitext_.source_;
	  const sentence_type& target = bitext.bitext_.target_;
	  
	  bitext.alignment_source_target_.clear();
	  bitext.alignment_target_source_.clear();
	  
	  if (! source.empty() && ! target.empty()) {
	    hmm_source_target_.forward(source, target, theta_source_target_, bitext.alignment_source_target_);
	    hmm_target_source_.forward(target, source, theta_target_source_, bitext.alignment_target_source_);
	    
	    log_likelihood_source_target_
		+= hmm_source_target_.backward(source,
					       target,
					       theta_source_target_,
					       gradient_source_target_batch_,
					       generator_);
	    
	    log_likelihood_target_source_
	      += hmm_target_source_.backward(target,
					     source,
					     theta_target_source_,
					     gradient_target_source_batch_,
					     generator_);
	    
	    ++ batch_learn;
	  }
	  
	  bitext_reducer_.push_swap(bitext);
	} else
	  learn_finished = true;
	
	if (batch_learn == batch_size_ || (learn_finished && batch_learn)) {
	  embedding_source_target_.assign(gradient_source_target_batch_, theta_source_target_);
	  embedding_target_source_.assign(gradient_target_source_batch_, theta_target_source_);
	  
	  learner_source_target_(theta_source_target_, gradient_source_target_batch_, embedding_target_source_);
	  learner_target_source_(theta_target_source_, gradient_target_source_batch_, embedding_source_target_);
	  
	  buffer_.clear();
	  
	  {
	    boost::iostreams::filtering_ostream os;
	    os.push(codec::lz4_compressor());
	    os.push(boost::iostreams::back_insert_device<buffer_type>(buffer_));
	    
	    os << gradient_source_target_batch_;
	    os << gradient_target_source_batch_;
	  }
	  
	  gradient_source_target_batch_.clear();
	  gradient_target_source_batch_.clear();
	  
	  gradient_mapper_.push(encoded_type(buffer_.begin(), buffer_.end()));
	  
	  batch_learn = 0;
	}
	
	if (learn_finished) {
	  gradient_mapper_.push(encoded_type());
	  
	  if (rank_)
	    bitext_reducer_.push(bitext_alignment_type());
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

    gradient_source_target_.clear();
    gradient_target_source_.clear();
    gradient_source_target_batch_.clear();
    gradient_target_source_batch_.clear();
  }
  
  Learner   learner_source_target_;
  Learner   learner_target_source_;
  Embedding embedding_source_target_;
  Embedding embedding_target_source_;
  
  model_type&             theta_source_target_;
  model_type&             theta_target_source_;

  queue_bitext_type& bitext_mapper_;
  queue_bitext_type& bitext_reducer_;
  queue_gradient_type& gradient_mapper_;
  queue_gradient_type& gradient_reducer_;
  
  hmm_type hmm_source_target_;
  hmm_type hmm_target_source_;

  gradient_type gradient_source_target_;
  gradient_type gradient_target_source_;
  gradient_type gradient_source_target_batch_;
  gradient_type gradient_target_source_batch_;
  
  log_likelihood_type log_likelihood_source_target_;
  log_likelihood_type log_likelihood_target_source_;

  size_type      batch_size_;
  int            rank_;
  boost::mt19937 generator_;
  
  buffer_type buffer_;
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

struct less_bitexts
{
  typedef size_t size_type;

  less_bitexts(const bitext_set_type& bitexts)
    : bitexts_(bitexts) {}
  
  bool operator()(const size_type& x, const size_type& y) const
  {
    return ((bitexts_[x].source_.size() + bitexts_[x].target_.size())
	    < (bitexts_[y].source_.size() + bitexts_[y].target_.size()));
  }

  const bitext_set_type& bitexts_;
};

template <typename Learner>
void learn_online_root(const Learner& learner,
		       const bitext_set_type& bitexts,
		       const dictionary_type& dict_source_target,
		       const dictionary_type& dict_target_source,
		       model_type& theta_source_target,
		       model_type& theta_target_source)
{
  typedef TaskAccumulate<Learner> task_type;

  typedef typename task_type::size_type size_type;

  typedef MapReduce        map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  typedef typename task_type::queue_bitext_type   queue_bitext_type;
  typedef typename task_type::queue_gradient_type queue_gradient_type;

  typedef typename task_type::bitext_alignment_type bitext_alignment_type;
  
  typedef typename task_type::log_likelihood_type log_likelihood_type;
  
  typedef typename task_type::encoded_type buffer_type;

  typedef boost::shared_ptr<buffer_type> buffer_ptr_type;
  typedef std::deque<buffer_ptr_type, std::allocator<buffer_ptr_type> >  buffer_set_type;
  typedef std::vector<buffer_set_type, std::allocator<buffer_set_type> > buffer_map_type;
  
  typedef utils::mpi_ostream        bitext_ostream_type;
  typedef utils::mpi_istream_simple bitext_istream_type;
  
  typedef utils::mpi_ostream_simple gradient_ostream_type;
  typedef utils::mpi_istream_simple gradient_istream_type;

  typedef boost::shared_ptr<bitext_ostream_type> bitext_ostream_ptr_type;
  typedef boost::shared_ptr<bitext_istream_type> bitext_istream_ptr_type;

  typedef boost::shared_ptr<gradient_ostream_type> gradient_ostream_ptr_type;
  typedef boost::shared_ptr<gradient_istream_type> gradient_istream_ptr_type;

  typedef std::vector<bitext_ostream_ptr_type, std::allocator<bitext_ostream_ptr_type> > bitext_ostream_ptr_set_type;
  typedef std::vector<bitext_istream_ptr_type, std::allocator<bitext_istream_ptr_type> > bitext_istream_ptr_set_type;
  
  typedef std::vector<gradient_ostream_ptr_type, std::allocator<gradient_ostream_ptr_type> > gradient_ostream_ptr_set_type;
  typedef std::vector<gradient_istream_ptr_type, std::allocator<gradient_istream_ptr_type> > gradient_istream_ptr_set_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  id_set_type ids(bitexts.size());
  for (size_type id = 0; id != bitexts.size(); ++ id)
    ids[id] = id;
  
  queue_bitext_type bitext_mapper(1);
  queue_bitext_type bitext_reducer;
  queue_gradient_type gradient_mapper;
  queue_gradient_type gradient_reducer;

  task_type task(learner,
		 dict_source_target,
		 dict_target_source,
		 theta_source_target,
		 theta_target_source,
		 batch_size, 
		 bitext_mapper,
		 bitext_reducer,
		 gradient_mapper,
		 gradient_reducer,
		 mpi_rank);
    
  std::string           line;
  bitext_alignment_type bitext;
  
  buffer_type          buffer;
  buffer_map_type      buffers(mpi_size);
  
  bitext_ostream_ptr_set_type bitext_ostream(mpi_size);
  bitext_istream_ptr_set_type bitext_istream(mpi_size);

  gradient_ostream_ptr_set_type gradient_ostream(mpi_size);
  gradient_istream_ptr_set_type gradient_istream(mpi_size);
  
  typename map_reduce_type::codec_type codec;

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
      typename id_set_type::iterator biter     = ids.begin();
      typename id_set_type::iterator biter_end = ids.end();
      
      while (biter < biter_end) {
	typename id_set_type::iterator iter_end = std::min(biter + (batch_size << 5), biter_end);
	
	std::sort(biter, iter_end, less_bitexts(bitexts));
	biter = iter_end;
      }
    }

    MPI::COMM_WORLD.Barrier();
    
    // prepare iostreams...
    for (int rank = 0; rank != mpi_size; ++ rank)
      if (rank != mpi_rank) {
	bitext_ostream[rank].reset(new bitext_ostream_type(rank, bitext_tag));
	bitext_istream[rank].reset(new bitext_istream_type(rank, bitext_tag));
	
	gradient_ostream[rank].reset(new gradient_ostream_type(rank, gradient_tag));
	gradient_istream[rank].reset(new gradient_istream_type(rank, gradient_tag));
      }
    
    std::auto_ptr<boost::progress_display> progress(debug
						    ? new boost::progress_display(bitexts.size(), std::cerr, "", "", "")
						    : 0);

    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);
    
    boost::thread output(output_alignment_type(! alignment_source_target_file.empty() && dump_mode
					       ? add_suffix(alignment_source_target_file, iter_tag)
					       : path_type(),
					       ! alignment_target_source_file.empty() && dump_mode
					       ? add_suffix(alignment_target_source_file, iter_tag)
					       : path_type(),
					       bitext_reducer));

    // create thread!
    boost::thread worker(boost::ref(task));
    
    utils::resource start;
    
    bool gradient_mapper_finished = false; // for gradients
    bool gradient_reducer_finished = false; // for gradients
    bool bitext_finished   = false; // for bitext
    
    size_type id = 0;
    
    int non_found_iter = 0;
    for (;;) {
      bool found = false;
      
      // mapping of bitexts...
      if (! bitext_finished)
	for (int rank = 1; rank != mpi_size && id != bitexts.size(); ++ rank)
	  if (bitext_ostream[rank]->test()) {
	    bitext.id_     = id;
	    bitext.bitext_ = bitexts[id];
	    bitext.alignment_source_target_.clear();
	    bitext.alignment_target_source_.clear();
	    
	    codec.encode(bitext, line);
	    bitext_ostream[rank]->write(line);
	    
	    if (progress.get())
	      ++ (*progress);
	    
	    ++ id;
	    found = true;
	  }
      
      if (! bitext_finished)
	if (bitext_mapper.empty() && id != bitexts.size()) {
	  bitext.id_     = id;
	  bitext.bitext_ = bitexts[id];
	  bitext.alignment_source_target_.clear();
	  bitext.alignment_target_source_.clear();

	  //std::cerr << "rank: " << mpi_rank << " bitext: " << bitext.id_ << std::endl;
	  
	  bitext_mapper.push_swap(bitext);
	  
	  if (progress.get())
	    ++ (*progress);
	  
	  ++ id;
	  found = true;
	}
      
      // finished bitext mapping
      if (! bitext_finished && id == bitexts.size()) {
	bitext_mapper.push(bitext_alignment_type());
	bitext_finished = true;
      }
      
      // terminate bitext mapping
      if (bitext_finished)
	for (int rank = 1; rank != mpi_size; ++ rank)
	  if (bitext_ostream[rank] && bitext_ostream[rank]->test()) {
	    if (! bitext_ostream[rank]->terminated())
	      bitext_ostream[rank]->terminate();
	    else
	      bitext_ostream[rank].reset();
	    
	    found = true;
	  }
      
      // reduce bitexts from others...
      for (int rank = 1; rank != mpi_size; ++ rank)
	if (bitext_istream[rank] && bitext_istream[rank]->test()) {
	  if (bitext_istream[rank]->read(line)) {
	    codec.decode(bitext, line);
	    
	    //std::cerr << "reduced: " << rank << " bitext: " << bitext.id_ << std::endl;

	    bitext_reducer.push_swap(bitext);
	  } else
	    bitext_istream[rank].reset();
	  
	  found = true;
	}
      
      // reduce gradients
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && gradient_istream[rank] && gradient_istream[rank]->test()) {
	  if (gradient_istream[rank]->read(buffer))
	    gradient_reducer.push_swap(buffer);
	  else
	    gradient_istream[rank].reset();
	  
	  buffer.clear();
	  found = true;
	}
      
      // check termination...
      if (! gradient_reducer_finished
	  && std::count(gradient_istream.begin(), gradient_istream.end(), gradient_istream_ptr_type()) == mpi_size) {
	gradient_reducer.push(buffer_type());
	gradient_reducer_finished = true;
      }
      
      // bcast...
      // first, get the encoded buffer from mapper
      if (! gradient_mapper_finished && gradient_mapper.pop_swap(buffer, true)) {
	buffer_ptr_type buffer_ptr;
	
	if (! buffer.empty()) {
	  buffer_ptr.reset(new buffer_type());
	  buffer_ptr->swap(buffer);
	  buffer.clear();
	} else
	  gradient_mapper_finished = true;
	
	for (int rank = 0; rank != mpi_size; ++ rank) 
	  if (rank != mpi_rank)
	    buffers[rank].push_back(buffer_ptr);
	
	found = true;
      }
      
      // second, bcast...
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && gradient_ostream[rank] && gradient_ostream[rank]->test() && ! buffers[rank].empty()) {
	  if (! buffers[rank].front()) {
	    // termination!
	    if (! gradient_ostream[rank]->terminated())
	      gradient_ostream[rank]->terminate();
	    else {
	      gradient_ostream[rank].reset();
	      buffers[rank].erase(buffers[rank].begin());
	    }
	  } else {
	    gradient_ostream[rank]->write(*(buffers[rank].front()));
	    buffers[rank].erase(buffers[rank].begin());
	  }
	  
	  found = true;
	}
      
      // termination condition
      if (bitext_finished && gradient_reducer_finished && gradient_mapper_finished
	  && std::count(bitext_istream.begin(), bitext_istream.end(), bitext_istream_ptr_type()) == mpi_size
	  && std::count(bitext_ostream.begin(), bitext_ostream.end(), bitext_ostream_ptr_type()) == mpi_size
	  && std::count(gradient_istream.begin(), gradient_istream.end(), gradient_istream_ptr_type()) == mpi_size
	  && std::count(gradient_ostream.begin(), gradient_ostream.end(), gradient_ostream_ptr_type()) == mpi_size) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
    
    worker.join();
    
    utils::resource end;

    bitext_reducer.push(bitext_alignment_type());
    
    log_likelihood_type log_likelihood_source_target = task.log_likelihood_source_target_;
    log_likelihood_type log_likelihood_target_source = task.log_likelihood_target_source_;

    for (int rank = 1; rank != mpi_size; ++ rank) {
      log_likelihood_type llst;
      log_likelihood_type llts;
      
      boost::iostreams::filtering_istream is;
      is.push(utils::mpi_device_source(rank, log_likelihood_tag, 4096));
      is.read((char*) &llst, sizeof(log_likelihood_type));
      is.read((char*) &llts, sizeof(log_likelihood_type));
      
      log_likelihood_source_target += llst;
      log_likelihood_target_source += llts;
    }
    
    if (debug)
      std::cerr << "log-likelihood P(target | source): " << static_cast<double>(log_likelihood_source_target) << std::endl
		<< "perplexity     P(target | source): " << std::exp(- static_cast<double>(log_likelihood_source_target)) << std::endl
		<< "log-likelihood P(source | target): " << static_cast<double>(log_likelihood_target_source) << std::endl
		<< "perplexity     P(source | target): " << std::exp(- static_cast<double>(log_likelihood_target_source)) << std::endl;

    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle bitexts!
    {
      typename id_set_type::iterator biter     = ids.begin();
      typename id_set_type::iterator biter_end = ids.end();
      
      while (biter < biter_end) {
	typename id_set_type::iterator iter_end = std::min(biter + (batch_size << 5), biter_end);
	
	std::random_shuffle(biter, iter_end);
	biter = iter_end;
      }
    }

    // mixing
    bcast_model(theta_source_target);
    bcast_model(theta_target_source);
    
    output.join();
  }
}

template <typename Learner>
void learn_online_others(const Learner& learner,
			 const bitext_set_type& bitexts,
			 const dictionary_type& dict_source_target,
			 const dictionary_type& dict_target_source,
			 model_type& theta_source_target,
			 model_type& theta_target_source)
{
  typedef TaskAccumulate<Learner> task_type;

  typedef typename task_type::size_type size_type;

  typedef MapReduce        map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  typedef typename task_type::queue_bitext_type   queue_bitext_type;
  typedef typename task_type::queue_gradient_type queue_gradient_type;

  typedef typename task_type::bitext_alignment_type bitext_alignment_type;
  
  typedef typename task_type::log_likelihood_type log_likelihood_type;
  
  typedef typename task_type::encoded_type buffer_type;

  typedef boost::shared_ptr<buffer_type> buffer_ptr_type;
  typedef std::deque<buffer_ptr_type, std::allocator<buffer_ptr_type> >  buffer_set_type;
  typedef std::vector<buffer_set_type, std::allocator<buffer_set_type> > buffer_map_type;
    
  typedef utils::mpi_ostream_simple bitext_ostream_type;
  typedef utils::mpi_istream        bitext_istream_type;
  
  typedef utils::mpi_ostream_simple gradient_ostream_type;
  typedef utils::mpi_istream_simple gradient_istream_type;

  typedef boost::shared_ptr<bitext_ostream_type> bitext_ostream_ptr_type;
  typedef boost::shared_ptr<bitext_istream_type> bitext_istream_ptr_type;

  typedef boost::shared_ptr<gradient_ostream_type> gradient_ostream_ptr_type;
  typedef boost::shared_ptr<gradient_istream_type> gradient_istream_ptr_type;

  typedef std::vector<bitext_ostream_ptr_type, std::allocator<bitext_ostream_ptr_type> > bitext_ostream_ptr_set_type;
  typedef std::vector<bitext_istream_ptr_type, std::allocator<bitext_istream_ptr_type> > bitext_istream_ptr_set_type;
  
  typedef std::vector<gradient_ostream_ptr_type, std::allocator<gradient_ostream_ptr_type> > gradient_ostream_ptr_set_type;
  typedef std::vector<gradient_istream_ptr_type, std::allocator<gradient_istream_ptr_type> > gradient_istream_ptr_set_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  queue_bitext_type bitext_mapper(1);
  queue_bitext_type bitext_reducer;
  queue_gradient_type gradient_mapper;
  queue_gradient_type gradient_reducer;

  task_type task(learner,
		 dict_source_target,
		 dict_target_source,
		 theta_source_target,
		 theta_target_source,
		 batch_size, 
		 bitext_mapper,
		 bitext_reducer,
		 gradient_mapper,
		 gradient_reducer,
		 mpi_rank);
  
  std::string           line;
  bitext_alignment_type bitext;
  
  buffer_type          buffer;
  buffer_map_type      buffers(mpi_size);
  
  gradient_ostream_ptr_set_type gradient_ostream(mpi_size);
  gradient_istream_ptr_set_type gradient_istream(mpi_size);
  
  typename map_reduce_type::codec_type codec;

  for (int t = 0; t < iteration; ++ t) {
    MPI::COMM_WORLD.Barrier();
    
    bitext_istream_ptr_type bitext_istream(new bitext_istream_type(0, bitext_tag));
    bitext_ostream_ptr_type bitext_ostream(new bitext_ostream_type(0, bitext_tag));
    
    // prepare iostreams...
    for (int rank = 0; rank != mpi_size; ++ rank)
      if (rank != mpi_rank) {
	gradient_ostream[rank].reset(new gradient_ostream_type(rank, gradient_tag));
	gradient_istream[rank].reset(new gradient_istream_type(rank, gradient_tag));
      }

    // create thread!
    boost::thread worker(boost::ref(task));
    
    bool gradient_mapper_finished = false; // for gradients
    bool gradient_reducer_finished = false; // for gradients
    bool bitext_finished   = false; // for bitext
    
    int non_found_iter = 0;
    for (;;) {
      bool found = false;
      
      // read bitexts mapped from root
      if (bitext_istream && bitext_istream->test() && bitext_mapper.empty()) {
	if (bitext_istream->read(line)) {
	  codec.decode(bitext, line);
	  
	  //std::cerr << "rank: " << mpi_rank << " bitext: " << bitext.id_ << std::endl;

	  bitext_mapper.push_swap(bitext);
	} else {
	  bitext_mapper.push(bitext_alignment_type());
	  bitext_istream.reset();
	}
	
	found = true;
      }
      
      // reduce derivations to root
      if (! bitext_finished)
	if (bitext_ostream && bitext_ostream->test() && bitext_reducer.pop_swap(bitext, true)) {
	  if (bitext.id_ == size_type(-1))
	    bitext_finished = true;
	  else {
	    codec.encode(bitext, line);
	    
	    //std::cerr << "rank: " << mpi_rank << " bitext: " << bitext.id_ << std::endl;

	    bitext_ostream->write(line);
	  }
	  
	  found = true;
	}
      
      if (bitext_finished && bitext_ostream && bitext_ostream->test()) {
	if (! bitext_ostream->terminated())
	  bitext_ostream->terminate();
	else
	  bitext_ostream.reset();
	
	found = true;
      }
      
      // reduce gradients
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && gradient_istream[rank] && gradient_istream[rank]->test()) {
	  if (gradient_istream[rank]->read(buffer))
	    gradient_reducer.push_swap(buffer);
	  else
	    gradient_istream[rank].reset();
	  
	  buffer.clear();
	  found = true;
	}
      
      // check termination...
      if (! gradient_reducer_finished
	  && std::count(gradient_istream.begin(), gradient_istream.end(), gradient_istream_ptr_type()) == mpi_size) {
	gradient_reducer.push(buffer_type());
	gradient_reducer_finished = true;
      }
      
      // bcast...
      // first, get the encoded buffer from mapper
      if (! gradient_mapper_finished && gradient_mapper.pop_swap(buffer, true)) {
	buffer_ptr_type buffer_ptr;
	
	if (! buffer.empty()) {
	  buffer_ptr.reset(new buffer_type());
	  buffer_ptr->swap(buffer);
	  buffer.clear();
	} else
	  gradient_mapper_finished = true;
	
	for (int rank = 0; rank != mpi_size; ++ rank) 
	  if (rank != mpi_rank)
	    buffers[rank].push_back(buffer_ptr);
	
	found = true;
      }
      
      // second, bcast...
      for (int rank = 0; rank != mpi_size; ++ rank)
	if (rank != mpi_rank && gradient_ostream[rank] && gradient_ostream[rank]->test() && ! buffers[rank].empty()) {
	  if (! buffers[rank].front()) {
	    // termination!
	    if (! gradient_ostream[rank]->terminated())
	      gradient_ostream[rank]->terminate();
	    else {
	      gradient_ostream[rank].reset();
	      buffers[rank].erase(buffers[rank].begin());
	    }
	  } else {
	    gradient_ostream[rank]->write(*(buffers[rank].front()));
	    buffers[rank].erase(buffers[rank].begin());
	  }
	  
	  found = true;
	}

      // termination condition
      if (! bitext_istream && ! bitext_ostream
	  && gradient_reducer_finished
	  && gradient_mapper_finished
	  && std::count(gradient_istream.begin(), gradient_istream.end(), gradient_istream_ptr_type()) == mpi_size
	  && std::count(gradient_ostream.begin(), gradient_ostream.end(), gradient_ostream_ptr_type()) == mpi_size) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }

    worker.join();
    
    // reduce loss and # of parsed
    {
      boost::iostreams::filtering_ostream os;
      os.push(utils::mpi_device_sink(0, log_likelihood_tag, 4096));
      os.write((char*) &task.log_likelihood_source_target_, sizeof(log_likelihood_type));
      os.write((char*) &task.log_likelihood_target_source_, sizeof(log_likelihood_type));
    }

    // mixing
    bcast_model(theta_source_target);
    bcast_model(theta_target_source);
  }
}

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta_source_target,
		  model_type& theta_target_source)
{
  if (MPI::COMM_WORLD.Get_rank() == 0)
    learn_online_root(learner, bitexts, dict_source_target, dict_target_source, theta_source_target, theta_target_source);
  else
    learn_online_others(learner, bitexts, dict_source_target, dict_target_source, theta_source_target, theta_target_source);
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

  typedef MapReduce map_reduce_type;

  typedef map_reduce_type::queue_type queue_type;
  typedef map_reduce_type::value_type bitext_alignment_type;
  
  TaskViterbi(const dictionary_type& dict_source_target,
	      const dictionary_type& dict_target_source,
	      const model_type& theta_source_target,
	      const model_type& theta_target_source,
	      queue_type& mapper,
	      queue_type& reducer,
	      const int rank)
    : theta_source_target_(theta_source_target),
      theta_target_source_(theta_target_source),
      mapper_(mapper),
      reducer_(reducer),
      hmm_source_target_(dict_source_target, sample_size, beam_size),
      hmm_target_source_(dict_target_source, sample_size, beam_size),
      rank_(rank) {}

  void operator()()
  {
    bitext_alignment_type bitext;
    
    for (;;) {
      mapper_.pop_swap(bitext);
      
      if (bitext.id_ == size_type(-1)) break;

      const sentence_type& source = bitext.bitext_.source_;
      const sentence_type& target = bitext.bitext_.target_;

      bitext.alignment_source_target_.clear();
      bitext.alignment_target_source_.clear();

      if (! source.empty() && ! target.empty()) {
	hmm_source_target_.viterbi(source, target, theta_source_target_, bitext.alignment_source_target_);
	
	hmm_target_source_.viterbi(target, source, theta_target_source_, bitext.alignment_target_source_);
      }
      
      // reduce alignment
      reducer_.push_swap(bitext);
    }
    
    if (rank_)
      reducer_.push(bitext_alignment_type());
  }

  void clear()
  {
  }

  const model_type& theta_source_target_;
  const model_type& theta_target_source_;
  
  queue_type&           mapper_;
  queue_type&           reducer_;
  
  hmm_type hmm_source_target_;
  hmm_type hmm_target_source_;

  int rank_;
};

void viterbi_root(const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  const model_type& theta_source_target,
		  const model_type& theta_target_source)
{
  typedef TaskViterbi task_type;

  typedef task_type::size_type size_type;

  typedef MapReduce        map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  typedef map_reduce_type::bitext_alignment_type bitext_alignment_type;

  typedef utils::mpi_ostream        ostream_type;
  typedef utils::mpi_istream_simple istream_type;
  
  typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
  typedef boost::shared_ptr<istream_type> istream_ptr_type;
  
  typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;
  typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  task_type::queue_type mapper(1);
  task_type::queue_type reducer;
  
  boost::thread worker(task_type(dict_source_target,
				 dict_target_source,
				 theta_source_target,
				 theta_target_source,
				 mapper,
				 reducer,
				 mpi_rank));
  
  boost::thread output(output_alignment_type(! alignment_source_target_file.empty()
					     ? alignment_source_target_file
					     : path_type(),
					     ! alignment_target_source_file.empty()
					     ? alignment_target_source_file
					     : path_type(),
					     reducer));
  
  ostream_ptr_set_type ostream(mpi_size);
  istream_ptr_set_type istream(mpi_size);
  
  for (int rank = 1; rank != mpi_size; ++ rank) {
    ostream[rank].reset(new ostream_type(rank, bitext_tag));
    istream[rank].reset(new istream_type(rank, bitext_tag));
  }
  
  std::string           line;
  bitext_alignment_type bitext;

  map_reduce_type::codec_type codec;
  
  if (debug)
    std::cerr << "Viterbi alignment" << std::endl;

  std::auto_ptr<boost::progress_display> progress(debug
						  ? new boost::progress_display(bitexts.size(), std::cerr, "", "", "")
						  : 0);
  
  utils::resource start;
  
  int non_found_iter = 0;
  
  size_type id = 0;
  while (id != bitexts.size()) {
    bool found = false;
    
    // mapping of bitexts...
    for (int rank = 1; rank != mpi_size && id != bitexts.size(); ++ rank)
      if (ostream[rank]->test()) {
	bitext.id_     = id;
	bitext.bitext_ = bitexts[id];
	bitext.alignment_source_target_.clear();
	bitext.alignment_target_source_.clear();
	
	codec.encode(bitext, line);
	ostream[rank]->write(line);

	if (progress.get())
	  ++ (*progress);
	
	++ id;
	found = true;
      }
    
    if (mapper.empty() && id != bitexts.size()) {
      bitext.id_     = id;
      bitext.bitext_ = bitexts[id];
      bitext.alignment_source_target_.clear();
      bitext.alignment_target_source_.clear();
      
      mapper.push_swap(bitext);
      
      if (progress.get())
	++ (*progress);
      
      ++ id;
      found = true;
    }
    
    // reduce from others...
    for (int rank = 1; rank != mpi_size; ++ rank)
      if (istream[rank] && istream[rank]->test()) {
	if (istream[rank]->read(line)) {
	  codec.decode(bitext, line);
	  reducer.push_swap(bitext);
	} else
	  istream[rank].reset();
	  
	found = true;
      }
      
    non_found_iter = loop_sleep(found, non_found_iter);
  }

  bool terminated = false;
  
  for (;;) {
    bool found = false;
    
    // termination
    if (! terminated && mapper.push(bitext_alignment_type(), true)) {
      terminated = true;
      found = true;
    }
    
    // termination!
    for (int rank = 1; rank != mpi_size; ++ rank)
      if (ostream[rank] && ostream[rank]->test()) {
	if (! ostream[rank]->terminated())
	  ostream[rank]->terminate();
	else
	  ostream[rank].reset();
	  
	found = true;
      }
      
    // reduce from others...
    for (int rank = 1; rank != mpi_size; ++ rank)
      if (istream[rank] && istream[rank]->test()) {
	if (istream[rank]->read(line)) {
	  codec.decode(bitext, line);
	  reducer.push_swap(bitext);
	} else
	  istream[rank].reset();
	
	found = true;
      }
      
    // termination condition!
    if (std::count(istream.begin(), istream.end(), istream_ptr_type()) == mpi_size
	&& std::count(ostream.begin(), ostream.end(), ostream_ptr_type()) == mpi_size)
      break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  utils::resource end;
  
  if (debug)
    std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
	      << "user time:   " << end.user_time() - start.user_time() << std::endl;
  
  reducer.push(bitext_alignment_type());
  
  output.join();
  worker.join();
}

void viterbi_others(const bitext_set_type& bitexts,
		    const dictionary_type& dict_source_target,
		    const dictionary_type& dict_target_source,
		    const model_type& theta_source_target,
		    const model_type& theta_target_source)
{
  typedef TaskViterbi task_type;

  typedef task_type::size_type size_type;

  typedef MapReduce        map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  typedef map_reduce_type::bitext_alignment_type bitext_alignment_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  task_type::queue_type mapper(1);
  task_type::queue_type reducer;
  
  boost::thread worker(task_type(dict_source_target,
				 dict_target_source,
				 theta_source_target,
				 theta_target_source,
				 mapper,
				 reducer,
				 mpi_rank));

  boost::shared_ptr<utils::mpi_istream>        is(new utils::mpi_istream(0, bitext_tag));
  boost::shared_ptr<utils::mpi_ostream_simple> os(new utils::mpi_ostream_simple(0, bitext_tag));

  std::string           line;
  bitext_alignment_type bitext;
  
  map_reduce_type::codec_type codec;
  
  bool terminated = false;
  
  int non_found_iter = 0;
  for (;;) {
    bool found = false;
    
    // read bitexts mapped from root
    if (is && is->test() && mapper.empty()) {
      if (is->read(line)) {
	codec.decode(bitext, line);
	mapper.push_swap(bitext);
      } else {
	mapper.push(bitext_alignment_type());
	is.reset();
      }
	
      found = true;
    }
      
    // reduce derivations to root
    if (! terminated) {
      if (os && os->test() && reducer.pop_swap(bitext, true)) {
	if (bitext.id_ == size_type(-1))
	  terminated = true;
	else {
	  codec.encode(bitext, line);
	  os->write(line);
	}
	
	found = true;
      }
    } else {
      if (os && os->test()) {
	if (! os->terminated())
	  os->terminate();
	else
	  os.reset();
	
	found = true;
      }
    }
    
    if (! is && ! os) break;
    
    non_found_iter = loop_sleep(found, non_found_iter);
  }
  
  worker.join();
}

void viterbi(const bitext_set_type& bitexts,
	     const dictionary_type& dict_source_target,
	     const dictionary_type& dict_target_source,
	     const model_type& theta_source_target,
	     const model_type& theta_target_source)
{
  if (MPI::COMM_WORLD.Get_rank() == 0)
    viterbi_root(bitexts, dict_source_target, dict_target_source, theta_source_target, theta_target_source);
  else
    viterbi_others(bitexts, dict_source_target, dict_target_source, theta_source_target, theta_target_source);
}

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

void bcast_dict(dictionary_type& dict)
{
  typedef dictionary_type::word_type  word_type;
  typedef dictionary_type::count_type count_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == 0) {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    karma::uint_generator<count_type> generate_count;
    
    boost::iostreams::filtering_ostream os;
    os.push(codec::lz4_compressor());
    os.push(utils::mpi_device_bcast_sink(0, 1024 * 1024));
    
    std::ostream_iterator<char> iter(os);
    
    for (word_type::id_type id = 0; id != dict.dicts_.size(); ++ id)
      if (dict.dicts_.exists(id)) {
	const word_type source(id);

	typedef dictionary_type::dict_type::count_set_type count_set_type;
	
	count_set_type::const_iterator titer_end = dict.dicts_[id].counts_.end();
	for (count_set_type::const_iterator titer = dict.dicts_[id].counts_.begin(); titer != titer_end; ++ titer)
	  karma::generate(iter,
			  standard::string << karma::lit(' ') << standard::string << karma::lit(' ') << generate_count
			  << karma::lit('\n'),
			  source, titer->first, titer->second);
      }
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
    is.push(utils::mpi_device_bcast_source(0, 1024 * 1024));
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
    
    std::string source;
    std::string target;
    count_type  count;
    
    while (iter != iter_end) {
      source.clear();
      target.clear();
      
      if (! qi::phrase_parse(iter, iter_end,
			     parse_word >> parse_word >> parse_count >> (qi::eol | qi::eoi),
			     standard::blank, source, target, count))
	if (iter != iter_end)
	  throw std::runtime_error("parsing failed");

      dict[source][target] = count;
    }
  }
}

void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       dictionary_type& dict_source_target,
	       dictionary_type& dict_target_source)
{
  typedef cicada::Symbol word_type;
  typedef cicada::Vocab  vocab_type;
  typedef bitext_type::sentence_type sentence_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  bitexts.clear();
  dict_source_target.clear();
  dict_target_source.clear();
  
  if (mpi_rank == 0) {
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
  }
  
  // bcast dict_source_target and dict_target_source
  bcast_dict(dict_source_target);
  bcast_dict(dict_target_source);
  
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
