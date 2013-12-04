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

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <set>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <boost/algorithm/string/trim.hpp>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"
#include "cicada/alignment.hpp"

#include "utils/alloc_vector.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bichart.hpp"
#include "utils/bithack.hpp"
#include "utils/compact_map.hpp"
#include "utils/compact_set.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/mathop.hpp"
#include "utils/unordered_map.hpp"
#include "utils/repository.hpp"
#include "utils/program_options.hpp"
#include "utils/random_seed.hpp"
#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/vector2.hpp"
#include "utils/sampler.hpp"
#include "utils/resource.hpp"

#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>

struct Average
{
  Average() : average_(0), count_(0) {}
  
  Average& operator+=(const double& x)
  {
    average_ += (x - average_) / (++ count_);
    return *this;
  }
  
  Average& operator+=(const Average& x)
  {
    const uint64_t total = count_ + x.count_;
    
    average_ = average_ * (double(count_) / total) + x.average_ * (double(x.count_) / total);
    count_ = total;
    
    return *this;
  }
  
  operator const double&() const { return average_; }
  
  double   average_;
  uint64_t count_;
};


struct Gradient
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef utils::unordered_map<word_type, tensor_type,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, tensor_type> > >::type embedding_type;
  
  Gradient() : dimension_embedding_(0), dimension_hidden_(0), order_(0) {}
  Gradient(const size_type& dimension_embedding,
	   const size_type& dimension_hidden,
	   const int order) 
    : dimension_embedding_(dimension_embedding),
      dimension_hidden_(dimension_hidden),
      order_(order)
  { initialize(dimension_embedding, dimension_hidden, order); }
  
  Gradient& operator-=(const Gradient& x)
  {
    embedding_type::const_iterator iiter_end = x.embedding_input_.end();
    for (embedding_type::const_iterator iiter = x.embedding_input_.begin(); iiter != iiter_end; ++ iiter) {
      tensor_type& embedding = embedding_input_[iiter->first];
      
      if (! embedding.rows())
	embedding = - iiter->second;
      else
	embedding -= iiter->second;
    }
    
    embedding_type::const_iterator oiter_end = x.embedding_output_.end();
    for (embedding_type::const_iterator oiter = x.embedding_output_.begin(); oiter != oiter_end; ++ oiter) {
      tensor_type& embedding = embedding_output_[oiter->first];
      
      if (! embedding.rows())
	embedding = - oiter->second;
      else
	embedding -= oiter->second;
    }
    
    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());
    
    if (! Wh_.rows())
      Wh_ = tensor_type::Zero(x.Wh_.rows(), x.Wh_.cols());
    if (! bh_.rows())
      bh_ = tensor_type::Zero(x.bh_.rows(), x.bh_.cols());
    
    Wc_ -= x.Wc_;
    bc_ -= x.bc_;

    Wh_ -= x.Wh_;
    bh_ -= x.bh_;

    count_ -= x.count_;
    
    return *this;
  }
  
  Gradient& operator+=(const Gradient& x)
  {
    embedding_type::const_iterator iiter_end = x.embedding_input_.end();
    for (embedding_type::const_iterator iiter = x.embedding_input_.begin(); iiter != iiter_end; ++ iiter) {
      tensor_type& embedding = embedding_input_[iiter->first];
      
      if (! embedding.rows())
	embedding = iiter->second;
      else
	embedding += iiter->second;
    }

    embedding_type::const_iterator oiter_end = x.embedding_output_.end();
    for (embedding_type::const_iterator oiter = x.embedding_output_.begin(); oiter != oiter_end; ++ oiter) {
      tensor_type& embedding = embedding_output_[oiter->first];
      
      if (! embedding.rows())
	embedding = oiter->second;
      else
	embedding += oiter->second;
    }

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    if (! Wh_.rows())
      Wh_ = tensor_type::Zero(x.Wh_.rows(), x.Wh_.cols());
    if (! bh_.rows())
      bh_ = tensor_type::Zero(x.bh_.rows(), x.bh_.cols());
    
    Wc_ += x.Wc_;
    bc_ += x.bc_;
    
    Wh_ += x.Wh_;
    bh_ += x.bh_;

    count_ += x.count_;

    return *this;
  }
  
  void clear()
  {
    // embedding
    embedding_input_.clear();
    embedding_output_.clear();
    
    // matrix for context
    Wc_.setZero();
    bc_.setZero();
    
    // matrix for hidden layer
    Wh_.setZero();
    bh_.setZero();

    count_ = 0;
  }

  
  tensor_type& embedding_input(const word_type& word)
  {
    tensor_type& embedding = embedding_input_[word];
    if (! embedding.cols())
      embedding = tensor_type::Zero(dimension_embedding_, 1);
    
    return embedding;
  }

  tensor_type& embedding_output(const word_type& word)
  {
    tensor_type& embedding = embedding_output_[word];
    if (! embedding.cols())
      embedding = tensor_type::Zero(dimension_embedding_ + 1, 1);
    
    return embedding;
  }

  
  void initialize(const size_type dimension_embedding, const size_type dimension_hidden, const int order)
  {
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (order <= 0)
      throw std::runtime_error("invalid order");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    
    order_ = order;
    
    clear();
    
    // initialize...
    Wc_ = tensor_type::Zero(dimension_hidden_, dimension_embedding_ * (order - 1));
    bc_ = tensor_type::Zero(dimension_hidden_, 1).array();
    
    Wh_ = tensor_type::Zero(dimension_embedding_, dimension_hidden_);
    bh_ = tensor_type::Zero(dimension_embedding_, 1);

    count_ = 0;
  }
  
  // dimension...
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type order_;
  
  // embedding
  embedding_type embedding_input_;
  embedding_type embedding_output_;
  
  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;
  
  // Wh and bh for hidden layer
  tensor_type Wh_;
  tensor_type bh_;

  size_type count_;
};

struct Model
{
  // model parameters
  
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef Gradient gradient_type;

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef boost::filesystem::path path_type;

  typedef std::vector<bool, std::allocator<bool> > unique_set_type;
  typedef std::vector<word_type, std::allocator<word_type> > word_set_type;
  
  Model() : dimension_embedding_(0), dimension_hidden_(0), order_(0), scale_(1) {}
  template <typename Unigram, typename Gen>
  Model(const size_type& dimension_embedding,
	const size_type& dimension_hidden,
	const int order,
	const Unigram& unigram,
	Gen& gen) 
    : dimension_embedding_(dimension_embedding),
      dimension_hidden_(dimension_hidden),
      order_(order),
      scale_(1)
  { initialize(dimension_embedding, dimension_hidden, order, unigram, gen); }
  
  
  void clear()
  {
    // embedding
    embedding_input_.setZero();
    embedding_output_.setZero();
    
    // matrix for context
    Wc_.setZero();
    bc_.setZero();
    
    // matrix for hidden layer
    Wh_.setZero();
    bh_.setZero();

    scale_ = 1.0;
  }
  
  
  template <typename Gen>
  struct randomize
  {
    randomize(Gen& gen, const double range=0.01) : gen_(gen), range_(range) {}
    
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return boost::random::uniform_real_distribution<Tp>(-range_, range_)(const_cast<Gen&>(gen_));
    }
    
    Gen& gen_;
    double range_;
  };
  
  template <typename Unigram, typename Gen>
  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const int order,
		  const Unigram& unigram,
		  Gen& gen)
  {
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (order <= 0)
      throw std::runtime_error("invalid order");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    
    order_ = order;
    
    clear();
    
    const size_type vocabulary_size = word_type::allocated();

    const double range_e = std::sqrt(6.0 / (dimension_embedding_ + 1));
    
    embedding_input_  = tensor_type::Zero(dimension_embedding_,     vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    embedding_output_ = tensor_type::Zero(dimension_embedding_ + 1, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    embedding_output_.row(dimension_embedding_).setZero();
    
    uniques_ = unique_set_type(vocabulary_size, false);
    
    // assign id... which is used to compute the final word embeddings...
    words_.clear();
    words_.push_back(vocab_type::BOS);
    words_.push_back(vocab_type::EOS);
    words_.push_back(vocab_type::EPSILON);
    words_.push_back(vocab_type::UNK);
    
    for (size_type pos = 0; pos != unigram.words_.size(); ++ pos)
      if (unigram.words_[pos] != vocab_type::EOS) {
	uniques_[unigram.words_[pos].id()] = true;

	if (unigram.words_[pos] != vocab_type::BOS
	    && unigram.words_[pos] != vocab_type::EPSILON
	    && unigram.words_[pos] != vocab_type::UNK)
	  words_.push_back(unigram.words_[pos]);
      }
    
    word_set_type(words_).swap(words_);
    
    uniques_[vocab_type::BOS.id()] = false;
    uniques_[vocab_type::EOS.id()] = false;

    const double range_c = std::sqrt(6.0 / (dimension_hidden_ + dimension_embedding_ * (order - 1)));
    const double range_h = std::sqrt(6.0 / (dimension_embedding_ + dimension_hidden_));
    
    Wc_ = tensor_type::Zero(dimension_hidden_, dimension_embedding_ * (order - 1)).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Zero(dimension_hidden_, 1);
    
    Wh_ = tensor_type::Zero(dimension_embedding_, dimension_hidden_).array().unaryExpr(randomize<Gen>(gen, range_h));
    bh_ = tensor_type::Zero(dimension_embedding_, 1);

    scale_ = 1.0;
  }

  void finalize()
  {
    // clear unused entries
    embedding_input_.col(vocab_type::EOS.id())  = tensor_type::Zero(dimension_embedding_, 1);
    embedding_output_.col(vocab_type::BOS.id()) = tensor_type::Zero(dimension_embedding_ + 1, 1);
    embedding_output_.col(vocab_type::EPSILON.id()) = tensor_type::Zero(dimension_embedding_ + 1, 1);
    
    const double factor = 1.0 / std::accumulate(uniques_.begin(), uniques_.end(), size_type(0));
    
    tensor_type average_input  = tensor_type::Zero(dimension_embedding_, 1);
    tensor_type average_output = tensor_type::Zero(dimension_embedding_ + 1, 1);
    
    bool has_unk = false;
    
    for (size_type pos = 0; pos != uniques_.size(); ++ pos) 
      if (uniques_[pos]) {
	average_input  += embedding_input_.col(pos) * factor;
	average_output += embedding_output_.col(pos) * factor;
	
	has_unk |= (word_type(pos) == vocab_type::UNK);
      }
    
    if (! has_unk) {
      embedding_input_.col(vocab_type::UNK.id())  = average_input;
      embedding_output_.col(vocab_type::UNK.id()) = average_output;
    }
    
    if (scale_ != 1.0) {
      embedding_input_.array() *= scale_;
      embedding_output_.block(0, 0, dimension_embedding_, embedding_output_.cols()).array() *= scale_;
      
      scale_ = 1.0;
    }
  }
  
  struct real_policy : boost::spirit::karma::real_policies<parameter_type>
  {
    static unsigned int precision(parameter_type)
    {
      return 10;
    }
  };
  
  boost::spirit::karma::real_generator<parameter_type, real_policy> float10;
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(dimension_embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(dimension_hidden_);    
    rep["order"]     = utils::lexical_cast<std::string>(order_);
    rep["size"]      = utils::lexical_cast<std::string>(words_.size());
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("input.gz"),  rep.path("input.bin"), embedding_input_);
    write_embedding(rep.path("output.gz"), rep.path("output.bin"), embedding_output_);
    
    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
    
    write(rep.path("Wh.txt.gz"), rep.path("Wh.bin"), Wh_);
    write(rep.path("bh.txt.gz"), rep.path("bh.bin"), bh_);
    
    // vocabulary...
    vocab_type vocab;
    
    vocab.open(rep.path("vocab"), words_.size() >> 1);

    word_set_type::const_iterator witer_end = words_.end();
    for (word_set_type::const_iterator witer = words_.begin(); witer != witer_end; ++ witer)
      vocab.insert(*witer);
    
    vocab.close();
  }
  
  void write_embedding(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    const tensor_type::Index rows = matrix.rows();
    const tensor_type::Index cols = matrix.cols();
    
    utils::compress_ostream os_txt(path_text, 1024 * 1024);
    utils::compress_ostream os_bin(path_binary, 1024 * 1024);
    std::ostream_iterator<char> iter(os_txt);
    
    word_set_type::const_iterator witer_end = words_.end();
    for (word_set_type::const_iterator witer = words_.begin(); witer != witer_end; ++ witer) {
      karma::generate(iter, standard::string, *witer);
      
      for (difference_type j = 0; j != rows; ++ j)
	karma::generate(iter, karma::lit(' ') << float10, matrix(j, witer->id()));
      
      karma::generate(iter, karma::lit('\n'));

      os_bin.write((char*) matrix.col(witer->id()).data(), sizeof(tensor_type::Scalar) * rows);
    }
  }

  void write(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix) const
  {
    {
      utils::compress_ostream os(path_text, 1024 * 1024);
      os.precision(10);
      os << matrix;
    }
    
    {
      utils::compress_ostream os(path_binary, 1024 * 1024);
      
      const tensor_type::Index rows = matrix.rows();
      const tensor_type::Index cols = matrix.cols();
      
      os.write((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
    }
  }
  
  // dimension...
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type order_;
  
  // embedding
  tensor_type embedding_input_;
  tensor_type embedding_output_;

  unique_set_type uniques_;
  word_set_type   words_;
  
  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;
  
  // Wh and bh for hidden layer
  tensor_type Wh_;
  tensor_type bh_;

  // scale for embedding...
  double scale_;
};

struct Unigram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef uint64_t count_type;

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef std::vector<double, std::allocator<double> >         logprob_set_type;
  typedef std::vector<count_type, std::allocator<count_type> > count_set_type;
  typedef std::vector<word_type, std::allocator<word_type> >   word_map_type;
  
  typedef boost::random::discrete_distribution<> distribution_type;
  
  template <typename Tp>
  struct compare_pair
  {
    bool operator()(const Tp& x, const Tp& y) const
    {
      return (x.second > y.second
	      || (x.second == y.second
		  && static_cast<const std::string&>(x.first) < static_cast<const std::string&>(y.first)));
    }
  };
  
  template <typename Iterator>
  Unigram(Iterator first, Iterator last)
  {
    typedef std::pair<word_type, count_type> word_count_type;
    typedef std::vector<word_count_type, std::allocator<word_count_type> > word_count_set_type;
    
    word_count_set_type word_counts(first, last);
    std::sort(word_counts.begin(), word_counts.end(), compare_pair<word_count_type>());

    words_.reserve(word_counts.size());
    counts_.reserve(word_counts.size());
    logprobs_.reserve(word_type::allocated());
    
    word_count_set_type::const_iterator witer_end = word_counts.end();
    for (word_count_set_type::const_iterator witer = word_counts.begin(); witer != witer_end; ++ witer) {
      words_.push_back(witer->first);
      counts_.push_back(witer->second);
    }
    
    // initialize logprobs and words
    const double norm = 1.0 / std::accumulate(counts_.begin(), counts_.end(), double(0));
    for (word_type::id_type id = 0; id != counts_.size(); ++ id)
      logprobs_[words_[id].id()] = std::log(norm * counts_[id]);
    
    // initialize distribution
    distribution_ = distribution_type(counts_.begin(), counts_.end());
  }
  
  double logprob(const word_type& word) const
  {
    return logprobs_[word.id()];
  }
  
  template <typename Gen>
  word_type draw(Gen& gen) const
  {
    return words_[distribution_(gen)];
  }
  
  logprob_set_type  logprobs_;
  count_set_type    counts_;
  word_map_type     words_;
  distribution_type distribution_;
};

struct NGram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;
  typedef Unigram  unigram_type;

  typedef Average log_likelihood_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type    tensor_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef std::vector<tensor_type*, std::allocator<tensor_type*> > gradient_embedding_type;

  NGram(const unigram_type& unigram,
	const size_type samples)
    : unigram_(unigram), samples_(samples), log_samples_(std::log(double(samples))) {}

  const unigram_type& unigram_;
  size_type           samples_;
  double              log_samples_;
  
  tensor_type layer_input_;
  tensor_type layer_input_back_;
  tensor_type layer_context_;
  tensor_type layer_hidden_;

  tensor_type delta_input_;
  tensor_type delta_context_;
  tensor_type delta_hidden_;

  gradient_embedding_type gradient_embedding_;
  gradient_embedding_type gradient_embedding_back_;
  
  struct hinge
  {
    // 50 for numerical stability...
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return std::min(std::max(x, Tp(0)), Tp(50));
    }
  };
  
  struct dhinge
  {
    // 50 for numerical stability...
    
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return Tp(0) < x && x < Tp(50);
    }
  };
  
  template <typename Gen>
  log_likelihood_type learn(const sentence_type& sentence,
			    const model_type& theta,
			    gradient_type& gradient,
			    Gen& gen)
  {
    const size_type dimension = theta.dimension_embedding_;
    const size_type order     = theta.order_;
    
    gradient_embedding_.resize(order - 1);
    
    if (! layer_input_.cols())
      layer_input_ = tensor_type::Zero(dimension * (order - 1), 1);
    
    tensor_type& gradient_embedding_bos = gradient.embedding_input(vocab_type::BOS);
    tensor_type& gradient_embedding_eps = gradient.embedding_input(vocab_type::EPSILON);

    size_type eps_size = (order > 2 ? order - 2 : size_type(0));
    
    if (order > 2)
      for (size_type i = 0; i < order - 2; ++ i) {
	layer_input_.block(dimension * i, 0, dimension, 1) = theta.embedding_input_.col(vocab_type::EPSILON.id()) * theta.scale_;
	
	gradient_embedding_[i] = &gradient_embedding_eps;
      }

    layer_input_.block(dimension * (order - 2), 0, dimension, 1) = theta.embedding_input_.col(vocab_type::BOS.id()) * theta.scale_;
    gradient_embedding_[order - 2] = &gradient_embedding_bos;
    
    log_likelihood_type log_likelihood;
    sentence_type::const_iterator siter_begin = sentence.begin();
    sentence_type::const_iterator siter_end   = sentence.end();
    for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
      log_likelihood += learn(*siter, theta, gradient, gen);
      
#if 1
      // shift layer_input...
      if (order > 2)
	for (size_type i = 0; i < order - 2; ++ i)
	  layer_input_.block(dimension * i, 0, dimension, 1) = layer_input_.block(dimension * (i + 1), 0, dimension, 1);
      layer_input_.block(dimension * (order - 2), 0, dimension, 1) = theta.embedding_input_.col(siter->id()) * theta.scale_;
      
      // shift context
      std::copy(gradient_embedding_.begin() + 1, gradient_embedding_.end(), gradient_embedding_.begin());
      gradient_embedding_.back() = &gradient.embedding_input(*siter);
#endif
      
#if 0
      layer_input_back_        = layer_input_;
      gradient_embedding_back_ = gradient_embedding_;
      
      // compute lower-order ngram language model
      for (size_type i = eps_size; i != order - 1; ++ i) {
	layer_input_.block(dimension * i, 0, dimension, 1) = theta.embedding_input_.col(vocab_type::EPSILON.id()) * theta.scale_;
	gradient_embedding_[i] = &gradient_embedding_eps;
	
	learn(*siter, theta, gradient, gen);
      }
      
      // shift layer_input...
      layer_input_.block(0, 0, dimension * (order - 2), 1) = layer_input_back_.block(dimension, 0, dimension * (order - 2), 1);
      layer_input_.block(dimension * (order - 2), 0, dimension, 1) = theta.embedding_input_.col(siter->id()) * theta.scale_;
      
      // shift context
      std::copy(gradient_embedding_back_.begin() + 1, gradient_embedding_back_.end(), gradient_embedding_.begin());
      gradient_embedding_.back() = &gradient.embedding_input(*siter);      
#endif
      
      
      eps_size -= (eps_size > 0);
    }
    
    // correct scoring
    log_likelihood += learn(vocab_type::EOS, theta, gradient, gen);
    
#if 0
    // compute lower-order ngram language model
    for (size_type i = eps_size; i != order - 1; ++ i) {
      layer_input_.block(dimension * i, 0, dimension, 1) = theta.embedding_input_.col(vocab_type::EPSILON.id()) * theta.scale_;
      gradient_embedding_[i] = &gradient_embedding_eps;
      
      learn(vocab_type::EOS, theta, gradient, gen);
    }
#endif
    
    return log_likelihood;
  }
  
  template <typename Gen>
  double learn(const word_type& word,
	       const model_type& theta,
	       gradient_type& gradient,
	       Gen& gen)
  {
    const size_type dimension = theta.dimension_embedding_;
    const size_type order     = theta.order_;
    
    layer_context_           = (theta.Wc_ * layer_input_  + theta.bc_).array().unaryExpr(hinge());
    layer_hidden_            = (theta.Wh_ * layer_context_ + theta.bh_).array().unaryExpr(hinge());
    
    const double score = (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1).transpose() * layer_hidden_ * theta.scale_
			  + theta.embedding_output_.col(word.id()).block(dimension, 0, 1, 1))(0, 0);
    const double score_noise = log_samples_ + unigram_.logprob(word);
    const double z = utils::mathop::logsum(score, score_noise);
    const double logprob = score - z;
    const double logprob_noise = score_noise - z;

    // we suffer loss...
    const double loss = - 1.0 + std::exp(logprob);
    double log_likelihood = logprob;
          
    // updated output embedding...
    tensor_type& dembedding = gradient.embedding_output(word);
    
    dembedding.block(0, 0, dimension, 1).array() += loss * layer_hidden_.array();
    dembedding.block(dimension, 0, 1, 1).array() += loss;
    
    // backward...
    delta_hidden_ = (layer_hidden_.array().unaryExpr(dhinge())
		     * (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1) * loss * theta.scale_).array());
    
    for (size_type k = 0; k != samples_; ++ k) {
      const word_type word = unigram_.draw(gen);

      const double score = (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1).transpose() * layer_hidden_ * theta.scale_
			    + theta.embedding_output_.col(word.id()).block(dimension, 0, 1, 1))(0, 0);
      const double score_noise = log_samples_ + unigram_.logprob(word);
      const double z = utils::mathop::logsum(score, score_noise);
      const double logprob = score - z;
      const double logprob_noise = score_noise - z;
      
      log_likelihood += logprob_noise;
      
      // we suffer loss...
      const double loss = std::exp(logprob);
            
      tensor_type& dembedding = gradient.embedding_output(word);
      
      dembedding.block(0, 0, dimension, 1).array() += loss * layer_hidden_.array();
      dembedding.block(dimension, 0, 1, 1).array() += loss;
      
      delta_hidden_.array() += (layer_hidden_.array().unaryExpr(dhinge())
				* (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1) * loss * theta.scale_).array());
    }
    
    gradient.Wh_.noalias() += delta_hidden_ * layer_context_.transpose();
    gradient.bh_.noalias() += delta_hidden_;
    
    delta_context_ = (layer_context_.array().unaryExpr(dhinge()) * (theta.Wh_.transpose() * delta_hidden_).array());
    
    gradient.Wc_.noalias() += delta_context_ * layer_input_.transpose();
    gradient.bc_.noalias() += delta_context_;
    
    delta_input_.noalias() = theta.Wc_.transpose() * delta_context_;
    
    // finally, input embedding...
    for (int i = 0; i != order - 1; ++ i)
      *gradient_embedding_[i] += delta_input_.block(dimension * i, 0, dimension, 1);

    // increment
    ++ gradient.count_;

    return log_likelihood;
  }
};

struct Learn
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef cicada::Symbol   word_type;
  
  typedef model_type::tensor_type tensor_type;
  
  typedef std::pair<model_type*, const gradient_type*> update_type;
  
  typedef utils::lockfree_list_queue<update_type, std::allocator<update_type> > queue_type;

  struct Counter
  {
    Counter() : counter(0) {}
  
    void increment()
    {
      utils::atomicop::fetch_and_add(counter, size_type(1));
    }
  
    void wait(size_type target)
    {
      for (;;) {
	for (int i = 0; i < 64; ++ i) {
	  if (counter == target)
	    return;
	  else
	    boost::thread::yield();
	}
	
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      }
    }
    
    void clear() { counter = 0; }
    
    volatile size_type counter;
  };
  typedef Counter counter_type;
};

struct LearnAdaGrad : public Learn
{
  
  struct Updator
  {
    Updator(LearnAdaGrad& learner,
	    queue_type& queue,
	    counter_type& counter,
	    size_type shard,
	    size_type size)
      : learner_(learner), queue_(queue), counter_(counter), shard_(shard), size_(size) {}

    void operator()()
    {
      update_type update;

      for (;;) {
	queue_.pop(update);
	
	if (! update.first) break;
	
	model_type& theta = *update.first;
	const gradient_type& gradient = *update.second;
	
	typedef gradient_type::embedding_type embedding_type;

	const double scale = 1.0 / gradient.count_;
	
	embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
	for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
	  if (iiter->first.id() % size_ == shard_)
	    learner_.update(iiter->first,
			    theta.embedding_input_,
			    learner_.embedding_input_,
			    iiter->second,
			    scale,
			    learner_.lambda_ != 0.0,
			    false);
	
	embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
	for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
	  if (oiter->first.id() % size_ == shard_)
	    learner_.update(oiter->first,
			    theta.embedding_output_,
			    learner_.embedding_output_,
			    oiter->second,
			    scale,
			    learner_.lambda_ != 0.0,
			    true);

	counter_.increment();
      }
    }
    
    LearnAdaGrad& learner_;
    queue_type&   queue_;
    counter_type& counter_;

    size_type shard_;
    size_type size_;
  };

  LearnAdaGrad(const size_type& dimension_embedding,
	       const size_type& dimension_hidden,
	       const int order,
	       const double& lambda,
	       const double& eta0,
	       const int threads)
    : dimension_embedding_(dimension_embedding),
      dimension_hidden_(dimension_hidden),
      order_(order),
      lambda_(lambda),
      eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    
    embedding_input_  = tensor_type::Zero(dimension_embedding_,     vocabulary_size);
    embedding_output_ = tensor_type::Zero(dimension_embedding_ + 1, vocabulary_size);
    
    // initialize...
    Wc_ = tensor_type::Zero(dimension_hidden_, dimension_embedding_ * (order - 1));
    bc_ = tensor_type::Zero(dimension_hidden_, 1);
    
    Wh_ = tensor_type::Zero(dimension_embedding_, dimension_hidden_);
    bh_ = tensor_type::Zero(dimension_embedding_, 1);
    
    for (int i = 0; i != threads; ++ i)
      workers_.add_thread(new boost::thread(Updator(*this, queue_, counter_, i, threads)));
  }
  
  ~LearnAdaGrad()
  {
    for (size_type i = 0; i != workers_.size(); ++ i)
      queue_.push(update_type(0, 0));
    
    workers_.join_all();
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    // parallelize here...
    const_cast<counter_type&>(counter_).clear();
    
    for (size_type i = 0; i != workers_.size(); ++ i)
      const_cast<queue_type&>(queue_).push(update_type(&theta, &gradient));
    
#if 0
    embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
    for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
      update(iiter->first,
	     theta.embedding_input_,
	     const_cast<tensor_type&>(embedding_input_),
	     iiter->second,
	     1.0 / gradient.count_,
	     lambda_ != 0.0,
	     false);

    embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
    for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
      update(oiter->first,
	     theta.embedding_output_,
	     const_cast<tensor_type&>(embedding_output_),
	     oiter->second,
	     1.0 / gradient.count_,
	     lambda_ != 0.0,
	     true);
#endif
    
    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, 1.0 / gradient.count_, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, 1.0 / gradient.count_, false);

    update(theta.Wh_, const_cast<tensor_type&>(Wh_), gradient.Wh_, 1.0 / gradient.count_, lambda_ != 0.0);
    update(theta.bh_, const_cast<tensor_type&>(bh_), gradient.bh_, 1.0 / gradient.count_, false);

    // wait...
    const_cast<counter_type&>(counter_).wait(workers_.size());
  }

  template <typename Theta, typename GradVar, typename Grad>
  struct update_visitor_regularize
  {
    update_visitor_regularize(Eigen::MatrixBase<Theta>& theta,
			      Eigen::MatrixBase<GradVar>& G,
			      const Eigen::MatrixBase<Grad>& g,
			      const double& scale,
			      const double& lambda,
			      const double& eta0)
      : theta_(theta), G_(G), g_(g), scale_(scale), lambda_(lambda), eta0_(eta0) {}
    
    void init(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      operator()(value, i, j);
    }
    
    void operator()(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      if (g_(i, j) == 0) return;
      
      G_(i, j) += g_(i, j) * g_(i, j) * scale_ * scale_;
      
      const double rate = eta0_ / std::sqrt(double(1.0) + G_(i, j));
      const double f = theta_(i, j) - rate * scale_ * g_(i, j);

      theta_(i, j) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
    }
    
    Eigen::MatrixBase<Theta>&      theta_;
    Eigen::MatrixBase<GradVar>&    G_;
    const Eigen::MatrixBase<Grad>& g_;
    
    const double scale_;
    const double lambda_;
    const double eta0_;
  };

  struct learning_rate
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return (x == 0.0 ? 0.0 : 1.0 / std::sqrt(double(1.0) + x));
    }
  };
  
  template <typename Theta, typename GradVar, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool regularize=true) const
  {
    if (regularize) {
      update_visitor_regularize<Theta, GradVar, Grad> visitor(theta, G, g, scale, lambda_, eta0_);
      
      theta.visit(visitor);
    } else {
      G.array() += g.array().square() * scale * scale;
      theta.array() -= eta0_ * scale * g.array() * G.array().unaryExpr(learning_rate());
    }
  }

  template <typename Theta, typename GradVar, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool regularize=true,
	      const bool bias_last=false) const
  {
    if (regularize) {
      for (int row = 0; row != g.rows() - bias_last; ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0) * scale * scale;
	  
	  const double rate = eta0_ / std::sqrt(double(1.0) + G(row, word.id()));
	  const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	  
	  theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
      
      if (bias_last) {
	const int row = g.rows() - 1;
	
	if (g(row, 0) != 0) {
	  G(row, word.id()) += g(row, 0) * g(row, 0) * scale * scale;
	  theta(row, word.id()) -= eta0_ * scale * g(row, 0) / std::sqrt(double(1.0) + G(row, word.id()));
	}
      }
    } else {
      G.col(word.id()).array() += g.array().square() * scale * scale;
      theta.col(word.id()).array() -= eta0_ * scale * g.array() * G.col(word.id()).array().unaryExpr(learning_rate());
    }
  }
  
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type order_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type embedding_input_;
  tensor_type embedding_output_;

  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;
  
  // Wh and bh for hidden layer
  tensor_type Wh_;
  tensor_type bh_;  

  queue_type   queue_;
  counter_type counter_;
  boost::thread_group workers_;
};

struct LearnSGD : public Learn
{  

  struct Updator
  {
    Updator(LearnSGD& learner,
	    queue_type& queue,
	    counter_type& counter,
	    size_type shard,
	    size_type size)
      : learner_(learner), queue_(queue), counter_(counter), shard_(shard), size_(size) {}

    void operator()()
    {
      update_type update;

      for (;;) {
	queue_.pop(update);
	
	if (! update.first) break;
	
	model_type& theta = *update.first;
	const gradient_type& gradient = *update.second;
	
	typedef gradient_type::embedding_type embedding_type;
	
	const double scale = 1.0 / gradient.count_;
	
	embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
	for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
	  if (iiter->first.id() % size_ == shard_)
	    learner_.update(iiter->first,
			    theta.embedding_input_,
			    iiter->second,
			    scale,
			    theta.scale_,
			    false);
	
	embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
	for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
	  if (oiter->first.id() % size_ == shard_)
	    learner_.update(oiter->first,
			    theta.embedding_output_,
			    oiter->second,
			    scale,
			    theta.scale_,
			    true);

	counter_.increment();
      }
    }
    
    LearnSGD&     learner_;
    queue_type&   queue_;
    counter_type& counter_;

    size_type shard_;
    size_type size_;
  };

  LearnSGD(const double& lambda,
	   const double& eta0,
	   const int threads)
    : lambda_(lambda),
      eta0_(eta0),
      epoch_(0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    for (int i = 0; i != threads; ++ i)
      workers_.add_thread(new boost::thread(Updator(*this, queue_, counter_, i, threads)));
  }
  
  ~LearnSGD()
  {
    for (size_type i = 0; i != workers_.size(); ++ i)
      queue_.push(update_type(0, 0));
    
    workers_.join_all();
  }

  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    //++ const_cast<size_type&>(epoch_);
    
    const double eta = eta0_ / (epoch_ + 1);
    
    if (lambda_ != 0.0)
      theta.scale_ *= 1.0 - eta * lambda_;

    // parallelize here...
    const_cast<counter_type&>(counter_).clear();
    
    for (size_type i = 0; i != workers_.size(); ++ i)
      const_cast<queue_type&>(queue_).push(update_type(&theta, &gradient));
    
#if 0
    embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
    for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
      update(iiter->first,
	     theta.embedding_input_,
	     iiter->second,
	     1.0 / gradient.count_,
	     theta.scale_,
	     false);    
    
    embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
    for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
      update(oiter->first,
	     theta.embedding_output_,
	     oiter->second,
	     1.0 / gradient.count_,
	     theta.scale_,
	     true);
#endif
    
    update(theta.Wc_, gradient.Wc_, 1.0 / gradient.count_, lambda_ != 0.0);
    update(theta.bc_, gradient.bc_, 1.0 / gradient.count_, false);
    
    update(theta.Wh_, gradient.Wh_, 1.0 / gradient.count_, lambda_ != 0.0);
    update(theta.bh_, gradient.bh_, 1.0 / gradient.count_, false);

    // wait...
    const_cast<counter_type&>(counter_).wait(workers_.size());
  }
  
  template <typename Theta, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool regularize=true) const
  {
    const double eta = eta0_ / (epoch_ + 1);

    if (regularize)
      theta *= 1.0 - eta * lambda_;
    
    theta -= eta * scale * g;
  }
  
  template <typename Theta, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const double theta_scale,
	      const bool bias_last=false) const
  {
    const double eta = eta0_ / (epoch_ + 1);

    if (bias_last) {
      const size_type rows = g.rows();
      
      theta.col(word.id()).block(0, 0, rows - 1, 1) -= (eta * scale /  theta_scale) * g.block(0, 0, rows - 1, 1);
      theta.col(word.id()).block(rows - 1, 0, 1, 1) -= eta * scale * g.block(rows - 1, 0, 1, 1);
    } else
      theta.col(word.id()) -= (eta * scale /  theta_scale) * g;
  }

  double lambda_;
  double eta0_;
  
  size_type epoch_;

  queue_type   queue_;
  counter_type counter_;
  boost::thread_group workers_;
};

typedef boost::filesystem::path path_type;

typedef cicada::Sentence sentence_type;
typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;

typedef Model    model_type;
typedef Unigram  unigram_type;

typedef uint64_t count_type;
typedef utils::unordered_map<unigram_type::word_type, count_type,
			     boost::hash<unigram_type::word_type>, std::equal_to<unigram_type::word_type>,
			     std::allocator<std::pair<const unigram_type::word_type, count_type> > >::type word_set_type;


static const size_t DEBUG_DOT  = 100000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type input_file;
path_type list_file;
path_type embedding_file;
path_type output_model_file;

int dimension_embedding = 32;
int dimension_hidden = 256;
int order = 5;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 1024;
int samples = 100;
int cutoff = 3;
double lambda = 0;
double eta0 = 0.1;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const sentence_set_type& sentences,
		  const unigram_type& unigram,
		  model_type& theta);
void read_data(const path_type& input_file,
	       const path_type& list_file,
	       sentence_set_type& sentences,
	       word_set_type& words);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
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
    
    sentence_set_type sentences;
    word_set_type     words;
    
    read_data(input_file, list_file, sentences, words);

    if (debug)
      std::cerr << "# of sentences: " << sentences.size() << std::endl
		<< "vocabulary: " << (words.size() - 1) << std::endl;
    
    unigram_type unigram(words.begin(), words.end());
    
    model_type theta(dimension_embedding, dimension_hidden, order, unigram, generator);
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, order, lambda, eta0, threads), sentences, unigram, theta);
      else
	learn_online(LearnSGD(lambda, eta0, threads), sentences, unigram, theta);
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
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_type;
  
  struct Counter
  {
    Counter() : counter(0) {}
  
    void increment()
    {
      utils::atomicop::fetch_and_add(counter, size_type(1));
    }
  
    void wait(size_type target)
    {
      for (;;) {
	for (int i = 0; i < 64; ++ i) {
	  if (counter == target)
	    return;
	  else
	    boost::thread::yield();
	}
      
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      }
    }

    void clear() { counter = 0; }
  
    volatile size_type counter;
  };
  typedef Counter counter_type;

  TaskAccumulate(const sentence_set_type& sentences,
		 const unigram_type& unigram,
		 const size_type& samples,
		 const model_type& theta,
		 queue_type& queue,
		 counter_type& counter)
    : sentences_(sentences),
      theta_(theta),
      queue_(queue),
      counter_(counter),
      ngram_(unigram, samples),
      gradient_(theta.dimension_embedding_, theta.dimension_hidden_, theta.order_),
      log_likelihood_() {}

  void operator()()
  {
    clear();

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    size_type sentence_id;
    for (;;) {
      queue_.pop(sentence_id);
      
      if (sentence_id == size_type(-1)) break;
      
      const sentence_type& sentence = sentences_[sentence_id];
      
      if (! sentence.empty())
	log_likelihood_ += ngram_.learn(sentence, theta_, gradient_, generator);
      
      counter_.increment();
    }
  }

  void clear()
  {
    gradient_.clear();
    log_likelihood_ = log_likelihood_type();
  }

  const sentence_set_type& sentences_;
  const model_type& theta_;
  
  queue_type&            queue_;
  counter_type&          counter_;
  
  ngram_type ngram_;

  gradient_type       gradient_;
  log_likelihood_type log_likelihood_;
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
		  const sentence_set_type& sentences,
		  const unigram_type& unigram,
		  model_type& theta)
{
  typedef TaskAccumulate task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  task_type::queue_type   mapper(256 * threads);
  task_type::counter_type reducer;
  
  task_set_type tasks(threads, task_type(sentences,
					 unigram,
					 samples,
					 theta,
					 mapper,
					 reducer));
  
  id_set_type ids(sentences.size());
  for (size_type i = 0; i != ids.size(); ++ i)
    ids[i] = i;
  
  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;

    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);

    id_set_type::const_iterator biter     = ids.begin();
    id_set_type::const_iterator biter_end = ids.end();

    task_type::log_likelihood_type log_likelihood;
    size_type num_text = 0;

    utils::resource start;
    
    while (biter < biter_end) {
      // clear gradients...
      for (size_type i = 0; i != tasks.size(); ++ i)
	tasks[i].clear();
      
      // clear reducer
      reducer.clear();
      
      // map bitexts
      id_set_type::const_iterator iter_end = std::min(biter + batch_size, biter_end);
      for (id_set_type::const_iterator iter = biter; iter != iter_end; ++ iter) {
	mapper.push(*iter);
	
	++ num_text;
	if (debug) {
	  if (num_text % DEBUG_DOT == 0)
	    std::cerr << '.';
	  if (num_text % DEBUG_LINE == 0)
	    std::cerr << '\n';
	}
      }
      
      // wait...
      reducer.wait(iter_end - biter);
      biter = iter_end;
      
      // merge gradients
      log_likelihood += tasks.front().log_likelihood_;
      for (size_type i = 1; i != tasks.size(); ++ i) {
	tasks.front().gradient_ += tasks[i].gradient_;
	log_likelihood += tasks[i].log_likelihood_;
      }
      
      // update model parameters
      learner(theta, tasks.front().gradient_);
    }

    utils::resource end;
    
    if (debug && ((num_text / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;

    if (debug)
      std::cerr << "log-likelihood: " << static_cast<double>(log_likelihood) << std::endl
		<< "perplexity: " << std::exp(- static_cast<double>(log_likelihood)) << std::endl;

    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle bitexts!
    std::random_shuffle(ids.begin(), ids.end());
  }

  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    mapper.push(size_type(-1));

  workers.join_all();

  // finalize model...
  theta.finalize();
}

struct Reader
{
  void operator()(const sentence_type& sentence)
  {
    typedef cicada::Vocab vocab_type;
    
    if (sentence.empty()) return;
    
    sentences_.push_back(sentence);
    
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter)
      ++ words_[*siter];
    
    ++ words_[vocab_type::EOS];
  }
  
  sentence_set_type sentences_;
  word_set_type     words_;
};

struct ReaderFile : public Reader
{
  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
  
  ReaderFile(queue_type& queue) : queue_(queue) {}
  
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
  
  ReaderLines(queue_type& queue) : queue_(queue) {}
  
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
	       sentence_set_type& sentences,
	       word_set_type& words)
{
  typedef cicada::Vocab vocab_type;

  sentences.clear();
  words.clear();
  
  if (! input_file.empty()) {
    if (input_file != "-" && ! boost::filesystem::exists(input_file))
      throw std::runtime_error("no input file? " + input_file.string());

    ReaderLines::queue_type queue;

    std::vector<ReaderLines> tasks(threads, ReaderLines(queue));
    
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
    
    // join data...
    size_t data_size = sentences.size();
    for (size_t i = 0; i != tasks.size(); ++ i)
      data_size += tasks[i].sentences_.size();
    
    sentences.reserve(data_size);

    for (size_t i = 0; i != tasks.size(); ++ i) {
      sentences.insert(sentences.end(), tasks[i].sentences_.begin(), tasks[i].sentences_.end());
      
      word_set_type::const_iterator witer_end = tasks[i].words_.end();
      for (word_set_type::const_iterator witer = tasks[i].words_.begin(); witer != witer_end; ++ witer)
	words[witer->first] += witer->second;
    }
  }
  
  if (! list_file.empty()) {
    if (list_file != "-" && ! boost::filesystem::exists(list_file))
	throw std::runtime_error("no list file? " + list_file.string());
    
    ReaderFile::queue_type queue;
    
    std::vector<ReaderFile> tasks(threads, ReaderFile(queue));
    
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
    
    // join data...
    size_t data_size = sentences.size();
    for (size_t i = 0; i != tasks.size(); ++ i)
      data_size += tasks[i].sentences_.size();
    
    sentences.reserve(data_size);
    
    for (size_t i = 0; i != tasks.size(); ++ i) {
      sentences.insert(sentences.end(), tasks[i].sentences_.begin(), tasks[i].sentences_.end());
      
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
    
    // enumerate sentences and replace by UNK
    
    sentence_set_type::iterator siter_end = sentences.end();
    for (sentence_set_type::iterator siter = sentences.begin(); siter != siter_end; ++ siter) {
      sentence_type::iterator iter_end = siter->end();
      for (sentence_type::iterator iter = siter->begin(); iter != iter_end; ++ iter)
	if (words.find(*iter) == words.end())
	  *iter = vocab_type::UNK;
    }
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
