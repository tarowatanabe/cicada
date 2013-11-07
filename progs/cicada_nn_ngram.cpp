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

#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>

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
  
  Model() : dimension_embedding_(0), dimension_hidden_(0), order_(0) {}  
  template <typename Unigram, typename Gen>
  Model(const size_type& dimension_embedding,
	const size_type& dimension_hidden,
	const int order,
	const Unigram& unigram,
	Gen& gen) 
    : dimension_embedding_(dimension_embedding),
      dimension_hidden_(dimension_hidden),
      order_(order)
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
  }
  
  
  template <typename Gen>
  struct randomize
  {
    randomize(Gen& gen) : gen_(gen) {}
    
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return boost::random::uniform_real_distribution<Tp>(-0.01, 0.01)(const_cast<Gen&>(gen_));
    }
    
    Gen& gen_;
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
    
    embedding_input_  = tensor_type::Zero(dimension_embedding_,     vocabulary_size).array().unaryExpr(randomize<Gen>(gen));
    embedding_output_ = tensor_type::Zero(dimension_embedding_ + 1, vocabulary_size).array().unaryExpr(randomize<Gen>(gen));
    
    words_ = unique_set_type(vocabulary_size, false);
    
    for (size_type pos = 0; pos != unigram.words_.size(); ++ pos)
      if (unigram.words_[pos] != vocab_type::EOS)
	words_[unigram.words_[pos].id()] = true;
    
    words_[vocab_type::BOS.id()] = false;
    words_[vocab_type::EOS.id()] = false;

    Wc_ = tensor_type::Zero(dimension_hidden_, dimension_embedding_ * (order - 1)).array().unaryExpr(randomize<Gen>(gen));
    bc_ = tensor_type::Zero(dimension_hidden_, 1).array().unaryExpr(randomize<Gen>(gen));
    
    Wh_ = tensor_type::Zero(dimension_embedding_, dimension_hidden_).array().unaryExpr(randomize<Gen>(gen));
    bh_ = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen));
  }

  void finalize()
  {
    embedding_input_.col(vocab_type::EPSILON.id()) = tensor_type::Zero(dimension_embedding_, 1);
    embedding_output_.col(vocab_type::EPSILON.id()) = tensor_type::Zero(dimension_embedding_ + 1, 1);

    const double factor = 1.0 / std::accumulate(words_.begin(), words_.end(), size_type(0));
    
    for (size_type pos = 0; pos != words_.size(); ++ pos) {
      embedding_input_.col(vocab_type::EPSILON.id())  += embedding_input_.col(pos) * factor;
      embedding_output_.col(vocab_type::EPSILON.id()) += embedding_output_.col(pos) * factor;
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
    
    write_embedding(rep.path("input.gz"),  rep.path("input.bin"), embedding_input_, vocab_type::BOS);
    write_embedding(rep.path("output.gz"), rep.path("output.bin"), embedding_output_, vocab_type::EOS);
    
    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
    
    write(rep.path("Wh.txt.gz"), rep.path("Wh.bin"), Wh_);
    write(rep.path("bh.txt.gz"), rep.path("bh.bin"), bh_);
    
    // vocabulary...
    word_type::write(rep.path("vocab"));
  }
  
  void write_embedding(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix, const word_type& add) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    {
      utils::compress_ostream os(path_text, 1024 * 1024);
      std::ostream_iterator<char> iter(os);

      const word_type::id_type id_max = matrix.cols();
      for (word_type::id_type i = 0; i != id_max; ++ i) {
	const word_type word(i);
	
	if (! word.empty() && (words_[word.id()] || word == vocab_type::EPSILON || word == add)) {
	  karma::generate(iter, standard::string, word);
	  
	  for (difference_type j = 0; j != matrix.rows(); ++ j)
	    karma::generate(iter, karma::lit(' ') << float10, matrix(j, i));
	  
	  karma::generate(iter, karma::lit('\n'));
	}
      }
    }
    
    {
      utils::compress_ostream os(path_binary, 1024 * 1024);
      
      const tensor_type::Index rows = matrix.rows();
      const tensor_type::Index cols = matrix.cols();
      
      os.write((char*) &rows, sizeof(tensor_type::Index));
      os.write((char*) &cols, sizeof(tensor_type::Index));
      os.write((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
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
      
      os.write((char*) &rows, sizeof(tensor_type::Index));
      os.write((char*) &cols, sizeof(tensor_type::Index));
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

  unique_set_type words_;
  
  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;
  
  // Wh and bh for hidden layer
  tensor_type Wh_;
  tensor_type bh_;
};

struct Unigram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef std::vector<double, std::allocator<double> > logprob_set_type;
  typedef std::vector<word_type, std::allocator<word_type> > word_map_type;

  typedef boost::random::discrete_distribution<> distribution_type;

  template <typename Iterator>
  Unigram(Iterator first, Iterator last)
  {
    typedef utils::indexed_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> > word_set_type;
    typedef std::vector<size_type, std::allocator<size_type> > count_set_type;
    
    word_set_type  uniques;

    word_set_type::iterator iter = uniques.insert(vocab_type::EOS).first;
    const int index_eos = (iter - uniques.begin());
    
    count_set_type counts(index_eos + 1);
    
    for (/**/; first != last; ++ first) {
      accumulate(first->begin(), first->end(), uniques, counts);
      
      ++ counts[index_eos];
    }
    
    logprobs_.reserve(word_type::allocated());
    logprobs_.resize(word_type::allocated());
    words_.reserve(uniques.size());
    words_.resize(uniques.size());
    
    // initialize logprobs and words
    const double norm = 1.0 / std::accumulate(counts.begin(), counts.end(), double(0));
    for (word_type::id_type id = 0; id != counts.size(); ++ id) {
      words_[id] = uniques[id];
      logprobs_[words_[id].id()] = std::log(norm * counts[id]);
    }
    
    // initialize distribution
    distribution_ = distribution_type(counts.begin(), counts.end());
  }
  
  template <typename Iterator, typename Uniques, typename Counts>
  void accumulate(Iterator first, Iterator last, Uniques& uniques, Counts& counts)
  {
    for (/**/; first != last; ++ first) {
      typename Uniques::iterator iter = uniques.insert(*first).first;
      const int index = iter - uniques.begin();
      
      if (index >= counts.size())
	counts.resize(index + 1);
      
      ++ counts[index];
    }
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

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type    tensor_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef std::vector<tensor_type*, std::allocator<tensor_type*> > gradient_embedding_type;

  NGram(const unigram_type& unigram,
	const size_type samples)
    : unigram_(unigram), samples_(samples), log_samples_(std::log(samples)) {}

  const unigram_type& unigram_;
  size_type           samples_;
  double              log_samples_;
  
  tensor_type layer_input_;
  tensor_type layer_context_;
  tensor_type layer_hidden_;

  tensor_type delta_input_;
  tensor_type delta_context_;
  tensor_type delta_hidden_;

  gradient_embedding_type gradient_embedding_;
  
  struct hinge
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return std::max(x, Tp(0));
    }
  };
  
  struct dhinge
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return x > Tp(0);
    }
  };
  
  template <typename Gen>
  double learn(const sentence_type& sentence,
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
    
    for (size_type i = 0; i < order - 1; ++ i) {
      layer_input_.block(dimension * i, 0, dimension, 1) = theta.embedding_input_.col(vocab_type::BOS.id());
      
      gradient_embedding_[i] = &gradient_embedding_bos;
    }
    
    double log_likelihood = 0.0;
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
      log_likelihood += learn(*siter, theta, gradient, gen);
      
      // shift layer_input...
      layer_input_.block(0, 0, dimension * (order - 2), 1) = layer_input_.block(dimension, 0, dimension * (order - 2), 1).eval();
      layer_input_.block(dimension * (order - 2), 0, dimension, 1) = theta.embedding_input_.col(siter->id());
      
      // shift context
      std::copy(gradient_embedding_.begin() + 1, gradient_embedding_.end(), gradient_embedding_.begin());
      gradient_embedding_.back() = &gradient.embedding_input(*siter);
    }
    
    // correct scoring
    log_likelihood += learn(vocab_type::EOS, theta, gradient, gen);
    
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
    
    const double score = (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1).transpose() * layer_hidden_
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
		     * (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1) * loss).array());
    
    for (size_type k = 0; k != samples_; ++ k) {
      const word_type word = unigram_.draw(gen);

      const double score = (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1).transpose() * layer_hidden_
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
				* (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1) * loss).array());
    }
    
    gradient.Wh_ += delta_hidden_ * layer_context_.transpose();
    gradient.bh_ += delta_hidden_;
    
    delta_context_ = (layer_context_.array().unaryExpr(dhinge()) * (theta.Wh_.transpose() * delta_hidden_).array());
    
    gradient.Wc_ += delta_context_ * layer_input_.transpose();
    gradient.bc_ += delta_context_;
    
    delta_input_ = theta.Wc_.transpose() * delta_context_;
    
    // finally, input embedding...
    for (int i = 0; i != order - 1; ++ i)
      *gradient_embedding_[i] += delta_input_.block(dimension * i, 0, dimension, 1);

    return log_likelihood;
  }
};

struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef cicada::Symbol   word_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& dimension_embedding,
	       const size_type& dimension_hidden,
	       const int order,
	       const double& lambda,
	       const double& eta0)
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
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
    for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
      update(iiter->first,
	     theta.embedding_input_,
	     const_cast<tensor_type&>(embedding_input_),
	     iiter->second,
	     lambda_ != 0.0,
	     false);

    embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
    for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
      update(oiter->first,
	     theta.embedding_output_,
	     const_cast<tensor_type&>(embedding_output_),
	     oiter->second,
	     lambda_ != 0.0,
	     true);
    
    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, false);

    update(theta.Wh_, const_cast<tensor_type&>(Wh_), gradient.Wh_, lambda_ != 0.0);
    update(theta.bh_, const_cast<tensor_type&>(bh_), gradient.bh_, false);
  }

  template <typename Theta, typename GradVar, typename Grad>
  struct update_visitor_regularize
  {
    update_visitor_regularize(Eigen::MatrixBase<Theta>& theta,
			      Eigen::MatrixBase<GradVar>& G,
			      const Eigen::MatrixBase<Grad>& g,
			      const double& lambda,
			      const double& eta0)
      : theta_(theta), G_(G), g_(g), lambda_(lambda), eta0_(eta0) {}
    
    void init(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      operator()(value, i, j);
    }
    
    void operator()(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      G_(i, j) = std::max(std::min(G_(i, j) + g_(i, j) * g_(i, j), tensor_type::Scalar(1e+30)), tensor_type::Scalar(1e-30));
      
      const double rate = eta0_ / std::sqrt(G_(i, j));
      const double f = theta_(i, j) - rate * g_(i, j);

      theta_(i, j) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
    }
    
    Eigen::MatrixBase<Theta>&      theta_;
    Eigen::MatrixBase<GradVar>&    G_;
    const Eigen::MatrixBase<Grad>& g_;
    
    const double lambda_;
    const double eta0_;
  };
  
  template <typename Theta, typename GradVar, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const bool regularize=true) const
  {
    if (regularize) {
      update_visitor_regularize<Theta, GradVar, Grad> visitor(theta, G, g, lambda_, eta0_);
      
      theta.visit(visitor);
    } else {
      G.array() = (G.array() + g.array().square()).min(1e+30).max(1e-30);
      theta.array() -= eta0_ * g.array() / G.array().sqrt();
    }
  }

  template <typename Theta, typename GradVar, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const bool regularize=true,
	      const bool bias_last=false) const
  {
    if (regularize) {
      for (int row = 0; row != g.rows() - bias_last; ++ row) {
	G(row, word.id()) = std::max(std::min(G(row, word.id()) + g(row, 0) * g(row, 0), tensor_type::Scalar(1e+30)), tensor_type::Scalar(1e-30));
	
	const double rate = eta0_ / std::sqrt(G(row, word.id()));
	const double f = theta(row, word.id()) - rate * g(row, 0);
	
	theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
      }
      
      if (bias_last) {
	const int row = g.rows() - 1;
	
	G(row, word.id()) = std::max(std::min(G(row, word.id()) + g(row, 0) * g(row, 0), tensor_type::Scalar(1e+30)), tensor_type::Scalar(1e-30));
	theta(row, word.id()) -= eta0_ * g(row, 0) / std::sqrt(G(row, word.id()));
      }
    } else {
      G.col(word.id()).array() = (G.col(word.id()).array() + g.array().square()).min(1e+30).max(1e-30);
      theta.col(word.id()).array() -= eta0_ * g.array() / G.col(word.id()).array().sqrt();
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
};


typedef boost::filesystem::path path_type;

typedef cicada::Sentence sentence_type;
typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;

typedef Model    model_type;
typedef Unigram  unigram_type;

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type input_file;
path_type embedding_file;
path_type output_model_file;

int dimension_embedding = 16;
int dimension_hidden = 128;
int order = 5;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 1024;
int samples = 100;
double lambda = 1e-5;
double eta0 = 1;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const sentence_set_type& sentences,
		  const unigram_type& unigram,
		  model_type& theta);
void read_data(const path_type& input_file,
	       sentence_set_type& sentences);

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
      optimize_adagrad = true;
    
    threads = utils::bithack::max(threads, 1);
    
    // srand is used in Eigen
    std::srand(utils::random_seed());
  
    // this is optional, but safe to set this
    ::srandom(utils::random_seed());
        
    boost::mt19937 generator;
    generator.seed(utils::random_seed());

    if (input_file.empty())
      throw std::runtime_error("no data?");
    
    sentence_set_type sentences;
    
    read_data(input_file, sentences);
    
    unigram_type unigram(sentences.begin(), sentences.end());
    
    model_type theta(dimension_embedding, dimension_hidden, order, unigram, generator);
    
    if (iteration > 0)
      learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, order, lambda, eta0), sentences, unigram, theta);
    
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
      log_likelihood_(0),
      samples_(0),
      words_(0) {}

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
      
      if (! sentence.empty()) {
	log_likelihood_ += ngram_.learn(sentence, theta_, gradient_, generator);
	++ samples_;
	words_ += sentence.size() + 1;
      }
      
      counter_.increment();
    }
  }

  void clear()
  {
    gradient_.clear();
    log_likelihood_ = 0;
    samples_ = 0;
    words_ = 0;
  }

  const sentence_set_type& sentences_;
  const model_type& theta_;
  
  queue_type&            queue_;
  counter_type&          counter_;
  
  ngram_type ngram_;

  gradient_type gradient_;
  double        log_likelihood_;
  size_type     samples_;
  size_type     words_;
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

    double log_likelihood = 0.0;
    size_type samples = 0;
    size_type words = 0;
    size_type num_text = 0;
    
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
      samples        += tasks.front().samples_;
      words          += tasks.front().words_;
      for (size_type i = 1; i != tasks.size(); ++ i) {
	tasks.front().gradient_ += tasks[i].gradient_;
	log_likelihood += tasks[i].log_likelihood_;
	samples        += tasks[i].samples_;
	words          += tasks[i].words_;
      }
      
      // update model parameters
      learner(theta, tasks.front().gradient_);
    }
    
    if (debug && ((num_text / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;
    
    if (debug)
      std::cerr << "log_likelihood (per sentence): " << (log_likelihood / samples) << std::endl
		<< "log_likelihood (per word): "     << (log_likelihood / words) << std::endl
		<< "# of sentences: " << samples << std::endl
		<< "# of words: " << words << std::endl;
    
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


void read_data(const path_type& input_file,
	       sentence_set_type& sentences)
{
  sentences.clear();
  
  utils::compress_istream is(input_file, 1024 * 1024);
  
  sentence_type sentence;
  
  while (is >> sentence)
    sentences.push_back(sentence);

  sentence_set_type(sentences).swap(sentences);
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("input",    po::value<path_type>(&input_file),    "input file")
    
    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("order",     po::value<int>(&order)->default_value(order),         "context order size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("samples",           po::value<int>(&samples)->default_value(samples),       "# of NCE samples")
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
