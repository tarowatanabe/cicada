//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// an implementation for neural network alignment model
//
// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by AdaGrad (default)
//
// an implementation of NCE estimate for hiden markov alignment model...
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
#include "cicada/bitext.hpp"

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

struct Gradient
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;
  
  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;
  
  typedef utils::unordered_map<word_type, tensor_type,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, tensor_type> > >::type embedding_type;
  
  Gradient() : embedding_(0), hidden_(0), alignment_(0) {}
  Gradient(const size_type& embedding,
	   const size_type& hidden,
	   const int alignemnt) 
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment)
  { initialize(embedding, hidden, alignment); }
  
  Gradient& operator-=(const Gradient& x)
  {
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      tensor_type& embedding = source_[siter->first];
      
      if (! embedding.rows())
	embedding = -siter->second;
      else
	embedding -= siter->second;
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      tensor_type& embedding = target_[titer->first];
      
      if (! embedding.rows())
	embedding = -titer->second;
      else
	embedding -= titer->second;
    }

    if (! Wt_.rows())
      Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    if (! bt_.rows())
      bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());

    if (! Wa_.rows())
      Wa_ = tensor_type::Zero(x.Wa_.rows(), x.Wa_.cols());
    if (! ba_.rows())
      ba_ = tensor_type::Zero(x.ba_.rows(), x.ba_.cols());

    if (! Wn_.rows())
      Wn_ = tensor_type::Zero(x.Wn_.rows(), x.Wn_.cols());
    if (! ba_.rows())
      bn_ = tensor_type::Zero(x.bn_.rows(), x.bn_.cols());

    if (! Wi_.rows())
      Wi_ = tensor_type::Zero(x.Wi_.rows(), x.Wi_.cols());
    
    Wt_ -= x.Wt_;
    bt_ -= x.bt_;
    
    Wa_ -= x.Wa_;
    ba_ -= x.ba_;

    Wn_ -= x.Wn_;
    bn_ -= x.bn_;
    
    Wi_ -= x.Wi_;

    return *this;
  }
  
  Gradient& operator+=(const Gradient& x)
  {
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      tensor_type& embedding = source_[siter->first];
      
      if (! embedding.rows())
	embedding = siter->second;
      else
	embedding += siter->second;
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      tensor_type& embedding = target_[titer->first];
      
      if (! embedding.rows())
	embedding = titer->second;
      else
	embedding += titer->second;
    }

    if (! Wt_.rows())
      Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    if (! bt_.rows())
      bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());

    if (! Wa_.rows())
      Wa_ = tensor_type::Zero(x.Wa_.rows(), x.Wa_.cols());
    if (! ba_.rows())
      ba_ = tensor_type::Zero(x.ba_.rows(), x.ba_.cols());

    if (! Wn_.rows())
      Wn_ = tensor_type::Zero(x.Wn_.rows(), x.Wn_.cols());
    if (! bn_.rows())
      bn_ = tensor_type::Zero(x.bn_.rows(), x.bn_.cols()); 
   
    if (! Wi_.rows())
      Wi_ = tensor_type::Zero(x.Wi_.rows(), x.Wi_.cols());

    Wt_ += x.Wt_;
    bt_ += x.bt_;
    
    Wa_ += x.Wa_;
    ba_ += x.ba_;

    Wn_ += x.Wn_;
    bn_ += x.n_;
    
    Wi_ += x.Wi_;

    return *this;
  }
  
  void clear()
  {
    // embedding
    source_.clear();
    target_.clear();
    
    Wt_.setZero();
    bt_.setZero();
    
    Wa_.setZero();
    ba_.setZero();

    Wn_.setZero();
    bn_.setZero();
    
    Wi_.setZero();
  }

  
  tensor_type& source(const word_type& word)
  {
    tensor_type& embedding = source_[word];
    if (! embedding.cols())
      embedding = tensor_type::Zero(embedding_, 1);
    
    return embedding;
  }

  tensor_type& target(const word_type& word)
  {
    tensor_type& embedding = target_[word];
    if (! embedding.cols())
      embedding = tensor_type::Zero(embedding_ + 1, 1);
    
    return embedding;
  }

  
  void initialize(const size_type embedding, const size_type hidden, const int alignment)
  {
    if (hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (alignment <= 0)
      throw std::runtime_error("invalid alignment");
    
    embedding_ = embedding;
    hidden_    = hidden;
    alignment_ = alignment;
    
    clear();
    
    // initialize...
    Wt_ = tensor_type::Zero(embedding_, hidden_ + embedding_);
    bt_ = tensor_type::Zero(embedding_, 1);
    
    Wa_ = tensor_type::Zero(hidden_ * (alignemnt * 2 + 1), hidden_ + embedding_);
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1);

    Wn_ = tensor_type::Zero(hidden_, hidden_ + embedding_);
    bn_ = tensor_type::Zero(hidden_, 1);
    
    Wi_ = tensor_type::Zero(hidden_, 1);
  }
  
  // dimension...
  size_type embedding_;
  size_type hidden_;
  int alignment_;
  
  // embedding
  embedding_type source_;
  embedding_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;
  
  tensor_type Wa_;
  tensor_type ba_;

  tensor_type Wn_;
  tensor_type bn_;

  tensor_type Wi_;
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

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;
  
  typedef boost::filesystem::path path_type;

  
  Model() : embedding_(0), hidden_(0), alignment_(0) {}  
  template <typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	const int alignment,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment)
  { initialize(embedding, hidden, alignment, gen); }
  
  void clear()
  {
    // embedding
    source_.setZero();
    target_.setZero();
    
    Wt_.setZero();
    bt_.setZero();
    
    Wa_.setZero();
    ba_.setZero();

    Wn_.setZero();
    bn_.setZero();
    
    Wi_.setZero();
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
  
  template <typename Gen>
  void initialize(const size_type embedding,
		  const size_type hidden,
		  const int alignment,
		  Gen& gen)
  {
    if (hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (alignment <= 0)
      throw std::runtime_error("invalid alignment");
    
    embedding_ = embedding;
    hidden_    = hidden;
    alignment_ = alignment;
    
    clear();
    
    const size_type vocabulary_size = word_type::allocated();
    
    source_ = tensor_type::Zero(embedding_,     vocabulary_size).array().unaryExpr(randomize<Gen>(gen));
    target_ = tensor_type::Zero(embedding_ + 1, vocabulary_size).array().unaryExpr(randomize<Gen>(gen));
    
    Wt_ = tensor_type::Zero(embedding_, hidden_ + embedding_).array().unaryExpr(randomize<Gen>(gen));
    bt_ = tensor_type::Zero(embedding_, 1).array().unaryExpr(randomize<Gen>(gen));
    
    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), hidden_ + embedding_).array().unaryExpr(randomize<Gen>(gen));
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1).array().unaryExpr(randomize<Gen>(gen));

    Wn_ = tensor_type::Zero(hidden_, hidden_ + embedding_).array().unaryExpr(randomize<Gen>(gen));
    bn_ = tensor_type::Zero(hidden_, 1).array().unaryExpr(randomize<Gen>(gen));
    
    Wi_ = tensor_type::Zero(hidden_, 1).array().unaryExpr(randomize<Gen>(gen));
  }

  void finalize()
  {
    
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
    
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);    
    rep["alignment"] = utils::lexical_cast<std::string>(alignment_);    
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), target_);
    
    write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);
    
    write(rep.path("Wa.txt.gz"), rep.path("Wa.bin"), Wa_);
    write(rep.path("ba.txt.gz"), rep.path("ba.bin"), ba_);

    write(rep.path("Wn.txt.gz"), rep.path("Wn.bin"), Wn_);
    write(rep.path("bn.txt.gz"), rep.path("bn.bin"), bn_);
    
    write(rep.path("Wi.txt.gz"), rep.path("Wi.bin"), Wi_);
    
    // vocabulary...
    word_type::write(rep.path("vocab"));
  }
  
  void write_embedding(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    const word_type::id_type rows = matrix.rows();
    const word_type::id_type cols = matrix.cols();
    
    utils::compress_ostream os_txt(path_text, 1024 * 1024);
    utils::compress_ostream os_bin(path_binary, 1024 * 1024);
    std::ostream_iterator<char> iter(os_txt);
    
    for (word_type::id_type id = 0; id != cols; ++ id)  {
      const word_type word(id);
      
      if (word != word_type()) {
	karma::generate(iter, standard::string, word);
	
	for (difference_type j = 0; j != rows; ++ j)
	  karma::generate(iter, karma::lit(' ') << float10, matrix(j, id));
	
	karma::generate(iter, karma::lit('\n'));
      }
    }
    
    os_bin.write((char*) matrix.data(), sizeof(tensor_type::Scalar) * cols * rows);
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
  size_type embedding_;
  size_type hidden_;
  int       alignemnt_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;
  
  tensor_type Wa_;
  tensor_type ba_;

  tensor_type Wn_;
  tensor_type bn_;

  tensor_type Wi_;
};

struct Unigram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef uint64_t count_type;

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;

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

struct HMM
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;
  typedef Unigram  unigram_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type    tensor_type;

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;
  
  typedef std::vector<tensor_type*, std::allocator<tensor_type*> > gradient_embedding_type;
  
  
  HMM(const unigram_type& source,
      const unigram_type& target,
      const size_type samples)
    : source_(source), target_(target), samples_(samples), log_samples_(std::log(samples)) {}
  
  const unigram_type& source_;
  const unigram_type& target_;
  size_type           samples_;
  double              log_samples_;
  
  tensor_type alpha_;
  tensor_type trans_;
  back_type   back_;
  tensor_type accu_;
  tensor_type score_;

  tensor_type layer_alpha_;
  tensor_type layer_alpha0_;
  tensor_type layer_trans_;

  tensor_type delta_beta_;
  tensor_type delta_alpha_;
  tensor_type delta_trans_;
  
  gradient_embedding_type gradient_embedding_;
  
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const model_type& theta)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    const size_type state_size = theta.embedding_ + theta.hidden_;

    const parameter_type zero = - std::numeric_limits<parameter_type>::infinity();
    
    alpha_.resize((source_size + 2) * 2 * state_size, target_size + 2);
    trans_.resize((source_size + 2) * 2 * theta.embedding_, target_size + 2);
    
    score_.resize((source_size + 2) * 2, target_size + 2);
    score_.setConstant(zero);

    accu_.resize((source_size + 2) * 2, target_size + 2);
    accu_.setConstant(zero);
    
    back_.resize((source_size + 2) * 2, target_size + 2);
    back_.setConstant(- 1);
    
    alpha_.block(0, 0, theta.embedding_, 1) = theta.source_.col(vocab_type::BOS.id());
    alpha_.block(theta.embedding_, 0, theta.hidden_, 1) = theta.Wi_;
    accu_(0, 0) = 0.0;

    layer_alpha0_.block(0, 0, theta.embedding_, 1) = theta_.source_.col(vocab_type::NONE.id());
    
    for (size_type trg = 1; trg != target_size + 2; ++ trg) {
      const word_type target_next = (trg == 0
				     ? vocab_type::BOS
				     : (trg == target_size + 1
					? vocab_type::EOS
					: target[trg - 1]));
      
      // for EOS, we need to match with EOS, otherwise, do not match..
      const size_type next_first = utils::bithack::branch(trg == target_size + 1, source_size + 1, size_type(1));
      const size_type next_last  = utils::bithack::branch(trg == target_size + 1, source_size + 2, source_size + 1);
      
      const size_type prev1_first = utils::bithack::branch(trg == 1, 0, 1);
      const size_type prev1_last  = utils::bithack::branch(trg == 1, size_type(1), source_size + 2 - 1);
      
      const size_type prev2_first = source_size + 2 + (source_size + 2 - 1) * (trg == 1);
      const size_type prev2_last  = soruce_size + 2 + source_size + 2 - 1;
      
      const size_type none_first = utils::bithack::branch(trg == target_size + 1, source_size + 2 - 1, 0);
      const size_type none_last  = source_size + 2 - 1;
      
      // alignment into aligned...
      for (size_type next = next_first; next != next_last; ++ next) {
	const word_type source_next = (next == 0
				       ? vocab_type::BOS
				       : (next == source_size + 1
					  ? vocab_type::EOS
					  : source[next - 1]));
	
	layer_alpha_.block(0, 0, theta.embedding_, 1) = theta_.source_.col(source_next.id());

	for (size_type prev1 = prev1_first; prev1 != prev1_last; ++ prev1) {
	  const size_type shift = utils::bithack::min(utils::bithack::max((difference_type(next)
									   - difference_type(prev1))
									  + theta.alignment_,
									  difference_type(0)),
						      difference_type(theta.alignment_ * 2));
	  
	  layer_alpha_.block(theta.embedding_, 0, theta.hidden_, 1)
	    = (theta.Wa_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size)
	       * alpha_.block(state_size * prev1, trg - 1, state_size, 1)
	       + theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)).array().unaryExpr(hinge());
	  
	  layer_trans_ = (theta.Wt_ * layer_alpha_ + theta.bt_).array().unaryExpr(hinge());
	  
	  const parameter_type score
	    = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_
	       + theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  if (score + accu_(prev1, trg - 1) > accu_(next, trg)) {
	    alpha_.block(state_size * next, trg, state_size, 1) = layer_alpha_;
	    trans_.block(theta.embedding_ * next, trg, theta.embedding_, 1) = layer_trans_;
	    back_(next, trg) = prev1;
	    accu_(next, trg) = score + accu_(prev1, trg - 1);
	    score_(next, trg) = score;
	  }
	}
	
	for (size_type prev2 = prev2_first; prev2 != prev2_last; ++ prev2) {
	  const size_type shift = utils::bithack::min(utils::bithack::max((difference_type(next) + source_size + 2
									   - difference_type(prev1))
									  + theta.alignment_,
									  difference_type(0)),
						      difference_type(theta.alignment_ * 2));
	  
	  layer_alpha_.block(theta.embedding_, 0, theta.hidden_, 1)
	    = (theta.Wa_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size)
	       * alpha_.block(state_size * prev1, trg - 1, state_size, 1)
	       + theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)).array().unaryExpr(hinge());
	  
	  layer_trans_ = (theta.Wt_ * layer_alpha_ + theta.bt_).array().unaryExpr(hinge());
	  
	  const parameter_type score
	    = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_
	       + theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  if (score + accu_(prev2, trg - 1) > accu_(next, trg)) {
	    alpha_.block(state_size * next, trg, state_size, 1) = layer_alpha_;
	    trans_.block(theta.embedding_ * next, trg, theta.embedding_, 1) = layer_trans_;
	    back_(next, trg) = prev2;
	    accu_(next, trg) = score + accu_(prev2, trg - 1);
	    score_(next, trg) = score;
	  }
	}
      }
      
      // alingment into none..
      for (size_type none = none_first; none != none_last; ++ none) {
	const size_type next  = none + source_size + 2;
	const size_type prev1 = none;
	const size_type prev2 = none + source_size + 2;
	
	if (accu_(prev1, trg - 1) > zero) {
	  layer_alpha0_.block(theta.embedding_, 0, theta.hidden_, 1)
	    = (theta.Wn_ * alpha_.block(state_size * prev1, trg - 1, state_size, 1) + theta.bn_).array().unaryExpr(hinge());
	  
	  layer_trans_ = (theta.Wt_ * layer_alpha0_ + theta.bt_).array().unaryExpr(hinge());
	  
	  const parameter_type score
	    = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_
	       + theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  if (score + accu_(prev1, trg - 1)> accu_(next, trg)) {
	    alpha_.block(state_size * next, trg, state_size, 1) = layer_alpha0_;
	    trans_.block(theta.embedding_ * next, trg, theta.embedding_, 1) = layer_trans_;
	    back_(next, trg) = prev1;
	    accu_(next, trg) = score + accu_(prev1, trg - 1);
	    score_(next, trg) = score;
	  }
	}
	
	if (accu_(prev2, trg - 1) > zero) {
	  layer_alpha0_.block(theta.embedding_, 0, theta.hidden_, 1)
	    = (theta.Wn_ * alpha_.block(state_size * prev2, trg - 1, state_size, 1) + theta.bn_).array().unaryExpr(hinge());
	  
	  layer_trans_ = (theta.Wt_ * layer_alpha0_ + theta.bt_).array().unaryExpr(hinge());
	  
	  const parameter_type score
	    = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_
	       + theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  if (score + accu_(prev2, trg - 1) > accu_(next, trg)) {
	    alpha_.block(state_size * next, trg, state_size, 1) = layer_alpha0_;
	    trans_.block(theta.embedding_ * next, trg, theta.embedding_, 1) = layer_trans_;
	    back_(next, trg) = prev2;
	    accu_(next, trg) = score + accu_(prev2, trg - 1);
	    score_(next, trg) = score;
	  }
	}
      }
    }
  }

  template <typename Gen>
  double backward(const sentence_type& source,
		  const sentence_type& target,
		  const model_type& theta,
		  gradient_type& gradient,
		  alignment_type& alignment,
		  Gen& gen)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    double log_likelihood = 0.0;
    
    // one-to-many alignment...
    alignment.clear();
    
    size_type src = source_size + 1;
    for (difference_type trg = target_size + 1; trg > 0; -- trg) {
      const size_type src_prev = back_(src, trg);
      
      log_likelihood += backward(source, target, theta, gradient, gen,
				 trg,
				 src,
				 src_prev);
      
      if (trg > 0 && trg <= target_size && src < source_size + 2)
	alignment.push_back(std::make_pair(src - 1, trg - 1));
      
      src = src_prev;
    }
    
    // propabate delta_beta_ to BOS and Wi_;
    gradient.source(vocab_type::BOS) += delta_beta_.block(0, 0, theta.embedding_, 1);
    Wi_ += delta_beta_.block(theta.embedding_, 0, theta.hidden_, 1);
    
    return log_likelihood;
  }
  
  template <typename Gen>
  double backward(const sentence_type& source,
		  const sentence_type& target,
		  const model_type& theta,
		  gradient_type& gradient,
		  alignment_type& alignment,
		  Gen& gen,
		  const size_type target_pos,
		  const size_type source_pos,
		  const size_type source_prev)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type state_size = theta.embedding_ + theta.hidden_;
    
    const bool is_none = (souorce_pos >= source_size + 2);
    
    // compute the error at the EOS...
    const word_type word_source(is_none
				? vocab_type::EPSILON
				: (source_pos == 0
				   ? vocab_type::BOS
				   : (source_pos == source_size + 1
				      ? vocab_type::EOS
				      : source[source_pos - 1])));
    const word_type word_target(target_pos == 0
				? vocab_type::BOS
				: (target_pos == target_size + 1
				   ? vocab_type::EOS
				   : target[target_pos - 1]));
    
    const double score = score_(source_pos, target_pos);
    const double score_noise = log_samples_ + target_.logprob(word_target);
    const double z = utils::mathop::logsum(score, score_noise);
    const double logprob = score - z;
    const double logprob_noise = score_noise - z;
    
    // we suffer loss...
    const double loss = - 1.0 + std::exp(logprob);
    double log_likelihood = logprob;
    
    tensor_type& dembedding = gradient.target(word_target);
    
    dembedding.block(0, 0, theta.embedding_, 1).array() += loss * trans_.block(theta.embedding_ * source_pos,
									       target_pos,
									       theta.embedding_,
									       1).array();
    dembedding.block(theta.embedding_, 0, 1, 1).array() += loss;
    
    delta_trans_ = (trans_.block(theta.embedding_ * source_pos,
				 target_pos,
				 theta.embedding_,
				 1).array().unaryExpr(dhinge())
		    * (theta.target_.col(word_target.id()).block(0, 0, theta.embedding_, 1) * loss).array());
    
    for (size_type k = 0; k != samples; ++ k) {
      const word_type word_target = target_.draw(gen);
      
      const double score = (theta.target_.col(word_target.id()).block(0, 0, theta.embedding_, 1).transpose()
			    * trans_.block(theta.embedding_ * source_pos,
					   target_pos,
					   theta.embedding_,
					   1)
			    + theta.target_.col(word_target.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
      const double score_noise = log_samples_ + target_.logprob(word_target);
      const double z = utils::mathop::logsum(score, score_noise);
      const double logprob = score - z;
      const double logprob_noise = score_noise - z;
      
      const double loss = std::exp(logprob);
      log_likelihood += logprob_noise;
      
      tensor_type& dembedding = gradient.target(word_target);
      
      dembedding.block(0, 0, theta.embedding_, 1).array() += loss * trans_.block(theta.embedding_ * source_pos,
										 target_pos,
										 theta.embedding_,
										 1).array();
      dembedding.block(theta.embedding_, 0, 1, 1).array() += loss;
      
      delta_trans_.array() += (trans_.block(theta.embedding_ * source_pos,
					    target_pos,
					    theta.embedding_,
					    1).array().unaryExpr(dhinge())
			       * (theta.target_.col(word_target.id()).block(0, 0, theta.embedding_, 1) * loss).array());
    }
    
    gradient.Wt_ += delta_trans_ * alpha_.block(state_size * source_pos, target_pos, state_size, 1).transpose();
    gradient.bt_ += delta_trans_;
    
    delta_alpha_ = (alpha_.block(state_size * source_pos, target_pos, state_size, 1).array().unaryExpr(dhinge())
		    * (delta_beta_ // propagate backward!
		       + theta.Wt_.transpose() * delta_trans_).array());
    
    if (is_none) {
      // word embedding...
      gradient.source(vocab_type::NONE) += delta_alpha_.block(0, 0, theta.embedding_, 1);

      gradient.Wn_ += delta_alpha_ * alpha_.block(state_size * source_prev, target_pos - 1, state_size, 1).transpose();
      gradient.bn_ += delta_alpha_;
      
      // this will be propagated backward!
      delta_beta_ = theta.Wn_.transpose() * delta_alpha_;
    } else {
      // accumulate for word embedding
      gradient.source(word_source) += delta_alpha_.block(0, 0, theta.embedding_, 1);
      
      const size_type shift = utils::bithack::min(utils::bithack::max((difference_type(source_pos)
								       - (difference_type(source_prev)
									  - utils::bithack::branch(source_prev >= source_size + 2,
												   source_size + 2,
												   size_type(0))))
								      + theta.alignment_,
								      difference_type(0)),
						  difference_type(theta.alignment_ * 2));
      
      egradient.Wa_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size)
	+= delta_alpha_ * alpha_.block(state_size * source_prev, target_pos - 1, state_size, 1).transpose();
      gradient.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
	+= delta_alpha_;
      
      // this will be propagated backward!
      delta_beta_ = theta.Wn_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size).transpose() * delta_alpha_;
    }
    
    return log_likelihoo;
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
  
  LearnAdaGrad(const size_type& embedding,
	       const size_type& hidden,
	       const int alignment,
	       const double& lambda,
	       const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment),
      lambda_(lambda),
      eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    
    source_ = tensor_type::Zero(embedding_,     vocabulary_size);
    target_ = tensor_type::Zero(embedding_ + 1, vocabulary_size);
    
    Wt_ = tensor_type::Zero(embedding_, hidden_ + embedding_);
    bt_ = tensor_type::Zero(embedding_, 1);
    
    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), hidden_ + embedding_);
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1);
    
    Wn_ = tensor_type::Zero(hidden_, hidden_ + embedding_);
    bn_ = tensor_type::Zero(hidden_, 1);
    
    Wi_ = tensor_type::Zero(hidden_, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     const_cast<tensor_type&>(source_),
	     siter->second,
	     lambda_ != 0.0,
	     false);
    
    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     const_cast<tensor_type&>(target_),
	     titer->second,
	     lambda_ != 0.0,
	     true);
    
    update(theta.Wt_, const_cast<tensor_type&>(Wt_), gradient.Wt_, lambda_ != 0.0);
    update(theta.bt_, const_cast<tensor_type&>(bt_), gradient.bt_, false);
    
    update(theta.Wa_, const_cast<tensor_type&>(Wa_), gradient.Wa_, lambda_ != 0.0);
    update(theta.ba_, const_cast<tensor_type&>(ba_), gradient.ba_, false);

    update(theta.Wn_, const_cast<tensor_type&>(Wn_), gradient.Wn_, lambda_ != 0.0);
    update(theta.bn_, const_cast<tensor_type&>(bn_), gradient.bn_, false);
    
    update(theta.Wi_, const_cast<tensor_type&>(Wi_), gradient.Wi_, lambda_ != 0.0);
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
      if (g_(i, j) == 0) return;
      
      G_(i, j) += g_(i, j) * g_(i, j);
      
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

  struct learning_rate
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return (x == 0.0 ? 0.0 : 1.0 / std::sqrt(x));
    }
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
      G.array() += g.array().square();
      theta.array() -= eta0_ * g.array() * G.array().unaryExpr(learning_rate());
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
      for (int row = 0; row != g.rows() - bias_last; ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0);
	  
	  const double rate = eta0_ / std::sqrt(G(row, word.id()));
	  const double f = theta(row, word.id()) - rate * g(row, 0);
	  
	  theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
      
      if (bias_last) {
	const int row = g.rows() - 1;
	
	if (g(row, 0) != 0) {
	  G(row, word.id()) += g(row, 0) * g(row, 0);
	  theta(row, word.id()) -= eta0_ * g(row, 0) / std::sqrt(G(row, word.id()));
	}
      }
    } else {
      G.col(word.id()).array() += g.array().square();
      theta.col(word.id()).array() -= eta0_ * g.array() * G.col(word.id()).array().unaryExpr(learning_rate());
    }
  }
  
  size_type embedding_;
  size_type hidden_;
  int       alignment_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;
  
  tensor_type Wa_;
  tensor_type ba_;

  tensor_type Wn_;
  tensor_type bn_;

  tensor_type Wi_;
};


typedef boost::filesystem::path path_type;

typedef cicada::Sentence sentence_type;
typedef cicada::Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;

typedef Model    model_type;
typedef Unigram  unigram_type;

typedef uint64_t count_type;
typedef utils::unordered_map<unigram_type::word_type, count_type,
			     boost::hash<unigram_type::word_type>, std::equal_to<unigram_type::word_type>,
			     std::allocator<std::pair<const unigram_type::word_type, count_type> > >::type word_set_type;


static const size_t DEBUG_DOT  = 100000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type source_file;
path_type target_file;
path_type output_source_target_file;
path_type output_target_source_file;
path_type alignment_source_target_file;
path_type alignment_target_source_file;

int dimension_embedding = 16;
int dimension_hidden = 128;
int alignment = 5;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 128;
int samples = 100;
int cutoff = 2;
double lambda = 1e-5;
double eta0 = 1;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const unigram_type& sources,
		  const unigram_type& targets,
		  model_type& theta_source_target,
		  model_type& theta_target_source);
void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       unigram_type& sources,
	       unigram_type& targets);

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

    if (source_file.empty() || target_file.empty())
      throw std::runtime_error("no data?");

    if (source_file != "-" && ! boost::filesystem::exists(source_file))
      throw std::runtime_error("no source file? " + source_file.string());
    
    if (target_file != "-" && ! boost::filesystem::exists(target_file))
      throw std::runtime_error("no target file? " + target_file.string());
    
    bitext_set_type bitexts;
    word_set_type   sources;
    word_set_type   targets;
    
    read_data(source_file, target_file, bitexts, sources, targets);

    if (debug)
      std::cerr << "# of sentences: " << bitexts.size() << std::endl
		<< "source vocabulary: " << sources.size() << std::endl
		<< "target vocabulary: " << targets.size() << std::endl;
    
    unigram_type unigram_sources(sources.begin(), sources.end());
    unigram_type unigram_targets(targets.begin(), targets.end());
    
    model_type theta_source_target(dimension_embedding, dimension_hidden, alignment, generator);
    model_type theta_target_source(dimension_embedding, dimension_hidden, alignment, generator);
    
    if (iteration > 0)
      learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, alignment, lambda, eta0),
		   bitexts,
		   unigram_sources,
		   unigram_targets,
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
  typedef std::set<value_type, compare_value, std::allocator<value_type> > bitext_set_type;
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
      bitext_set_type bitexts;
      bitext_alignment_type bitext;
      size_type id = 0;
      
      utils::compress_ostream os(path_, 1024 * 1024);
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;

	if (bitext.id_ == id) {
	  os << bitext.alignemnt_ << '\n';
	  ++ id;
	} else
	  bitexts.insert(bitext);

	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  os << bitexts.begin()->alignment_ << '\n';
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }

      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	os << bitexts.begin()->alignment_ << '\n';
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while writing alignment output?");
    }
  }

  path_type   path_;
  queue_type& queue_;
};


struct TaskAccumulate
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef HMM hmm_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef OutputMapReduce output_map_reduce_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_type;

  typedef output_map_reduce_type::queu_type queue_alignment_type;
  typedef output_map_reduce_type::value_type bitext_alignment_type;
  
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

  TaskAccumulate(const bitext_set_type& bitexts,
		 const unigram_type& sources,
		 const unigram_type& targets,
		 const size_type& samples,
		 const model_type& theta_source_target,
		 const model_type& theta_target_source,
		 queue_type& queue,
		 queue_type& queue_source_target,
		 queue_type& queue_target_source,
		 counter_type& counter)
    : bitexts_(bitexts),
      theta_source_target_(theta_source_target),
      theta_target_source_(theta_target_source),
      queue_(queue),
      queue_source_target_(queue_source_target),
      queue_target_source_(queue_target_source),
      counter_(counter),
      hmm_source_target_(sources, targets, samples),
      hmm_target_source_(targets, sources, samples),
      gradient_source_target_(theta_source_target.dimension_embedding_,
			      theta_source_target.dimension_hidden_,
			      theta_source_target.alignment_),
      gradient_target_source_(theta_target_source.dimension_embedding_,
			      theta_target_source.dimension_hidden_,
			      theta_target_source.alignment_),
      log_likelihood_source_target_(0),
      log_likelihood_target_source_(0) {}

  void operator()()
  {
    clear();

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    bitext_alignment_type bitext_source_target;
    bitext_alignment_type bitext_target_source;
    
    size_type sentence_id;
    for (;;) {
      queue_.pop(sentence_id);
      
      if (sentence_id == size_type(-1)) break;
      
      const bitext_type& bitext = bitexts_[sentence_id];
      
      bitext_source_target.bitext_.source_ = bitext.source_;
      bitext_source_target.bitext_.target_ = bitext.target_;
      bitext_source_target.alignemnt_.clear();
      
      bitext_target_source.bitext_.source_ = bitext.target_;
      bitext_target_source.bitext_.target_ = bitext.source_;
      bitext_target_source.alignment_.clear();

      if (! bitext.source_.empty() && ! bitext.target_.empty()) {
	hmm_source_target_.forward(bitext.source_, bitext.target_, theta_source_target_);
	hmm_target_source_.forward(bitext.target_, bitext.source_, theta_target_source_);
	
	log_likelihood_source_target_
	  += hmm_source_target_.backward(bitext.source_,
					 bitext.target_,
					 theta_source_target_,
					 gradient_source_target_,
					 bitext_source_target.alignment_,
					 generator);
	
	log_likelihood_target_source_
	  += hmm_target_source_.backward(bitext.target_,
					 bitext.source_,
					 theta_target_source_,
					 gradient_target_source_,
					 bitext_target_source.alignment_,
					 generator);
      }
      
      // reduce alignment
      queue_source_target_.push_swap(bitext_source_target);
      queue_target_source_.push_swap(bitext_target_source);
      
      // increment counter for synchronization
      counter_.increment();
    }
  }

  void clear()
  {
    gradient_.clear();
    log_likelihood_source_target_ = 0;
    log_likelihood_target_source_ = 0;
  }

  const bitext_set_type& bitexts_;
  const model_type& theta_source_target_;
  const model_type& theta_target_source_;
  
  queue_type&           queue_;
  queue_alignment_type& queue_source_target_;
  queue_alignment_type& queue_target_source_;
  counter_type&         counter_;
  
  hmm_type hmm_source_target_;
  hmm_type hmm_target_source_;
  
  gradient_type gradient_source_target_;
  gradient_type gradient_target_source_;
  double        log_likelihood_source_target_;
  double        log_likelihood_target_source_;
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
		  const bitext_set_type& bitexts,
		  const unigram_type& sources,
		  const unigram_type& targets,
		  model_type& theta_source_target,
		  model_type& theta_target_source)
{
  typedef TaskAccumulate task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputAlignment  output_alignment_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  task_type::queue_type   mapper(256 * threads);
  task_type::counter_type reducer;
  
  output_map_reduce_type::queue_type queue_source_target;
  output_map_reduce_type::queue_type queue_target_source;
  
  task_set_type tasks(threads, task_type(bitexts,
					 sources,
					 targets,
					 samples,
					 theta_source_target,
					 theta_target_source,
					 mapper,
					 queue_source_taget,
					 queue_target_source,
					 reducer));
  
  id_set_type ids(bitexts.size());
  for (size_type i = 0; i != ids.size(); ++ i)
    ids[i] = i;
  
  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;

    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);
    
    boost::thread output_source_target(output_alignment_type(! alignment_source_target_file.empty() && dump_mode
							     ? add_suffix(alignment_source_target_file, iter_tag)
							     : path_type(),
							     queue_source_target));
    boost::thread output_target_source(output_alignment_type(! alignment_target_source_file.empty() && dump_mode
							     ? add_suffix(alignment_target_source_file, iter_tag)
							     : path_type(),
							     queue_target_source));
    
    id_set_type::const_iterator biter     = ids.begin();
    id_set_type::const_iterator biter_end = ids.end();

    double log_likelihood_source_target = 0.0;
    double log_likelihood_target_source = 0.0;
    size_type samples = 0;
    size_type words_source = 0;
    size_type words_target = 0;
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
	if (! bitexts[*iter].source_.empty() && ! bitexts[*iter].target_.empty()) {
	  ++ samples;
	  words_source += bitexts[*iter].source_.size();
	  words_target += bitexts[*iter].target_.size();
	}
	
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
      log_likelihood_source_target += tasks.front().log_likelihood_source_target_;
      log_likelihood_target_source += tasks.front().log_likelihood_target_source_;
      for (size_type i = 1; i != tasks.size(); ++ i) {
	tasks.front().gradient_source_target_ += tasks[i].gradient_source_target_;
	tasks.front().gradient_target_source_ += tasks[i].gradient_target_source_;
	
	log_likelihood_source_target += tasks[i].log_likelihood_source_target_;
	log_likelihood_target_source += tasks[i].log_likelihood_target_source_;
      }
      
      // update model parameters
      learner(theta_source_target, tasks.front().gradient_source_target_);
      learner(theta_target_source, tasks.front().gradient_target_source_);
    }
    
    utils::resource end;
    
    queue_source_target.push(output_map_reduce_type::value_type());
    queue_target_source.push(output_map_reduce_type::value_type());
    
    if (debug && ((num_text / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;
    
    if (debug)
      std::cerr << "log_likelihood P(target | source) (per sentence): " << (log_likelihood_source_target / samples) << std::endl
		<< "log_likelihood P(target | source) (per word): "     << (log_likelihood_source_target / words_target) << std::endl
		<< "log_likelihood P(source | target) (per sentence): " << (log_likelihood_source_target / samples) << std::endl
		<< "log_likelihood P(source | target) (per word): "     << (log_likelihood_target_source / words_source) << std::endl;

    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle bitexts!
    std::random_shuffle(ids.begin(), ids.end());
    
    output_source_target.join();
    output_target_source.join();
  }

  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    mapper.push(size_type(-1));

  workers.join_all();

  // finalize model...
  theta.finalize();
}


void read_data(const path_type& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts,
	       word_set_type& sources,
	       word_set_type& targets)
{
  bitexts.clear();
  sources.clear();
  targets.clear();

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
    
    sentence_type::const_iterator siter_end = source.end();
    for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
      ++ sources[*siter];
    
    sentence_type::const_iterator titer_end = target.end();
    for (sentence_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
      ++ targets[*titer];    
  }
  
  if (src || trg)
    throw std::runtime_error("# of sentences does not match");
  
  if (cutoff > 1) {
    word_set_type sources_new;
    word_set_type targets_new;
    count_type source_unk = 0;
    count_type target_unk = 0;

    word_set_type::const_iterator siter_end = sources.end();
    for (word_set_type::const_iterator siter = sources.begin(); siter != siter_end; ++ siter)
      if (siter->second >= cutoff)
	sources_new.insert(*siter);
      else
	source_unk += siter->second;

    word_set_type::const_iterator titer_end = targets.end();
    for (word_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
      if (titer->second >= cutoff)
	targets_new.insert(*titer);
      else
	target_unk += titer->second;
    
    sources_new[vocab_type::UNK] = source_unk;
    targets_new[vocab_type::UNK] = target_unk;
    
    sources_new.swap(sources);
    targets_new.swap(targets);
    
    // enumerate sentences and replace by UNK

    bitext_set_type::iterator biter_end = bitexts.end();
    for (bitext_set_type::iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {
      
      sentence_type::const_iterator siter_end = biter->source_.end();
      for (sentence_type::const_iterator siter = biter->source_.begin(); siter != siter_end; ++ siter)
	if (sources.find(*suter) == sources.end())
	  *siter = vocab_type::UNK;
      
      sentence_type::const_iterator titer_end = biter->target_.end();
      for (sentence_type::const_iterator titer = biter->target_.begin(); titer != titer_end; ++ titer)
	if (targets.find(*suter) == targets.end())
	  *titer = vocab_type::UNK;
    }
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    
    ("output-source-target", po::value<path_type>(&output_source_target_file), "output model parameter for P(target | source)")
    ("output-target-source", po::value<path_type>(&output_target_source_file), "output model parameter for P(source | target)")

    ("alignment-source-target", po::value<path_type>(&alignment_source_target_file), "output alignment for P(target | source)")
    ("alignment-target-source", po::value<path_type>(&alignment_target_source_file), "output alignment for P(source | target)")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("alignment",           po::value<int>(&alignment)->default_value(alignemnt),                     "alignment model size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
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
