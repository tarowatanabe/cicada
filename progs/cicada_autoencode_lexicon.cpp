//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// 
// Lexical autoencoder
//
//  reconstruction
//    internal ------ sampled loss
//   embedding
//
// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by Pegasos
// SGD with L2 regularizer inspired by AdaGrad (default)
// SGD with L2/L2 regularizer from RDA (TODO)
//
// + batch algorithm using LBFGS (TODO)
//

#include <cstdlib>
#include <cmath>
#include <climits>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <deque>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"
#include "cicada/alignment.hpp"
#include "cicada/bitext.hpp"

#include "utils/alloc_vector.hpp"
#include "utils/lexical_cast.hpp"
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
#include "utils/sampler.hpp"
#include "utils/resource.hpp"

#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/progress.hpp>

struct Average
{
  Average() : average_(0), count_(0) {}
  Average(const double& x) : average_(x), count_(1) {}
  
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
  // model parameters
  
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef cicada::Bitext bitext_type;
  typedef cicada::Vocab vocab_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  
  typedef utils::unordered_map<word_type, tensor_type,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, tensor_type> > >::type embedding_type;

  typedef boost::filesystem::path path_type;
  
  Gradient() : dimension_embedding_(0), dimension_hidden_(0), window_(0), count_(0), shared_(0) {}
  Gradient(const size_type& dimension_embedding, const size_type& dimension_hidden, const size_type& window) 
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden),
      window_(window), count_(0), shared_(0) { initialize(dimension_embedding, dimension_hidden, window); }

  Gradient& operator-=(const Gradient& x)
  {
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      tensor_type& embedding = source_[siter->first];

      if (! embedding.rows())
	embedding = - siter->second;
      else
	embedding -= siter->second;
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      tensor_type& embedding = target_[titer->first];
      
      if (! embedding.rows())
	embedding = - titer->second;
      else
	embedding -= titer->second;
    }
        
    Wt1_ -= x.Wt1_;
    bt1_ -= x.bt1_;
    Wt2_ -= x.Wt2_;
    bt2_ -= x.bt2_;
    
    Wc_ -= x.Wc_;
    bc_ -= x.bc_;

    count_ -= x.count_;
    
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
    
    Wt1_ += x.Wt1_;
    bt1_ += x.bt1_;
    Wt2_ += x.Wt2_;
    bt2_ += x.bt2_;
    
    Wc_ += x.Wc_;
    bc_ += x.bc_;

    count_ += x.count_;
    
    return *this;
  }
  
  void clear()
  {
    // embedding
    source_.clear();
    target_.clear();
    
    // matrices
    Wt1_.setZero();
    bt1_.setZero();
    Wt2_.setZero();
    bt2_.setZero();    

    Wc_.setZero();
    bc_.setZero();

    count_ = 0;
    shared_ = 0;
  }
    
  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const size_type window)
  {
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    window_              = window;
    
    // embedding
    source_.clear();
    target_.clear();
        
    // lexicon
    Wt1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1));
    bt1_ = tensor_type::Zero(dimension_hidden, 1);
    
    // lexicon reconstruction
    Wt2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden);
    bt2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);

    // classification
    Wc_ = tensor_type::Zero(1, dimension_hidden);
    bc_ = tensor_type::Zero(1, 1);

    count_ = 0;
    shared_ = 0;
  }
  
public:
  void increment()
  {
    utils::atomicop::add_and_fetch(shared_, size_type(1));
  }

  size_type shared() const
  {
    const size_type ret = shared_;
    utils::atomicop::memory_barrier();
    return ret;
  }
  
  // dimension...
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type window_;
  
  // Embedding
  embedding_type source_;
  embedding_type target_;
  
  // Wt1 and bt1 for encoding
  tensor_type Wt1_;
  tensor_type bt1_;
  
  // Wt2 and bt2 for reconstruction
  tensor_type Wt2_;
  tensor_type bt2_;  

  // Wc and bc for classification
  tensor_type Wc_;
  tensor_type bc_;

  size_type count_;
  size_type shared_;
};

struct Model
{
  // model parameters
  
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef cicada::Bitext bitext_type;
  typedef cicada::Vocab vocab_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;

  typedef std::vector<bool, std::allocator<bool> > unique_set_type;

  typedef boost::filesystem::path path_type;
  
  Model() : dimension_embedding_(0), dimension_hidden_(0), window_(0), alpha_(0), beta_(0), scale_(1) {}
  template <typename Words, typename Gen>
  Model(const size_type& dimension_embedding, const size_type& dimension_hidden, const size_type& window,
	const double& alpha, const double& beta,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden),
      window_(window), alpha_(alpha), beta_(beta), scale_(1)
  {
    initialize(dimension_embedding, dimension_hidden, window, words_source, words_target, gen);
  }

  
  void clear()
  {
    // embedding
    source_.setZero();
    target_.setZero();
    
    // matrices
    Wt1_.setZero();
    bt1_.setZero();
    Wt2_.setZero();
    bt2_.setZero();    

    Wc_.setZero();
    bc_.setZero();

    scale_ = 1;
  }

  void finalize()
  {
    if (scale_ == 1) return;
    
    source_ *= scale_;
    target_ *= scale_;

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

  template <typename Words, typename Gen>
  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const size_type window,
		  Words& words_source,
		  Words& words_target,
		  Gen& gen)
  {
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    window_              = window;

    const size_type vocabulary_size = word_type::allocated();

    const double range_e = std::sqrt(6.0 / (dimension_embedding + 1));
    const double range_l = std::sqrt(6.0 / (dimension_hidden + dimension_embedding * 2 * (window * 2 + 1)));
    const double range_c = std::sqrt(6.0 / (dimension_hidden + 1));

    // intialize randomly...
    
    // embedding
    source_ = tensor_type::Zero(dimension_embedding, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_ = tensor_type::Zero(dimension_embedding, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    
    // lexicon
    Wt1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1)).array().unaryExpr(randomize<Gen>(gen, range_l));
    bt1_ = tensor_type::Zero(dimension_hidden, 1);
    
    // lexicon reconstruction
    Wt2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden).array().unaryExpr(randomize<Gen>(gen, range_l));
    bt2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);

    // classification
    Wc_ = tensor_type::Zero(1, dimension_hidden).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Ones(1, 1).array();

    scale_ = 1;

    
    words_source_.clear();
    words_target_.clear();
    words_source_.resize(vocabulary_size, false);
    words_target_.resize(vocabulary_size, false);

    words_source_[vocab_type::EPSILON.id()] = true;
    words_target_[vocab_type::EPSILON.id()] = true;
    words_source_[vocab_type::BOS.id()] = true;
    words_target_[vocab_type::BOS.id()] = true;
    words_source_[vocab_type::EOS.id()] = true;
    words_target_[vocab_type::EOS.id()] = true;

    for (typename Words::const_iterator siter = words_source.begin(); siter != words_source.end(); ++ siter)
      words_source_[siter->id()] = true;
    for (typename Words::const_iterator titer = words_target.begin(); titer != words_target.end(); ++ titer)
      words_target_[titer->id()] = true;
  }
  
  
  struct real_policy : boost::spirit::karma::real_policies<parameter_type>
  {
    static unsigned int precision(parameter_type)
    {
      return 10;
    }
  };
  
  void read_embedding(const path_type& source_file, const path_type& target_file)
  {
    // we will overwrite embedding... thus we will not clear embedding_
    
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;
    typedef boost::fusion::tuple<std::string, parameter_set_type > embedding_parsed_type;
    typedef boost::spirit::istream_iterator iterator_type;
    
    qi::rule<iterator_type, std::string(), standard::blank_type>           word;
    qi::rule<iterator_type, embedding_parsed_type(), standard::blank_type> parser; 
    
    word   %= qi::lexeme[+(standard::char_ - standard::space)];
    parser %= word >> *qi::double_ >> (qi::eol | qi::eoi);

    if (! source_file.empty()) {
      
      if (source_file != "-" && ! boost::filesystem::exists(source_file))
	throw std::runtime_error("no embedding: " + source_file.string());
      
      utils::compress_istream is(source_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
    
      iterator_type iter(is);
      iterator_type iter_end;
      
      embedding_parsed_type parsed;
      
      while (iter != iter_end) {
	boost::fusion::get<0>(parsed).clear();
	boost::fusion::get<1>(parsed).clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, parsed))
	  if (iter != iter_end)
	    throw std::runtime_error("embedding parsing failed");
	
	if (boost::fusion::get<1>(parsed).size() != dimension_embedding_)
	  throw std::runtime_error("invalid embedding size");

	const word_type word = boost::fusion::get<0>(parsed);
	
	if (word.id() >= source_.cols())
	  source_.conservativeResize(Eigen::NoChange, word.id() + 1);
	if (word.id() >= words_source_.size())
	  words_source_.resize(word.id() + 1, false);

	source_.col(word.id()).block(0, 0, dimension_embedding_, 1)
	  = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), dimension_embedding_, 1);

	words_source_[word.id()] = true;
      }
    }

    if (! target_file.empty()) {
      
      if (target_file != "-" && ! boost::filesystem::exists(target_file))
	throw std::runtime_error("no embedding: " + target_file.string());
      
      utils::compress_istream is(target_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
    
      iterator_type iter(is);
      iterator_type iter_end;
      
      embedding_parsed_type parsed;
      
      while (iter != iter_end) {
	boost::fusion::get<0>(parsed).clear();
	boost::fusion::get<1>(parsed).clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, parsed))
	  if (iter != iter_end)
	    throw std::runtime_error("embedding parsing failed");
	
	if (boost::fusion::get<1>(parsed).size() != dimension_embedding_)
	  throw std::runtime_error("invalid embedding size");
	
	const word_type word = boost::fusion::get<0>(parsed);

	if (word.id() >= target_.cols())
	  target_.conservativeResize(Eigen::NoChange, word.id() + 1);
	if (word.id() >= words_target_.size())
	  words_target_.resize(word.id() + 1, false);
	
	target_.col(word.id()).block(0, 0, dimension_embedding_, 1)
	  = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), dimension_embedding_, 1);

	words_target_[word.id()] = true;
      }
    }
  }

  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    karma::real_generator<double, real_policy> float10;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(dimension_embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(dimension_hidden_);
    rep["window"]    = utils::lexical_cast<std::string>(window_);
    rep["alpha"]     = utils::lexical_cast<std::string>(alpha_);
    rep["beta"]      = utils::lexical_cast<std::string>(beta_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), rep.path("vocab-source"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), rep.path("vocab-target"), target_, words_target_);
    
    // dump matrices...
    write(rep.path("Wt1.txt.gz"), rep.path("Wt1.bin"), Wt1_);
    write(rep.path("bt1.txt.gz"), rep.path("bt1.bin"), bt1_);

    write(rep.path("Wt2.txt.gz"), rep.path("Wt2.bin"), Wt2_);
    write(rep.path("bt2.txt.gz"), rep.path("bt2.bin"), bt2_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
  }

  void write_embedding(const path_type& path_text,
		       const path_type& path_binary,
		       const path_type& path_vocab,
		       const tensor_type& matrix,
		       const unique_set_type& words) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    karma::real_generator<double, real_policy> float10;
    
    const word_type::id_type rows = matrix.rows();
    const word_type::id_type cols = std::min(static_cast<size_type>(matrix.cols()), words.size());
    
    utils::compress_ostream os_txt(path_text, 1024 * 1024);
    utils::compress_ostream os_bin(path_binary, 1024 * 1024);
    std::ostream_iterator<char> iter(os_txt);
    
    vocab_type vocab;
    vocab.open(path_vocab, words.size());
    
    for (word_type::id_type id = 0; id != cols; ++ id)  
      if (words[id]) {
	const word_type word(id);
	
	karma::generate(iter, standard::string, word);
	
	for (difference_type j = 0; j != rows; ++ j)
	  karma::generate(iter, karma::lit(' ') << float10, matrix(j, id));
	
	karma::generate(iter, karma::lit('\n'));
	
	os_bin.write((char*) matrix.col(id).data(), sizeof(tensor_type::Scalar) * rows);
	
	vocab.insert(word);
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
  size_type window_;
  
  // hyperparameter
  double alpha_;
  double beta_;
  
  // Embedding
  tensor_type source_;
  tensor_type target_;

  unique_set_type words_source_;
  unique_set_type words_target_;
  
  // Wt1 and bt1 for encoding
  tensor_type Wt1_;
  tensor_type bt1_;
  
  // Wt2 and bt2 for reconstruction
  tensor_type Wt2_;
  tensor_type bt2_;  

  // Wc and bc for classification
  tensor_type Wc_;
  tensor_type bc_;

  size_type scale_;
};

struct Dictionary
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;
  
  typedef uint64_t count_type;
  
  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;
  
  struct Dict
  {
    typedef utils::compact_map<word_type, double,
			       utils::unassigned<word_type>, utils::unassigned<word_type>,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, double> > > logprob_set_type;
    typedef utils::compact_map<word_type, count_type,
			       utils::unassigned<word_type>, utils::unassigned<word_type>,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, count_type> > > count_set_type;
    
    typedef std::vector<word_type, std::allocator<word_type> >   word_set_type;
    typedef boost::random::discrete_distribution<>               distribution_type;
    
    Dict() {}

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
    
    void initialize()
    {
      typedef std::pair<word_type, double> word_prob_type;
      typedef std::vector<word_prob_type, std::allocator<word_prob_type> > word_prob_set_type;
      typedef std::vector<double, std::allocator<double> > prob_set_type;
      
      logprobs_.clear();
      words_.clear();
      
      word_prob_set_type word_probs(counts_.begin(), counts_.end());
      std::sort(word_probs.begin(), word_probs.end(), compare_pair<word_prob_type>());
      
      prob_set_type probs;
      words_.reserve(word_probs.size());
      probs.reserve(word_probs.size());

      word_prob_set_type::const_iterator witer_end = word_probs.end();
      for (word_prob_set_type::const_iterator witer = word_probs.begin(); witer != witer_end; ++ witer) {
	words_.push_back(witer->first);
	probs.push_back(witer->second);
      }
      
      const double norm = 1.0 / std::accumulate(probs.begin(), probs.end(), double(0));
      for (word_prob_set_type::const_iterator witer = word_probs.begin(); witer != witer_end; ++ witer)
	logprobs_[witer->first] = std::log(witer->second * norm);
      
      // initialize distribution
      distribution_ = distribution_type(probs.begin(), probs.end());
    }
    
    double logprob(const word_type& word) const
    {
      logprob_set_type::const_iterator liter = logprobs_.find(word);
      if (liter != logprobs_.end())
	return liter->second;
      else
	return - std::numeric_limits<double>::infinity();
    }
    
    template <typename Gen>
    word_type draw(Gen& gen) const
    {
      return words_[distribution_(gen)];
    }

    count_type& operator[](const word_type& word) { return counts_[word]; }
    
    count_set_type    counts_;
    logprob_set_type  logprobs_;
    word_set_type     words_;
    distribution_type distribution_;
  };

  typedef Dict dict_type;
  typedef utils::alloc_vector<dict_type, std::allocator<dict_type> > dict_set_type;
  
  Dictionary() {}

  dict_type& operator[](const word_type& word) { return dicts_[word.id()]; }

  void swap(Dictionary& x) { dicts_.swap(x.dicts_); }

  void initialize()
  {
    for (size_type i = 0; i != dicts_.size(); ++ i)
      if (dicts_.exists(i))
	dicts_[i].initialize();
  }
  
  void clear()
  {
    dicts_.clear();
  }
  
  double logprob(const word_type& source, const word_type& target) const
  {
    if (dicts_.exists(source.id()))
      return dicts_[source.id()].logprob(target);
    else
      return dicts_[vocab_type::UNK.id()].logprob(target);
  }
  
  template <typename Gen>
  word_type draw(const word_type& source, Gen& gen) const
  {
    if (dicts_.exists(source.id()))
      return dicts_[source.id()].draw(gen);
    else
      return dicts_[vocab_type::UNK.id()].draw(gen);
  }

  size_type size(const word_type& source) const
  {
    return (dicts_.exists(source.id()) ? dicts_[source.id()].words_.size() : size_type(0));
  }

  dict_set_type dicts_;
};

struct Lexicon
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model      model_type;
  typedef Gradient   gradient_type;
  typedef Dictionary dictionary_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type tensor_type;

  typedef Average loss_type;
  
  typedef cicada::Bitext bitext_type;
  typedef cicada::Vocab vocab_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  
  Lexicon(const dictionary_type& dict_source_target,
	  const dictionary_type& dict_target_source)
    : dict_source_target_(dict_source_target),
      dict_target_source_(dict_target_source) {}

  template <typename Function, typename Derivative, typename Gen>
  std::pair<loss_type, loss_type> operator()(const sentence_type& source,
					     const sentence_type& target,
					     const model_type& theta,
					     gradient_type& gradient,
					     Function   func,
					     Derivative deriv,
					     Gen& gen)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
    typedef gradient_type::embedding_type embedding_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type dimension_embedding = theta.dimension_embedding_;
    const size_type dimension_hidden    = theta.dimension_hidden_;
    const size_type window    = theta.window_;
    
    loss_type error;
    loss_type error_classification;

    const size_type input_size = dimension_embedding * 2 * (window * 2 + 1);
    
    tensor_type input(input_size, 1);
    tensor_type input_sampled(input_size, 1);
    tensor_type y_minus_c(input_size, 1);
    
    boost::random::uniform_int_distribution<> uniform_source(0, source_size - 1);
    boost::random::uniform_int_distribution<> uniform_target(0, target_size - 1);
    
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	++ gradient.count_;
	
	// forward...
	if (! src) {
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    input.block(dimension_embedding * i, 0, dimension_embedding, 1)
	      = theta.source_.col(vocab_type::EPSILON.id()) * theta.scale_;
	} else {
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& embedding_source = (src + shift <= 0
						 ? vocab_type::BOS
						 : (src + shift > source_size
						    ? vocab_type::EOS
						    : source[src + shift - 1]));
	    
	    input.block(dimension_embedding * i, 0, dimension_embedding, 1)
	      = theta.source_.col(embedding_source.id()) * theta.scale_;
	  }
	}
	
	if (! trg) {
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    input.block(dimension_embedding * i + offset, 0, dimension_embedding, 1)
	      = theta.target_.col(vocab_type::EPSILON.id()) * theta.scale_;
	} else {
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& embedding_target = (trg + shift <= 0
						 ? vocab_type::BOS
						 : (trg + shift > target_size
						    ? vocab_type::EOS
						    : target[trg + shift - 1]));
	    
	    input.block(dimension_embedding * i + offset, 0, dimension_embedding, 1)
	      = theta.target_.col(embedding_target.id()) * theta.scale_;
	  }
	}
	
	// compute sampled input...
	
	input_sampled = input;
	
	word_type source_sampled = vocab_type::EPSILON;
	word_type target_sampled = vocab_type::EPSILON;
	
	if (src) {
	  if (trg && dict_target_source_.size(target[trg - 1]) > 1) {
	    source_sampled = dict_target_source_.draw(target[trg - 1], gen);
	    while (source_sampled == source[src - 1])
	      source_sampled = dict_target_source_.draw(target[trg - 1], gen);
	  } else {
	    source_sampled = dict_target_source_.draw(target[uniform_target(gen)], gen);
	    while (source_sampled == source[src - 1])
	      source_sampled = dict_target_source_.draw(target[uniform_target(gen)], gen);
	  }
	  
	  input_sampled.block(dimension_embedding * window, 0, dimension_embedding, 1)
	    = theta.source_.col(source_sampled.id()) * theta.scale_;
	}
	
	if (trg) {
	  if (src && dict_source_target_.size(source[src - 1]) > 1) {
	    target_sampled = dict_source_target_.draw(source[src - 1], gen);
	    while (target_sampled == target[trg - 1])
	      target_sampled = dict_source_target_.draw(source[src - 1], gen);
	  } else {
	    target_sampled = dict_source_target_.draw(source[uniform_source(gen)], gen);
	    while (target_sampled == target[trg - 1])
	      target_sampled = dict_source_target_.draw(source[uniform_source(gen)], gen);
	  }
	  
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  input_sampled.block(offset + dimension_embedding * window, 0, dimension_embedding, 1)
	    = theta.target_.col(target_sampled.id()) * theta.scale_;
	}

#if 0
	std::cerr << "input:   rows: " << input.rows() << " cols: " << input.cols() << std::endl
		  << "sampled: rows: " << input_sampled.rows() << " cols: " << input_sampled.cols() << std::endl;
#endif
	
	const tensor_type p = (theta.Wt1_ * input + theta.bt1_).array().unaryExpr(func);
	const tensor_type p_norm = p.normalized();
	
	const tensor_type y = (theta.Wt2_ * p_norm + theta.bt2_).array().unaryExpr(func);
	
	for (size_type i = 0; i != 2 * (window * 2 + 1); ++ i)
	  y_minus_c.block(i * dimension_embedding, 0, dimension_embedding, 1)
	    = (y.block(i * dimension_embedding, 0, dimension_embedding, 1).normalized()
	       - input.block(i * dimension_embedding, 0, dimension_embedding, 1));
	
	const tensor_type p_sampled = (theta.Wt1_ * input_sampled + theta.bt1_).array().unaryExpr(func);
	const tensor_type p_sampled_norm = p_sampled.normalized();
	
	const double e = theta.alpha_ * 0.5 * y_minus_c.squaredNorm();

	//std::cerr << "error: " << e << std::endl;
	
	const tensor_type reconstruction       = y_minus_c.array() * theta.alpha_;
	const tensor_type delta_reconstruction = y.array().unaryExpr(deriv) * y_minus_c.array() * theta.alpha_;
	
	const double y_p = (theta.Wc_ * p_norm + theta.bc_)(0,0);
	const double y_m = (theta.Wc_ * p_sampled_norm + theta.bc_)(0,0);
	
	const double e_classification = std::max(1.0 - (y_p - y_m), 0.0) * theta.beta_;
	
	//std::cerr << "classification: " << e_classification << std::endl;
	
	const double delta_classification_p = - (e_classification > 0.0) * theta.beta_;
	const double delta_classification_m =   (e_classification > 0.0) * theta.beta_;

	// update error...
	error                += e;
	error_classification += e_classification;

	// backward...
	
	const tensor_type delta = (p.array().unaryExpr(deriv)
				   * (theta.Wt2_.transpose() * delta_reconstruction
				      + theta.Wc_.transpose() * delta_classification_p).array());
	
	const tensor_type delta_sampled = (p_sampled.array().unaryExpr(deriv)
					   * (theta.Wc_.transpose() * delta_classification_m).array());
	
#if 0
	std::cerr << "delta:   rows: " << delta.rows() << " cols: " << delta.cols() << std::endl
		  << "sampled: rows: " << delta_sampled.rows() << " cols: " << delta_sampled.cols() << std::endl;
#endif
	
	gradient.Wt1_ += delta * input.transpose();
	gradient.bt1_ += delta;
	
	gradient.Wt1_ += delta_sampled * input_sampled.transpose();
	gradient.bt1_ += delta_sampled;
	
	gradient.Wt2_ += delta_reconstruction * p_norm.transpose();
	gradient.bt2_ += delta_reconstruction;

	//std::cerr << "update classification" << std::endl;
	
	gradient.Wc_         += delta_classification_p * p_norm.transpose();
	gradient.bc_.array() += delta_classification_p;
	
	gradient.Wc_         += delta_classification_m * p_sampled_norm.transpose();
	gradient.bc_.array() += delta_classification_m;
	
	//std::cerr << "update embedding" << std::endl;
	
	// update embedding
	const tensor_type delta_embedding_p = theta.Wt1_.transpose() * delta - reconstruction;
	const tensor_type delta_embedding_m = theta.Wt1_.transpose() * delta_sampled;
	
	if (! src) {
	  tensor_type& dsource = gradient.source_[vocab_type::EPSILON];
	  
	  if (! dsource.cols() || ! dsource.rows())
	    dsource = tensor_type::Zero(dimension_embedding, 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    dsource += delta_embedding_p.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    dsource += delta_embedding_m.block(dimension_embedding * i, 0, dimension_embedding, 1);
	  }
	} else {
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& word = (src + shift <= 0
				     ? vocab_type::BOS
				     : (src + shift > source_size
					? vocab_type::EOS
					: source[src + shift - 1]));
	    
	    tensor_type& dsource = gradient.source_[word];
	    
	    if (! dsource.cols() || ! dsource.rows())
	      dsource = tensor_type::Zero(dimension_embedding, 1);
	    
	    dsource += delta_embedding_p.block(dimension_embedding * i, 0, dimension_embedding, 1);

	    if (shift != 0)
	      dsource += delta_embedding_m.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    else {
	      tensor_type& dsource = gradient.source_[source_sampled];
	      
	      if (! dsource.cols() || ! dsource.rows())
		dsource = tensor_type::Zero(dimension_embedding, 1);
	      
	      dsource += delta_embedding_m.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    }
	  }
	}

	if (! trg) {
	  tensor_type& dtarget = gradient.target_[vocab_type::EPSILON];
	  
	  if (! dtarget.cols() || ! dtarget.rows())
	    dtarget = tensor_type::Zero(dimension_embedding, 1);
	  
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    dtarget += delta_embedding_p.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    dtarget += delta_embedding_m.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	  }
	} else {
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& word = (trg + shift <= 0
				     ? vocab_type::BOS
				     : (trg + shift > target_size
					? vocab_type::EOS
					: target[trg + shift - 1]));
	    
	    tensor_type& dtarget = gradient.target_[word];
	    
	    if (! dtarget.cols() || ! dtarget.rows())
	      dtarget = tensor_type::Zero(dimension_embedding, 1);
	    
	    dtarget += delta_embedding_p.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);

	    if (shift != 0)
	      dtarget += delta_embedding_m.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    else {
	      tensor_type& dtarget = gradient.target_[target_sampled];
	      
	      if (! dtarget.cols() || ! dtarget.rows())
		dtarget = tensor_type::Zero(dimension_embedding, 1);
	      
	      dtarget += delta_embedding_m.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    }
	  }
	}
      }
    
    return std::make_pair(error, error_classification);
  }

  const dictionary_type& dict_source_target_;
  const dictionary_type& dict_target_source_;
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
	       const size_type& window,
	       const double& lambda,
	       const double& eta0)
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden),
      window_(window), lambda_(lambda), eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    source_ = tensor_type::Zero(dimension_embedding, vocabulary_size);
    target_ = tensor_type::Zero(dimension_embedding, vocabulary_size);
    
    // initialize...
    Wt1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1));
    bt1_ = tensor_type::Zero(dimension_hidden, 1);
    
    Wt2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden);
    bt2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);    
    
    Wc_ = tensor_type::Zero(1, dimension_hidden);
    bc_ = tensor_type::Zero(1, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    const double scale = 1.0 / gradient.count_;

    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first, theta.source_, const_cast<tensor_type&>(source_), siter->second, scale, lambda_ != 0.0);
    
    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first, theta.target_, const_cast<tensor_type&>(target_), titer->second, scale, lambda_ != 0.0);
    
    update(theta.Wt1_, const_cast<tensor_type&>(Wt1_), gradient.Wt1_, scale, lambda_ != 0.0);
    update(theta.bt1_, const_cast<tensor_type&>(bt1_), gradient.bt1_, scale, false);

    update(theta.Wt2_, const_cast<tensor_type&>(Wt2_), gradient.Wt2_, scale, lambda_ != 0.0);
    update(theta.bt2_, const_cast<tensor_type&>(bt2_), gradient.bt2_, scale, false);

    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, scale, false);
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
	      const bool regularize=true) const
  {
    if (regularize) {
      for (int row = 0; row != g.rows(); ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) += g(row, 0) * g(row, 0) * scale * scale;
	  
	  const double rate = eta0_ / std::sqrt(double(1.0) + G(row, word.id()));
	  const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	  
	  theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
    } else {
      G.col(word.id()).array() += g.array().square() * scale * scale;
      theta.col(word.id()).array() -= eta0_ * scale * g.array() * G.col(word.id()).array().unaryExpr(learning_rate());
    }
  }
  
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type window_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  // Wt1 and bt1 for encoding
  tensor_type Wt1_;
  tensor_type bt1_;
  
  // Wt2 and bt2 for reconstruction
  tensor_type Wt2_;
  tensor_type bt2_;

  // Wc and bc for classification
  tensor_type Wc_;
  tensor_type bc_;
};


struct LearnSGD
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef cicada::Symbol   word_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnSGD(const double& lambda,
	   const double& eta0)
    : lambda_(lambda), eta0_(eta0), epoch_(0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    //++ const_cast<size_type&>(epoch_);
    
#if 0
    // this regularization seems to be problematic...why?
    if (lambda_ != 0.0) {
      const double eta = eta0_ / (epoch_ + 1);
      
      theta.scale_ *= 1.0 - eta * lambda_;
    }
#endif

    const double scale = 1.0 / gradient.count_;
    
    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first, theta.source_, siter->second, scale, theta.scale_);
    
    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first, theta.target_, titer->second, scale, theta.scale_);
    
    update(theta.Wt1_, gradient.Wt1_, scale, lambda_ != 0.0);
    update(theta.bt1_, gradient.bt1_, scale, false);

    update(theta.Wt2_, gradient.Wt2_, scale, lambda_ != 0.0);
    update(theta.bt2_, gradient.bt2_, scale, false);

    update(theta.Wc_, gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, gradient.bc_, scale, false);
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
    
    theta.noalias() -= (eta * scale) * g;
  }

  template <typename Theta, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const double theta_scale) const
  {
    const double eta = eta0_ / (epoch_ + 1);
    
    theta.col(word.id()) -= (eta * scale / theta_scale) * g;
  }
  
  double lambda_;
  double eta0_;

  size_type epoch_;
};

typedef boost::filesystem::path path_type;

typedef cicada::Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;

typedef Model model_type;
typedef Dictionary dictionary_type;

path_type source_file;
path_type target_file;

path_type embedding_source_file;
path_type embedding_target_file;

path_type output_model_file;

double alpha = 0.99;
double beta = 0.01;
int dimension_embedding = 32;
int dimension_hidden = 128;
int window = 0;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 4;
double lambda = 0;
double eta0 = 0.1;
int cutoff = 3;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta);
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
    if (window < 0)
      throw std::runtime_error("window size should be >= 0");
    
    if (alpha < 0.0)
      throw std::runtime_error("alpha should be >= 0.0");
    if (beta < 0.0)
      throw std::runtime_error("beta should be >= 0.0");
    
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
    
    model_type theta(dimension_embedding, dimension_hidden, window, alpha, beta, sources, targets, generator);

    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());

      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta.read_embedding(embedding_source_file, embedding_target_file);
    }
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, window, lambda, eta0),
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

  typedef Lexicon lexicon_type;

  typedef lexicon_type::loss_type loss_type;

  typedef cicada::Bitext bitext_type;
  typedef cicada::Vocab  vocab_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  
  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_mapper_type;
  typedef utils::lockfree_list_queue<gradient_type*, std::allocator<gradient_type*> > queue_merger_type;
  typedef std::vector<queue_merger_type, std::allocator<queue_merger_type> > queue_merger_set_type;
  
  typedef std::deque<gradient_type, std::allocator<gradient_type> > gradient_set_type;

  TaskAccumulate(const Learner& learner,
		 const bitext_set_type& bitexts,
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 const model_type& theta,
		 const size_type batch_size,
		 queue_mapper_type& mapper,
		 queue_merger_set_type& mergers)
    : learner_(learner),
      bitexts_(bitexts),
      theta_(theta),
      mapper_(mapper),
      mergers_(mergers),
      lexicon_(dict_source_target, dict_target_source),
      batch_size_(batch_size)
  {
    generator_.seed(utils::random_seed());
  }
  
  struct sigmoid
  {
    double operator()(const double& x) const
    {
      return 1.0 / (std::exp(- x) + 1.0);
    }
  };

  struct dsigmoid
  {
    double operator()(const double& x) const
    {
      const double m = 1.0 / (std::exp(- x) + 1.0);
      
      return m * (1.0 - m);
    }
  };

  struct tanh
  {
    double operator()(const double& x) const
    {
      return std::tanh(x);
    }
  };
  
  struct dtanh
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return Tp(1) - x * x;
    }
  };
  
  struct htanh
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return std::min(std::max(x, Tp(- 1)), Tp(1));
    }
  };
  
  struct dhtanh
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return Tp(- 1) < x && x < Tp(1);
    }
  };

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

  void operator()()
  {
    clear();
    
    const size_type shard_size = mergers_.size();
    
    size_type      batch = 0;
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
	    if (gradients_[j].shared() == shard_size) {
	      grad = &gradients_[j];
	      break;
	    }
	  
	  if (! grad) {
	    gradients_.push_back(gradient_type(theta_.dimension_embedding_, theta_.dimension_hidden_, theta_.window_));
	    grad = &gradients_.back();
	  }
	  
	  grad->clear();
	  
	  const size_type first = batch * batch_size_;
	  const size_type last  = utils::bithack::min(first + batch_size_, bitexts_.size());
	  
	  for (size_type id = first; id != last; ++ id) {
	    const sentence_type& source = bitexts_[id].source_;
	    const sentence_type& target = bitexts_[id].target_;
	    
	    if (! source.empty() && ! target.empty()) {
	      std::pair<loss_type, loss_type> errors = lexicon_(source, target, theta_, *grad, htanh(), dhtanh(),
								generator_);
	      
	      error_          += errors.first;
	      classification_ += errors.second;
	    }
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
    error_ = loss_type();
    classification_ = loss_type();
  }

  Learner                learner_;
  const bitext_set_type& bitexts_;
  model_type             theta_;

  queue_mapper_type&     mapper_;
  queue_merger_set_type& mergers_;
  
  lexicon_type lexicon_;
  
  gradient_set_type gradients_;
  
  loss_type error_;
  loss_type classification_;
  
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
  
  typedef typename task_type::queue_mapper_type     queue_mapper_type;
  typedef typename task_type::queue_merger_set_type queue_merger_set_type;

  typedef typename task_type::loss_type loss_type;

  typedef std::vector<size_type, std::allocator<size_type> > batch_set_type;
  
  const size_type batches_size = (bitexts.size() + batch_size - 1) / batch_size;
  
  batch_set_type batches(batches_size);
  for (size_type batch = 0; batch != batches_size; ++ batch)
    batches[batch] = batch;
  
  queue_mapper_type     mapper(threads);
  queue_merger_set_type mergers(threads);
  
  task_set_type tasks(threads, task_type(learner,
					 bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta,
					 batch_size,
					 mapper,
					 mergers));
  
  // assign shard id
  for (size_type shard = 0; shard != tasks.size(); ++ shard)
    tasks[shard].shard_ = shard;
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;
    
    std::unique_ptr<boost::progress_display> progress(debug
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

    loss_type error;
    loss_type classification;
    for (size_type i = 0; i != tasks.size(); ++ i) {
      error          += tasks[i].error_;
      classification += tasks[i].classification_;
    }
    
    if (debug)
      std::cerr << "reconstruction error: " << static_cast<double>(error) << std::endl
		<< "classification error: " << static_cast<double>(classification) << std::endl;
    
    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle bitexts!
    boost::random_number_generator<boost::mt19937> gen(tasks.front().generator_);
    std::random_shuffle(batches.begin(), batches.end(), gen);
  }
  
  theta = tasks.front().theta_;
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
    
    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("alpha",     po::value<double>(&alpha)->default_value(alpha),      "parameter for reconstruction error")
    ("beta",      po::value<double>(&beta)->default_value(beta),        "parameter for classificaiton error")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("window",              po::value<int>(&window)->default_value(window),                           "context window size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
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
