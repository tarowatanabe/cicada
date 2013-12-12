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

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <set>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <boost/algorithm/string/trim.hpp>

#include <boost/functional/hash/hash.hpp>

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
#include "utils/compress_stream.hpp"
#include "utils/vector2.hpp"
#include "utils/sampler.hpp"
#include "utils/resource.hpp"
#include "utils/simple_vector.hpp"

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
  
  Gradient() : embedding_(0), hidden_(0), window_(0), count_(0), shared_(0) {}
  Gradient(const size_type& embedding,
	   const size_type& hidden,
	   const size_type& window) 
    : embedding_(embedding),
      hidden_(hidden),
      window_(window),
      count_(0),
      shared_(0)
  { initialize(embedding, hidden, window); }
  
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

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());
    
    Wt_ -= x.Wt_;
    bt_ -= x.bt_;

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

    if (! Wt_.rows())
      Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    if (! bt_.rows())
      bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    Wt_ += x.Wt_;
    bt_ += x.bt_;

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
    
    Wt_.setZero();
    bt_.setZero();

    Wc_.setZero();
    bc_.setZero();
    
    count_ = 0;
    shared_ = 0;
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
      embedding = tensor_type::Zero(embedding_, 1);
    
    return embedding;
  }

  
  void initialize(const size_type embedding, const size_type hidden, const size_type window)
  {
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    
    embedding_ = embedding;
    hidden_    = hidden;
    window_    = window;
    
    clear();
    
    const size_type state_size = embedding * (window * 2 + 1) * 2;
    
    // initialize...
    Wt_ = tensor_type::Zero(hidden_, state_size);
    bt_ = tensor_type::Zero(hidden_, 1);
    
    Wc_ = tensor_type::Zero(1, hidden_);
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
  size_type embedding_;
  size_type hidden_;
  size_type window_;
  
  // embedding
  embedding_type source_;
  embedding_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;

  tensor_type Wc_;
  tensor_type bc_;
  
  size_type count_;
  size_type shared_;
};

struct Embedding
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef cicada::Symbol    word_type;
  
  typedef Gradient gradient_type;

  Embedding(const size_type& embedding) : embedding_(embedding) { initialize(embedding); }

  void initialize(const size_type& embedding)
  {
    embedding_ = embedding;
    
    const size_type vocabulary_size = word_type::allocated();
    
    source_ = tensor_type::Zero(embedding_, vocabulary_size);
    target_ = tensor_type::Zero(embedding_, vocabulary_size);
  }

  template <typename Model>
  void assign(const gradient_type& x, const Model& theta)
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    gradient_embedding_type::const_iterator siter_end = x.source_.end();
    for (gradient_embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter)
      source_.col(siter->first.id()) = theta.source_.col(siter->first.id()) * theta.scale_;
    
    gradient_embedding_type::const_iterator titer_end = x.target_.end();
    for (gradient_embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer)
      target_.col(titer->first.id()) = theta.target_.col(titer->first.id()) * theta.scale_;
  }
  
  size_type embedding_;
  
  tensor_type source_;
  tensor_type target_;
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

  typedef std::vector<bool, std::allocator<bool> > word_unique_type;
  
  Model() : embedding_(0), hidden_(0), window_(0), scale_(1) {}  
  template <typename Words, typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	const size_type& window,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      window_(window),
      scale_(1)
  { initialize(embedding, hidden, window, words_source, words_target, gen); }
  
  void clear()
  {
    // embedding
    source_.setZero();
    target_.setZero();
    
    Wt_.setZero();
    bt_.setZero();

    Wc_.setZero();
    bc_.setZero();
    
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
  void initialize(const size_type embedding,
		  const size_type hidden,
		  const size_type window,
		  Words& words_source,
		  Words& words_target,
		  Gen& gen)
  {
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    
    embedding_ = embedding;
    hidden_    = hidden;
    window_    = window;
    
    clear();
    
    const size_type vocabulary_size = word_type::allocated();
    const size_type state_size = embedding_ * (window * 2 + 1) * 2;
    
    const double range_e = std::sqrt(6.0 / (embedding_ + 1));
    const double range_t = std::sqrt(6.0 / (hidden_ + state_size));
    const double range_c = std::sqrt(6.0 / (1 + hidden_));

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
 
    source_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    
    Wt_ = tensor_type::Zero(hidden_, state_size).array().unaryExpr(randomize<Gen>(gen, range_t));
    bt_ = tensor_type::Zero(hidden_, 1);

    Wc_ = tensor_type::Zero(1, hidden_).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Zero(1, 1);
    
    scale_ = 1.0;
  }
  
  void finalize()
  {
    if (scale_ == 1.0) return;
    
    source_ *= scale_;
    target_.block(0, 0, embedding_, target_.cols()) *= scale_;

    scale_ = 1.0;    
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
	
	if (boost::fusion::get<1>(parsed).size() != embedding_)
	  throw std::runtime_error("invalid embedding size");

	const word_type word = boost::fusion::get<0>(parsed);
	
	if (word.id() < source_.cols())
	  source_.col(word.id()).block(0, 0, embedding_, 1)
	    = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), embedding_, 1);
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
	
	if (boost::fusion::get<1>(parsed).size() != embedding_)
	  throw std::runtime_error("invalid embedding size");
	
	const word_type word = boost::fusion::get<0>(parsed);

	if (word.id() < target_.cols())
	  target_.col(word.id()).block(0, 0, embedding_, 1)
	    = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), embedding_, 1);
      }
    }
  }
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);
    rep["window"]    = utils::lexical_cast<std::string>(window_);
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), target_, words_target_);
    
    write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
    
    // vocabulary...
    vocab_type vocab;

    const word_type::id_type vocabulary_size = utils::bithack::max(words_source_.size(), words_target_.size());
    
    vocab.open(rep.path("vocab"), vocabulary_size >> 1);
    
    for (word_type::id_type id = 0; id != vocabulary_size; ++ id)
      if ((id < words_source_.size() && words_source_[id]) || (id < words_target_.size() && words_target_[id])) {
	const word_type word(id);
	
	vocab.insert(word);
      }
    
    vocab.close();
  }
  
  void write_embedding(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix, const word_unique_type& words) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    karma::real_generator<double, real_policy> float10;
    
    const word_type::id_type rows = matrix.rows();
    const word_type::id_type cols = utils::bithack::min(static_cast<size_type>(matrix.cols()), words.size());
    
    utils::compress_ostream os_txt(path_text, 1024 * 1024);
    utils::compress_ostream os_bin(path_binary, 1024 * 1024);
    std::ostream_iterator<char> iter(os_txt);
    
    for (word_type::id_type id = 0; id != cols; ++ id)  
      if (words[id]) {
	const word_type word(id);
	
	karma::generate(iter, standard::string, word);
	
	for (difference_type j = 0; j != rows; ++ j)
	  karma::generate(iter, karma::lit(' ') << float10, matrix(j, id));
	
	karma::generate(iter, karma::lit('\n'));
	
	os_bin.write((char*) matrix.col(id).data(), sizeof(tensor_type::Scalar) * rows);
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
  size_type embedding_;
  size_type hidden_;
  size_type window_;

  word_unique_type words_source_;
  word_unique_type words_target_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;

  tensor_type Wc_;
  tensor_type bc_;
  
  double scale_;
};


struct Dictionary
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  
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
	logprobs_[witer->first] = std::log(witer->second);
      }
      
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
  typedef model_type::tensor_type    tensor_type;

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;

  typedef std::vector<parameter_type, std::allocator<parameter_type> > layer_type;
  typedef std::vector<double, std::allocator<double> > score_set_type;

  typedef Average loss_type;
  
  Lexicon(const dictionary_type& dict,
	  const size_type samples)
    : dict_(dict), samples_(samples), log_samples_(std::log(double(samples))) {}
  
  const dictionary_type& dict_;

  size_type samples_;
  double    log_samples_;

  layer_type layer_input_;
  layer_type layer_hidden_;
  layer_type delta_hidden_;
  
  layer_type layer_input_sampled_;
  layer_type layer_hidden_sampled_;
  layer_type delta_hidden_sampled_;
  
  score_set_type scores_;
  score_set_type errors_;
  
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
  
  
  template <typename Embedding>
  void copy_embedding(const sentence_type& source,
		      const sentence_type& target,
		      const difference_type src,
		      const difference_type trg,
		      const model_type& theta,
		      Eigen::MatrixBase<Embedding>& embedding)
  {
    const difference_type source_size    = source.size();
    const difference_type target_size    = target.size();
    const difference_type window_size    = theta.window_;
    const difference_type embedding_size = theta.embedding_;

    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_size * (window_size * 2 + 1);

    if (src == 0) {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	embedding.block(offset_source + embedding_size * i, 0, embedding_size, 1)
	  = theta.source_.col(vocab_type::EPSILON.id()) * theta.scale_;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (src + shift <= 0
				 ? vocab_type::BOS
				 : (src + shift > source_size
				    ? vocab_type::EOS
				    : source[src + shift - 1]));
	
	embedding.block(offset_source + embedding_size * i, 0, embedding_size, 1)
	  = theta.source_.col(word.id()) * theta.scale_;
      }
    }
    
    if (trg == 0) {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	embedding.block(offset_target + embedding_size * i, 0, embedding_size, 1)
	  = theta.target_.col(vocab_type::EPSILON.id()) * theta.scale_;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (trg + shift <= 0
				 ? vocab_type::BOS
				 : (trg + shift > target_size
				    ? vocab_type::EOS
				    : target[trg + shift - 1]));
	
	embedding.block(offset_target + embedding_size * i, 0, embedding_size, 1)
	  = theta.target_.col(word.id()) * theta.scale_;
      }
    }
  }
  
  template <typename Embedding>
  void propagate_embedding(const sentence_type& source,
			   const sentence_type& target,
			   const difference_type src,
			   const difference_type trg,
			   const model_type& theta,
			   gradient_type& gradient,
			   const Eigen::MatrixBase<Embedding>& delta)
  {
    const difference_type source_size    = source.size();
    const difference_type target_size    = target.size();
    
    const difference_type embedding_size = theta.embedding_;
    const difference_type hidden_size    = theta.hidden_;
    const difference_type window_size    = theta.window_;

    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_size * (window_size * 2 + 1);
    
    if (src == 0) {
      tensor_type& dembedding = gradient.source(vocab_type::EPSILON.id());
      
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	dembedding
	  += theta.Wt_.block(0, offset_source + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (src + shift <= 0
				 ? vocab_type::BOS
				 : (src + shift > source_size
				    ? vocab_type::EOS
				    : source[src + shift - 1]));
	gradient.source(word.id())
	  += theta.Wt_.block(0, offset_source + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
      }
    }
    
    if (trg == 0) {
      tensor_type& dembedding = gradient.target(vocab_type::EPSILON.id());
      
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	dembedding
	  += theta.Wt_.block(0, offset_target + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (trg + shift <= 0
				 ? vocab_type::BOS
				 : (trg + shift > target_size
				    ? vocab_type::EOS
				    : target[trg + shift - 1]));
	gradient.target(word.id())
	  += theta.Wt_.block(0, offset_target + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
      }
    }
  }

  template <typename Embedding>
  void propagate_embedding(const sentence_type& source,
			   const sentence_type& target,
			   const difference_type src,
			   const difference_type trg,
			   const word_type& sampled,
			   const model_type& theta,
			   gradient_type& gradient,
			   const Eigen::MatrixBase<Embedding>& delta)
  {
    const difference_type source_size    = source.size();
    const difference_type target_size    = target.size();
    
    const difference_type embedding_size = theta.embedding_;
    const difference_type hidden_size    = theta.hidden_;
    const difference_type window_size    = theta.window_;

    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_size * (window_size * 2 + 1);
    
    if (src == 0) {
      tensor_type& dembedding = gradient.source(vocab_type::EPSILON.id());
      
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	dembedding
	  += theta.Wt_.block(0, offset_source + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (src + shift <= 0
				 ? vocab_type::BOS
				 : (src + shift > source_size
				    ? vocab_type::EOS
				    : source[src + shift - 1]));
	gradient.source(word.id())
	  += theta.Wt_.block(0, offset_source + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
      }
    }
    
    if (trg == 0) {
      tensor_type& dembedding = gradient.target(vocab_type::EPSILON.id());
      
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	dembedding
	  += theta.Wt_.block(0, offset_target + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (i == window_size
				 ? sampled
				 : (trg + shift <= 0
				    ? vocab_type::BOS
				    : (trg + shift > target_size
				       ? vocab_type::EOS
				       : target[trg + shift - 1])));
	gradient.target(word.id())
	  += theta.Wt_.block(0, offset_target + i * embedding_size, hidden_size, embedding_size).transpose() * delta;
      }
    }
  }
  
  template <typename Generator>
  loss_type learn(const sentence_type& source,
		  const sentence_type& target,
		  const model_type& theta,
		  gradient_type& gradient,
		  alignment_type& alignment,
		  Generator& gen)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type embedding_size = theta.embedding_;
    const size_type hidden_size    = theta.hidden_;
    const size_type window_size    = theta.window_;
    
    const size_type state_size = embedding_size * (window_size * 2 + 1) * 2;
    
    const size_type offset_source = 0;
    const size_type offset_target = embedding_size * (window_size * 2 + 1);

    loss_type loss;
    
    alignment.clear();
    
    // construct source-side input and hidden layers
    layer_input_.resize((source_size + 1) * state_size);
    layer_hidden_.resize((source_size + 1) * hidden_size);
    delta_hidden_.resize((source_size + 1) * hidden_size);

    layer_input_sampled_.resize((source_size + 1) * state_size);
    layer_hidden_sampled_.resize((source_size + 1) * hidden_size);
    delta_hidden_sampled_.resize((source_size + 1) * hidden_size);
    
    boost::random::uniform_int_distribution<> uniform_source(0, source_size - 1);
    
    // enumerate target side
    for (size_type trg = 1; trg <= target_size; ++ trg) {
      const word_type word_target(target[trg - 1]);
      
      size_type align_best = 0;
      double    score_best = - std::numeric_limits<double>::infinity();
      
      scores_.clear();

      errors_.clear();
      errors_.resize(source_size + 1, 0.0);
      
      for (size_type src = 0; src <= source_size; ++ src) {
	const word_type word_source(src == 0 ? vocab_type::EPSILON : source[src - 1]);
	
	matrix_type embedding(&(*layer_input_.begin()) + state_size * src, state_size, 1);
	matrix_type hidden(&(*layer_hidden_.begin()) + hidden_size * src, hidden_size, 1);
	
	copy_embedding(source, target, src, trg, theta, embedding);
	
	hidden = (theta.Wt_ * embedding + theta.bt_).array().unaryExpr(hinge());
	
	const double score = (theta.Wc_ * hidden + theta.bc_)(0, 0);
	
	scores_.push_back(score);
	
	if (score > score_best) {
	  align_best = src;
	  score_best = score;
	}
      }
      
      if (align_best)
	alignment.push_back(std::make_pair(align_best - 1, trg - 1));
      
      const double loss_factor = 1.0 / ((source_size + 1) * samples_);
      
      double loss_target = 0.0;
      
      // copy input layer
      layer_input_sampled_ = layer_input_;
      
      // perform sampling
      for (size_type k = 0; k != samples_; ++ k) {
	// first, choose one source-word-position
	// then, sample a target word
	
	word_type sampled_target = dict_.draw(source[uniform_source(gen)], gen);
	while (sampled_target == word_target)
	  sampled_target = dict_.draw(source[uniform_source(gen)], gen);
	
	for (size_type src_sampled = 0; src_sampled <= source_size; ++ src_sampled) {
	  const word_type sampled_source(src_sampled == 0 ? vocab_type::EPSILON : source[src_sampled - 1]);
	  
	  matrix_type embedding_sampled(&(*layer_input_sampled_.begin()) + state_size * src_sampled, state_size, 1);
	  matrix_type hidden_sampled(&(*layer_hidden_sampled_.begin()) + hidden_size * src_sampled, hidden_size, 1);
	  matrix_type delta_sampled(&(*delta_hidden_sampled_.begin()) + hidden_size * src_sampled, hidden_size, 1);
	  
	  embedding_sampled.block(offset_target + embedding_size * window_size, 0, embedding_size, 1)
	    = theta.target_.col(sampled_target.id()) * theta.scale_;
	  
	  hidden_sampled = (theta.Wt_ * embedding_sampled + theta.bt_).array().unaryExpr(hinge());
	  
	  const double score = (theta.Wc_ * hidden_sampled + theta.bc_)(0, 0);
	  
	  double errors = 0.0;
	  
	  for (size_type src = 0; src <= source_size; ++ src) {
	    const word_type word_source(src == 0 ? vocab_type::EPSILON : source[src - 1]);
	    
	    const double error = std::max(1.0 - (scores_[src] - score), 0.0) * loss_factor;
	    
	    if (error == 0.0) continue;
	    
	    errors       += error;
	    errors_[src] -= error;
	  }
	  
	  loss_target += errors;
	    
	  // we will propagate - error for correct instance, and error for sampled instance
	  
	  gradient.Wc_.array() += errors * hidden_sampled.transpose().array();
	  gradient.bc_.array() += errors;
	  
	  delta_sampled = (hidden_sampled.array().unaryExpr(dhinge()) * (theta.Wc_.transpose() * errors).array());
	  
	  gradient.Wt_.noalias() += delta_sampled * embedding_sampled.transpose();
	  gradient.bt_.noalias() += delta_sampled;
	  
	  propagate_embedding(source, target, src_sampled, trg, sampled_target, theta, gradient, delta_sampled);
	}
      }

      for (size_type src = 0; src <= source_size; ++ src) {
	const matrix_type embedding(&(*layer_input_.begin()) + state_size * src, state_size, 1);
	const matrix_type hidden(&(*layer_hidden_.begin()) + hidden_size * src, hidden_size, 1);
	/**/  matrix_type delta(&(*delta_hidden_.begin()) + hidden_size * src, hidden_size, 1);
	
	gradient.Wc_.array() += errors_[src] * hidden.transpose().array();
	gradient.bc_.array() += errors_[src];
	
	delta = (hidden.array().unaryExpr(dhinge()) * (theta.Wc_.transpose() * errors_[src]).array());
	
	gradient.Wt_.noalias() += delta * embedding.transpose();
	gradient.bt_.noalias() += delta;
	
	propagate_embedding(source, target, src, trg, theta, gradient, delta);
      }
      
      ++ gradient.count_;
      loss += loss_target;
    }
    
    return loss;
  }
  
    
  double viterbi(const sentence_type& source,
		 const sentence_type& target,
		 const model_type& theta,
		 alignment_type& alignment)
  {
    typedef Eigen::Map<tensor_type> matrix_type;

    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type embedding_size = theta.embedding_;
    const size_type hidden_size    = theta.hidden_;
    const size_type window_size    = theta.window_;
    
    const size_type state_size = embedding_size * (window_size * 2 + 1) * 2;
    
    double total = 0.0;
    
    alignment.clear();
    
    // construct source-side input and hidden layers
    layer_input_.resize((source_size + 1) * state_size);
    layer_hidden_.resize((source_size + 1) * hidden_size);
    
    // enumerate target side
    for (size_type trg = 1; trg <= target_size; ++ trg) {
      const word_type word(target[trg - 1]);
      
      size_type align_best = 0;
      double    score_best = - std::numeric_limits<double>::infinity();
      
      for (size_type src = 0; src <= source_size; ++ src) {
	matrix_type embedding(&(*layer_input_.begin()) + state_size * src, state_size, 1);
	
	copy_embedding(source, target, src, trg, theta, embedding);
	
	const double score = (theta.Wc_ * (theta.Wt_ * embedding + theta.bt_).array().unaryExpr(hinge()).matrix()
			      + theta.bc_)(0, 0);
	
	if (score > score_best) {
	  align_best = src;
	  score_best = score;
	}
      }
      
      if (align_best)
	alignment.push_back(std::make_pair(align_best - 1, trg - 1));
      
      total += score_best;
    }

    return total;
  }
};


struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model     model_type;
  typedef Gradient  gradient_type;
  typedef Embedding embedding_type;

  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& embedding,
	       const size_type& hidden,
	       const size_type& window,
	       const double& lambda,
	       const double& lambda2,
	       const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      window_(window),
      lambda_(lambda),
      lambda2_(lambda2),
      eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");

    if (lambda2_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    
    source_ = tensor_type::Zero(embedding_, vocabulary_size);
    target_ = tensor_type::Zero(embedding_, vocabulary_size);

    const size_type state_size = embedding * (window * 2 + 1) * 2;
    
    Wt_ = tensor_type::Zero(hidden_, state_size);
    bt_ = tensor_type::Zero(hidden_, 1);

    Wc_ = tensor_type::Zero(1, hidden_);
    bc_ = tensor_type::Zero(1, 1);
  }

  
  void operator()(model_type& theta,
		  const gradient_type& gradient,
		  const embedding_type& embedding) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;
    
    const double scale = 1.0 / gradient.count_;

    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     const_cast<tensor_type&>(source_),
	     siter->second,
	     embedding.target_.col(siter->first.id()),
	     scale);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     const_cast<tensor_type&>(target_),
	     titer->second,
	     embedding.source_.col(titer->first.id()),
	     scale);
    
    update(theta.Wt_, const_cast<tensor_type&>(Wt_), gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, const_cast<tensor_type&>(bt_), gradient.bt_, scale, false);

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
    learning_rate(const double& eta0) : eta0_(eta0) {}

    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return (x == 0.0 ? 0.0 : eta0_ / std::sqrt(double(1.0) + x));
    }

    const double& eta0_;
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
      theta.array() -= scale * g.array() * G.array().unaryExpr(learning_rate(eta0_));
    }
  }

  template <typename Theta, typename GradVar, typename Grad, typename GradCross>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const Eigen::MatrixBase<GradCross>& c,
	      const double scale) const
  {
    if (lambda2_ > 0.0) {
      for (int row = 0; row != g.rows(); ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0) * scale * scale;
	  
	  const double rate = eta0_ / std::sqrt(double(1.0) + G(row, word.id()));
	  const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	  const double value = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	  
	  const double shared = c(row, 0);
	  const double diff = value - shared;
	  
	  theta(row, word.id()) = utils::mathop::sgn(diff) * std::max(0.0, std::fabs(diff) - rate * lambda2_) + shared;
	}
    } else {
      for (int row = 0; row != g.rows(); ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0) * scale * scale;
	  
	  const double rate = eta0_ / std::sqrt(double(1.0) + G(row, word.id()));
	  const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	  
	  theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
    }    
  }
  
  size_type embedding_;
  size_type hidden_;
  size_type window_;
  
  double lambda_;
  double lambda2_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;

  tensor_type Wc_;
  tensor_type bc_;
};

struct LearnSGD
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model     model_type;
  typedef Gradient  gradient_type;
  typedef Embedding embedding_type;

  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnSGD(const double& lambda,
	   const double& lambda2,
	   const double& eta0)
    : lambda_(lambda),
      lambda2_(lambda2),
      eta0_(eta0),
      epoch_(0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");

    if (lambda2_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
  }

  void operator()(model_type& theta,
		  const gradient_type& gradient,
		  const embedding_type& embedding) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    //++ const_cast<size_type&>(epoch_);
    
    const double eta = eta0_ / (epoch_ + 1);
    const double scale = 1.0 / gradient.count_;
    
    if (lambda_ != 0.0)
      theta.scale_ *= 1.0 - eta * lambda_;
    
    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     siter->second,
	     embedding.target_.col(siter->first.id()),
	     scale,
	     theta.scale_);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     titer->second,
	     embedding.source_.col(titer->first.id()),
	     scale,
	     theta.scale_);
    
    update(theta.Wt_, gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, gradient.bt_, scale, false);

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
  
  template <typename Theta, typename Grad, typename GradCross>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      const Eigen::MatrixBase<Grad>& g,
	      const Eigen::MatrixBase<GradCross>& c,
	      const double scale,
	      const double theta_scale) const
  {
    const double eta = eta0_ / (epoch_ + 1);
    
    if (lambda2_ != 0.0)
      theta.col(word.id()) -= eta * lambda2_ * (theta.col(word.id()) - c / theta_scale);

    theta.col(word.id()) -= (eta * scale / theta_scale) * g;
  }

  double lambda_;
  double lambda2_;
  double eta0_;
  
  size_type epoch_;
};


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

path_type model_source_target_file;
path_type model_target_source_file;

path_type output_source_target_file;
path_type output_target_source_file;
path_type alignment_source_target_file;
path_type alignment_target_source_file;

int dimension_embedding = 32;
int dimension_hidden = 128;
int window = 0;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 4;
int sample_size = 10;
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
    if (window < 0)
      throw std::runtime_error("window size should be positive");

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

    if (debug)
      std::cerr << "# of sentences: " << bitexts.size() << std::endl;
    
    const dictionary_type::dict_type::word_set_type& sources = dict_target_source[cicada::Vocab::EPSILON].words_;
    const dictionary_type::dict_type::word_set_type& targets = dict_source_target[cicada::Vocab::EPSILON].words_;

    model_type theta_source_target(dimension_embedding, dimension_hidden, window, sources, targets, generator);
    model_type theta_target_source(dimension_embedding, dimension_hidden, window, targets, sources, generator);

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
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, window, lambda, lambda2, eta0),
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

  typedef Lexicon lexicon_type;
  
  typedef lexicon_type::loss_type loss_type;
  
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
      lexicon_source_target_(dict_source_target, sample_size),
      lexicon_target_source_(dict_target_source, sample_size),
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
    const size_type window_size    = theta_source_target_.window_;
    
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
	    gradients_.push_back(gradient_type(embedding_size, hidden_size, window_size));
	    grad_source_target = &gradients_.back();
	  }
	  
	  if (! grad_target_source) {
	    gradients_.push_back(gradient_type(embedding_size, hidden_size, window_size));
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
	      loss_source_target_
		+= lexicon_source_target_.learn(bitext.source_,
						bitext.target_,
						theta_source_target_,
						*grad_source_target,
						bitext_source_target.alignment_,
						generator_);
	      
	      loss_target_source_
		+= lexicon_target_source_.learn(bitext.target_,
						bitext.source_,
						theta_target_source_,
						*grad_target_source,
						bitext_target_source.alignment_,
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
    loss_source_target_ = loss_type();
    loss_target_source_ = loss_type();
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
  
  lexicon_type lexicon_source_target_;
  lexicon_type lexicon_target_source_;
  
  gradient_set_type   gradients_;
  loss_type loss_source_target_;
  loss_type loss_target_source_;

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
  
  typedef typename task_type::loss_type loss_type;

  typedef std::vector<size_type, std::allocator<size_type> > batch_set_type;
  
  const size_type batches_size = (bitexts.size() + batch_size - 1) / batch_size;
  
  batch_set_type batches(batches_size);
  for (size_type batch = 0; batch != batches_size; ++ batch)
    batches[batch] = batch;
  
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
    
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;

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
    
    loss_type loss_source_target;
    loss_type loss_target_source;
    
    for (size_type i = 0; i != tasks.size(); ++ i) {
      loss_source_target += tasks[i].loss_source_target_;
      loss_target_source += tasks[i].loss_target_source_;
    }
    
    if (debug)
      std::cerr << "loss P(target | source): " << static_cast<double>(loss_source_target) << std::endl
		<< "loss P(source | target): " << static_cast<double>(loss_target_source) << std::endl;

    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;
    
    // shuffle bitexts!
    {
      typename batch_set_type::iterator biter     = batches.begin();
      typename batch_set_type::iterator biter_end = batches.end();
      
      while (biter < biter_end) {
	typename batch_set_type::iterator iter_end = std::min(biter + (batch_size << 5), biter_end);
	
	std::random_shuffle(biter, iter_end);
	biter = iter_end;
      }
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

  typedef Lexicon lexicon_type;
  
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
      lexicon_source_target_(dict_source_target, sample_size),
      lexicon_target_source_(dict_target_source, sample_size) {}

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
	lexicon_source_target_.viterbi(bitext.source_, bitext.target_, theta_source_target_, bitext_source_target.alignment_);
	
	lexicon_target_source_.viterbi(bitext.target_, bitext.source_, theta_target_source_, bitext_target_source.alignment_);
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
  
  lexicon_type lexicon_source_target_;
  lexicon_type lexicon_target_source_;
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
  size_type num_text = 0;
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

    ("model-source-target", po::value<path_type>(&model_source_target_file), "model parameter for P(target | source)")
    ("model-target-source", po::value<path_type>(&model_target_source_file), "model parameter for P(source | target)")
    
    ("output-source-target", po::value<path_type>(&output_source_target_file), "output model parameter for P(target | source)")
    ("output-target-source", po::value<path_type>(&output_target_source_file), "output model parameter for P(source | target)")

    ("alignment-source-target", po::value<path_type>(&alignment_source_target_file), "output alignment for P(target | source)")
    ("alignment-target-source", po::value<path_type>(&alignment_target_source_file), "output alignment for P(source | target)")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("window",              po::value<int>(&window)->default_value(window),                           "context window size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),     "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size),   "mini-batch size")
    ("sample",            po::value<int>(&sample_size)->default_value(sample_size), "sample size")
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
