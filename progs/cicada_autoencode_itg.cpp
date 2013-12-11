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

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <set>
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
#include "utils/vector2.hpp"
#include "utils/sampler.hpp"
#include "utils/resource.hpp"

#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>
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
  
  Gradient() : embedding_(0), hidden_(0), count_(0) {}
  Gradient(const size_type& embedding,
	   const size_type& hidden) 
    : embedding_(embedding),
      hidden_(hidden),
      count_(0),
      shared_(0)
  { initialize(embedding, hidden); }

  Gradient& operator-=(const Gradient& x)
  {
    decrement(source_, x.source_);
    decrement(target_, x.target_);

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    if (! Ws1_.rows())
      Ws1_ = tensor_type::Zero(x.Ws1_.rows(), x.Ws1_.cols());
    if (! bs1_.rows())
      bs1_ = tensor_type::Zero(x.bs1_.rows(), x.bs1_.cols());
    if (! Ws2_.rows())
      Ws2_ = tensor_type::Zero(x.Ws2_.rows(), x.Ws2_.cols());
    if (! bs2_.rows())
      bs2_ = tensor_type::Zero(x.bs2_.rows(), x.bs2_.cols());

    if (! Wi1_.rows())
      Wi1_ = tensor_type::Zero(x.Wi1_.rows(), x.Wi1_.cols());
    if (! bi1_.rows())
      bi1_ = tensor_type::Zero(x.bi1_.rows(), x.bi1_.cols());
    if (! Wi2_.rows())
      Wi2_ = tensor_type::Zero(x.Wi2_.rows(), x.Wi2_.cols());
    if (! bi2_.rows())
      bi2_ = tensor_type::Zero(x.bi2_.rows(), x.bi2_.cols());

    if (! Wt1_.rows())
      Wt1_ = tensor_type::Zero(x.Wt1_.rows(), x.Wt1_.cols());
    if (! bt1_.rows())
      bt1_ = tensor_type::Zero(x.bt1_.rows(), x.bt1_.cols());
    if (! Wt2_.rows())
      Wt2_ = tensor_type::Zero(x.Wt2_.rows(), x.Wt2_.cols());
    if (! bt2_.rows())
      bt2_ = tensor_type::Zero(x.bt2_.rows(), x.bt2_.cols());
    
    Wc_ -= x.Wc_;
    bc_ -= x.bc_;

    Ws1_ -= x.Ws1_;
    bs1_ -= x.bs1_;
    Ws2_ -= x.Ws2_;
    bs2_ -= x.bs2_;
    
    Wi1_ -= x.Wi1_;
    bi1_ -= x.bi1_;
    Wi2_ -= x.Wi2_;
    bi2_ -= x.bi2_;
    
    Wt1_ -= x.Wt1_;
    bt1_ -= x.bt1_;    
    Wt2_ -= x.Wt2_;
    bt2_ -= x.bt2_;    
    
    count_ -= x.count_;

    return *this;
  }
  
  Gradient& operator+=(const Gradient& x)
  {
    increment(source_, x.source_);
    increment(target_, x.target_);

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    if (! Ws1_.rows())
      Ws1_ = tensor_type::Zero(x.Ws1_.rows(), x.Ws1_.cols());
    if (! bs1_.rows())
      bs1_ = tensor_type::Zero(x.bs1_.rows(), x.bs1_.cols());
    if (! Ws2_.rows())
      Ws2_ = tensor_type::Zero(x.Ws2_.rows(), x.Ws2_.cols());
    if (! bs2_.rows())
      bs2_ = tensor_type::Zero(x.bs2_.rows(), x.bs2_.cols());

    if (! Wi1_.rows())
      Wi1_ = tensor_type::Zero(x.Wi1_.rows(), x.Wi1_.cols());
    if (! bi1_.rows())
      bi1_ = tensor_type::Zero(x.bi1_.rows(), x.bi1_.cols());
    if (! Wi2_.rows())
      Wi2_ = tensor_type::Zero(x.Wi2_.rows(), x.Wi2_.cols());
    if (! bi2_.rows())
      bi2_ = tensor_type::Zero(x.bi2_.rows(), x.bi2_.cols());

    if (! Wt1_.rows())
      Wt1_ = tensor_type::Zero(x.Wt1_.rows(), x.Wt1_.cols());
    if (! bt1_.rows())
      bt1_ = tensor_type::Zero(x.bt1_.rows(), x.bt1_.cols());
    if (! Wt2_.rows())
      Wt2_ = tensor_type::Zero(x.Wt2_.rows(), x.Wt2_.cols());
    if (! bt2_.rows())
      bt2_ = tensor_type::Zero(x.bt2_.rows(), x.bt2_.cols());


    Wc_ += x.Wc_;
    bc_ += x.bc_;

    Ws1_ += x.Ws1_;
    bs1_ += x.bs1_;
    Ws2_ += x.Ws2_;
    bs2_ += x.bs2_;
    
    Wi1_ += x.Wi1_;
    bi1_ += x.bi1_;
    Wi2_ += x.Wi2_;
    bi2_ += x.bi2_;
    
    Wt1_ += x.Wt1_;
    bt1_ += x.bt1_;    
    Wt2_ += x.Wt2_;
    bt2_ += x.bt2_;    
    
    count_ += x.count_;
    
    return *this;
  }
  
  void clear()
  {
    // embedding
    source_.clear();
    target_.clear();

    Wc_.setZero();
    bc_.setZero();

    Ws1_.setZero();
    bs1_.setZero();
    Ws2_.setZero();
    bs2_.setZero();
    
    Wi1_.setZero();
    bi1_.setZero();
    Wi2_.setZero();
    bi2_.setZero();

    Wt1_.setZero();
    bt1_.setZero();
    Wt2_.setZero();
    bt2_.setZero();
    
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
  
  void initialize(const size_type embedding, const size_type hidden)
  {
    if (hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    
    embedding_ = embedding;
    hidden_    = hidden;
    
    clear();
    
    Wc_ = tensor_type::Zero(embedding, hidden_);
    bc_ = tensor_type::Zero(embedding, 1);
    
    Ws1_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    bs1_ = tensor_type::Zero(hidden_, 1);

    Ws2_ = tensor_type::Zero(hidden_ + hidden_, hidden_);
    bs2_ = tensor_type::Zero(hidden_ + hidden_, 1);

    Wi1_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    bi1_ = tensor_type::Zero(hidden_, 1);
    
    Wi2_ = tensor_type::Zero(hidden_ + hidden_,  hidden_);
    bi2_ = tensor_type::Zero(hidden_ + hidden_, 1);
    
    Wt1_ = tensor_type::Zero(hidden_, embedding_ + embedding_);
    bt1_ = tensor_type::Zero(hidden_, 1);

    Wt2_ = tensor_type::Zero(embedding_ + embedding_, hidden_);
    bt2_ = tensor_type::Zero(embedding_ + embedding_, 1);
    
    count_ = 0;
    shared_ = 0;
  }
  
private:
  void increment(embedding_type& embedding, const embedding_type& x)
  {
    embedding_type::const_iterator siter_end = x.end();
    for (embedding_type::const_iterator siter = x.begin(); siter != siter_end; ++ siter) {
      tensor_type& matrix = embedding[siter->first];
      
      if (! matrix.rows())
	matrix = siter->second;
      else
	matrix += siter->second;
    }
  }
  
  void decrement(embedding_type& embedding, const embedding_type& x)
  {
    embedding_type::const_iterator siter_end = x.end();
    for (embedding_type::const_iterator siter = x.begin(); siter != siter_end; ++ siter) {
      tensor_type& matrix = embedding[siter->first];
      
      if (! matrix.rows())
	matrix = - siter->second;
      else
	matrix += siter->second;
    }
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
  
public:
  // dimension...
  size_type embedding_;
  size_type hidden_;
  
  // embedding
  embedding_type source_;
  embedding_type target_;

  // classification
  tensor_type Wc_;
  tensor_type bc_;

  // straight
  tensor_type Ws1_;
  tensor_type bs1_;
  tensor_type Ws2_;
  tensor_type bs2_;

  // inversion
  tensor_type Wi1_;
  tensor_type bi1_;  
  tensor_type Wi2_;
  tensor_type bi2_;  
  
  // terminal
  tensor_type Wt1_;
  tensor_type bt1_;
  tensor_type Wt2_;
  tensor_type bt2_;
  
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
  
  typedef Gradient gradient_type;
  
  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;
  
  typedef boost::filesystem::path path_type;

  typedef std::vector<bool, std::allocator<bool> > word_unique_type;
  
  Model() : embedding_(0), hidden_(0), scale_(1) {}  
  template <typename Words, typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      scale_(1)
  { initialize(embedding, hidden, words_source, words_target, gen); }
  
  void clear()
  {
    // embedding
    source_.setZero();
    target_.setZero();

    Wc_.setZero();
    bc_.setZero();

    Ws1_.setZero();
    bs1_.setZero();
    Ws2_.setZero();
    bs2_.setZero();
    
    Wi1_.setZero();
    bi1_.setZero();
    Wi2_.setZero();
    bi2_.setZero();

    Wt1_.setZero();
    bt1_.setZero();
    Wt2_.setZero();
    bt2_.setZero();
    
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
		  Words& words_source,
		  Words& words_target,
		  Gen& gen)
  {
    if (hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    
    embedding_ = embedding;
    hidden_    = hidden;
    
    clear();
    
    const size_type vocabulary_size = word_type::allocated();
    
    const double range_e = std::sqrt(6.0 / (embedding_ + 1));
    const double range_c = std::sqrt(6.0 / (embedding_ + hidden_));
    const double range_s = std::sqrt(6.0 / (hidden_ + hidden_ + hidden_));
    const double range_i = std::sqrt(6.0 / (hidden_ + hidden_ + hidden_));
    const double range_t = std::sqrt(6.0 / (hidden + embedding_ + embedding_));

    source_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));

    Wc_ = tensor_type::Zero(embedding, hidden_).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Zero(embedding, 1);
    
    Ws1_ = tensor_type::Zero(hidden_, hidden_ + hidden_).array().unaryExpr(randomize<Gen>(gen, range_s));
    bs1_ = tensor_type::Zero(hidden_, 1);
    Ws2_ = tensor_type::Zero(hidden_ + hidden_,  hidden_).array().unaryExpr(randomize<Gen>(gen, range_s));
    bs2_ = tensor_type::Zero(hidden_ + hidden_, 1);
    
    Wi1_ = tensor_type::Zero(hidden_, hidden_ + hidden_).array().unaryExpr(randomize<Gen>(gen, range_i));
    bi1_ = tensor_type::Zero(hidden_, 1);
    Wi2_ = tensor_type::Zero(hidden_ + hidden_, hidden_).array().unaryExpr(randomize<Gen>(gen, range_i));
    bi2_ = tensor_type::Zero(hidden_ + hidden_, 1);
    
    Wt1_ = tensor_type::Zero(hidden_, embedding_ + embedding_).array().unaryExpr(randomize<Gen>(gen, range_t));
    bt1_ = tensor_type::Zero(hidden_, 1);
    Wt2_ = tensor_type::Zero(embedding_ + embedding_, hidden_).array().unaryExpr(randomize<Gen>(gen, range_t));
    bt2_ = tensor_type::Zero(embedding_ + embedding_, 1);
    
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
    
    scale_ = 1.0;
  }

  void finalize()
  {
    if (scale_ == 1.0) return;
    
    source_ *= scale_;
    target_ *= scale_;

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
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), target_, words_target_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);

    write(rep.path("Ws1.txt.gz"), rep.path("Ws1.bin"), Ws1_);
    write(rep.path("bs1.txt.gz"), rep.path("bs1.bin"), bs1_);
    write(rep.path("Ws2.txt.gz"), rep.path("Ws2.bin"), Ws2_);
    write(rep.path("bs2.txt.gz"), rep.path("bs2.bin"), bs2_);
    
    write(rep.path("Wi1.txt.gz"), rep.path("Wi1.bin"), Wi1_);
    write(rep.path("bi1.txt.gz"), rep.path("bi1.bin"), bi1_);
    write(rep.path("Wi2.txt.gz"), rep.path("Wi2.bin"), Wi2_);
    write(rep.path("bi2.txt.gz"), rep.path("bi2.bin"), bi2_);
    
    write(rep.path("Wt1.txt.gz"), rep.path("Wt1.bin"), Wt1_);
    write(rep.path("bt1.txt.gz"), rep.path("bt1.bin"), bt1_);
    write(rep.path("Wt2.txt.gz"), rep.path("Wt2.bin"), Wt2_);
    write(rep.path("bt2.txt.gz"), rep.path("bt2.bin"), bt2_);
    
    // vocabulary...
    vocab_type vocab;

    const word_type::id_type vocabulary_size = utils::bithack::min(words_source_.size(), words_target_.size());
    
    vocab.open(rep.path("vocab"), vocabulary_size >> 1);
    
    for (word_type::id_type id = 0; id != vocabulary_size; ++ id)
      if (words_source_[id] || words_target_[id]) {
	const word_type word(id);
	
	vocab.insert(word);
      }
    
    vocab.close();
  }
  
  void write_embedding(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix, const word_unique_type& words) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    karma::real_generator<parameter_type, real_policy> float10;

    const word_type::id_type rows = matrix.rows();
    const word_type::id_type cols = std::min(static_cast<size_type>(matrix.cols()), words.size());
    
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

  word_unique_type words_source_;
  word_unique_type words_target_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;

  tensor_type Wc_;
  tensor_type bc_;

  tensor_type Ws1_;
  tensor_type bs1_;
  tensor_type Ws2_;
  tensor_type bs2_;

  tensor_type Wi1_;
  tensor_type bi1_;
  tensor_type Wi2_;
  tensor_type bi2_;
  
  tensor_type Wt1_;
  tensor_type bt1_;
  tensor_type Wt2_;
  tensor_type bt2_;

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

  dict_set_type dicts_;
};

struct ITG
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef Dictionary dictionary_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type tensor_type;

  typedef cicada::Bitext    bitext_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Vocab     vocab_type;
  
  typedef bitext_type::word_type     word_type;
  typedef bitext_type::sentence_type sentence_type;
  
  struct Span
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef uint32_t index_type;
    typedef index_type value_type;
    
    Span(const index_type& first, const index_type& last) : first_(first), last_(last) {}
    template <typename _Index>
    Span(const std::pair<_Index, _Index>& x) : first_(x.first), last_(x.second) {}
    Span() : first_(0), last_(0) {}
    
    bool empty() const { return first_ == last_; }
    size_type size() const { return last_ - first_; }

    friend
    bool operator==(const Span& x, const Span& y)
    {
      return x.first_ == y.first_ && x.last_ == y.last_;
    }

    friend
    bool operator!=(const Span& x, const Span& y)
    {
      return x.first_ != y.first_ || x.last_ != y.last_;
    }

    friend
    std::ostream& operator<<(std::ostream& os, const Span& x)
    {
      os << x.first_ << ".." << x.last_;
      return os;
    }
    
    index_type first_;
    index_type last_;
  };
  
  typedef Span span_type;
  
  struct SpanPair
  {
    typedef Span span_type;
    typedef span_type value_type;
    
    typedef span_type::index_type index_type;
    
    SpanPair(const span_type& source, const span_type& target)
      : source_(source), target_(target) {}
    SpanPair(const index_type& w, const index_type& x, const index_type& y, const index_type& z)
      : source_(w, x), target_(y, z) {}
    template <typename _Span>
    SpanPair(const std::pair<_Span, _Span>& x)
      : source_(x.first), target_(x.second) {}
    SpanPair() : source_(), target_() {}

    bool empty() const { return source_.empty() && target_.empty(); }
    size_type size() const { return source_.size() + target_.size(); }

    friend
    bool operator==(const SpanPair& x, const SpanPair& y)
    {
      return x.source_ == y.source_ && x.target_ == y.target_;
    }

    friend
    bool operator!=(const SpanPair& x, const SpanPair& y)
    {
      return x.source_ != y.source_ || x.target_ != y.target_;
    }

    friend
    std::ostream& operator<<(std::ostream& os, const SpanPair& x)
    {
      os << x.source_ << "-" << x.target_;
      return os;
    }
    
    span_type source_;
    span_type target_;
  };
  
  typedef SpanPair span_pair_type;

  struct Hyperedge
  {
    span_pair_type span_;
    span_pair_type left_;
    span_pair_type right_;

    Hyperedge()
      : span_(), left_(), right_() {}
    Hyperedge(const span_pair_type& span)
      : span_(span), left_(), right_() {}
    Hyperedge(const span_pair_type& span, const span_pair_type& left, const span_pair_type& right)
      : span_(span), left_(left), right_(right) {}
    
    bool aligned() const { return ! span_.source_.empty() && ! span_.target_.empty(); }
    bool terminal() const { return left_.empty() && right_.empty(); }
    bool straight() const { return ! terminal() && left_.target_.last_  == right_.target_.first_; }
    bool inverted() const { return ! terminal() && left_.target_.first_ == right_.target_.last_; }
  };

  typedef Hyperedge hyperedge_type;

  struct State
  {
    typedef State state_type;
    
    State() : span_(), left_(0), right_(0) {}
    State(const span_pair_type& span, const state_type* left=0, const state_type* right=0)
      : span_(span), left_(left), right_(right) {}
    
    bool aligned()  const { return ! span_.source_.empty() && ! span_.target_.empty(); }
    bool terminal() const { return ! left_  && ! right_; }
    bool straight() const { return ! terminal() && left_->span_.target_.last_ == right_->span_.target_.first_; }
    bool inverted() const { return ! terminal() && left_->span_.target_.first_ == right_->span_.target_.last_; }
    
    span_pair_type    span_;
    const state_type* left_;
    const state_type* right_;
    
    // other state related data
    double loss_;
    double cost_;
    
    tensor_type layer_;
    tensor_type layer_norm_;
    tensor_type reconstruction_;
    tensor_type delta_reconstruction_;
    tensor_type delta_;
  };

  typedef State state_type;
  
  struct heap_compare
  {
    // sort by greater item so that we can pop from less items
    bool operator()(const state_type* x, const state_type* y) const
    {
      return x->cost_ > y->cost_;
    }
  };

  typedef utils::chunk_vector<state_type, 1024 * 16 / sizeof(state_type), std::allocator<state_type> > state_set_type;
  typedef utils::vector2<state_type, std::allocator<state_type> > terminal_set_type;

  typedef std::vector<const state_type*, std::allocator<const state_type*> > heap_type;
  typedef std::vector<heap_type, std::allocator<heap_type> > agenda_type;
  typedef utils::bichart<heap_type, std::allocator<heap_type> > chart_type;

  typedef std::vector<hyperedge_type, std::allocator<hyperedge_type> > derivation_type;
  typedef std::vector<const state_type*, std::allocator<const state_type*> > stack_type;

  typedef std::pair<const state_type*, const state_type*> state_pair_type;
  typedef utils::compact_set<state_pair_type,
			     utils::unassigned<state_pair_type>, utils::unassigned<state_pair_type>,
			     utils::hashmurmur3<size_t>, std::equal_to<state_pair_type>,
			     std::allocator<state_pair_type> > state_pair_unique_type;

  struct RestCost
  {
    double cost_;
    double alpha_;
    double beta_;
    
    RestCost()
      : cost_(std::numeric_limits<double>::infinity()),
	alpha_(std::numeric_limits<double>::infinity()),
	beta_(std::numeric_limits<double>::infinity()) {}
  };
  
  typedef RestCost rest_cost_type;
  typedef std::vector<rest_cost_type, std::allocator<rest_cost_type> > rest_cost_set_type;
  
  
  ITG(const dictionary_type& dict_source_target,
      const dictionary_type& dict_target_source,
      const size_type beam)
    : dict_source_target_(dict_source_target),
      dict_target_source_(dict_target_source),
      beam_(beam) {}
  
  const dictionary_type& dict_source_target_;
  const dictionary_type& dict_target_source_;

  size_type beam_;
  
  chart_type     chart_;
  agenda_type    agenda_;
  state_set_type states_;

  terminal_set_type terminals_;

  rest_cost_set_type costs_source_;
  rest_cost_set_type costs_target_;
  
  state_pair_unique_type uniques_;
  stack_type stack_;

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

  void clear()
  {
    chart_.clear();
    agenda_.clear();
    states_.clear();

    terminals_.clear();
    costs_source_.clear();
    costs_target_.clear();
    
    uniques_.clear();
    stack_.clear();
  }
  
  double forward(const sentence_type& source,
		 const sentence_type& target,
		 const model_type& theta)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
#if 0
    std::cerr << "forward source: " << source << std::endl
	      << "forward target: " << target << std::endl;
#endif


    clear();
    
    chart_.reserve(source_size + 1, target_size + 1);
    chart_.resize(source_size + 1, target_size + 1);
    
    agenda_.reserve(source_size + target_size + 1);
    agenda_.resize(source_size + target_size + 1);
    
    terminals_.reserve(source_size + 1, target_size + 1);
    terminals_.resize(source_size + 1, target_size + 1);
    
    costs_source_.reserve(source_size + 1);
    costs_target_.reserve(target_size + 1);

    costs_source_.resize(source_size + 1);
    costs_target_.resize(target_size + 1);

    // initialize terminals...
    forward_terminals(source, target, theta);
    
    // initialize chart...
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = 0; trg <= target_size; ++ trg)
	if (src < source_size || trg < target_size) {
	  // epsilon at target
	  if (src < source_size)
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg)), theta);
	  
	  // epsilon at source
	  if (trg < target_size)
	    forward(source, target, span_pair_type(span_type(src, src), span_type(trg, trg + 1)), theta);
	  
	  // word-pair
	  if (src < source_size && trg < target_size)
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg + 1)), theta);
	}

#if 0
    std::cerr << "biparsing source: " << source << std::endl
	      << "biparsing target: " << target << std::endl;
#endif

    // forward actual forward path
    const size_type length_max = source_size + target_size;
    
    for (size_type length = 1; length != length_max; ++ length) 
      if (! agenda_[length].empty()) {
	
	//std::cerr << "length: " << length << std::endl;

	heap_type& heap = agenda_[length];

	heap_type::iterator hiter_begin = heap.begin();
	heap_type::iterator hiter       = heap.end();
	heap_type::iterator hiter_end   = heap.end();
	
	if (length > 2 && std::distance(hiter_begin, hiter_end) > beam_) {
	  std::make_heap(hiter_begin, hiter_end, heap_compare());
	  
	  for (/**/; hiter_begin != hiter && std::distance(hiter, hiter_end) != beam_; -- hiter)
	    std::pop_heap(hiter_begin, hiter, heap_compare());
	} else
	  hiter = hiter_begin;
	
	// clear uniques
	uniques_.clear();

	// first, put into chart
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const state_type* state = *iter;
	  const span_pair_type& span = state->span_;
	  
	  if (! state->terminal())
	    chart_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(state);
	}
	
	// then, enumerate new hypotheses
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const state_type* state = *iter;
	  
	  const span_pair_type& span = state->span_;
	  
	  // we borrow the notation...
	  const difference_type l = length;
	  const difference_type s = span.source_.first_;
	  const difference_type t = span.source_.last_;
	  const difference_type u = span.target_.first_;
	  const difference_type v = span.target_.last_;
	  
	  const difference_type T = source_size;
	  const difference_type V = target_size;
	    
	  // remember, we have processed only upto length l. thus, do not try to combine with spans greather than l!
	  // also, keep used pair of spans in order to remove duplicates.
	    
	  for (difference_type S = utils::bithack::max(s - l, difference_type(0)); S <= s; ++ S) {
	    const difference_type L = l - (s - S);
	      
	    // straight
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == s); ++ U) {
	      // parent span: StUv
	      // span1: SsUu
	      // span2: stuv
	      
	      const heap_type& states = chart_(S, s, U, u);
	      
	      heap_type::const_iterator siter_end = states.end();
	      for (heap_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(*siter, state)).second)
		  forward(source, target, span_pair_type(S, t, U, v), *siter, state, theta, true);
	    }

	    // inversion
	    for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: StuU
	      // span1: SsvU
	      // span2: stuv

	      const heap_type& states = chart_(S, s, v, U);
		
	      heap_type::const_iterator siter_end = states.end();
	      for (heap_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(*siter, state)).second)
		  forward(source, target, span_pair_type(S, t, u, U), *siter, state, theta, false);
	    }
	  }
	  
	  for (difference_type S = t; S <= utils::bithack::min(t + l, T); ++ S) {
	    const difference_type L = l - (S - t);
	    
	    // inversion
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == t); ++ U) {
	      // parent span: sSUv
	      // span1: stuv
	      // span2: tSUu
	      
	      const heap_type& states = chart_(t, S, U, u);
	      
	      heap_type::const_iterator siter_end = states.end();
	      for (heap_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(state, *siter)).second)
		  forward(source, target, span_pair_type(s, S, U, v), state, *siter, theta, false);
	    }
	    
	    // straight
	    for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: sSuU
	      // span1: stuv
	      // span2: tSvU
		
	      const heap_type& states = chart_(t, S, v, U);
	      
	      heap_type::const_iterator siter_end = states.end();
	      for (heap_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(state, *siter)).second)
		  forward(source, target, span_pair_type(s, S, u, U), state, *siter, theta, true);
	    }
	  }
	}
      }
    
    if (agenda_[length_max].empty())
      return std::numeric_limits<double>::infinity();
    else {
      heap_type& heap = agenda_[length_max];
      
      std::make_heap(heap.begin(), heap.end(), heap_compare());
      
      return heap.front()->loss_;
    }
  }
  
  void forward_terminals(const sentence_type& source,
			 const sentence_type& target,
			 const model_type& theta)
  {
#if 0
    std::cerr << "terminals source: " << source << std::endl
	      << "terminals target: " << target << std::endl;
#endif

    const size_type hidden_size    = theta.hidden_;
    const size_type embedding_size = theta.embedding_;

    const size_type offset_source = 0;
    const size_type offset_target = embedding_size;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	state_type& state = terminals_(src, trg);
	
	const word_type& word_source = (src == 0 ? vocab_type::EPSILON : source[src - 1]);
	const word_type& word_target = (trg == 0 ? vocab_type::EPSILON : target[trg - 1]);

	state.layer_ = (theta.bt1_
			+ (theta.Wt1_.block(0, offset_source, hidden_size, embedding_size)
			   * theta.source_.col(word_source.id())
			   * theta.scale_)
			+ (theta.Wt1_.block(0, offset_target, hidden_size, embedding_size)
			   * theta.target_.col(word_target.id())
			   * theta.scale_)).array().unaryExpr(htanh());
	
	state.layer_norm_ = state.layer_.normalized();
	
	const tensor_type y = (theta.Wt2_ * state.layer_norm_ + theta.bt2_).array().unaryExpr(htanh());
	
	tensor_type y_minus_c = tensor_type::Zero(embedding_size * 2, 1);
	y_minus_c.block(offset_source, 0, embedding_size, 1) = (y.block(offset_source, 0, embedding_size, 1).normalized()
								- theta.source_.col(word_source.id()) * theta.scale_);
	y_minus_c.block(offset_target, 0, embedding_size, 1) = (y.block(offset_target, 0, embedding_size, 1).normalized()
								- theta.target_.col(word_target.id()) * theta.scale_);
	
	state.loss_                 = 0.5 * y_minus_c.squaredNorm();
	state.reconstruction_       = y_minus_c;
	state.delta_reconstruction_ = y.array().unaryExpr(dhtanh()) * y_minus_c.array();
	
	if (src)
	  costs_source_[src - 1].cost_ = std::min(costs_source_[src - 1].cost_, state.loss_);
	if (trg)
	  costs_target_[trg - 1].cost_ = std::min(costs_target_[trg - 1].cost_, state.loss_);
      }

    // estimate rest-costs
    forward_backward(costs_source_);
    forward_backward(costs_target_);
  }

  void forward_backward(rest_cost_set_type& costs)
  {
    const size_type sentence_size = costs.size() - 1;
    
    // forward...
    costs[0].alpha_ = 0;
    for (size_type last = 1; last <= sentence_size; ++ last) {
      const size_type first = last - 1;
      
      costs[last].alpha_ = std::min(costs[last].alpha_, costs[first].alpha_ + costs[first].cost_);
    }
    
    // backward...
    costs[sentence_size].beta_ = 0;
    for (difference_type first = sentence_size - 1; first >= 0; -- first) {
      const size_type last = first + 1;
      
      costs[first].beta_ = std::min(costs[first].beta_, costs[first].cost_ + costs[last].beta_);
    }
  }
  
  // terminal rules
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const span_pair_type& span,
	       const model_type& theta)
  {
    state_type& state = allocate(span);
    
    state = terminals_(span.source_.empty() ? 0 : span.source_.first_ + 1, span.target_.empty() ? 0 : span.target_.first_ + 1);
    
    state.span_ = span;
    state.cost_ = state.loss_ + std::max(costs_source_[span.source_.first_].alpha_ + costs_source_[span.source_.last_].beta_,
					 costs_target_[span.target_.first_].alpha_ + costs_target_[span.target_.last_].beta_);
    
    chart_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(&state);
    agenda_[span.size()].push_back(&state);
  }
  
  // binary rules
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const span_pair_type& span,
	       const state_type* state1,
	       const state_type* state2,
	       const model_type& theta,
	       const bool straight)
  {
    const size_type hidden_size    = theta.hidden_;
    const size_type embedding_size = theta.embedding_;

    const size_type offset_left  = 0;
    const size_type offset_right = hidden_size;

    state_type& state = allocate(span, state1, state2);
    
    const span_pair_type& span1 = state1->span_;
    const span_pair_type& span2 = state2->span_;

    //std::cerr << "span: " << span << " left: " << span1 << " right: " << span2 << std::endl;

    const tensor_type& Wr1 = (straight ? theta.Ws1_ : theta.Wi1_);
    const tensor_type& br1 = (straight ? theta.bs1_ : theta.bi1_);
    const tensor_type& Wr2 = (straight ? theta.Ws2_ : theta.Wi2_);
    const tensor_type& br2 = (straight ? theta.bs2_ : theta.bi2_);
    
    state.layer_ = (br1
		    + Wr1.block(0, offset_left, hidden_size, hidden_size) * state1->layer_
		    + Wr1.block(0, offset_right, hidden_size, hidden_size) * state2->layer_).array().unaryExpr(htanh());
    
    state.layer_norm_ = state.layer_.normalized();
    
    const tensor_type y = (Wr2 * state.layer_norm_ + br2).array().unaryExpr(htanh());
    
    tensor_type y_minus_c = tensor_type::Zero(hidden_size * 2, 1);
    y_minus_c.block(offset_left, 0, hidden_size, 1)  = (y.block(offset_left, 0, hidden_size, 1).normalized()
							- state1->layer_);
    y_minus_c.block(offset_right, 0, hidden_size, 1) = (y.block(offset_right, 0, hidden_size, 1).normalized()
							- state2->layer_);
    
    state.loss_                 = 0.5 * y_minus_c.squaredNorm() + state1->loss_ + state2->loss_;
    state.cost_                 = state.loss_ + std::max(costs_source_[span.source_.first_].alpha_
							 + costs_source_[span.source_.last_].beta_,
							 costs_target_[span.target_.first_].alpha_
							 + costs_target_[span.target_.last_].beta_);
    state.reconstruction_       = y_minus_c;
    state.delta_reconstruction_ = y.array().unaryExpr(dhtanh()) * y_minus_c.array();
    
    //chart_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(&state);
    agenda_[span.size()].push_back(&state);
  }
  
  template <typename Gen>
  void backward(const sentence_type& source,
		const sentence_type& target,
		const model_type& theta,
		gradient_type& gradient,
		Gen& gen)
  {
    const heap_type& heap = agenda_[source.size() + target.size()];
    
#if 0
    std::cerr << "backward source: " << source << std::endl
	      << "backward target: " << target << std::endl;
#endif

    if (heap.empty()) return;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type hidden_size    = theta.hidden_;
    const size_type embedding_size = theta.embedding_;
    
    const size_type offset_source = 0;
    const size_type offset_target = embedding_size;
    
    const size_type offset_left  = 0;
    const size_type offset_right = hidden_size;

    ++ gradient.count_;

    const_cast<tensor_type&>(heap.front()->delta_) = tensor_type::Zero(hidden_size, 1);
    
    stack_.clear();
    stack_.push_back(heap.front());
    
    while (! stack_.empty()) {
      const state_type* state = stack_.back();
      stack_.pop_back();

      if (state->terminal()) {
	const span_pair_type& span = state->span_;

	const word_type& word_source = (span.source_.empty() ? vocab_type::EPSILON : source[span.source_.first_]);
	const word_type& word_target = (span.target_.empty() ? vocab_type::EPSILON : target[span.target_.first_]);
	
	// increment delta from reconstruction
	const_cast<tensor_type&>(state->delta_).array() += (state->layer_.array().unaryExpr(dhtanh())
							    * (theta.Wt2_.transpose() * state->delta_reconstruction_).array());
	
	gradient.Wt1_.block(0, offset_source, hidden_size, embedding_size) += (state->delta_
									       * theta.source_.col(word_source.id()).transpose()
									       * theta.scale_);
	gradient.Wt1_.block(0, offset_target, hidden_size, embedding_size) += (state->delta_
									       * theta.target_.col(word_target.id()).transpose()
									       * theta.scale_);
	gradient.bt1_ += state->delta_;
	
	gradient.Wt2_ += state->delta_reconstruction_ * state->layer_norm_.transpose();
	gradient.bt2_ += state->delta_reconstruction_;
	
	gradient.source(word_source) += (theta.Wt1_.block(0, offset_source, hidden_size, embedding_size).transpose()
					 * state->delta_
					 - state->reconstruction_.block(offset_source, 0, embedding_size, 1));
	gradient.target(word_target) += (theta.Wt1_.block(0, offset_target, hidden_size, embedding_size).transpose()
					 * state->delta_
					 - state->reconstruction_.block(offset_target, 0, embedding_size, 1));
      } else {
	stack_.push_back(state->left_);
	stack_.push_back(state->right_);
	
	const bool straight = state->straight();

	const tensor_type& Wr1 = (straight ? theta.Ws1_ : theta.Wi1_);
	const tensor_type& Wr2 = (straight ? theta.Ws2_ : theta.Wi2_);
	
	tensor_type& dWr1 = (straight ? gradient.Ws1_ : gradient.Wi1_);
	tensor_type& dbr1 = (straight ? gradient.bs1_ : gradient.bi1_);
	tensor_type& dWr2 = (straight ? gradient.Ws2_ : gradient.Wi2_);
	tensor_type& dbr2 = (straight ? gradient.bs2_ : gradient.bi2_);

	// increment delta from reconstruction
	const_cast<tensor_type&>(state->delta_).array() += (state->layer_.array().unaryExpr(dhtanh())
							    * (Wr2.transpose() * state->delta_reconstruction_).array());

	
	
	dWr1.block(0, offset_left,  hidden_size, hidden_size) += state->delta_ * state->left_->layer_.transpose();
	dWr1.block(0, offset_right, hidden_size, hidden_size) += state->delta_ * state->right_->layer_.transpose();
	dbr1 += state->delta_;
	
	dWr2 += state->delta_reconstruction_ * state->layer_norm_.transpose();
	dbr2 += state->delta_reconstruction_;
	
	tensor_type& delta_left  = const_cast<tensor_type&>(state->left_->delta_);
	tensor_type& delta_right = const_cast<tensor_type&>(state->right_->delta_);
	
	delta_left  = (state->left_->layer_.array().unaryExpr(dhtanh())
		       * (Wr1.block(0, offset_left, hidden_size, hidden_size).transpose() * state->delta_
			  - state->reconstruction_.block(offset_left, 0, hidden_size, 1)).array());
	delta_right = (state->right_->layer_.array().unaryExpr(dhtanh())
		       * (Wr1.block(0, offset_right, hidden_size, hidden_size).transpose() * state->delta_
			  - state->reconstruction_.block(offset_right, 0, hidden_size, 1)).array());
      }
    }
  }
  
  void derivation(const sentence_type& source,
		  const sentence_type& target,
		  derivation_type& d)
  {
    d.clear();
    
    const heap_type& heap = agenda_[source.size() + target.size()];
    
    if (heap.empty()) return;
    
    stack_.clear();
    stack_.push_back(heap.front());

    while (! stack_.empty()) {
      const state_type* state = stack_.back();
      stack_.pop_back();
      
      if (state->terminal())
	d.push_back(hyperedge_type(state->span_));
      else {
	d.push_back(hyperedge_type(state->span_, state->left_->span_, state->right_->span_));
	
	// we will push in right-to-left order...!
	stack_.push_back(state->right_);
	stack_.push_back(state->left_);
      }
    }
  }

  state_type& allocate(const span_pair_type& span,
		       const state_type* state1 = 0,
		       const state_type* state2 = 0)
  {
    states_.push_back(state_type(span, state1, state2));
    
    return states_.back();
  }
};

struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model     model_type;
  typedef Gradient  gradient_type;

  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& embedding,
	       const size_type& hidden,
	       const double& lambda,
	       const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      lambda_(lambda),
      eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    
    source_ = tensor_type::Zero(embedding_, vocabulary_size);
    target_ = tensor_type::Zero(embedding_, vocabulary_size);
    
    Wc_ = tensor_type::Zero(embedding, hidden_);
    bc_ = tensor_type::Zero(embedding, 1);
    
    Ws1_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    bs1_ = tensor_type::Zero(hidden_, 1);

    Ws2_ = tensor_type::Zero(hidden_ + hidden_, hidden_);
    bs2_ = tensor_type::Zero(hidden_ + hidden_, 1);

    Wi1_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    bi1_ = tensor_type::Zero(hidden_, 1);
    
    Wi2_ = tensor_type::Zero(hidden_ + hidden_,  hidden_);
    bi2_ = tensor_type::Zero(hidden_ + hidden_, 1);
    
    Wt1_ = tensor_type::Zero(hidden_, embedding_ + embedding_);
    bt1_ = tensor_type::Zero(hidden_, 1);

    Wt2_ = tensor_type::Zero(embedding_ + embedding_, hidden_);
    bt2_ = tensor_type::Zero(embedding_ + embedding_, 1);
  }

  
  void operator()(model_type& theta,
		  const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;
    
    const double scale = 1.0 / gradient.count_;

    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     const_cast<tensor_type&>(source_),
	     siter->second,
	     scale);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     const_cast<tensor_type&>(target_),
	     titer->second,
	     scale);

    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, scale, false);
    
    update(theta.Ws1_, const_cast<tensor_type&>(Ws1_), gradient.Ws1_, scale, lambda_ != 0.0);
    update(theta.bs1_, const_cast<tensor_type&>(bs1_), gradient.bs1_, scale, false);
    update(theta.Ws2_, const_cast<tensor_type&>(Ws2_), gradient.Ws2_, scale, lambda_ != 0.0);
    update(theta.bs2_, const_cast<tensor_type&>(bs2_), gradient.bs2_, scale, false);

    update(theta.Wi1_, const_cast<tensor_type&>(Wi1_), gradient.Wi1_, scale, lambda_ != 0.0);
    update(theta.bi1_, const_cast<tensor_type&>(bi1_), gradient.bi1_, scale, false);
    update(theta.Wi2_, const_cast<tensor_type&>(Wi2_), gradient.Wi2_, scale, lambda_ != 0.0);
    update(theta.bi2_, const_cast<tensor_type&>(bi2_), gradient.bi2_, scale, false);
    
    update(theta.Wt1_, const_cast<tensor_type&>(Wt1_), gradient.Wt1_, scale, lambda_ != 0.0);
    update(theta.bt1_, const_cast<tensor_type&>(bt1_), gradient.bt1_, scale, false);
    update(theta.Wt2_, const_cast<tensor_type&>(Wt2_), gradient.Wt2_, scale, lambda_ != 0.0);
    update(theta.bt2_, const_cast<tensor_type&>(bt2_), gradient.bt2_, scale, false);
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

  template <typename Theta, typename GradVar, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale) const
  {
    for (int row = 0; row != g.rows(); ++ row) 
      if (g(row, 0) != 0) {
	G(row, word.id()) +=  g(row, 0) * g(row, 0) * scale * scale;
	
	const double rate = eta0_ / std::sqrt(double(1.0) + G(row, word.id()));
	const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	
	theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
      }
  }
  
  size_type embedding_;
  size_type hidden_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;

  // classification
  tensor_type Wc_;
  tensor_type bc_;

  // straight
  tensor_type Ws1_;
  tensor_type bs1_;
  tensor_type Ws2_;
  tensor_type bs2_;

  // inversion
  tensor_type Wi1_;
  tensor_type bi1_;  
  tensor_type Wi2_;
  tensor_type bi2_;  
  
  // terminal
  tensor_type Wt1_;
  tensor_type bt1_;
  tensor_type Wt2_;
  tensor_type bt2_;
};

struct LearnSGD
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model     model_type;
  typedef Gradient  gradient_type;

  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef model_type::tensor_type tensor_type;

  LearnSGD(const double& lambda,
	   const double& eta0)
    : lambda_(lambda),
      eta0_(eta0),
      epoch_(0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
  }
  
  void operator()(model_type& theta,
		  const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    //++ const_cast<size_type&>(epoch_);

    const double scale = 1.0 / gradient.count_;
    
    if (lambda_ != 0.0) {
      const double eta = eta0_ / (epoch_ + 1);
      
      theta.scale_ *= 1.0 - eta * lambda_;
    }
    
    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     siter->second,
	     scale,
	     theta.scale_);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     titer->second,
	     scale,
	     theta.scale_);

    update(theta.Wc_, gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, gradient.bc_, scale, false);

    update(theta.Ws1_, gradient.Ws1_, scale, lambda_ != 0.0);
    update(theta.bs1_, gradient.bs1_, scale, false);
    update(theta.Ws2_, gradient.Ws2_, scale, lambda_ != 0.0);
    update(theta.bs2_, gradient.bs2_, scale, false);

    update(theta.Wi1_, gradient.Wi1_, scale, lambda_ != 0.0);
    update(theta.bi1_, gradient.bi1_, scale, false);
    update(theta.Wi2_, gradient.Wi2_, scale, lambda_ != 0.0);
    update(theta.bi2_, gradient.bi2_, scale, false);
    
    update(theta.Wt1_, gradient.Wt1_, scale, lambda_ != 0.0);
    update(theta.bt1_, gradient.bt1_, scale, false);
    update(theta.Wt2_, gradient.Wt2_, scale, lambda_ != 0.0);
    update(theta.bt2_, gradient.bt2_, scale, false);
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
    
    theta.col(word.id()).noalias() -= (eta * scale / theta_scale) * g;
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

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

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

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
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
    
    model_type theta(dimension_embedding, dimension_hidden, sources, targets, generator);
    
    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());
      
      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta.read_embedding(embedding_source_file, embedding_target_file);
    }
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, lambda, eta0),
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
    bitext_type     bitext_;
    derivation_type derivation_;
    
    bitext_derivation_type() : id_(size_type(-1)), bitext_(), derivation_() {}
    bitext_derivation_type(const size_type& id,
			   const bitext_type& bitext,
			   const derivation_type& derivation)
      : id_(id), bitext_(bitext), derivation_(derivation) {}
    
    void swap(bitext_derivation_type& x)
    {
      std::swap(id_, x.id_);
      bitext_.swap(x.bitext_);
      derivation_.swap(x.derivation_);
    }

    void clear()
    {
      id_ = size_type(-1);
      bitext_.clear();
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

struct OutputDerivation : OutputMapReduce
{
  typedef std::vector<std::string, std::allocator<std::string> > stack_type;
	  
  OutputDerivation(const path_type& path,
		   queue_type& queue)
    : path_(path), queue_(queue) {}
  
  void operator()()
  {
    if (path_.empty()) {
      bitext_derivation_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_set_type bitexts;
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
				   ? bitext.bitext_.source_[diter->span_.source_.first_]
				   : vocab_type::EPSILON);
	const word_type& target = (! diter->span_.target_.empty()
				   ? bitext.bitext_.target_[diter->span_.target_.first_]
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
  
  path_type   path_;
  queue_type& queue_;
  
  stack_type stack_;
};

struct OutputAlignment : OutputMapReduce
{
  typedef cicada::Alignment alignment_type;

  OutputAlignment(const path_type& path_source_target,
		  const path_type& path_target_source,
		  queue_type& queue)
    : path_source_target_(path_source_target),
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
      bitext_set_type bitexts;
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
	output(*os_source_target, bitext.id_, bitext.bitext_.source_, bitext.bitext_.target_, alignment_);
    }
    
    if (os_target_source) {
      alignment_.inverse();
      
      std::sort(alignment_.begin(), alignment_.end());

      if (moses_mode)
	*os_target_source << alignment_ << '\n';
      else
	output(*os_target_source, bitext.id_, bitext.bitext_.target_, bitext.bitext_.source_, alignment_);
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

  path_type  path_source_target_;
  path_type  path_target_source_;
  queue_type& queue_;
  
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
	    gradients_.push_back(gradient_type(theta_.embedding_, theta_.hidden_));
	    grad = &gradients_.back();
	  }
	  
	  grad->clear();
	  
	  const size_type first = batch * batch_size_;
	  const size_type last  = utils::bithack::min(first + batch_size_, bitexts_.size());
	  
	  for (size_type id = first; id != last; ++ id) {
	    const sentence_type& source = bitexts_[id].source_;
	    const sentence_type& target = bitexts_[id].target_;
	    
	    bitext_derivation.id_     = id;
	    bitext_derivation.bitext_ = bitexts_[id];
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
  for (size_type batch = 0; batch != batches_size; ++ batch)
    batches[batch] = batch;
  
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
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;

    std::auto_ptr<boost::progress_display> progress(debug
						    ? new boost::progress_display(batches_size, std::cerr, "", "", "")
						    : 0);
    
    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);

    boost::thread output_derivation(output_derivation_type(! derivation_file.empty() && dump_mode
							   ? add_suffix(derivation_file, iter_tag)
							   : path_type(),
							   queue_derivation));
    boost::thread output_alignment(output_alignment_type(! alignment_source_target_file.empty() && dump_mode
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
	typename batch_set_type::iterator iter_end = std::min(biter + (batch_size << 5), biter_end);
	
	std::random_shuffle(biter, iter_end);
	biter = iter_end;
      }
    }
    
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
      
      bitext_derivation.id_     = bitext_id;
      bitext_derivation.bitext_ = bitexts_[bitext_id];
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

  task_type::queue_type   mapper(256 * threads);
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
  workers_dump.add_thread(new boost::thread(output_derivation_type(! derivation_file.empty()
								   ? derivation_file
								   : path_type(),
								   queue_derivation)));
  workers_dump.add_thread(new boost::thread(output_alignment_type(! alignment_source_target_file.empty()
								  ? alignment_source_target_file
								  : path_type(),
								  ! alignment_target_source_file.empty()
								  ? alignment_target_source_file
								  : path_type(),
								  queue_alignment)));

  if (debug)
    std::cerr << "max derivation" << std::endl;

  utils::resource start;
  
  for (size_type i = 0; i != bitexts.size(); ++ i)
    mapper.push(i);
  
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
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
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
