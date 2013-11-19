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
#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/vector2.hpp"
#include "utils/sampler.hpp"
#include "utils/resource.hpp"
#include "utils/simple_vector.hpp"

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
	   const int alignment) 
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
    bn_ += x.bn_;
    
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
    
    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), hidden_ + embedding_);
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
    
    source_ = tensor_type::Zero(embedding_,     vocabulary_size);
    target_ = tensor_type::Zero(embedding_ + 1, vocabulary_size);
  }

  void assign(const gradient_type& x)
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    gradient_embedding_type::const_iterator siter_end = x.source_.end();
    for (gradient_embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter)
      source_.col(siter->first.id()) = siter->second;
    
    gradient_embedding_type::const_iterator titer_end = x.target_.end();
    for (gradient_embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer)
      target_.col(titer->first.id()) = titer->second;
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
  
  Model() : embedding_(0), hidden_(0), alignment_(0) {}  
  template <typename Words, typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	const int alignment,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment)
  { initialize(embedding, hidden, alignment, words_source, words_target, gen); }
  
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
		  const int alignment,
		  Words& words_source,
		  Words& words_target,
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
    
    const double range_e = std::sqrt(6.0 / (embedding_ + 1));
    const double range_t = std::sqrt(6.0 / (hidden_ + embedding_ + embedding_));
    const double range_a = std::sqrt(6.0 / (hidden_ + hidden_ + embedding_));
    const double range_n = std::sqrt(6.0 / (hidden_ + hidden_ + embedding_));
    const double range_i = std::sqrt(6.0 / (hidden_ + 1));

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
 
    source_ = tensor_type::Zero(embedding_,     vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_ = tensor_type::Zero(embedding_ + 1, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_.row(embedding_).setZero();
    
    Wt_ = tensor_type::Zero(embedding_, hidden_ + embedding_).array().unaryExpr(randomize<Gen>(gen, range_t));
    bt_ = tensor_type::Zero(embedding_, 1);
    
    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), hidden_ + embedding_).array().unaryExpr(randomize<Gen>(gen, range_a));
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1);

    Wn_ = tensor_type::Zero(hidden_, hidden_ + embedding_).array().unaryExpr(randomize<Gen>(gen, range_n));
    bn_ = tensor_type::Zero(hidden_, 1);
    
    Wi_ = tensor_type::Zero(hidden_, 1).array().unaryExpr(randomize<Gen>(gen, range_i));
  }

  size_type shift(const difference_type source_size,
		  const difference_type target_size,
		  const difference_type prev,
		  const difference_type next) const
  {
    const difference_type prev_adjusted = utils::bithack::branch(prev >= source_size + 2, prev - source_size - 2, prev);
    const difference_type next_adjusted = utils::bithack::branch(next >= source_size + 2, next - source_size - 2, next);
    
    return utils::bithack::min(utils::bithack::max(next_adjusted - prev_adjusted + difference_type(alignment_),
						   difference_type(0)),
			       difference_type(alignment_ * 2));
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
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), target_, words_target_);
    
    write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);
    
    write(rep.path("Wa.txt.gz"), rep.path("Wa.bin"), Wa_);
    write(rep.path("ba.txt.gz"), rep.path("ba.bin"), ba_);

    write(rep.path("Wn.txt.gz"), rep.path("Wn.bin"), Wn_);
    write(rep.path("bn.txt.gz"), rep.path("bn.bin"), bn_);
    
    write(rep.path("Wi.txt.gz"), rep.path("Wi.bin"), Wi_);

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
  int       alignment_;

  word_unique_type words_source_;
  word_unique_type words_target_;
  
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

class State
{
public:
  typedef char*  pointer;
  typedef int    index_type;
  typedef float  parameter_type;
  typedef State  state_type;
  
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;

  static const size_type offset_prev   = 0;
  static const size_type offset_index  = offset_prev  + sizeof(pointer);
  static const size_type offset_score  = offset_index + sizeof(index_type);
  static const size_type offset_matrix = (offset_score + sizeof(parameter_type) + 15) & (~15);
  
public:
  State(pointer __base) : base_(__base) {}
  State() : base_(0) {}
  
public:
  static inline
  size_type size(const size_type state_size) { return offset_matrix + sizeof(parameter_type) * state_size; }
  
  bool empty() const { return ! base_; }

  inline       state_type& prev()       { return *reinterpret_cast<state_type*>(base_ + offset_prev); }
  inline const state_type& prev() const { return *reinterpret_cast<const state_type*>(base_ + offset_prev); }
  
  inline       index_type& index()       { return *reinterpret_cast<index_type*>(base_ + offset_index); }
  inline const index_type& index() const { return *reinterpret_cast<const index_type*>(base_ + offset_index); }

  inline       parameter_type& score()       { return *reinterpret_cast<parameter_type*>(base_ + offset_score); }
  inline const parameter_type& score() const { return *reinterpret_cast<const parameter_type*>(base_ + offset_score); }
  
  inline parameter_type* matrix()       { return reinterpret_cast<parameter_type*>(base_ + offset_matrix); }
  inline parameter_type* matrix() const { return reinterpret_cast<parameter_type*>(const_cast<char*>(base_) + offset_matrix); }

public:
  friend
  size_t hash_value(State const& x)
  {
    return boost::hash<pointer>()(x.base_);
  }

  friend
  bool operator==(const State& x, const State& y)
  {
    return x.base_ == y.base_;
  }

  friend
  bool operator!=(const State& x, const State& y)
  {
    return x.base_ != y.base_;
  }
  
public:
  pointer base_;
};

class StateAllocator : public std::allocator<char>
{
public:
  typedef State      state_type;
  typedef char*      pointer;
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;
  
  typedef std::allocator<char> allocator_type;
  typedef utils::simple_vector<pointer, std::allocator<pointer> > state_set_type;
  
public:
  static const size_type chunk_size  = 1024 * 4;
  static const size_type chunk_mask  = chunk_size - 1;
  
public:

  StateAllocator()
    : states(), state_iterator(0),
      cache(0),
      state_size(0),
      state_alloc_size(0),
      state_chunk_size(0) {}
  
  StateAllocator(size_type __state_size)
    : states(), state_iterator(0),
      cache(0),
      state_size(__state_size),
      state_alloc_size(0),
      state_chunk_size(0)
  {
    if (state_size != 0) {
      // sizeof(char*) aligned size..
      const size_type pointer_size = sizeof(pointer);
      const size_type pointer_mask = ~(pointer_size - 1);
      
      state_alloc_size = (state_size + pointer_size - 1) & pointer_mask;
      state_chunk_size = state_alloc_size * chunk_size;
    }
  }
  
  StateAllocator(const StateAllocator& x) 
    : allocator_type(static_cast<const allocator_type&>(x)),
      states(), state_iterator(0),
      cache(0),
      state_size(x.state_size),
      state_alloc_size(x.state_alloc_size),
      state_chunk_size(x.state_chunk_size) {}
  
  StateAllocator& operator=(const StateAllocator& x)
  { 
    clear();
      
    static_cast<allocator_type&>(*this) = static_cast<const allocator_type&>(x);
      
    state_size = x.state_size;
    state_alloc_size = x.state_alloc_size;
    state_chunk_size = x.state_chunk_size;

    return *this;
  }
  
  ~StateAllocator() { clear(); }
  
public:  
  state_type allocate()
  {
    if (state_size == 0) return 0;
    
    if (cache) {
      pointer state = cache;
      cache = *reinterpret_cast<pointer*>(state);
      
      // clear buffer...
      std::fill(state, state + state_alloc_size, 0);
      
      return state_type(state);
    }
    
    const size_type chunk_pos = state_iterator & chunk_mask;
    
    if (chunk_pos == 0) {
      states.push_back(allocator_type::allocate(state_chunk_size));
      std::uninitialized_fill(states.back(), states.back() + state_chunk_size, 0);
      
      //std::cerr << "allocated: " << ((void*) states.back()) << " size: " << states.size() << std::endl;
    }
    
    ++ state_iterator;
    
    return state_type(states.back() + chunk_pos * state_alloc_size);
  }
  
  void deallocate(const state_type& state)
  {
    if (state.empty() || state_size == 0 || states.empty()) return;
    
    *reinterpret_cast<pointer*>(const_cast<state_type&>(state).base_) = cache;
    cache = state.base_;
  }
  
  state_type clone(const state_type& state)
  {
    if (state_size == 0)
      return state_type();
      
    state_type state_new = allocate();
      
    std::copy(state.base_, state.base_ + state_alloc_size, state_new.base_);
      
    return state_new;
  }

  void assign(size_type __state_size)
  {
    if (state_size != __state_size)
      *this = StateAllocator(__state_size);
    else
      clear();
  }
  
  void clear()
  {
    //std::cerr << "allocated states size: " << states.size() << std::endl;

    state_set_type::iterator siter_end = states.end();
    for (state_set_type::iterator siter = states.begin(); siter != siter_end; ++ siter) {
      //std::cerr << "de-allocated: " << ((void*) *siter) << std::endl;
      allocator_type::deallocate(*siter, state_chunk_size);
    }
    
    states.clear();
    state_iterator = 0;
    cache = 0;
  }
  
private:  
  state_set_type states;
  size_type state_iterator;
  pointer   cache;

  size_type state_size;
  size_type state_alloc_size;
  size_type state_chunk_size;
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

struct HMM
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model      model_type;
  typedef Gradient   gradient_type;
  typedef Dictionary dictionary_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type    tensor_type;

  typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> back_type;
  typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> visited_type;

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Bitext    bitext_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;

  typedef State          state_type;
  typedef StateAllocator state_allocator_type;

  typedef std::vector<state_type, std::allocator<state_type> > heap_type;
  typedef std::vector<heap_type, std::allocator<heap_type> > heap_set_type;
  
  typedef utils::unordered_map<state_type, state_type, boost::hash<state_type>, std::equal_to<state_type>,
			       std::allocator<std::pair<const state_type, state_type> > >::type state_set_type;
  typedef std::vector<state_set_type, std::allocator<state_set_type> > state_map_type;

  struct heap_compare
  {
    // comparison by less, so that we can pop in greater order
    bool operator()(const state_type& x, const state_type& y) const
    {
      return x.score() < y.score();
    }
  };
  
  HMM(const dictionary_type& dict,
      const size_type& beam)
    : dict_(dict), beam_(beam) {}
  
  const dictionary_type& dict_;

  size_type beam_;
  
  tensor_type alpha_;
  tensor_type beta_;
  tensor_type trans_;
  tensor_type visited_;
  
  tensor_type forw_;
  back_type   back_;
  tensor_type score_;

  tensor_type layer_trans_;
  
  tensor_type delta_beta_;
  tensor_type delta_alpha_;
  tensor_type delta_trans_;
  
  state_allocator_type state_allocator_;
  heap_set_type  heaps_;
  state_map_type states_;

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

  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const model_type& theta,
	       alignment_type& alignment)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    const size_type offset_word   = 0;
    const size_type offset_matrix = theta.embedding_;
    
    state_allocator_.assign(state_type::size(theta.hidden_ + theta.embedding_));
    
    alignment.clear();
    
    heaps_.clear();
    heaps_.resize(target_size + 2);
    
    state_type state_bos = state_allocator_.allocate();
    
    state_bos.prev() = state_type();
    state_bos.index() = 0;
    state_bos.score() = 0.0;
    
    matrix_type(state_bos.matrix(), theta.hidden_, 1) = theta.Wi_;
    
    heaps_[0].push_back(state_bos);
    
    for (size_type trg = 1; trg != target_size + 2; ++ trg) {
      const word_type target_next = (trg == 0
				     ? vocab_type::BOS
				     : (trg == target_size + 1
					? vocab_type::EOS
					: target[trg - 1]));
      
      heap_type& heap      = heaps_[trg - 1];
      heap_type& heap_next = heaps_[trg];
      
      const size_type next_first = utils::bithack::branch(trg == target_size + 1, source_size + 1, size_type(1));
      const size_type next_last  = utils::bithack::branch(trg == target_size + 1, source_size + 2, source_size + 1);
      
      // perform pruning
      heap_type::iterator hiter_begin = heap.begin();
      heap_type::iterator hiter       = heap.end();
      heap_type::iterator hiter_end   = heap.end();
      
      if (heap.size() > beam_) {
	for (/**/; hiter_begin != hiter && hiter_end - hiter != beam_; -- hiter)
	  std::pop_heap(hiter_begin, hiter, heap_compare());
	
	// deallocate unused states
	for (/**/; hiter_begin != hiter; ++ hiter_begin)
	  state_allocator_.deallocate(*hiter_begin);
	
      } else
	hiter = hiter_begin;
      
      for (/**/; hiter != hiter_end; ++ hiter) {
	const state_type& state = *hiter;
	
	const size_type prev = state.index();
	const word_type source_aligned_prev = (prev >= source_size + 2
					       ? (prev - source_size - 2 == 0
						  ? vocab_type::BOS
						  : (prev - source_size - 2 == source_size + 1
						     ? vocab_type::EOS
						     : source[prev - source_size - 2 - 1]))
					       : (prev == 0
						  ? vocab_type::BOS
						  : (prev == source_size + 1
						     ? vocab_type::EOS
						     : source[prev - 1])));
	const word_type source_none_prev = (prev >= source_size + 2
					    ? vocab_type::EPSILON
					    : (prev == 0
					       ? vocab_type::BOS
					       : (prev == source_size + 1
						  ? vocab_type::EOS
						  : source[prev - 1])));

	// alignment into aligned...
	for (size_type next = next_first; next != next_last; ++ next) {
	  const word_type source_next = (next == 0
					 ? vocab_type::BOS
					 : (next == source_size + 1
					    ? vocab_type::EOS
					    : source[next - 1]));
	  
	  const size_type shift = theta.shift(source_size, target_size, prev, next);

	  state_type state_next = state_allocator_.allocate();
	  state_next.prev() = state;
	  state_next.index() = next;
	  
	  matrix_type layer_alpha(state_next.matrix(), theta.hidden_, 1);
	  matrix_type layer_trans(state_next.matrix() + theta.hidden_, theta.embedding_, 1);
	  
	  layer_alpha
	    = ((theta.Wa_.block(theta.hidden_ * shift, offset_word, theta.hidden_, theta.embedding_)
		* theta.source_.col(source_aligned_prev.id()))
	       
	       + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
		  * matrix_type(state.matrix(), theta.hidden_, 1))
	       
	       + theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)).array().unaryExpr(htanh());

	  layer_trans = ((theta.Wt_.block(0, offset_word, theta.embedding_, theta.embedding_)
			  * theta.source_.col(source_next.id()))
			 
			 + (theta.Wt_.block(0, offset_matrix, theta.embedding_, theta.hidden_)
			    * layer_alpha)
			 
			 + theta.bt_).array().unaryExpr(htanh());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans
				+ theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  state_next.score() = score + state.score();

	  heap_next.push_back(state_next);
	  std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	}
	
	// none...
	if (trg != target_size + 1) {
	  const size_type next = utils::bithack::branch(prev >= source_size + 2, prev, prev + source_size + 2);
	  
	  state_type state_next = state_allocator_.allocate();
	  state_next.prev() = state;
	  state_next.index() = next;
	  
	  matrix_type layer_alpha(state_next.matrix(), theta.hidden_, 1);
	  matrix_type layer_trans(state_next.matrix() + theta.hidden_, theta.embedding_, 1);
	  
	  layer_alpha
	    = ((theta.Wn_.block(0, 0, theta.hidden_, theta.embedding_)
		* theta.source_.col(source_none_prev.id()))
	       
	       + (theta.Wn_.block(0, theta.embedding_, theta.hidden_, theta.hidden_)
		  * matrix_type(state.matrix(), theta.hidden_, 1))
	       
	       + theta.bn_.block(0, 0, theta.hidden_, 1)).array().unaryExpr(htanh());
	  
	  layer_trans = ((theta.Wt_.block(0, 0, theta.embedding_, theta.embedding_)
			  * theta.source_.col(vocab_type::EPSILON.id()))
			 
			 + (theta.Wt_.block(0, theta.embedding_, theta.embedding_, theta.hidden_)
			    * layer_alpha)
			 
			 + theta.bt_).array().unaryExpr(htanh());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans
				+ theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  state_next.score() = score + state.score();
	  
	  heap_next.push_back(state_next);
	  std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	}
      }
    }
    
    heap_type& heap = heaps_[target_size + 1];
    
    // perform pruning!
    {
      heap_type::iterator hiter_begin = heap.begin();
      heap_type::iterator hiter       = heap.end();
      heap_type::iterator hiter_end   = heap.end();
      
      for (/**/; hiter_begin != hiter && hiter_end - hiter != beam_; -- hiter)
	std::pop_heap(hiter_begin, hiter, heap_compare());
      
      // deallocate unused states
      for (/**/; hiter_begin != hiter; ++ hiter_begin)
	state_allocator_.deallocate(*hiter_begin);
    }
    
    state_type state = heap.back();
    
    for (size_type trg = target_size + 1; trg >= 1; -- trg) {
      const size_type src = state.index();
      
      if (trg < target_size + 1 && src < source_size + 1)
	alignment.push_back(std::make_pair(src - 1, trg - 1));
      
      state = state.prev();
    }
  }
  
  template <typename Gen>
  double backward(const sentence_type& source,
		  const sentence_type& target,
		  const model_type& theta,
		  gradient_type& gradient,
		  Gen& gen)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    const size_type offset_word   = 0;
    const size_type offset_matrix = theta.embedding_;
    
    states_.clear();
    states_.resize(target_size + 2);

    {
      state_set_type& states = states_[target_size + 1];
      const heap_type& heap = heaps_[target_size + 1];
      
      heap_type::const_iterator hiter_end = heap.end();
      for (heap_type::const_iterator hiter = std::max(heap.begin(), hiter_end - beam_); hiter != hiter_end; ++ hiter) {
	state_type buffer = state_allocator_.allocate();
	
	matrix_type(buffer.matrix(), theta.hidden_, 1).setZero();
	
	states[*hiter] = buffer;
      }
    }

    double log_likelihood = 0.0;
    
    for (size_type trg = target_size + 1; trg > 0; -- trg) {
      const state_set_type& states_next = states_[trg];
      state_set_type& states_prev = states_[trg - 1];

      const word_type target_next(trg == 0
				  ? vocab_type::BOS
				  : (trg == target_size + 1
				     ? vocab_type::EOS
				     : target[trg - 1]));
      
      state_set_type::const_iterator siter_end = states_next.end();
      for (state_set_type::const_iterator siter = states_next.begin(); siter != siter_end; ++ siter) {
	const state_type state_next(siter->first);
	const state_type state_prev(state_next.prev());
	
	const matrix_type beta_next(siter->second.matrix(), theta.hidden_, 1);
	
	state_set_type::iterator piter = states_prev.find(state_prev);
	if (piter == states_prev.end()) {
	  state_type buffer = state_allocator_.allocate();
	  
	  matrix_type(buffer.matrix(), theta.hidden_, 1).setZero();
	  
	  piter = states_prev.insert(std::make_pair(state_prev, buffer)).first;
	}
	
	matrix_type beta_prev(piter->second.matrix(), theta.hidden_, 1);
	
	const size_type next = state_next.index();
	const size_type prev = state_prev.index();
	
	const word_type source_next(next >= source_size + 2
				    ? vocab_type::EPSILON
				    : (next == 0
				       ? vocab_type::BOS
				       : (next == source_size + 1
					  ? vocab_type::EOS
					  : source[next - 1])));
	const word_type source_aligned_prev = (prev >= source_size + 2
					       ? (prev - source_size - 2 == 0
						  ? vocab_type::BOS
						  : (prev - source_size - 2 == source_size + 1
						     ? vocab_type::EOS
						     : source[prev - source_size - 2 - 1]))
					       : (prev == 0
						  ? vocab_type::BOS
						  : (prev == source_size + 1
						     ? vocab_type::EOS
						     : source[prev - 1])));
	const word_type source_none_prev = (prev >= source_size + 2
					    ? vocab_type::EPSILON
					    : (prev == 0
					       ? vocab_type::BOS
					       : (prev == source_size + 1
						  ? vocab_type::EOS
						  : source[prev - 1])));
	
	const matrix_type layer_alpha(state_next.matrix(), theta.hidden_, 1);
	const matrix_type layer_trans(state_next.matrix() + theta.hidden_, theta.embedding_, 1);
	
	const word_type target_sampled = dict_.draw(source_none_prev, gen);
	
	const double score_c = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans
				+ theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	const double score_m = (theta.target_.col(target_sampled.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans
				+ theta.target_.col(target_sampled.id()).block(theta.embedding_, 0, 1, 1))(0, 0);

	const double score_noise_c = 0.0 + dict_.logprob(source_none_prev, target_next);
	const double score_noise_m = 0.0 + dict_.logprob(source_none_prev, target_sampled);
	
	const double z_c = utils::mathop::logsum(score_c, score_noise_c);
	const double z_m = utils::mathop::logsum(score_m, score_noise_m);
	
	const double logprob_c = score_c - z_c;
	const double logprob_m = score_m - z_m;
	
	const double logprob_noise_c = score_noise_c - z_c;
	const double logprob_noise_m = score_noise_m - z_m;
	
	log_likelihood += logprob_c + logprob_noise_m;
	
	const double loss_c = - 1.0 + std::exp(logprob_c);
	const double loss_m =         std::exp(logprob_m);
	
	tensor_type& dembedding_c = gradient.target(target_next);
	tensor_type& dembedding_m = gradient.target(target_sampled);
	
	dembedding_c.block(0, 0, theta.embedding_, 1).array() += loss_c * layer_trans.array();
	dembedding_c.block(theta.embedding_, 0, 1, 1).array() += loss_c;
	dembedding_m.block(0, 0, theta.embedding_, 1).array() += loss_m * layer_trans.array();
	dembedding_m.block(theta.embedding_, 0, 1, 1).array() += loss_m;
	
	delta_trans_ = (layer_trans.array().unaryExpr(dhtanh())
			* (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1) * loss_c
			   + theta.target_.col(target_sampled.id()).block(0, 0, theta.embedding_, 1) * loss_m).array());
	
	gradient.Wt_.block(0, offset_word, theta.embedding_, theta.embedding_)
	  += delta_trans_ * theta.source_.col(source_next.id()).transpose();
	gradient.Wt_.block(0, offset_matrix, theta.embedding_, theta.hidden_)
	  += delta_trans_ * layer_alpha.transpose();
	gradient.bt_ += delta_trans_;
	
	delta_beta_ = theta.Wt_.transpose() * delta_trans_;
	
	gradient.source(source_next) += delta_beta_.block(0, 0, theta.embedding_, 1);
	
	delta_alpha_ = (layer_alpha.array().unaryExpr(dtanh())
			* (delta_beta_.block(theta.embedding_, 0, theta.hidden_, 1) + beta_next).array());
	
	if (next >= source_size + 2) {
	  gradient.Wn_.block(0, offset_word, theta.hidden_, theta.embedding_)
	    += (delta_alpha_ * theta.source_.col(source_none_prev.id()).transpose());
	  
	  gradient.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
	    += (delta_alpha_ * matrix_type(state_prev.matrix(), theta.hidden_, 1).transpose());
	  
	  gradient.bn_ += delta_alpha_;

	  gradient.source(source_none_prev) += (theta.Wn_.block(0, offset_word, theta.hidden_, theta.embedding_).transpose()
						* delta_alpha_);
	  
	  beta_prev += (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_).transpose() * delta_alpha_);
	} else {
	  const size_type shift = theta.shift(source_size, target_size, prev, next);
	  
	  gradient.Wa_.block(theta.hidden_ * shift, offset_word, theta.hidden_, theta.embedding_)
	    += (delta_alpha_ * theta.source_.col(source_aligned_prev.id()).transpose());
	  
	  gradient.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
	    += (delta_alpha_ * matrix_type(state_prev.matrix(), theta.hidden_, 1).transpose());
	  
	  gradient.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
	    += delta_alpha_;
	  
	  gradient.source(source_aligned_prev) += (theta.Wa_.block(theta.hidden_ * shift, offset_word, theta.hidden_, theta.embedding_).transpose()
						   * delta_alpha_);
	  
	  beta_prev += (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_).transpose()
			* delta_alpha_);
	}
      }
    }
    
    // final propagation...
    if (states_[0].empty())
      throw std::runtime_error("no initial state?");
    if (states_[0].size() != 1)
      throw std::runtime_error("multiple initial states?");

    const matrix_type beta(states_[0].begin()->second.matrix(), theta.embedding_ + theta.hidden_, 1);
    
    gradient.source(vocab_type::BOS) += beta.block(0, 0, theta.embedding_, 1);
    gradient.Wi_ += beta.block(theta.embedding_, 0, theta.hidden_, 1);
    
    return log_likelihood;
  }
    
  double viterbi(const sentence_type& source,
		 const sentence_type& target,
		 const model_type& theta,
		 alignment_type& alignment)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type offset_word   = 0;
    const size_type offset_matrix = theta.embedding_;
    
    state_allocator_.assign(state_type::size(theta.hidden_));
    
    alignment.clear();
    
    heaps_.clear();
    heaps_.resize(target_size + 2);
    
    state_type state_bos = state_allocator_.allocate();
    
    state_bos.prev() = state_type();
    state_bos.index() = 0;
    state_bos.score() = 0.0;
    
    matrix_type(state_bos.matrix(), theta.hidden_, 1) = theta.Wi_;
    
    heaps_[0].push_back(state_bos);
    
    for (size_type trg = 1; trg != target_size + 2; ++ trg) {
      const word_type target_next = (trg == 0
				     ? vocab_type::BOS
				     : (trg == target_size + 1
					? vocab_type::EOS
					: target[trg - 1]));
      
      heap_type& heap      = heaps_[trg - 1];
      heap_type& heap_next = heaps_[trg];
      
      const size_type next_first = utils::bithack::branch(trg == target_size + 1, source_size + 1, size_type(1));
      const size_type next_last  = utils::bithack::branch(trg == target_size + 1, source_size + 2, source_size + 1);
      
      // perform pruning
      heap_type::iterator hiter_begin = heap.begin();
      heap_type::iterator hiter       = heap.end();
      heap_type::iterator hiter_end   = heap.end();
      
      if (heap.size() > beam_) {
	for (/**/; hiter_begin != hiter && hiter_end - hiter != beam_; -- hiter)
	  std::pop_heap(hiter_begin, hiter, heap_compare());
	
	// deallocate unused states
	for (/**/; hiter_begin != hiter; ++ hiter_begin)
	  state_allocator_.deallocate(*hiter_begin);
	
      } else
	hiter = hiter_begin;
      
      for (/**/; hiter != hiter_end; ++ hiter) {
	const state_type& state = *hiter;
	
	const size_type prev = state.index();
	const word_type source_aligned_prev = (prev >= source_size + 2
					       ? (prev - source_size - 2 == 0
						  ? vocab_type::BOS
						  : (prev - source_size - 2 == source_size + 1
						     ? vocab_type::EOS
						     : source[prev - source_size - 2 - 1]))
					       : (prev == 0
						  ? vocab_type::BOS
						  : (prev == source_size + 1
						     ? vocab_type::EOS
						     : source[prev - 1])));
	const word_type source_none_prev = (prev >= source_size + 2
					    ? vocab_type::EPSILON
					    : (prev == 0
					       ? vocab_type::BOS
					       : (prev == source_size + 1
						  ? vocab_type::EOS
						  : source[prev - 1])));

	// alignment into aligned...
	for (size_type next = next_first; next != next_last; ++ next) {
	  const word_type source_next = (next == 0
					 ? vocab_type::BOS
					 : (next == source_size + 1
					    ? vocab_type::EOS
					    : source[next - 1]));
	  
	  const size_type shift = theta.shift(source_size, target_size, prev, next);

	  state_type state_next = state_allocator_.allocate();
	  state_next.prev() = state;
	  state_next.index() = next;
	  
	  matrix_type layer_alpha(state_next.matrix(), theta.hidden_, 1);
	  
	  layer_alpha
	    = ((theta.Wa_.block(theta.hidden_ * shift, offset_word, theta.hidden_, theta.embedding_)
		* theta.source_.col(source_aligned_prev.id()))
	       
	       + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
		  * matrix_type(state.matrix(), theta.hidden_, 1))
	       
	       + theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)).array().unaryExpr(htanh());

	  layer_trans_ = ((theta.Wt_.block(0, offset_word, theta.embedding_, theta.embedding_)
			   * theta.source_.col(source_next.id()))
			  
			  + (theta.Wt_.block(0, offset_matrix, theta.embedding_, theta.hidden_)
			     * layer_alpha)
			  
			  + theta.bt_).array().unaryExpr(htanh());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_
				+ theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  state_next.score() = score + state.score();

	  heap_next.push_back(state_next);
	  std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	}
	
	// none...
	if (trg != target_size + 1) {
	  const size_type next = utils::bithack::branch(prev >= source_size + 2, prev, prev + source_size + 2);
	  
	  state_type state_next = state_allocator_.allocate();
	  state_next.prev() = state;
	  state_next.index() = next;
	  
	  matrix_type layer_alpha(state_next.matrix(), theta.hidden_, 1);
	  
	  layer_alpha
	    = ((theta.Wn_.block(0, offset_word, theta.hidden_, theta.embedding_)
		* theta.source_.col(source_none_prev.id()))
	       
	       + (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
		  * matrix_type(state.matrix(), theta.hidden_, 1))
	       
	       + theta.bn_.block(0, 0, theta.hidden_, 1)).array().unaryExpr(htanh());
	  
	  layer_trans_ = ((theta.Wt_.block(0, offset_word, theta.embedding_, theta.embedding_)
			   * theta.source_.col(vocab_type::EPSILON.id()))
			  
			  + (theta.Wt_.block(0, offset_matrix, theta.embedding_, theta.hidden_)
			     * layer_alpha)
			  
			  + theta.bt_).array().unaryExpr(htanh());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_
				+ theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	  
	  state_next.score() = score + state.score();
	  
	  heap_next.push_back(state_next);
	  std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	}
      }
    }
    
#if 0
    std::cerr << "source: " << source << std::endl
	      << "target: " << target << std::endl;
#endif

    heap_type& heap = heaps_[target_size + 1];
    
    std::pop_heap(heap.begin(), heap.end(), heap_compare());
    
    state_type state = heap.back();
    const double score = state.score();
    
    for (size_type trg = target_size + 1; trg >= 1; -- trg) {
      const size_type src = state.index();

      if (trg < target_size + 1 && src < source_size + 1)
	alignment.push_back(std::make_pair(src - 1, trg - 1));
      
      state = state.prev();
    }

    return score;
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
	       const int alignment,
	       const double& lambda,
	       const double& lambda2,
	       const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment),
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

  
  void operator()(model_type& theta,
		  const gradient_type& gradient,
		  const embedding_type& embedding) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     const_cast<tensor_type&>(source_),
	     siter->second,
	     embedding.target_.col(siter->first.id()),
	     false);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     const_cast<tensor_type&>(target_),
	     titer->second,
	     embedding.source_.col(titer->first.id()),
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
      
      const double rate = eta0_ / std::sqrt(double(G_(i, j)));
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
    learning_rate(const double& eta0) : eta0_(eta0) {}

    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return (x == 0.0 ? 0.0 : eta0_ / std::sqrt(double(x)));
    }

    const double& eta0_;
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
      theta.array() -= g.array() * G.array().unaryExpr(learning_rate(eta0_));
    }
  }

  template <typename Theta, typename GradVar, typename Grad, typename GradCross>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const Eigen::MatrixBase<GradCross>& c,
	      const bool bias_last=false) const
  {
    if (lambda2_ > 0.0) {
      for (int row = 0; row != g.rows() - bias_last; ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0);
	  
	  const double rate = eta0_ / std::sqrt(double(G(row, word.id())));
	  const double f = theta(row, word.id()) - rate * g(row, 0);
	  const double value = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	  
	  const double shared = c(row, 0);
	  const double diff = value - shared;
	  
	  theta(row, word.id()) = utils::mathop::sgn(diff) * std::max(0.0, std::fabs(diff) - rate * lambda2_) + shared;
	}
    } else {
      for (int row = 0; row != g.rows() - bias_last; ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0);
	  
	  const double rate = eta0_ / std::sqrt(double(G(row, word.id())));
	  const double f = theta(row, word.id()) - rate * g(row, 0);
	  
	  theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
    }
    
    if (bias_last) {
      const int row = g.rows() - 1;
      
      if (g(row, 0) != 0) {
	G(row, word.id()) += g(row, 0) * g(row, 0);
	theta(row, word.id()) -= eta0_ * g(row, 0) / std::sqrt(double(G(row, word.id())));
      }
    }
  }
  
  size_type embedding_;
  size_type hidden_;
  int       alignment_;
  
  double lambda_;
  double lambda2_;
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

typedef Model      model_type;
typedef Dictionary dictionary_type;

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

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
int alignment = 8;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 128;
int beam_size = 10;
int cutoff = 3;
double lambda = 1e-5;
double lambda2 = 0.01;
double eta0 = 1;

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

    if (batch_size <= 0)
      throw std::runtime_error("invalid batch size");

    if (int(giza_mode) + moses_mode > 1)
      throw std::runtime_error("either giza style output or moses style output");

    if (int(giza_mode) + moses_mode == 0)
      moses_mode = true;
        
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
        
    dictionary_type dict_source_target;
    dictionary_type dict_target_source;
    
    read_data(source_file, target_file, bitexts, dict_source_target, dict_target_source);

    if (debug)
      std::cerr << "# of sentences: " << bitexts.size() << std::endl;
    
    const dictionary_type::dict_type::word_set_type& sources = dict_target_source[cicada::Vocab::EPSILON].words_;
    const dictionary_type::dict_type::word_set_type& targets = dict_source_target[cicada::Vocab::EPSILON].words_;

    model_type theta_source_target(dimension_embedding, dimension_hidden, alignment, sources, targets, generator);
    model_type theta_target_source(dimension_embedding, dimension_hidden, alignment, targets, sources, generator);

    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());
      
      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta_source_target.read_embedding(embedding_source_file, embedding_target_file);
    }
    
    const size_t cols = utils::bithack::min(utils::bithack::min(theta_source_target.source_.cols(),
								theta_source_target.target_.cols()),
					    utils::bithack::min(theta_target_source.source_.cols(),
								theta_target_source.target_.cols()));
    
    theta_source_target.source_.block(0, 0, dimension_embedding, cols)
      = theta_target_source.target_.block(0, 0, dimension_embedding, cols);
    theta_source_target.target_.block(0, 0, dimension_embedding, cols)
      = theta_target_source.source_.block(0, 0, dimension_embedding, cols);
    
    if (iteration > 0)
      learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, alignment, lambda, lambda2, eta0),
		   bitexts,
		   dict_source_target,
		   dict_target_source,
		   theta_source_target,
		   theta_target_source);

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
    
    if (bitext.alignment_.empty() || bitext.bitext_.source_.empty() || bitext.bitext_.target_.empty()) {
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

  typedef output_map_reduce_type::queue_type queue_alignment_type;
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
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 const model_type& theta_source_target,
		 const model_type& theta_target_source,
		 queue_type& queue,
		 queue_alignment_type& queue_source_target,
		 queue_alignment_type& queue_target_source,
		 counter_type& counter)
    : bitexts_(bitexts),
      theta_source_target_(theta_source_target),
      theta_target_source_(theta_target_source),
      queue_(queue),
      queue_source_target_(queue_source_target),
      queue_target_source_(queue_target_source),
      counter_(counter),
      hmm_source_target_(dict_source_target, beam_size),
      hmm_target_source_(dict_target_source, beam_size),
      gradient_source_target_(theta_source_target.embedding_,
			      theta_source_target.hidden_,
			      theta_source_target.alignment_),
      gradient_target_source_(theta_target_source.embedding_,
			      theta_target_source.hidden_,
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
      
      bitext_source_target.id_ = sentence_id;
      bitext_source_target.bitext_.source_ = bitext.source_;
      bitext_source_target.bitext_.target_ = bitext.target_;
      bitext_source_target.alignment_.clear();
      
      bitext_target_source.id_ = sentence_id;
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
					 gradient_source_target_,
					 generator);
	
	log_likelihood_target_source_
	  += hmm_target_source_.backward(bitext.target_,
					 bitext.source_,
					 theta_target_source_,
					 gradient_target_source_,
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
    gradient_source_target_.clear();
    gradient_target_source_.clear();
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
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
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

  Learner learner_source_target = learner;
  Learner learner_target_source = learner;

  Embedding embedding_source_target(theta_source_target.embedding_);
  Embedding embedding_target_source(theta_target_source.embedding_);
  
  task_set_type tasks(threads, task_type(bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta_source_target,
					 theta_target_source,
					 mapper,
					 queue_source_target,
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

      embedding_source_target.assign(tasks.front().gradient_source_target_);
      embedding_target_source.assign(tasks.front().gradient_target_source_);
      
      // update model parameters
      learner_source_target(theta_source_target, tasks.front().gradient_source_target_, embedding_target_source);
      learner_target_source(theta_target_source, tasks.front().gradient_target_source_, embedding_source_target);
    }
    
    utils::resource end;
    
    queue_source_target.push(output_map_reduce_type::value_type());
    queue_target_source.push(output_map_reduce_type::value_type());
    
    if (debug && ((num_text / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;
    
    if (debug)
      std::cerr << "log-likelihood P(target | source) (per sentence): " << (log_likelihood_source_target / samples) << std::endl
		<< "log-likelihood P(target | source) (per word): "     << (log_likelihood_source_target / words_target) << std::endl
		<< "log-likelihood P(source | target) (per sentence): " << (log_likelihood_target_source / samples) << std::endl
		<< "log-likelihood P(source | target) (per word): "     << (log_likelihood_target_source / words_source) << std::endl;

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
  theta_source_target.finalize();
  theta_target_source.finalize();
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
      hmm_source_target_(dict_source_target, beam_size),
      hmm_target_source_(dict_target_source, beam_size),
      log_likelihood_source_target_(0),
      log_likelihood_target_source_(0) {}

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
	log_likelihood_source_target_
	  += hmm_source_target_.viterbi(bitext.source_, bitext.target_, theta_source_target_, bitext_source_target.alignment_);
	log_likelihood_target_source_
	  += hmm_target_source_.viterbi(bitext.target_, bitext.source_, theta_target_source_, bitext_target_source.alignment_);
      }
      
      // reduce alignment
      queue_source_target_.push_swap(bitext_source_target);
      queue_target_source_.push_swap(bitext_target_source);
    }
  }

  void clear()
  {
    log_likelihood_source_target_ = 0;
    log_likelihood_target_source_ = 0;
  }

  const bitext_set_type& bitexts_;
  const model_type& theta_source_target_;
  const model_type& theta_target_source_;
  
  queue_type&           queue_;
  queue_alignment_type& queue_source_target_;
  queue_alignment_type& queue_target_source_;
  
  hmm_type hmm_source_target_;
  hmm_type hmm_target_source_;
  
  double        log_likelihood_source_target_;
  double        log_likelihood_target_source_;

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

  task_type::queue_type   mapper(256 * threads);
  
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

  utils::resource start;
  
  // actually run...
  size_type num_text = 0;
  for (size_type id = 0; id != bitexts.size(); ++ id) {
    mapper.push(id);
    
    ++ num_text;
    if (debug) {
      if (num_text % DEBUG_DOT == 0)
	std::cerr << '.';
      if (num_text % DEBUG_LINE == 0)
	std::cerr << '\n';
    }
  }
  
  if (debug && ((num_text / DEBUG_DOT) % DEBUG_WRAP))
    std::cerr << std::endl;
  
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
    ("alignment",           po::value<int>(&alignment)->default_value(alignment),                     "alignment model size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("beam",              po::value<int>(&beam_size)->default_value(beam_size),   "histogram beam size")
    ("cutoff",            po::value<int>(&cutoff)->default_value(cutoff),         "cutoff count for vocabulary (<= 1 to keep all)")
    ("lambda",            po::value<double>(&lambda)->default_value(lambda),      "regularization constant")
    ("lambda2",           po::value<double>(&lambda2)->default_value(lambda2),    "regularization constant for bilingual agreement")
    ("eta0",              po::value<double>(&eta0)->default_value(eta0),          "\\eta_0 for decay")

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
