//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_RNN_ALIGNMENT_IMPL__HPP__
#define __CICADA_RNN_ALIGNMENT_IMPL__HPP__ 1

#include <cstdlib>
#include <cmath>
#include <climits>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>


#include <boost/functional/hash/hash.hpp>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"
#include "cicada/alignment.hpp"
#include "cicada/bitext.hpp"

#include "utils/alloc_vector.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/compact_map.hpp"
#include "utils/compact_set.hpp"
#include "utils/mathop.hpp"
#include "utils/unordered_map.hpp"
#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/simple_vector.hpp"

#include "codec/lz4.hpp"

#include <boost/random.hpp>

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
  
  Gradient() : embedding_(0), hidden_(0), alignment_(0), count_(0), shared_(0) {}
  Gradient(const size_type& embedding,
	   const size_type& hidden,
	   const size_type& alignment) 
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment),
      count_(0),
      shared_(0)
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
    if (! bn_.rows())
      bn_ = tensor_type::Zero(x.bn_.rows(), x.bn_.cols());

    if (! bi_.rows())
      bi_ = tensor_type::Zero(x.bi_.rows(), x.bi_.cols());
    
    Wt_ -= x.Wt_;
    bt_ -= x.bt_;
    
    Wa_ -= x.Wa_;
    ba_ -= x.ba_;

    Wn_ -= x.Wn_;
    bn_ -= x.bn_;
    
    bi_ -= x.bi_;

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

    if (! Wa_.rows())
      Wa_ = tensor_type::Zero(x.Wa_.rows(), x.Wa_.cols());
    if (! ba_.rows())
      ba_ = tensor_type::Zero(x.ba_.rows(), x.ba_.cols());

    if (! Wn_.rows())
      Wn_ = tensor_type::Zero(x.Wn_.rows(), x.Wn_.cols());
    if (! bn_.rows())
      bn_ = tensor_type::Zero(x.bn_.rows(), x.bn_.cols()); 

    if (! bi_.rows())
      bi_ = tensor_type::Zero(x.bi_.rows(), x.bi_.cols());

    Wt_ += x.Wt_;
    bt_ += x.bt_;
    
    Wa_ += x.Wa_;
    ba_ += x.ba_;

    Wn_ += x.Wn_;
    bn_ += x.bn_;

    bi_ += x.bi_;

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
    
    Wa_.setZero();
    ba_.setZero();

    Wn_.setZero();
    bn_.setZero();
    
    bi_.setZero();

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
      embedding = tensor_type::Zero(embedding_ + 1, 1);
    
    return embedding;
  }

  
  void initialize(const size_type embedding, const size_type hidden, const size_type alignment)
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
    
    bi_ = tensor_type::Zero(hidden_, 1);

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

public:
  typedef std::vector<char, std::allocator<char> > buffer_type;

  buffer_type buffer_;
  
  void encode(std::string& encoded) const
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    buffer.clear();
    
    {
      boost::iostreams::filtering_ostream os;
      os.push(codec::lz4_compressor());
      os.push(boost::iostreams::back_insert_device<buffer_type>(buffer));
      
      write(os);
    }

    encoded = std::string(buffer.begin(), buffer.end());
  }

  void decode(const std::string& encoded) 
  {
    boost::iostreams::filtering_istream is;
    is.push(codec::lz4_decompressor());
    is.push(boost::iostreams::array_source(&(*encoded.begin()), encoded.size()));
    
    read(is);
  }

  friend
  std::ostream& operator<<(std::ostream& os, const Gradient& x)
  {
    x.write(os);
    return os;
  }

  friend
  std::istream& operator>>(std::istream& is, Gradient& x)
  {
    x.read(is);
    return is;
  }

private:
  void write(std::ostream& os) const
  {
    os.write((char*) &embedding_, sizeof(size_type));
    os.write((char*) &hidden_,    sizeof(size_type));
    os.write((char*) &alignment_, sizeof(size_type));
    os.write((char*) &count_,     sizeof(size_type));
    
    write(os, Wt_);
    write(os, bt_);

    write(os, Wa_);
    write(os, ba_);

    write(os, Wn_);
    write(os, bn_);

    write(os, bi_);
    
    write(os, source_, false);
    write(os, target_, true);
  }
  
  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &alignment_, sizeof(size_type));
    is.read((char*) &count_,     sizeof(size_type));
    
    read(is, Wt_);
    read(is, bt_);
    
    read(is, Wa_);
    read(is, ba_);

    read(is, Wn_);
    read(is, bn_);

    read(is, bi_);

    read(is, source_, false);
    read(is, target_, true);
  }

  void write(std::ostream& os, const embedding_type& embedding, const bool bias_last) const
  {
    const size_type size = embedding.size();
    
    os.write((char*) &size, sizeof(size_type));
    
    embedding_type::const_iterator eiter_end = embedding.end();
    for (embedding_type::const_iterator eiter = embedding.begin(); eiter != eiter_end; ++ eiter) {
      const size_type word_size = eiter->first.size();
      
      os.write((char*) &word_size, sizeof(size_type));
      os.write((char*) &(*eiter->first.begin()), word_size);
      os.write((char*) eiter->second.data(), sizeof(tensor_type::Scalar) * eiter->second.rows());
    }
  }

  void read(std::istream& is, embedding_type& embedding, const bool bias_last)
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    
    embedding.clear();
    
    size_type size = 0;
    
    is.read((char*) &size, sizeof(size_type));
    
    for (size_type i = 0; i != size; ++ i) {
      size_type word_size = 0;
      is.read((char*) &word_size, sizeof(size_type));
      
      buffer.resize(word_size);
      is.read((char*) &(*buffer.begin()), word_size);
      
      tensor_type& matrix = embedding[word_type(buffer.begin(), buffer.end())];
      
      matrix.resize(embedding_ + bias_last, 1);
      
      is.read((char*) matrix.data(), sizeof(tensor_type::Scalar) * matrix.rows());
    }
  }

  void write(std::ostream& os, const tensor_type& matrix) const
  {
    const tensor_type::Index rows = matrix.rows();
    const tensor_type::Index cols = matrix.cols();

    os.write((char*) &rows, sizeof(tensor_type::Index));
    os.write((char*) &cols, sizeof(tensor_type::Index));
    
    os.write((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
  }
  
  void read(std::istream& is, tensor_type& matrix)
  {
    tensor_type::Index rows = 0;
    tensor_type::Index cols = 0;
    
    is.read((char*) &rows, sizeof(tensor_type::Index));
    is.read((char*) &cols, sizeof(tensor_type::Index));
    
    matrix.resize(rows, cols);
    
    is.read((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
  }

public:
  // dimension...
  size_type embedding_;
  size_type hidden_;
  size_type alignment_;
  
  // embedding
  embedding_type source_;
  embedding_type target_;
  
  tensor_type Wt_;
  tensor_type bt_;
  
  tensor_type Wa_;
  tensor_type ba_;

  tensor_type Wn_;
  tensor_type bn_;

  tensor_type bi_;
  
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
    
    source_ = tensor_type::Zero(embedding_,     vocabulary_size);
    target_ = tensor_type::Zero(embedding_ + 1, vocabulary_size);
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
  
  Model() : embedding_(0), hidden_(0), alignment_(0), scale_(1) {}  
  template <typename Words, typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	const size_type& alignment,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      alignment_(alignment),
      scale_(1)
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
    
    bi_.setZero();
    
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
		  const size_type alignment,
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

    bi_ = tensor_type::Zero(hidden_, 1).array().unaryExpr(randomize<Gen>(gen, range_i));

    scale_ = 1.0;
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
    if (scale_ == 1.0) return;
    
    source_ *= scale_;
    target_.block(0, 0, embedding_, target_.cols()) *= scale_;

    scale_ = 1.0;    
  }

  Model& operator+=(const Model& x)
  {
    if (scale_ != x.scale_)
      throw std::runtime_error("different scaling");

    source_ += x.source_;
    target_ += x.target_;
    
    Wt_ += x.Wt_;
    bt_ += x.bt_;

    Wa_ += x.Wa_;
    ba_ += x.ba_;

    Wn_ += x.Wn_;
    bn_ += x.bn_;

    bi_ += x.bi_;
    
    return *this;
  }

  Model& operator*=(const double& x)
  {
    source_ *= x;
    target_ *= x;
    
    Wt_ *= x;
    bt_ *= x;

    Wa_ *= x;
    ba_ *= x;

    Wn_ *= x;
    bn_ *= x;

    bi_ *= x;
    
    return *this;
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

	if (word.id() >= source_.cols())
	  source_.conservativeResize(Eigen::NoChange, word.id() + 1);
	if (word.id() >= words_source_.size())
	  words_source_.resize(word.id() + 1, false);
	
	source_.col(word.id()).block(0, 0, embedding_, 1)
	  = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), embedding_, 1);
	
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
	
	if (boost::fusion::get<1>(parsed).size() != embedding_)
	  throw std::runtime_error("invalid embedding size");
	
	const word_type word = boost::fusion::get<0>(parsed);

	if (word.id() >= target_.cols())
	  target_.conservativeResize(Eigen::NoChange, word.id() + 1);
	if (word.id() >= words_target_.size())
	  words_target_.resize(word.id() + 1, false);

	target_.col(word.id()).block(0, 0, embedding_, 1)
	  = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), embedding_, 1);
	
	words_target_[word.id()] = true;
      }
    }
  }

  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);    
    rep["alignment"] = utils::lexical_cast<std::string>(alignment_);    
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), rep.path("vocab-source"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), rep.path("vocab-target"), target_, words_target_);
    
    write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);
    
    write(rep.path("Wa.txt.gz"), rep.path("Wa.bin"), Wa_);
    write(rep.path("ba.txt.gz"), rep.path("ba.bin"), ba_);

    write(rep.path("Wn.txt.gz"), rep.path("Wn.bin"), Wn_);
    write(rep.path("bn.txt.gz"), rep.path("bn.bin"), bn_);
    
    write(rep.path("bi.txt.gz"), rep.path("bi.bin"), bi_);
  }

  void write_embedding(const path_type& path_text,
		       const path_type& path_binary,
		       const path_type& path_vocab,
		       const tensor_type& matrix,
		       const word_unique_type& words) const
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
  
public:
  typedef std::vector<char, std::allocator<char> > buffer_type;

  buffer_type buffer_;

  void encode(std::string& encoded) const
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    buffer.clear();
    
    {
      boost::iostreams::filtering_ostream os;
      os.push(codec::lz4_compressor());
      os.push(boost::iostreams::back_insert_device<buffer_type>(buffer));
      
      write(os);
    }
    
    encoded = std::string(buffer.begin(), buffer.end());
  }

  void decode(const std::string& encoded) 
  {
    boost::iostreams::filtering_istream is;
    is.push(codec::lz4_decompressor());
    is.push(boost::iostreams::array_source(&(*encoded.begin()), encoded.size()));
    
    read(is);
  }

  friend
  std::ostream& operator<<(std::ostream& os, const Model& x)
  {
    x.write(os);
    return os;
  }

  friend
  std::istream& operator>>(std::istream& is, Model& x)
  {
    x.read(is);
    return is;
  }

private:
  void write(std::ostream& os) const
  {
    os.write((char*) &embedding_, sizeof(size_type));
    os.write((char*) &hidden_,    sizeof(size_type));
    os.write((char*) &alignment_, sizeof(size_type));
    os.write((char*) &scale_,     sizeof(double));
    
    write(os, Wt_);
    write(os, bt_);

    write(os, Wa_);
    write(os, ba_);

    write(os, Wn_);
    write(os, bn_);
    
    write(os, bi_);
    
    write_embedding(os, source_, false);
    write_embedding(os, target_, true);
  }

  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &alignment_, sizeof(size_type));
    is.read((char*) &scale_,     sizeof(double));

    read(is, Wt_);
    read(is, bt_);

    read(is, Wa_);
    read(is, ba_);

    read(is, Wn_);
    read(is, bn_);

    read(is, bi_);
    
    read_embedding(is, source_, false);
    read_embedding(is, target_, true);

    // checking...
  }

  void write_embedding(std::ostream& os, const tensor_type& embedding, const bool bias_last) const
  {
    const size_type rows = embedding.rows();
    const size_type cols = embedding.cols();
    
    os.write((char*) &cols, sizeof(size_type));
    
    for (word_type::id_type id = 0; id != cols; ++ id) {
      const word_type word(id);
      
      const size_type word_size = word.size();
      
      os.write((char*) &word_size, sizeof(size_type));
      os.write((char*) &(*word.begin()), word_size);
      os.write((char*) embedding.col(id).data(), sizeof(tensor_type::Scalar) * rows);
    }
  }

  void read_embedding(std::istream& is, tensor_type& embedding, const bool bias_last)
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    
    size_type cols = 0;
    is.read((char*) &cols, sizeof(size_type));
    
    if (cols > embedding.cols())
      embedding.conservativeResize(embedding_ + bias_last, cols);
    
    for (size_type i = 0; i != cols; ++ i) {
      size_type word_size = 0;
      is.read((char*) &word_size, sizeof(size_type));
      
      buffer.resize(word_size);
      is.read((char*) &(*buffer.begin()), word_size);
      
      const word_type word(buffer.begin(), buffer.end());

      if (word.id() >= embedding.cols())
	embedding.conservativeResize(embedding_ + bias_last, word.id() + 1);
      
      is.read((char*) embedding.col(word.id()).data(), sizeof(tensor_type::Scalar) * (embedding_ + bias_last));
    }
  }

  void write(std::ostream& os, const tensor_type& matrix) const
  {
    const tensor_type::Index rows = matrix.rows();
    const tensor_type::Index cols = matrix.cols();

    os.write((char*) &rows, sizeof(tensor_type::Index));
    os.write((char*) &cols, sizeof(tensor_type::Index));
    
    os.write((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
  }

  void read(std::istream& is, tensor_type& matrix)
  {
    tensor_type::Index rows = 0;
    tensor_type::Index cols = 0;
    
    is.read((char*) &rows, sizeof(tensor_type::Index));
    is.read((char*) &cols, sizeof(tensor_type::Index));
    
    matrix.resize(rows, cols);
    
    is.read((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
  }

public:
  // dimension...
  size_type embedding_;
  size_type hidden_;
  size_type alignment_;

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

  tensor_type bi_;

  double scale_;
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

struct HMM
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

  typedef Average log_likelihood_type;

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

  typedef utils::compact_set<word_type,
			     utils::unassigned<word_type>, utils::unassigned<word_type>,
			     boost::hash<word_type>, std::equal_to<word_type>,
			     std::allocator<word_type> > word_set_type;
  
  HMM(const dictionary_type& dict,
      const size_type& samples,
      const size_type& beam)
    : dict_(dict), samples_(samples), log_samples_(std::log(double(samples))), beam_(beam) {}
  
  const dictionary_type& dict_;
  
  size_type samples_;
  double log_samples_;
  size_type beam_;
  
  tensor_type layer_trans_;
  
  tensor_type delta_beta_;
  tensor_type delta_alpha_;
  tensor_type delta_trans_;
  
  state_allocator_type state_allocator_;
  heap_set_type  heaps_;
  state_map_type states_;

  word_set_type sources_;
  word_set_type targets_;

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
    
#if 0
    std::cerr << "forward source: " << source << std::endl
	      << "forward target: " << target << std::endl;
#endif

    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type state_size = theta.embedding_ + theta.hidden_;

    const size_type offset_word   = 0;
    const size_type offset_matrix = theta.embedding_;
    
    state_allocator_.assign(state_type::size(state_size + theta.embedding_));
    
    alignment.clear();
    
    heaps_.clear();
    heaps_.resize(target_size + 2);
    
    state_type state_bos = state_allocator_.allocate();
    
    state_bos.prev() = state_type();
    state_bos.index() = 0;
    state_bos.score() = 0.0;
    
    matrix_type(state_bos.matrix(), theta.hidden_, 1) = theta.bi_;
    
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
	const word_type source_prev = (prev >= source_size + 2
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
	  
	  matrix_type layer_alpha(state_next.matrix(), state_size, 1);
	  matrix_type layer_trans(state_next.matrix() + state_size, theta.embedding_, 1);
	  
	  layer_alpha.block(offset_word, 0, theta.embedding_, 1) = theta.source_.col(source_next.id()) * theta.scale_;
	  
	  layer_alpha.block(offset_matrix, 0, theta.hidden_, 1)
	    = (theta.Wa_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size)
	       * matrix_type(state.matrix(), state_size, 1)
	       + theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)).array().unaryExpr(hinge());
	  
	  layer_trans = (theta.Wt_ * layer_alpha + theta.bt_).array().unaryExpr(hinge());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans * theta.scale_
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
	  
	  matrix_type layer_alpha(state_next.matrix(), state_size, 1);
	  matrix_type layer_trans(state_next.matrix() + state_size, theta.embedding_, 1);
	  
	  layer_alpha.block(offset_word, 0, theta.embedding_, 1) = theta.source_.col(vocab_type::EPSILON.id()) * theta.scale_;
	  
	  layer_alpha.block(offset_matrix, 0, theta.hidden_, 1)
	    = (theta.Wn_ * matrix_type(state.matrix(), state_size, 1) + theta.bn_).array().unaryExpr(hinge());
	  
	  layer_trans = (theta.Wt_ * layer_alpha + theta.bt_).array().unaryExpr(hinge());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans * theta.scale_
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
  log_likelihood_type backward(const sentence_type& source,
			       const sentence_type& target,
			       const model_type& theta,
			       gradient_type& gradient,
			       Gen& gen)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
    
#if 0
    std::cerr << "backward source: " << source << std::endl
	      << "backward target: " << target << std::endl;
#endif

    const size_type source_size = source.size();
    const size_type target_size = target.size();

    const size_type state_size = theta.embedding_ + theta.hidden_;

    const size_type offset_word   = 0;
    const size_type offset_matrix = theta.embedding_;

    sources_.clear();
    targets_.clear();
    
    sources_.insert(source.begin(), source.end());
    targets_.insert(target.begin(), target.end());
    
    states_.clear();
    states_.resize(target_size + 2);

    {
      state_set_type& states = states_[target_size + 1];
      const heap_type& heap = heaps_[target_size + 1];
      
      heap_type::const_iterator hiter_end = heap.end();
      
      double logsum = - std::numeric_limits<double>::infinity();
      for (heap_type::const_iterator hiter = std::max(heap.begin(), hiter_end - beam_); hiter != hiter_end; ++ hiter)
	logsum = utils::mathop::logsum(logsum, double(hiter->score()));

      for (heap_type::const_iterator hiter = std::max(heap.begin(), hiter_end - beam_); hiter != hiter_end; ++ hiter) {
	state_type buffer = state_allocator_.allocate();
	
	buffer.score() = std::exp(double(hiter->score()) - logsum);
	matrix_type(buffer.matrix(), state_size, 1).setZero();
	
	states[*hiter] = buffer;
      }
    }

    boost::random::uniform_int_distribution<> uniform_source(0, source_size - 1);

    log_likelihood_type log_likelihood;

    ++ gradient.count_;
    
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
	
	state_set_type::iterator piter = states_prev.find(state_prev);
	if (piter == states_prev.end()) {
	  state_type buffer = state_allocator_.allocate();
	  
	  buffer.score() = 0.0;
	  matrix_type(buffer.matrix(), state_size, 1).setZero();
	  
	  piter = states_prev.insert(std::make_pair(state_prev, buffer)).first;
	}
	
	const matrix_type beta_next(siter->second.matrix(), state_size, 1);
	/**/  matrix_type beta_prev(piter->second.matrix(), state_size, 1);
	
	const size_type next = state_next.index();
	const size_type prev = state_prev.index();

	const double weight = siter->second.score();
	piter->second.score() += siter->second.score();

	//std::cerr << "target: " << trg << " weight: " << weight << std::endl;
	
	const word_type source_next(next >= source_size + 2
				    ? vocab_type::EPSILON
				    : (next == 0
				       ? vocab_type::BOS
				       : (next == source_size + 1
					  ? vocab_type::EOS
					  : source[next - 1])));
	const word_type source_prev = (prev >= source_size + 2
				       ? vocab_type::EPSILON
				       : (prev == 0
					  ? vocab_type::BOS
					  : (prev == source_size + 1
					     ? vocab_type::EOS
					     : source[prev - 1])));
	
	const matrix_type alpha_next(state_next.matrix(), state_size, 1);
	const matrix_type trans_next(state_next.matrix() + state_size, theta.embedding_, 1);
	const matrix_type alpha_prev(state_prev.matrix(), state_size, 1);
	const matrix_type trans_prev(state_prev.matrix() + state_size, theta.embedding_, 1);

	word_type target_sampled = target_next;
	
	if (dict_.size(source_next) > 1) {
	  if (source_next == vocab_type::EPSILON) {
	    target_sampled = dict_.draw(source[uniform_source(gen)], gen);
	    
	    while (targets_.find(target_sampled) != targets_.end())
	      target_sampled = dict_.draw(source[uniform_source(gen)], gen);
	  } else {
	    target_sampled = dict_.draw(source_next, gen);
	    
	    while (target_sampled == target_next)
	      target_sampled = dict_.draw(source_next, gen);
	  } 
	}
	
	const double score_c = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * trans_next * theta.scale_
				+ theta.target_.col(target_next.id()).block(theta.embedding_, 0, 1, 1))(0, 0);
	const double score_m = (theta.target_.col(target_sampled.id()).block(0, 0, theta.embedding_, 1).transpose() * trans_next * theta.scale_
				+ theta.target_.col(target_sampled.id()).block(theta.embedding_, 0, 1, 1))(0, 0);

	const double score_noise_c = 0.0 + dict_.logprob(source_next, target_next);
	const double score_noise_m = 0.0 + dict_.logprob(source_next, target_sampled);
	
	const double z_c = utils::mathop::logsum(score_c, score_noise_c);
	const double z_m = utils::mathop::logsum(score_m, score_noise_m);
	
	const double logprob_c = score_c - z_c;
	const double logprob_m = score_m - z_m;
	
	const double logprob_noise_c = score_noise_c - z_c;
	const double logprob_noise_m = score_noise_m - z_m;
	
	log_likelihood += (logprob_c + logprob_noise_m) * weight;
	
	const double loss_c = weight * (- 1.0 + std::exp(logprob_c));
	const double loss_m = weight * (        std::exp(logprob_m));
	
	tensor_type& dembedding_c = gradient.target(target_next);
	tensor_type& dembedding_m = gradient.target(target_sampled);
	
	dembedding_c.block(0, 0, theta.embedding_, 1).array() += loss_c * trans_next.array();
	dembedding_c.block(theta.embedding_, 0, 1, 1).array() += loss_c;
	dembedding_m.block(0, 0, theta.embedding_, 1).array() += loss_m * trans_next.array();
	dembedding_m.block(theta.embedding_, 0, 1, 1).array() += loss_m;
	
	delta_trans_ = (trans_next.array().unaryExpr(dhinge())
			* (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1) * loss_c * theta.scale_
			   + theta.target_.col(target_sampled.id()).block(0, 0, theta.embedding_, 1) * loss_m * theta.scale_).array());
	
	gradient.Wt_ += delta_trans_ * alpha_next.transpose();
	gradient.bt_ += delta_trans_;

	delta_alpha_ = theta.Wt_.transpose() * delta_trans_;

	if (! delta_beta_.rows())
	  delta_beta_.resize(state_size, 1);
	
	delta_beta_.block(offset_word, 0, theta.embedding_, 1)
	  = (delta_alpha_.block(offset_word, 0, theta.embedding_, 1)
	     + beta_next.block(offset_word, 0, theta.embedding_, 1));
	
	delta_beta_.block(offset_matrix, 0, theta.hidden_, 1)
	  = (alpha_next.block(offset_matrix, 0, theta.hidden_, 1).array().unaryExpr(dhinge())
	      * (delta_alpha_.block(offset_matrix, 0, theta.hidden_, 1)
		 + beta_next.block(offset_matrix, 0, theta.hidden_, 1)).array());

	if (next >= source_size + 2) {
	  gradient.source(vocab_type::EPSILON.id()) += delta_beta_.block(offset_word, 0, theta.embedding_, 1);
	  
	  gradient.Wn_ += delta_beta_.block(offset_matrix, 0, theta.hidden_, 1) * alpha_prev.transpose();
	  gradient.bn_ += delta_beta_.block(offset_matrix, 0, theta.hidden_, 1);
	  
	  beta_prev += theta.Wn_.transpose() * delta_beta_.block(offset_matrix, 0, theta.hidden_, 1);
	} else {
	  const size_type shift = theta.shift(source_size, target_size, prev, next);
	  
	  gradient.source(source_next.id()) += delta_beta_.block(offset_word, 0, theta.embedding_, 1);
	  
	  gradient.Wa_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size)
	    += delta_beta_.block(offset_matrix, 0, theta.hidden_, 1) * alpha_prev.transpose();
	  gradient.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
	    += delta_beta_.block(offset_matrix, 0, theta.hidden_, 1);
	  
	  beta_prev += (theta.Wa_.block(theta.hidden_ * shift, 0, theta.hidden_, state_size).transpose()
			* delta_beta_.block(offset_matrix, 0, theta.hidden_, 1));
	}
      }
    }
    
    // final propagation...
    if (states_[0].empty())
      throw std::runtime_error("no initial state?");
    if (states_[0].size() != 1)
      throw std::runtime_error("multiple initial states?");
    
    matrix_type beta_bos(states_[0].begin()->second.matrix(), state_size, 1);
    
    gradient.source(vocab_type::BOS.id()) += beta_bos.block(offset_word, 0, theta.embedding_, 1);
    gradient.bi_ += beta_bos.block(offset_matrix, 0, theta.hidden_, 1);
    
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
    
    matrix_type(state_bos.matrix(), theta.hidden_, 1) = theta.bi_;
    
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
	const word_type source_prev = (prev >= source_size + 2
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
		* theta.source_.col(source_prev.id()) * theta.scale_)
	       
	       + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
		  * matrix_type(state.matrix(), theta.hidden_, 1))
	       
	       + theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)).array().unaryExpr(hinge());
	  
	  layer_trans_ = ((theta.Wt_.block(0, offset_word, theta.embedding_, theta.embedding_)
			   * theta.source_.col(source_next.id()) * theta.scale_)
			  
			  + (theta.Wt_.block(0, offset_matrix, theta.embedding_, theta.hidden_)
			     * layer_alpha)
			  
			  + theta.bt_).array().unaryExpr(hinge());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_ * theta.scale_
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
		* theta.source_.col(source_prev.id()) * theta.scale_)
	       
	       + (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
		  * matrix_type(state.matrix(), theta.hidden_, 1))
	       
	       + theta.bn_.block(0, 0, theta.hidden_, 1)).array().unaryExpr(hinge());
	  
	  layer_trans_ = ((theta.Wt_.block(0, offset_word, theta.embedding_, theta.embedding_)
			   * theta.source_.col(vocab_type::EPSILON.id()) * theta.scale_)
			  
			  + (theta.Wt_.block(0, offset_matrix, theta.embedding_, theta.hidden_)
			     * layer_alpha)
			  
			  + theta.bt_).array().unaryExpr(hinge());
	  
	  const double score = (theta.target_.col(target_next.id()).block(0, 0, theta.embedding_, 1).transpose() * layer_trans_ * theta.scale_
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

    bi_ = tensor_type::Zero(hidden_, 1);
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
	     scale,
	     false);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     const_cast<tensor_type&>(target_),
	     titer->second,
	     embedding.source_.col(titer->first.id()),
	     scale,
	     true);
    
    update(theta.Wt_, const_cast<tensor_type&>(Wt_), gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, const_cast<tensor_type&>(bt_), gradient.bt_, scale, false);
    
    update(theta.Wa_, const_cast<tensor_type&>(Wa_), gradient.Wa_, scale, lambda_ != 0.0);
    update(theta.ba_, const_cast<tensor_type&>(ba_), gradient.ba_, scale, false);

    update(theta.Wn_, const_cast<tensor_type&>(Wn_), gradient.Wn_, scale, lambda_ != 0.0);
    update(theta.bn_, const_cast<tensor_type&>(bn_), gradient.bn_, scale, false);
    
    update(theta.bi_, const_cast<tensor_type&>(bi_), gradient.bi_, scale, lambda_ != 0.0);
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
	      const double scale,
	      const bool bias_last=false) const
  {
    if (lambda2_ > 0.0) {
      for (int row = 0; row != g.rows() - bias_last; ++ row) 
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
      for (int row = 0; row != g.rows() - bias_last; ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) +=  g(row, 0) * g(row, 0) * scale * scale;
	  
	  const double rate = eta0_ / std::sqrt(double(1.0) + G(row, word.id()));
	  const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	  
	  theta(row, word.id()) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
    }
    
    if (bias_last) {
      const int row = g.rows() - 1;
      
      if (g(row, 0) != 0) {
	G(row, word.id()) += g(row, 0) * g(row, 0) * scale * scale;
	theta(row, word.id()) -= eta0_ * scale * g(row, 0) / std::sqrt(double(1.0) + G(row, word.id()));
      }
    }
  }

  template <typename Theta, typename GradVar, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool bias_last=false) const
  {
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
  }
  
  size_type embedding_;
  size_type hidden_;
  size_type alignment_;
  
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

  tensor_type bi_;
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
	     theta.scale_,
	     false);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     titer->second,
	     embedding.source_.col(titer->first.id()),
	     scale,
	     theta.scale_,
	     true);
    
    update(theta.Wt_, gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, gradient.bt_, scale, false);
    
    update(theta.Wa_, gradient.Wa_, scale, lambda_ != 0.0);
    update(theta.ba_, gradient.ba_, scale, false);

    update(theta.Wn_, gradient.Wn_, scale, lambda_ != 0.0);
    update(theta.bn_, gradient.bn_, scale, false);
    
    update(theta.bi_, gradient.bi_, scale, lambda_ != 0.0);
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
	      const double theta_scale,
	      const bool bias_last=false) const
  {
    const double eta = eta0_ / (epoch_ + 1);

    if (lambda2_ != 0.0) {
      const size_type rows = g.rows();
      
      if (bias_last) {
	theta.col(word.id()).block(0, 0, rows - 1, 1) -= eta * lambda2_ * (theta.col(word.id()).block(0, 0, rows - 1, 1)
									   - c.block(0, 0, rows - 1, 1) / theta_scale);
	theta.col(word.id()).block(0, 0, rows - 1, 1) -= (eta * scale / theta_scale) * g.block(0, 0, rows - 1, 1);
	theta.col(word.id()).block(rows - 1, 0, 1, 1) -= eta * scale * g.block(rows - 1, 0, 1, 1);
      } else {
	theta.col(word.id()) -= eta * lambda2_ * (theta.col(word.id()) - c.block(0, 0, rows, 1) / theta_scale);
	theta.col(word.id()) -= (eta * scale / theta_scale) * g;
      }
    } else {
      if (bias_last) {
	const size_type rows = g.rows();
	
	theta.col(word.id()).block(0, 0, rows - 1, 1) -= (eta * scale / theta_scale) * g.block(0, 0, rows - 1, 1);
	theta.col(word.id()).block(rows - 1, 0, 1, 1) -= eta * scale * g.block(rows - 1, 0, 1, 1);
      } else
	theta.col(word.id()) -= (eta * scale / theta_scale) * g;
    }
  }

  double lambda_;
  double lambda2_;
  double eta0_;
  
  size_type epoch_;
};

#endif
