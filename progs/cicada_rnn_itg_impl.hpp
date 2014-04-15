//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_RNN_ITG_IMPL__HPP__
#define __CICADA_RNN_ITG_IMPL__HPP__ 1

// rnn-ITG
//
// bi-phrase representation constructed by ITG
// then, classify the contexts out of the bi-phrases, inspired by the skip-gram
// for computing embeddings (Mikolov et al., 2013).
// we need input/output source/target embedding: source/target, output_source/target
// Wt/bt for terminal emission from hidden layers
// Ws/bs and Wi/bi for straight/inversion of ITG
// we perform negative sampling of Mikolov et al., (2013)

#include <cstdlib>
#include <cmath>
#include <climits>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"
#include "cicada/alignment.hpp"
#include "cicada/bitext.hpp"
#include "cicada/semiring/logprob.hpp"

#include "utils/alloc_vector.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bichart.hpp"
#include "utils/bithack.hpp"
#include "utils/chunk_vector.hpp"
#include "utils/compact_map.hpp"
#include "utils/compact_set.hpp"
#include "utils/mathop.hpp"
#include "utils/unordered_map.hpp"
#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/vector2.hpp"

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
  
  Gradient() : embedding_(0), hidden_(0), span_(0), count_(0), shared_(0) {}
  Gradient(const size_type& embedding,
	   const size_type& hidden,
	   const size_type& span) 
    : embedding_(embedding),
      hidden_(hidden),
      span_(span),
      count_(0),
      shared_(0)
  { initialize(embedding, hidden, span); }

  Gradient& operator-=(const Gradient& x)
  {
    decrement(source_, x.source_);
    decrement(target_, x.target_);

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    if (! Ws_.rows())
      Ws_ = tensor_type::Zero(x.Ws_.rows(), x.Ws_.cols());
    if (! bs_.rows())
      bs_ = tensor_type::Zero(x.bs_.rows(), x.bs_.cols());

    if (! Wi_.rows())
      Wi_ = tensor_type::Zero(x.Wi_.rows(), x.Wi_.cols());
    if (! bi_.rows())
      bi_ = tensor_type::Zero(x.bi_.rows(), x.bi_.cols());

    if (! Wt_.rows())
      Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    if (! bt_.rows())
      bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());

    Wc_ -= x.Wc_;
    bc_ -= x.bc_;
    
    Ws_ -= x.Ws_;
    bs_ -= x.bs_;
    
    Wi_ -= x.Wi_;
    bi_ -= x.bi_;
    
    Wt_ -= x.Wt_;
    bt_ -= x.bt_;    
    
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

    if (! Ws_.rows())
      Ws_ = tensor_type::Zero(x.Ws_.rows(), x.Ws_.cols());
    if (! bs_.rows())
      bs_ = tensor_type::Zero(x.bs_.rows(), x.bs_.cols());

    if (! Wi_.rows())
      Wi_ = tensor_type::Zero(x.Wi_.rows(), x.Wi_.cols());
    if (! bi_.rows())
      bi_ = tensor_type::Zero(x.bi_.rows(), x.bi_.cols());

    if (! Wt_.rows())
      Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    if (! bt_.rows())
      bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());

    Wc_ += x.Wc_;
    bc_ += x.bc_;

    Ws_ += x.Ws_;
    bs_ += x.bs_;
    
    Wi_ += x.Wi_;
    bi_ += x.bi_;
    
    Wt_ += x.Wt_;
    bt_ += x.bt_;    
    
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

    Ws_.setZero();
    bs_.setZero();
    
    Wi_.setZero();
    bi_.setZero();

    Wt_.setZero();
    bt_.setZero();
    
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
  
  void initialize(const size_type embedding, const size_type hidden, const size_type span)
  {
    if (hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    
    embedding_ = embedding;
    hidden_    = hidden;
    span_      = span;
    
    clear();

    Wc_ = tensor_type::Zero(1, hidden_);
    bc_ = tensor_type::Zero(1, 1);
    
    Ws_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    bs_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    Wi_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    bi_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    Wt_ = tensor_type::Zero(hidden_, embedding_ * 2);
    bt_ = tensor_type::Zero(hidden_, 1);
    
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
	matrix -= siter->second;
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
    os.write((char*) &span_,      sizeof(size_type));
    os.write((char*) &count_,     sizeof(size_type));

    write(os, Wc_);
    write(os, bc_);

    write(os, Ws_);
    write(os, bs_);

    write(os, Wi_);
    write(os, bi_);

    write(os, Wt_);
    write(os, bt_);
    
    write(os, source_, embedding_);
    write(os, target_, embedding_);
  }
  
  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &span_,      sizeof(size_type));
    is.read((char*) &count_,     sizeof(size_type));

    read(is, Wc_);
    read(is, bc_);
    
    read(is, Ws_);
    read(is, bs_);
    
    read(is, Wi_);
    read(is, bi_);
    
    read(is, Wt_);
    read(is, bt_);
    
    read(is, source_, embedding_);
    read(is, target_, embedding_);
  }

  void write(std::ostream& os, const embedding_type& embedding, const size_type& dimension) const
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

  void read(std::istream& is, embedding_type& embedding, const size_type& dimension)
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
      
      matrix.resize(dimension, 1);
      
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
  size_type span_;
  
  // embedding
  embedding_type source_;
  embedding_type target_;

  // classification
  tensor_type Wc_;
  tensor_type bc_;

  // straight
  tensor_type Ws_;
  tensor_type bs_;

  // inversion
  tensor_type Wi_;
  tensor_type bi_;  
  
  // terminal
  tensor_type Wt_;
  tensor_type bt_;
  
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
  
  Model() : embedding_(0), hidden_(0), span_(0), scale_(1) {}  
  template <typename Words, typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	const size_type& span,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      span_(span),
      scale_(1)
  { initialize(embedding, hidden, span, words_source, words_target, gen); }

  void clear()
  {
    // embedding
    source_.setZero();
    target_.setZero();

    Wc_.setZero();
    bc_.setZero();

    Ws_.setZero();
    bs_.setZero();
    
    Wi_.setZero();
    bi_.setZero();

    Wt_.setZero();
    bt_.setZero();
    
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
		  const size_type span,
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
    span_      = span;
    
    clear();
    
    const size_type vocabulary_size = word_type::allocated();
    
    const double range_e = std::sqrt(6.0 / (embedding_ + 1));
    const double range_c = std::sqrt(6.0 / (hidden_ + 1));
    const double range_s = std::sqrt(6.0 / (hidden_ + hidden_ + hidden_));
    const double range_i = std::sqrt(6.0 / (hidden_ + hidden_ + hidden_));
    const double range_t = std::sqrt(6.0 / (hidden + embedding_ * 2));

    source_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));

    Wc_ = tensor_type::Zero(1, hidden_).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Zero(1, 1);
    
    Ws_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_).array().unaryExpr(randomize<Gen>(gen, range_s));
    bs_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    Wi_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_).array().unaryExpr(randomize<Gen>(gen, range_i));
    bi_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    Wt_ = tensor_type::Zero(hidden_, embedding_ * 2).array().unaryExpr(randomize<Gen>(gen, range_t));
    bt_ = tensor_type::Zero(hidden_, 1);
    
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
    os.write((char*) &span_,      sizeof(size_type));
    os.write((char*) &scale_,     sizeof(double));

    write(os, Wc_);
    write(os, bc_);
    
    write(os, Ws_);
    write(os, bs_);

    write(os, Wi_);
    write(os, bi_);

    write(os, Wt_);
    write(os, bt_);
    
    write_embedding(os, source_, embedding_);
    write_embedding(os, target_, embedding_);
  }

  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &span_,      sizeof(size_type));
    is.read((char*) &scale_,     sizeof(double));

    read(is, Wc_);
    read(is, bc_);
    
    read(is, Ws_);
    read(is, bs_);

    read(is, Wi_);
    read(is, bi_);

    read(is, Wt_);
    read(is, bt_);
    
    read_embedding(is, source_, embedding_);
    read_embedding(is, target_, embedding_);
    
    // checking...
  }

  void write_embedding(std::ostream& os, const tensor_type& embedding, const size_type& dimension) const
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

  void read_embedding(std::istream& is, tensor_type& embedding, const size_type& dimension)
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    
    size_type cols = 0;
    is.read((char*) &cols, sizeof(size_type));
    
    if (cols > embedding.cols())
      embedding.conservativeResize(dimension, cols);
    
    for (size_type i = 0; i != cols; ++ i) {
      size_type word_size = 0;
      is.read((char*) &word_size, sizeof(size_type));
      
      buffer.resize(word_size);
      is.read((char*) &(*buffer.begin()), word_size);
      
      const word_type word(buffer.begin(), buffer.end());

      if (word.id() >= embedding.cols())
	embedding.conservativeResize(dimension, word.id() + 1);
      
      is.read((char*) embedding.col(word.id()).data(), sizeof(tensor_type::Scalar) * (dimension));
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
  
  struct real_policy : boost::spirit::karma::real_policies<parameter_type>
  {
    static unsigned int precision(parameter_type)
    {
      return 10;
    }
  };

public:
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
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);
    rep["span"]      = utils::lexical_cast<std::string>(span_);
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), rep.path("vocab-source"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), rep.path("vocab-target"), target_, words_target_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);

    write(rep.path("Ws.txt.gz"), rep.path("Ws.bin"), Ws_);
    write(rep.path("bs.txt.gz"), rep.path("bs.bin"), bs_);
    
    write(rep.path("Wi.txt.gz"), rep.path("Wi.bin"), Wi_);
    write(rep.path("bi.txt.gz"), rep.path("bi.bin"), bi_);
    
    write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);
  }
  
private:
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
  // dimension...
  size_type embedding_;
  size_type hidden_;
  size_type span_;

  word_unique_type words_source_;
  word_unique_type words_target_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;

  tensor_type Wc_;
  tensor_type bc_;

  tensor_type Ws_;
  tensor_type bs_;

  tensor_type Wi_;
  tensor_type bi_;
  
  tensor_type Wt_;
  tensor_type bt_;

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

  typedef std::vector<parameter_type, std::allocator<parameter_type> > buffer_type;
  
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
    
    State() : span_(),
	      left_(0),
	      right_(0),
	      score_(- std::numeric_limits<double>::infinity()),
	      error_(0),
	      loss_(0) {}
    State(const span_pair_type& span)
      : span_(span),
	left_(0),
	right_(0),
	score_(- std::numeric_limits<double>::infinity()),
	error_(0),
	loss_(0) {}
    State(const span_pair_type& span, const state_type* left, const state_type* right)
      : span_(span),
	left_(left),
	right_(right),
	score_(- std::numeric_limits<double>::infinity()),
	error_(0),
	loss_(0) {}
    
    bool terminal() const { return left_ == 0 && right_ == 0; }
    bool straight() const { return ! terminal() && left_->span_.target_.last_  == right_->span_.target_.first_; }
    bool inverted() const { return ! terminal() && left_->span_.target_.first_ == right_->span_.target_.last_; }
    
    span_pair_type span_;
    const state_type* left_;
    const state_type* right_;
    
    // score and # of errors
    double          score_;
    difference_type error_;
    double          loss_;
    
    // hidden layers
    tensor_type layer_;
    tensor_type delta_;
    
    // sampled word
    word_type source_;
    word_type target_;
  };
  
  typedef State state_type;
  
  typedef utils::chunk_vector<state_type, 4 * 1024 * 1024 / sizeof(state_type), std::allocator<state_type> > state_pool_type;
  
  typedef std::vector<const state_type*, std::allocator<const state_type*> > state_set_type;
  typedef utils::bichart<state_set_type, std::allocator<state_set_type> > chart_type;
  
  typedef std::pair<double, const state_type*> score_state_type;
  typedef std::vector<score_state_type, std::allocator<score_state_type> > heap_type;  
  
  struct heap_compare
  {
    // sort by less item so that we can pop from greater items
    bool operator()(const score_state_type& x, const score_state_type& y) const
    {
      return x.first < y.first;
    }
  };
  
  typedef std::vector<state_set_type, std::allocator<state_set_type> > agenda_type;

  typedef std::vector<char, std::allocator<char> > error_set_type;
  
  typedef std::vector<hyperedge_type, std::allocator<hyperedge_type> > derivation_type;
  typedef std::vector<const state_type*, std::allocator<const state_type*> > stack_type;

  typedef std::pair<const state_type*, const state_type*> state_pair_type;
  
  typedef utils::compact_set<state_pair_type,
			     utils::unassigned<state_pair_type>, utils::unassigned<state_pair_type>,
			     utils::hashmurmur3<size_t>, std::equal_to<state_pair_type>,
			     std::allocator<state_pair_type> > state_pair_unique_type;

  typedef utils::compact_set<const state_type*,
			     utils::unassigned<const state_type*>, utils::unassigned<const state_type*>,
			     utils::hashmurmur3<size_t>, std::equal_to<const state_type*>,
			     std::allocator<const state_type*> > state_unique_type;
  typedef std::vector<state_unique_type, std::allocator<state_unique_type> > state_unique_set_type;

  struct RestCost
  {
    double score_;
    double alpha_;
    double beta_;
    
    RestCost()
      : score_(- std::numeric_limits<double>::infinity()),
	alpha_(- std::numeric_limits<double>::infinity()),
	beta_(- std::numeric_limits<double>::infinity()) {}
  };
  
  typedef RestCost rest_cost_type;
  typedef std::vector<rest_cost_type, std::allocator<rest_cost_type> > rest_cost_set_type;
  

  typedef utils::compact_set<word_type,
			     utils::unassigned<word_type>, utils::unassigned<word_type>,
			     boost::hash<word_type>, std::equal_to<word_type>,
			     std::allocator<word_type> > word_set_type;
  
  
  ITG(const dictionary_type& dict_source_target,
      const dictionary_type& dict_target_source,
      const size_type& window,
      const size_type& samples,
      const size_type& beam)
    : dict_source_target_(dict_source_target),
      dict_target_source_(dict_target_source),
      window_(window),
      samples_(samples),
      beam_(beam) 
  {
  }
  
  const dictionary_type& dict_source_target_;
  const dictionary_type& dict_target_source_;

  size_type window_;
  size_type samples_;
  size_type beam_;
  
  chart_type  chart_;
  chart_type  chart_sampled_;
  agenda_type agenda_;
  heap_type   heap_;
  
  error_set_type errors_;
  
  rest_cost_set_type costs_source_;
  rest_cost_set_type costs_target_;
  
  state_pair_unique_type uniques_;
  stack_type stack_;

  state_unique_set_type backwards_;

  state_pool_type   states_;
  const state_type* linked_;

  word_set_type sources_;
  word_set_type targets_;
  word_set_type sampled_;

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
    chart_sampled_.clear();
    chart_.clear();
    agenda_.clear();
    heap_.clear();

    costs_source_.clear();
    costs_target_.clear();
    
    uniques_.clear();
    stack_.clear();

    backwards_.clear();

    states_.clear();
    linked_ = 0;
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
    
    chart_.clear();
    chart_.reserve(source_size + 1, target_size + 1);
    chart_.resize(source_size + 1, target_size + 1);
    
    agenda_.clear();
    agenda_.reserve(source_size + target_size + 1);
    agenda_.resize(source_size + target_size + 1);

    errors_.clear();
    errors_.resize(source_size + target_size + 1, false);
    
    costs_source_.clear();
    costs_target_.clear();

    costs_source_.reserve(source_size + 1);
    costs_target_.reserve(target_size + 1);

    costs_source_.resize(source_size + 1);
    costs_target_.resize(target_size + 1);
    
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
    
    // estimate rest-costs
    forward_backward(costs_source_);
    forward_backward(costs_target_);
    
#if 0
    std::cerr << "biparsing source: " << source << std::endl
	      << "biparsing target: " << target << std::endl;
#endif

    // forward actual forward path
    const size_type length_max = source_size + target_size;
    
    for (size_type length = 1; length != length_max; ++ length) 
      if (! agenda_[length].empty()) {
	
	//std::cerr << "length: " << length << std::endl;

	state_set_type& states = agenda_[length];
	
	heap_.clear();
	heap_.reserve(states.size());

	state_set_type::const_iterator siter_end = states.end();
	for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
	  const double score = ((*siter)->score_
				+ std::min(costs_source_[(*siter)->span_.source_.first_].alpha_
					   + costs_source_[(*siter)->span_.source_.last_].beta_,
					   costs_target_[(*siter)->span_.target_.first_].alpha_
					   + costs_target_[(*siter)->span_.target_.last_].beta_));
	  
	  heap_.push_back(score_state_type(score, *siter));
	  std::push_heap(heap_.begin(), heap_.end(), heap_compare());
	}
	
	heap_type::iterator hiter_begin = heap_.begin();
	heap_type::iterator hiter       = heap_.end();
	heap_type::iterator hiter_end   = heap_.end();
	
	if (length > 1 && std::distance(hiter_begin, hiter_end) > beam_) {
	  for (/**/; hiter_begin != hiter && std::distance(hiter, hiter_end) != beam_; -- hiter)
	    std::pop_heap(hiter_begin, hiter, heap_compare());
	  
	  for (heap_type::iterator iter = hiter_begin ; iter != hiter; ++ iter)
	    if (! iter->second->terminal())
	      deallocate(iter->second);
	} else
	  hiter = hiter_begin;
	
	// clear uniques
	uniques_.clear();

	// first, enumerate and put into chart...
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const state_type* state = iter->second;
	  const span_pair_type& span = state->span_;
	  
	  if (! state->terminal())
	    chart_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(state);
	}
	
	// then, enumerate new hypotheses
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  //std::cerr << "length: " << length << " loss: " << iter->first << std::endl;
	  
	  const state_type* state = iter->second;
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
	      
	      const state_set_type& states = chart_(S, s, U, u);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(*siter, state)).second)
		  forward(source, target, span_pair_type(S, t, U, v), *siter, state, theta, true);
	    }
	    
	    // inversion
	    for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: StuU
	      // span1: SsvU
	      // span2: stuv
	      
	      const state_set_type& states = chart_(S, s, v, U);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
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
	      
	      const state_set_type& states = chart_(t, S, U, u);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(*siter, state)).second)
		  forward(source, target, span_pair_type(s, S, U, v), state, *siter, theta, false);
	    }
	    
	    // straight
	    for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: sSuU
	      // span1: stuv
	      // span2: tSvU

	      const state_set_type& states = chart_(t, S, v, U);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(*siter, state)).second)
		  forward(source, target, span_pair_type(s, S, u, U), state, *siter, theta, true);
	    }
	  }
	}
      }
    
    // final enumeration...
    state_set_type& states = agenda_[length_max];
    
    if (states.empty())
      return - std::numeric_limits<double>::infinity();
    
    heap_.clear();
    heap_.reserve(states.size());
    
    state_set_type::const_iterator siter_end = states.end();
    for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
      heap_.push_back(score_state_type((*siter)->score_, *siter));
      std::push_heap(heap_.begin(), heap_.end(), heap_compare());
    }
    
    heap_type::iterator hiter_begin = heap_.begin();
    heap_type::iterator hiter       = heap_.end();
    heap_type::iterator hiter_end   = heap_.end();
    
    // perform sorting...
    for (/**/; hiter_begin != hiter && std::distance(hiter, hiter_end) != beam_; -- hiter)
      std::pop_heap(hiter_begin, hiter, heap_compare());
    
    for (heap_type::iterator iter = hiter_begin ; iter != hiter; ++ iter)
      if (! iter->second->terminal())
	deallocate(iter->second);
    
    // allocate new states...
    chart_(0, source_size, 0, target_size).clear();
    for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter)
      chart_(0, source_size, 0, target_size).push_back(iter->second);
    
    return chart_(0, source_size, 0, target_size).back()->score_;
  }

  template <typename Gen>
  double forward(const sentence_type& source,
		 const sentence_type& target,
		 const model_type& theta,
		 Gen& gen)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
#if 0
    std::cerr << "sampled forward source: " << source << std::endl
	      << "sampled forward target: " << target << std::endl;
#endif
    
    chart_sampled_.clear();
    chart_sampled_.reserve(source_size + 1, target_size + 1);
    chart_sampled_.resize(source_size + 1, target_size + 1);
    
    agenda_.clear();
    agenda_.reserve(source_size + target_size + 1);
    agenda_.resize(source_size + target_size + 1);

    errors_.clear();
    errors_.resize(source_size + target_size + 1, false);
    
    costs_source_.clear();
    costs_target_.clear();
    
    costs_source_.reserve(source_size + 1);
    costs_target_.reserve(target_size + 1);

    costs_source_.resize(source_size + 1);
    costs_target_.resize(target_size + 1);

    sources_.clear();
    targets_.clear();
    
    sources_.insert(source.begin(), source.end());
    targets_.insert(target.begin(), target.end());
    
    // initialize chart...
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = 0; trg <= target_size; ++ trg)
	if (src < source_size || trg < target_size) {
	  // epsilon at target
	  if (src < source_size)
	    forward(source, target, sources_, targets_, span_pair_type(span_type(src, src + 1), span_type(trg, trg)), theta, gen);
	  
	  // epsilon at source
	  if (trg < target_size)
	    forward(source, target, sources_, targets_, span_pair_type(span_type(src, src), span_type(trg, trg + 1)), theta, gen);
	  
	  // word-pair
	  if (src < source_size && trg < target_size)
	    forward(source, target, sources_, targets_, span_pair_type(span_type(src, src + 1), span_type(trg, trg + 1)), theta, gen);
	}
    
    // estimate rest-costs
    forward_backward(costs_source_);
    forward_backward(costs_target_);
    
#if 0
    std::cerr << "sampled biparsing source: " << source << std::endl
	      << "sampled biparsing target: " << target << std::endl;
#endif

    // forward actual forward path
    const size_type length_max = source_size + target_size;
    
    for (size_type length = 1; length != length_max; ++ length) 
      if (! agenda_[length].empty()) {
	
	//std::cerr << "length: " << length << std::endl;

	state_set_type& states = agenda_[length];
	
	heap_.clear();
	heap_.reserve(states.size());
	
	// here, we use loss-biased search...
	state_set_type::const_iterator siter_end = states.end();
	for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
	  const double score = ((*siter)->score_ + (*siter)->error_
				+ std::min(costs_source_[(*siter)->span_.source_.first_].alpha_
					   + costs_source_[(*siter)->span_.source_.last_].beta_,
					   costs_target_[(*siter)->span_.target_.first_].alpha_
					   + costs_target_[(*siter)->span_.target_.last_].beta_));
	  
	  heap_.push_back(score_state_type(score, *siter));
	  std::push_heap(heap_.begin(), heap_.end(), heap_compare());
	}
	
	heap_type::iterator hiter_begin = heap_.begin();
	heap_type::iterator hiter       = heap_.end();
	heap_type::iterator hiter_end   = heap_.end();
	
	if (length > 1 && std::distance(hiter_begin, hiter_end) > beam_) {
	  bool has_error = false;
	  
	  for (/**/; hiter_begin != hiter && std::distance(hiter, hiter_end) != beam_; -- hiter) {
	    std::pop_heap(hiter_begin, hiter, heap_compare());
	    
	    has_error |= ((hiter - 1)->second->error_ > 0);
	  }
	  
	  // search until we find an error...
	  if (! has_error && errors_[length])
	    for (/**/; hiter_begin != hiter && ! has_error; -- hiter) {
	      std::pop_heap(hiter_begin, hiter, heap_compare());
	      
	      has_error |= ((hiter - 1)->second->error_ > 0);
	    }
	  
	  for (heap_type::iterator iter = hiter_begin ; iter != hiter; ++ iter)
	    if (! iter->second->terminal())
	      deallocate(iter->second);
	} else
	  hiter = hiter_begin;
	
	// clear uniques
	uniques_.clear();

	// first, enumerate and put into chart...
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const state_type* state = iter->second;
	  const span_pair_type& span = state->span_;
	  
	  if (! state->terminal())
	    chart_sampled_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(state);
	}
	
	// then, enumerate new hypotheses
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  //std::cerr << "length: " << length << " loss: " << iter->first << std::endl;
	  
	  const state_type* state = iter->second;
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
	      
	      const state_set_type& states = chart_sampled_(S, s, U, u);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(*siter, state)).second)
		  forward(source, target, span_pair_type(S, t, U, v), *siter, state, theta, true);
	    }
	    
	    // inversion
	    for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: StuU
	      // span1: SsvU
	      // span2: stuv
	      
	      const state_set_type& states = chart_sampled_(S, s, v, U);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
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
	      
	      const state_set_type& states = chart_sampled_(t, S, U, u);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(state, *siter)).second)
		  forward(source, target, span_pair_type(s, S, U, v), state, *siter, theta, false);
	    }
	    
	    // straight
	    for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: sSuU
	      // span1: stuv
	      // span2: tSvU

	      const state_set_type& states = chart_sampled_(t, S, v, U);
	      
	      state_set_type::const_iterator siter_end = states.end();
	      for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter)
		if (uniques_.insert(std::make_pair(state, *siter)).second)
		  forward(source, target, span_pair_type(s, S, u, U), state, *siter, theta, true);
	    }
	  }
	}
      }
    
    // final enumeration...
    state_set_type& states = agenda_[length_max];
    
    if (states.empty())
      return - std::numeric_limits<double>::infinity();
    
    heap_.clear();
    heap_.reserve(states.size());
    
    // here, we use loss-biased search...
    state_set_type::const_iterator siter_end = states.end();
    for (state_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
      heap_.push_back(score_state_type((*siter)->score_ + (*siter)->error_, *siter));
      std::push_heap(heap_.begin(), heap_.end(), heap_compare());
    }
    
    heap_type::iterator hiter_begin = heap_.begin();
    heap_type::iterator hiter       = heap_.end();
    heap_type::iterator hiter_end   = heap_.end();
    
    // perform sorting...
    bool has_error = false;

    for (/**/; hiter_begin != hiter && std::distance(hiter, hiter_end) != beam_; -- hiter) {
      std::pop_heap(hiter_begin, hiter, heap_compare());
      
      has_error |= ((hiter - 1)->second->error_ > 0);
    }
    
    // search until we find an error...
    if (! has_error && errors_[length_max])
      for (/**/; hiter_begin != hiter && ! has_error; -- hiter) {
	std::pop_heap(hiter_begin, hiter, heap_compare());
	
	has_error |= ((hiter - 1)->second->error_ > 0);
      }
    
    for (heap_type::iterator iter = hiter_begin ; iter != hiter; ++ iter)
      if (! iter->second->terminal())
	deallocate(iter->second);
    
    // allocate new states...
    chart_sampled_(0, source_size, 0, target_size).clear();
    for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter)
      chart_sampled_(0, source_size, 0, target_size).push_back(iter->second);

    return chart_sampled_(0, source_size, 0, target_size).back()->score_;
  }

  void forward_backward(rest_cost_set_type& costs)
  {
    const size_type sentence_size = costs.size() - 1;
    
    // forward...
    costs[0].alpha_ = 0;
    for (size_type last = 1; last <= sentence_size; ++ last) {
      const size_type first = last - 1;
      
      costs[last].alpha_ = std::max(costs[last].alpha_, costs[first].alpha_ + costs[first].score_);
    }
    
    // backward...
    costs[sentence_size].beta_ = 0;
    for (difference_type first = sentence_size - 1; first >= 0; -- first) {
      const size_type last = first + 1;
      
      costs[first].beta_ = std::max(costs[first].beta_, costs[first].score_ + costs[last].beta_);
    }
  }
  
  // terminal rules
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const span_pair_type& span,
	       const model_type& theta)
  {
    typedef Eigen::Map<tensor_type> matrix_type;

    const difference_type embedding_size = theta.embedding_;
    const difference_type hidden_size    = theta.hidden_;
    const difference_type window_size    = window_;
    
    const difference_type source_size = source.size();
    const difference_type target_size = target.size();
    
    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_size;

    const word_type& word_source = (span.source_.empty() ? vocab_type::EPSILON : source[span.source_.first_]);
    const word_type& word_target = (span.target_.empty() ? vocab_type::EPSILON : target[span.target_.first_]);

    state_type& state = allocate();
    
    state.span_  = span;
    state.left_  = 0;
    state.right_ = 0;
    
    state.layer_ = (theta.bt_
		    + (theta.Wt_.block(0, offset_source, hidden_size, embedding_size)
		       * theta.source_.col(word_source.id())
		       * theta.scale_)
		    + (theta.Wt_.block(0, offset_target, hidden_size, embedding_size)
		       * theta.target_.col(word_target.id())
		       * theta.scale_)).array().unaryExpr(hinge());
    state.delta_.setZero();
    
    state.score_ = (theta.Wc_ * state.layer_ + theta.bc_)(0, 0);
    state.error_ = 0;
    state.loss_  = 0;

    state.source_ = word_source;
    state.target_ = word_target;
    
    if (! span.source_.empty())
      costs_source_[span.source_.first_].score_ = std::max(costs_source_[span.source_.first_].score_, state.score_);
    if (! span.target_.empty())
      costs_target_[span.target_.first_].score_ = std::max(costs_target_[span.target_.first_].score_, state.score_);
    
    // put into chart and agenda!
    chart_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(&state);
    agenda_[span.size()].push_back(&state);
    errors_[span.size()] |= (state.error_ > 0);
  }

  // perform sampling...
  template <typename Gen>
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const word_set_type& sources,
	       const word_set_type& targets,
	       const span_pair_type& span,
	       const model_type& theta,
	       Gen& gen)
  {
    typedef Eigen::Map<tensor_type> matrix_type;

    const difference_type embedding_size = theta.embedding_;
    const difference_type hidden_size    = theta.hidden_;
    const difference_type window_size    = window_;
    
    const difference_type source_size = source.size();
    const difference_type target_size = target.size();
    
    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_size;

    const word_type& word_source = (span.source_.empty() ? vocab_type::EPSILON : source[span.source_.first_]);
    const word_type& word_target = (span.target_.empty() ? vocab_type::EPSILON : target[span.target_.first_]);

    state_type& state = allocate();
    
    state.span_  = span;
    state.left_  = 0;
    state.right_ = 0;
    
    state.layer_ = (theta.bt_
		    + (theta.Wt_.block(0, offset_source, hidden_size, embedding_size)
		       * theta.source_.col(word_source.id())
		       * theta.scale_)
		    + (theta.Wt_.block(0, offset_target, hidden_size, embedding_size)
		       * theta.target_.col(word_target.id())
		       * theta.scale_)).array().unaryExpr(hinge());
    state.delta_.setZero();
    
    state.score_ = (theta.Wc_ * state.layer_ + theta.bc_)(0, 0);
    state.error_ = 0;
    state.loss_  = 0;
    
    state.source_ = word_source;
    state.target_ = word_target;
    
    if (! span.source_.empty())
      costs_source_[span.source_.first_].score_ = std::max(costs_source_[span.source_.first_].score_, state.score_ + state.error_);
    if (! span.target_.empty())
      costs_target_[span.target_.first_].score_ = std::max(costs_target_[span.target_.first_].score_, state.score_ + state.error_);
    
    // put into chart and agenda!
    chart_sampled_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(&state);
    agenda_[span.size()].push_back(&state);
    errors_[span.size()] |= (state.error_ > 0);
    
    boost::random::uniform_int_distribution<> uniform_source(0, source_size - 1);
    boost::random::uniform_int_distribution<> uniform_target(0, target_size - 1);
    
    if (! span.source_.empty() && dict_target_source_.size(word_target) > 1) {
      sampled_.clear();
      
      for (size_type sample = 0; sample != samples_; ++ sample) {
	word_type sampled_source = word_source;
	
	if (word_target == vocab_type::EPSILON) {
	  while (sampled_source == word_source)
	    sampled_source = dict_target_source_.draw(target[uniform_target(gen)], gen);
	} else {
	  while (sampled_source == word_source)
	    sampled_source = dict_target_source_.draw(word_target, gen);
	}
	
	// check if already sampled...
	if (sources.find(sampled_source) != sources.end() || ! sampled_.insert(sampled_source).second) continue;
	
	state_type& state = allocate();
	
	state.span_  = span;
	state.left_  = 0;
	state.right_ = 0;
	
	state.layer_ = (theta.bt_
			+ (theta.Wt_.block(0, offset_source, hidden_size, embedding_size)
			   * theta.source_.col(sampled_source.id())
			   * theta.scale_)
			+ (theta.Wt_.block(0, offset_target, hidden_size, embedding_size)
			   * theta.target_.col(word_target.id())
			   * theta.scale_)).array().unaryExpr(hinge());
	state.delta_.setZero();
	
	state.score_ = (theta.Wc_ * state.layer_ + theta.bc_)(0, 0);
	state.error_ = 1;
	state.loss_  = 0;
	
	state.source_ = sampled_source;
	state.target_ = word_target;
	
	// here, we use loss-biased search...
	if (! span.source_.empty())
	  costs_source_[span.source_.first_].score_ = std::max(costs_source_[span.source_.first_].score_, state.score_ + state.error_);
	if (! span.target_.empty())
	  costs_target_[span.target_.first_].score_ = std::max(costs_target_[span.target_.first_].score_, state.score_ + state.error_);
	
	// put into chart and agenda!
	chart_sampled_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(&state);
	agenda_[span.size()].push_back(&state);
	errors_[span.size()] |= true;
      }
    }
    
    if (! span.target_.empty() && dict_source_target_.size(word_source) > 1) {
      sampled_.clear();
      
      for (size_type sample = 0; sample != samples_; ++ sample) {
	word_type sampled_target = word_target;
	
	if (word_source == vocab_type::EPSILON) {
	  while (sampled_target == word_target)
	    sampled_target = dict_source_target_.draw(source[uniform_source(gen)], gen);
	} else {
	  while (sampled_target == word_target)
	    sampled_target = dict_source_target_.draw(word_source, gen);
	}
	
	// check if already sampled...
	if (targets.find(sampled_target) != targets.end() || ! sampled_.insert(sampled_target).second) continue;
	
	state_type& state = allocate();
	
	state.span_  = span;
	state.left_  = 0;
	state.right_ = 0;
	
	state.layer_ = (theta.bt_
			+ (theta.Wt_.block(0, offset_source, hidden_size, embedding_size)
			   * theta.source_.col(word_source.id())
			   * theta.scale_)
			+ (theta.Wt_.block(0, offset_target, hidden_size, embedding_size)
			   * theta.target_.col(sampled_target.id())
			   * theta.scale_)).array().unaryExpr(hinge());
	state.delta_.setZero();
	
	state.score_ = (theta.Wc_ * state.layer_ + theta.bc_)(0, 0);
	state.error_ = 1;
	state.loss_  = 0;
	
	state.source_ = word_source;
	state.target_ = sampled_target;
	
	// here, we use loss-biased search...
	if (! span.source_.empty())
	  costs_source_[span.source_.first_].score_ = std::max(costs_source_[span.source_.first_].score_, state.score_ + state.error_);
	if (! span.target_.empty())
	  costs_target_[span.target_.first_].score_ = std::max(costs_target_[span.target_.first_].score_, state.score_ + state.error_);
	
	// put into chart and agenda!
	chart_sampled_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_).push_back(&state);
	agenda_[span.size()].push_back(&state);
	errors_[span.size()] |= true;
      }
    }
  }

  size_type span_index(const span_pair_type& span, const model_type& theta)
  {
    // if span.size() is: (differentiate by the lengths of two)
    //  1 --  4 will be 0
    //  5 --  8 will be 1
    //  9 -- 12 will be 2
    // 13 -- 16 will be 3
    // 17 -- 20 will be 4
    // 21 -- 24 will be 5
    // 25 -- 28 will be 6
    // 29 -- 32 will be 7
    // 33 -- 36 will be 8
    
    return utils::bithack::min((span.size() - 1) / 4, theta.span_);
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
    typedef Eigen::Map<tensor_type> matrix_type;
    
    const difference_type embedding_size = theta.embedding_;
    const difference_type hidden_size    = theta.hidden_;
    const difference_type window_size    = window_;
    
    const difference_type source_size = source.size();
    const difference_type target_size = target.size();
    
    const difference_type offset_left  = 0;
    const difference_type offset_right = hidden_size;
    
    const span_pair_type& span1 = state1->span_;
    const span_pair_type& span2 = state2->span_;
    
    //std::cerr << "span: " << span << " left: " << span1 << " right: " << span2 << std::endl;

    const size_type offset_span = span_index(span, theta) * hidden_size;
    
    state_type& state = allocate();
    
    const tensor_type& Wr = (straight ? theta.Ws_ : theta.Wi_);
    const tensor_type& br = (straight ? theta.bs_ : theta.bi_);
    
    state.span_  = span;
    state.left_  = state1;
    state.right_ = state2;
    
    state.layer_ = (br.block(offset_span, 0, hidden_size, 1)
		    + Wr.block(offset_span, offset_left,  hidden_size, hidden_size) * state1->layer_
		    + Wr.block(offset_span, offset_right, hidden_size, hidden_size) * state2->layer_).array().unaryExpr(hinge());
    state.delta_.setZero();
    
    state.score_ = state1->score_ + state2->score_ + (theta.Wc_ * state.layer_ + theta.bc_)(0, 0);
    state.error_ = state1->error_ + state2->error_;
    state.loss_  = 0;
    
    // put into agenda, but not chart...
    agenda_[span.size()].push_back(&state);
    errors_[span.size()] |= (state.error_ > 0);
  }
  
  
  // we will perform sampling at the same time...
  template <typename Gen>
  double backward(const sentence_type& source,
		  const sentence_type& target,
		  const model_type& theta,
		  gradient_type& gradient,
		  Gen& gen)
  {
    typedef cicada::semiring::Logprob<double> weight_type;
    typedef Eigen::Map<tensor_type> matrix_type;
    
#if 0
    std::cerr << "backward source: " << source << std::endl
	      << "backward target: " << target << std::endl;
#endif

    const difference_type embedding_size = theta.embedding_;
    const difference_type hidden_size    = theta.hidden_;
    const difference_type window_size    = window_;
    
    const difference_type source_size = source.size();
    const difference_type target_size = target.size();

    const difference_type offset_left  = 0;
    const difference_type offset_right = hidden_size;
    
    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_size;
    
    const size_type length_max = source_size + target_size;
    
    backwards_.clear();
    backwards_.reserve(length_max + 1);
    backwards_.resize(length_max + 1);
    
    // compute hinge-loss    
    const state_set_type& states_correct = chart_(0, source_size, 0, target_size);
    const state_set_type& states_mistake = chart_sampled_(0, source_size, 0, target_size);

    state_set_type::const_iterator citer_begin = states_correct.begin();
    state_set_type::const_iterator citer_end   = states_correct.end();
    
    state_set_type::const_iterator miter_begin = states_mistake.begin();
    state_set_type::const_iterator miter_end   = states_mistake.end();
    
    // no errors!
    if (citer_begin == citer_end || miter_begin == miter_end) return 0.0;
    
    weight_type Z_correct;
    weight_type Z_mistake;
    size_type num_mistake = 0;
    
    for (state_set_type::const_iterator miter = miter_begin; miter != miter_end; ++ miter)
      if ((*miter)->error_ > 0) {
	Z_mistake += cicada::semiring::traits<weight_type>::exp((*miter)->score_);
	++ num_mistake;
      }
    
    if (! num_mistake) return 0.0;
    
    for (state_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) 
      Z_correct += cicada::semiring::traits<weight_type>::exp((*citer)->score_);

    bool found = false;
    double loss = 0.0;
    
    for (state_set_type::const_iterator miter = miter_begin; miter != miter_end; ++ miter)
      if ((*miter)->error_ > 0)
	for (state_set_type::const_iterator citer = citer_begin; citer != citer_end; ++ citer) {
	  const double error = std::max(double((*miter)->error_) - ((*citer)->score_ - (*miter)->score_), 0.0);
	  
	  if (error == 0.0) continue;
	  
	  const weight_type prob_correct = cicada::semiring::traits<weight_type>::exp((*citer)->score_) / Z_correct;
	  const weight_type prob_mistake = cicada::semiring::traits<weight_type>::exp((*miter)->score_) / Z_mistake;
	  
	  const double error_factor = prob_correct * prob_mistake;
	  
	  backwards_[(*citer)->span_.size()].insert(*citer);
	  backwards_[(*miter)->span_.size()].insert(*miter);
	  
	  const_cast<state_type&>(*(*citer)).loss_ -= error_factor;
	  const_cast<state_type&>(*(*miter)).loss_ += error_factor;
	  
	  if (! (*citer)->delta_.rows()) 
	    const_cast<state_type&>(*(*citer)).delta_ = tensor_type::Zero(hidden_size, 1);
	  if (! (*miter)->delta_.rows())
	    const_cast<state_type&>(*(*miter)).delta_ = tensor_type::Zero(hidden_size, 1);
	  
	  loss += error * error_factor;
	  found = true;
	}
    
    if (! found) return 0.0;
    
    //std::cerr << "loss: " << loss << std::endl;
    
    ++ gradient.count_;
    
    for (size_type length = length_max; length != 0; -- length) {
      //std::cerr << "backward: " << length << std::endl;

      const state_unique_type& states = backwards_[length];

      state_unique_type::const_iterator siter_end = states.end();
      for (state_unique_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
	state_type& state = const_cast<state_type&>(*(*siter));

#if 0
	std::cerr << "span: " << state.span_
		  << " score: " << state.score_
		  << " loss: " << state.loss_
		  << " error: " << state.error_
		  << std::endl;
#endif

	// classification
	gradient.Wc_         += state.loss_ * state.layer_.transpose();
	gradient.bc_.array() += state.loss_;
	
	// propagate to delta...
	state.delta_.array() += state.layer_.array().unaryExpr(dhinge()) * (theta.Wc_.transpose().array() * state.loss_);
	
	if (state.terminal()) {
	  gradient.Wt_.block(0, offset_source, hidden_size, embedding_size)
	    += state.delta_ * theta.source_.col(state.source_.id()).transpose() * theta.scale_;
	  gradient.Wt_.block(0, offset_target, hidden_size, embedding_size)
	    += state.delta_ * theta.target_.col(state.target_.id()).transpose() * theta.scale_;
	  gradient.bt_ += state.delta_;
	  
	  gradient.source(state.source_)
	    += theta.Wt_.block(0, offset_source, hidden_size, embedding_size).transpose() * state.delta_;
	  gradient.target(state.target_)
	    += theta.Wt_.block(0, offset_target, hidden_size, embedding_size).transpose() * state.delta_;
	} else {
	  backwards_[state.left_->span_.size()].insert(state.left_);
	  backwards_[state.right_->span_.size()].insert(state.right_);
	  
	  const bool straight = state.straight();
	  const size_type offset_span = span_index(state.span_, theta) * hidden_size;
	  
	  const tensor_type& Wr = (straight ? theta.Ws_ : theta.Wi_);
	  
	  tensor_type& dWr = (straight ? gradient.Ws_ : gradient.Wi_);
	  tensor_type& dbr = (straight ? gradient.bs_ : gradient.bi_);

	  state_type& state_left  = const_cast<state_type&>(*state.left_);
	  state_type& state_right = const_cast<state_type&>(*state.right_);
	  
	  // propagate loss
	  state_left.loss_  += state.loss_;
	  state_right.loss_ += state.loss_;
	  
	  dWr.block(offset_span, offset_left,  hidden_size, hidden_size) += state.delta_ * state_left.layer_.transpose();
	  dWr.block(offset_span, offset_right, hidden_size, hidden_size) += state.delta_ * state_right.layer_.transpose();
	  dbr.block(offset_span, 0, hidden_size, 1)                      += state.delta_;
	  
	  tensor_type& delta_left  = state_left.delta_;
	  tensor_type& delta_right = state_right.delta_;

	  if (! delta_left.rows())
	    delta_left = tensor_type::Zero(hidden_size, 1);
	  if (! delta_right.rows())
	    delta_right = tensor_type::Zero(hidden_size, 1);
	  
	  delta_left.array()  += (state_left.layer_.array().unaryExpr(dhinge())
				  * (Wr.block(offset_span, offset_left,  hidden_size, hidden_size).transpose()
				     * state.delta_).array());
	  delta_right.array() += (state_right.layer_.array().unaryExpr(dhinge())
				  * (Wr.block(offset_span, offset_right, hidden_size, hidden_size).transpose()
				     * state.delta_).array());
	}
      }
    }

    return loss;
  }
  
  void derivation(const sentence_type& source,
		  const sentence_type& target,
		  derivation_type& d)
  {
    d.clear();
    
    stack_.clear();
    stack_.push_back(chart_(0, source.size(), 0, target.size()).back());
    
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

  void deallocate(const state_type* state)
  {
    const_cast<state_type*>(state)->left_ = linked_;
    linked_ = state;
  }

  state_type& allocate()
  {
    if (linked_) {
      const state_type* state = linked_;
      linked_ = state->left_;
      return const_cast<state_type&>(*state);
    }

    states_.push_back(state_type());
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
	       const size_type& span,
	       const double& lambda,
	       const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      span_(span),
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
    
    Wc_ = tensor_type::Zero(1, hidden_);
    bc_ = tensor_type::Zero(1, 1);
    
    Ws_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    bs_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);

    Wi_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    bi_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    Wt_ = tensor_type::Zero(hidden_, embedding_ * 2);
    bt_ = tensor_type::Zero(hidden_, 1);
  }

  
  void operator()(model_type& theta,
		  const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    if (! gradient.count_) return;

    LearnAdaGrad& optimizer = const_cast<LearnAdaGrad&>(*this);
    
    const double scale = 1.0 / gradient.count_;

    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     optimizer.source_,
	     siter->second,
	     scale);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     optimizer.target_,
	     titer->second,
	     scale);

    update(theta.Wc_, optimizer.Wc_, gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, optimizer.bc_, gradient.bc_, scale, false);
    
    update(theta.Ws_, optimizer.Ws_, gradient.Ws_, scale, lambda_ != 0.0);
    update(theta.bs_, optimizer.bs_, gradient.bs_, scale, false);

    update(theta.Wi_, optimizer.Wi_, gradient.Wi_, scale, lambda_ != 0.0);
    update(theta.bi_, optimizer.bi_, gradient.bi_, scale, false);
    
    update(theta.Wt_, optimizer.Wt_, gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, optimizer.bt_, gradient.bt_, scale, false);
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
  size_type span_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  // classification
  tensor_type Wc_;
  tensor_type bc_;

  // straight
  tensor_type Ws_;
  tensor_type bs_;

  // inversion
  tensor_type Wi_;
  tensor_type bi_;  
  
  // terminal
  tensor_type Wt_;
  tensor_type bt_;
};

struct LearnAdaDelta
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model     model_type;
  typedef Gradient  gradient_type;

  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaDelta(const size_type& embedding,
		const size_type& hidden,
		const size_type& span,
		const double& lambda,
		const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      span_(span),
      lambda_(lambda),
      eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    
    gsource_ = tensor_type::Zero(embedding_, vocabulary_size);
    gtarget_ = tensor_type::Zero(embedding_, vocabulary_size);
    xsource_ = tensor_type::Zero(embedding_, vocabulary_size);
    xtarget_ = tensor_type::Zero(embedding_, vocabulary_size);
    
    gWc_ = tensor_type::Zero(1, hidden_);
    gbc_ = tensor_type::Zero(1, 1);
    xWc_ = tensor_type::Zero(1, hidden_);
    xbc_ = tensor_type::Zero(1, 1);
    
    gWs_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    gbs_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    xWs_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    xbs_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    gWi_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    gbi_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    xWi_ = tensor_type::Zero(hidden_ * (span_ + 1), hidden_ + hidden_);
    xbi_ = tensor_type::Zero(hidden_ * (span_ + 1), 1);
    
    gWt_ = tensor_type::Zero(hidden_, embedding_ * 2);
    gbt_ = tensor_type::Zero(hidden_, 1);
    xWt_ = tensor_type::Zero(hidden_, embedding_ * 2);
    xbt_ = tensor_type::Zero(hidden_, 1);
  }

  
  void operator()(model_type& theta,
		  const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type gradient_embedding_type;

    if (! gradient.count_) return;

    LearnAdaDelta& optimizer = const_cast<LearnAdaDelta&>(*this);
    
    const double scale = 1.0 / gradient.count_;

    gradient_embedding_type::const_iterator siter_end = gradient.source_.end();
    for (gradient_embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter)
      update(siter->first,
	     theta.source_,
	     optimizer.gsource_,
	     optimizer.xsource_,
	     siter->second,
	     scale);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     optimizer.gtarget_,
	     optimizer.xtarget_,
	     titer->second,
	     scale);

    update(theta.Wc_, optimizer.gWc_, optimizer.xWc_, gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, optimizer.gbc_, optimizer.xbc_, gradient.bc_, scale, false);
    
    update(theta.Ws_, optimizer.gWs_, optimizer.xWs_, gradient.Ws_, scale, lambda_ != 0.0);
    update(theta.bs_, optimizer.gbs_, optimizer.xbs_, gradient.bs_, scale, false);

    update(theta.Wi_, optimizer.gWi_, optimizer.xWi_, gradient.Wi_, scale, lambda_ != 0.0);
    update(theta.bi_, optimizer.gbi_, optimizer.xbi_, gradient.bi_, scale, false);
    
    update(theta.Wt_, optimizer.gWt_, optimizer.xWt_, gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, optimizer.gbt_, optimizer.xbt_, gradient.bt_, scale, false);
  }

  template <typename Theta, typename GradVar, typename XVar, typename Grad>
  struct update_visitor_regularize
  {
    update_visitor_regularize(Eigen::MatrixBase<Theta>& theta,
			      Eigen::MatrixBase<GradVar>& G,
			      Eigen::MatrixBase<XVar>&    X,
			      const Eigen::MatrixBase<Grad>& g,
			      const double& scale,
			      const double& lambda,
			      const double& eta0)
      : theta_(theta), G_(G), X_(X), g_(g), scale_(scale), lambda_(lambda), eta0_(eta0) {}
    
    void init(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      operator()(value, i, j);
    }
    
    void operator()(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      if (g_(i, j) == 0) return;
      
      G_(i, j) = G_(i, j) * 0.99 + g_(i, j) * g_(i, j) * scale_ * scale_;
      
      const double rate = std::sqrt(eta0_ + X_(i, j)) / std::sqrt(eta0_ + G_(i, j));
      const double f = theta_(i, j) - rate * scale_ * g_(i, j);
      const double x = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
      
      X_(i, j) = X_(i, j) * 0.99 + (x - theta_(i, j)) * (x - theta_(i, j));
      
      theta_(i, j) = x;
    }
    
    Eigen::MatrixBase<Theta>&      theta_;
    Eigen::MatrixBase<GradVar>&    G_;
    Eigen::MatrixBase<XVar>&       X_;
    const Eigen::MatrixBase<Grad>& g_;
    
    const double scale_;
    const double lambda_;
    const double eta0_;
  };

  template <typename Theta, typename GradVar, typename XVar, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      Eigen::MatrixBase<XVar>& X,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool regularize=true) const
  {
    update_visitor_regularize<Theta, GradVar, XVar, Grad> visitor(theta, G, X, g, scale, regularize ? lambda_ : 0.0, eta0_);
    
    theta.visit(visitor);
  }

  template <typename Theta, typename GradVar, typename XVar, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      Eigen::MatrixBase<XVar>& X,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale) const
  {
    for (int row = 0; row != g.rows(); ++ row) 
      if (g(row, 0) != 0) {
	G(row, word.id()) = G(row, word.id()) * 0.99 + g(row, 0) * g(row, 0) * scale * scale;
	
	const double rate = std::sqrt(eta0_ + X(row, word.id())) / std::sqrt(eta0_ + G(row, word.id()));
	const double f = theta(row, word.id()) - rate * scale * g(row, 0);
	const double x = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	
	X(row, word.id()) = X(row, word.id()) * 0.99 + (x - theta(row, word.id())) * (x - theta(row, word.id()));
	
	theta(row, word.id()) = x;
      }
  }
  
  size_type embedding_;
  size_type hidden_;
  size_type span_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type gsource_;
  tensor_type gtarget_;
  tensor_type xsource_;
  tensor_type xtarget_;
  
  // classification
  tensor_type gWc_;
  tensor_type gbc_;
  tensor_type xWc_;
  tensor_type xbc_;

  // straight
  tensor_type gWs_;
  tensor_type gbs_;
  tensor_type xWs_;
  tensor_type xbs_;

  // inversion
  tensor_type gWi_;
  tensor_type gbi_;  
  tensor_type xWi_;
  tensor_type xbi_;  
  
  // terminal
  tensor_type gWt_;
  tensor_type gbt_;
  tensor_type xWt_;
  tensor_type xbt_;

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

    if (! gradient.count_) return;

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

    update(theta.Ws_, gradient.Ws_, scale, lambda_ != 0.0);
    update(theta.bs_, gradient.bs_, scale, false);

    update(theta.Wi_, gradient.Wi_, scale, lambda_ != 0.0);
    update(theta.bi_, gradient.bi_, scale, false);
    
    update(theta.Wt_, gradient.Wt_, scale, lambda_ != 0.0);
    update(theta.bt_, gradient.bt_, scale, false);
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

#endif
