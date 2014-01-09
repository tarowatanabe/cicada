//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_RNN_JOINT_ALIGNMENT_IMPL__HPP__
#define __CICADA_RNN_JOINT_ALIGNMENT_IMPL__HPP__ 1

#include <cstdlib>
#include <cmath>
#include <climits>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <set>

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
  
  Gradient() : embedding_(0), hidden_(0), window_(), alignment_(0), count_(0), shared_(0) {}
  Gradient(const size_type& embedding,
	   const size_type& hidden,
	   const size_type& window,
	   const size_type& alignment) 
    : embedding_(embedding),
      hidden_(hidden),
      window_(window),
      alignment_(alignment),
      count_(0),
      shared_(0)
  { initialize(embedding, hidden, window, alignment); }
  
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

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    //if (! Wt_.rows())
    //  Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    //if (! bt_.rows())
    //  bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());
    
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

    Wc_ -= x.Wc_;
    bc_ -= x.bc_;
    
    //Wt_ -= x.Wt_;
    //bt_ -= x.bt_;
    
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

    if (! Wc_.rows())
      Wc_ = tensor_type::Zero(x.Wc_.rows(), x.Wc_.cols());
    if (! bc_.rows())
      bc_ = tensor_type::Zero(x.bc_.rows(), x.bc_.cols());

    //if (! Wt_.rows())
    //  Wt_ = tensor_type::Zero(x.Wt_.rows(), x.Wt_.cols());
    //if (! bt_.rows())
    //  bt_ = tensor_type::Zero(x.bt_.rows(), x.bt_.cols());

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

    Wc_ += x.Wc_;
    bc_ += x.bc_;

    //Wt_ += x.Wt_;
    //bt_ += x.bt_;
    
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

    Wc_.setZero();
    bc_.setZero();
    
    //Wt_.setZero();
    //bt_.setZero();
    
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
      embedding = tensor_type::Zero(embedding_, 1);
    
    return embedding;
  }

  
  void initialize(const size_type embedding, const size_type hidden, const size_type window, const size_type alignment)
  {
    if (hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (alignment <= 0)
      throw std::runtime_error("invalid alignment");
    
    embedding_ = embedding;
    hidden_    = hidden;
    window_    = window;
    alignment_ = alignment;
    
    clear();
    
    const size_type state_size = embedding_ * (window * 2 + 1) * 2 + hidden_;
    
    // initialize...
    Wc_ = tensor_type::Zero(1, hidden_);
    bc_ = tensor_type::Zero(1, 1);
    
    //Wt_ = tensor_type::Zero(hidden_, state_size);
    //bt_ = tensor_type::Zero(hidden_, 1);

    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), state_size);
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1);
    
    Wn_ = tensor_type::Zero(hidden_, state_size);
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
    os.write((char*) &window_,    sizeof(size_type));
    os.write((char*) &alignment_, sizeof(size_type));
    os.write((char*) &count_,     sizeof(size_type));
    
    write(os, Wc_);
    write(os, bc_);
    
    //write(os, Wt_);
    //write(os, bt_);

    write(os, Wa_);
    write(os, ba_);

    write(os, Wn_);
    write(os, bn_);

    write(os, bi_);
    
    write(os, source_);
    write(os, target_);
  }
  
  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &window_,    sizeof(size_type));
    is.read((char*) &alignment_, sizeof(size_type));
    is.read((char*) &count_,     sizeof(size_type));

    read(is, Wc_);
    read(is, bc_);
    
    //read(is, Wt_);
    //read(is, bt_);

    read(is, Wa_);
    read(is, ba_);

    read(is, Wn_);
    read(is, bn_);

    read(is, bi_);
    
    read(is, source_);
    read(is, target_);
  }

  void write(std::ostream& os, const embedding_type& embedding) const
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

  void read(std::istream& is, embedding_type& embedding)
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
      
      matrix.resize(embedding_, 1);
      
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
  size_type window_;
  size_type alignment_;
  
  // embedding
  embedding_type source_;
  embedding_type target_;

  tensor_type Wc_;
  tensor_type bc_;
  
  //tensor_type Wt_;
  //tensor_type bt_;
  
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
  
  Model() : embedding_(0), hidden_(0), window_(0), alignment_(0), scale_(1) {}  
  template <typename Words, typename Gen>
  Model(const size_type& embedding,
	const size_type& hidden,
	const size_type& window,
	const size_type& alignment,
	Words& words_source,
	Words& words_target,
	Gen& gen) 
    : embedding_(embedding),
      hidden_(hidden),
      window_(window),
      alignment_(alignment),
      scale_(1)
  { initialize(embedding, hidden, window, alignment, words_source, words_target, gen); }
  
  void clear()
  {
    // embedding
    source_.setZero();
    target_.setZero();

    Wc_.setZero();
    bc_.setZero();
    
    //Wt_.setZero();
    //bt_.setZero();
    
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
		  const size_type window,
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
    window_    = window;
    alignment_ = alignment;
    
    clear();

    const size_type state_size = embedding_ * (window * 2 + 1) * 2 + hidden_;
    
    const size_type vocabulary_size = word_type::allocated();
    
    const double range_e = std::sqrt(6.0 / (embedding_ + 1));
    const double range_c = std::sqrt(6.0 / (hidden_ + 1));
    const double range_t = std::sqrt(6.0 / (hidden_ + state_size));
    const double range_a = std::sqrt(6.0 / (hidden_ + state_size));
    const double range_n = std::sqrt(6.0 / (hidden_ + state_size));
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
 
    source_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    target_ = tensor_type::Zero(embedding_, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));

    Wc_ = tensor_type::Zero(1, hidden_).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Ones(1, 1);
    
    //Wt_ = tensor_type::Zero(hidden_, state_size).array().unaryExpr(randomize<Gen>(gen, range_t));
    //bt_ = tensor_type::Zero(hidden_, 1);
    
    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), state_size).array().unaryExpr(randomize<Gen>(gen, range_a));
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1);

    Wn_ = tensor_type::Zero(hidden_, state_size).array().unaryExpr(randomize<Gen>(gen, range_n));
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
    target_ *= scale_;

    scale_ = 1.0;
  }
  
  Model& operator+=(const Model& x)
  {
    if (scale_ != x.scale_)
      throw std::runtime_error("different scaling");

    source_ += x.source_;
    target_ += x.target_;

    Wc_ += x.Wc_;
    bc_ += x.bc_;
    
    //Wt_ += x.Wt_;
    //bt_ += x.bt_;

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

    Wc_ *= x;
    bc_ *= x;
    
    //Wt_ *= x;
    //bt_ *= x;

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
  
  template <typename Rep, typename Tp>
  inline
  bool key_value(const Rep& rep, const std::string& key, Tp& value)
  {
    typename Rep::const_iterator iter = rep.find(key);
    if (iter == rep.end())
      return false;
    value = utils::lexical_cast<Tp>(iter->second);
    return true;
  }

  void read(const path_type& path) 
  {
    // we use a repository structure...
    typedef utils::repository repository_type;

    if (! boost::filesystem::exists(path))
      throw std::runtime_error("no file? " + path.string());
    
    repository_type rep(path, repository_type::read);
    
    size_type embedding;
    size_type hidden;
    size_type window;
    size_type alignment;
    
    if (! key_value(rep, "embedding", embedding))
      throw std::runtime_error("no embedding?");
    if (! key_value(rep, "hidden", hidden))
      throw std::runtime_error("no hidden?");
    if (! key_value(rep, "window", window))
      throw std::runtime_error("no window?");

    if (embedding_ != embedding)
      throw std::runtime_error("embedding size differ");
    if (hidden_ != hidden)
      throw std::runtime_error("hidden size differ");
    if (window_ != window)
      throw std::runtime_error("window size differ");

    if (! boost::filesystem::exists(rep.path("source.gz")))
      throw std::runtime_error("no source embedding?");
    if (! boost::filesystem::exists(rep.path("target.gz")))
      throw std::runtime_error("no target embedding?");
    
    // read embedding...
    read_embedding(rep.path("source.gz"), rep.path("target.gz"));
    
    const bool has_alignment =  key_value(rep, "alignment", alignment);
    
    if (has_alignment) {
      // we are reading a complete data
      
      read(rep.path("Wc.bin"), Wc_);
      read(rep.path("bc.bin"), bc_);

      //read(rep.path("Wt.bin"), Wt_);
      //read(rep.path("bt.bin"), bt_);

      read(rep.path("Wa.bin"), Wa_);
      read(rep.path("ba.bin"), ba_);

      read(rep.path("Wn.bin"), Wn_);
      read(rep.path("bn.bin"), bn_);

      read(rep.path("bi.bin"), bi_);
    } else {
      // we are reading an incomplete data... we do not read Wa,ba,Wn,bn,bi
      // also, Wt and bt is assumed to be partially trained
      
      read(rep.path("Wc.bin"), Wc_);
      read(rep.path("bc.bin"), bc_);

      const size_type embedding_size = embedding_ * (window_ * 2 + 1) * 2;

      tensor_type Wt(hidden_, embedding_size);
      tensor_type bt(hidden_, 1);
      
      read(rep.path("Wt.bin"), Wt);
      read(rep.path("bt.bin"), bt);
      
      Wn_.block(0, 0, hidden_, embedding_size) = Wt;
      bn_.block(0, 0, hidden_, 1) = bt;
      
      for (size_type shift = 0; shift != window_ * 2 + 1; ++ shift) {
	Wa_.block(hidden_ * shift, 0, hidden_, embedding_size) = Wt;
	ba_.block(hidden_ * shift, 0, hidden_, 1) = bt;
      }
    }
  }
  
  
  void read(const path_type& path, tensor_type& matrix)
  {
    if (! boost::filesystem::exists(path))
      throw std::runtime_error("no file? " + path.string());
    
    const size_type file_size = boost::filesystem::file_size(path);
    
    if (file_size != sizeof(tensor_type::Scalar) * matrix.rows() * matrix.cols())
      throw std::runtime_error("read size differ!");
    
    utils::compress_istream is(path, 1024 * 1024);
    
    is.read((char*) matrix.data(), sizeof(tensor_type::Scalar) * matrix.rows() * matrix.cols());
  }
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);
    rep["window"]    = utils::lexical_cast<std::string>(window_);
    rep["alignment"] = utils::lexical_cast<std::string>(alignment_);
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), rep.path("vocab-source"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), rep.path("vocab-target"), target_, words_target_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
    
    //write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    //write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);
    
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
    os.write((char*) &window_,    sizeof(size_type));
    os.write((char*) &alignment_, sizeof(size_type));
    os.write((char*) &scale_,     sizeof(double));

    write(os, Wc_);
    write(os, bc_);
    
    //write(os, Wt_);
    //write(os, bt_);

    write(os, Wa_);
    write(os, ba_);
    
    write(os, Wn_);
    write(os, bn_);
    
    write(os, bi_);
    
    write_embedding(os, source_);
    write_embedding(os, target_);
  }

  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &window_,    sizeof(size_type));
    is.read((char*) &alignment_, sizeof(size_type));
    is.read((char*) &scale_,     sizeof(double));

    read(is, Wc_);
    read(is, bc_);

    //read(is, Wt_);
    //read(is, bt_);

    read(is, Wa_);
    read(is, ba_);

    read(is, Wn_);
    read(is, bn_);

    read(is, bi_);
    
    read_embedding(is, source_);
    read_embedding(is, target_);
    
    // checking...
  }

  void write_embedding(std::ostream& os, const tensor_type& embedding) const
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

  void read_embedding(std::istream& is, tensor_type& embedding)
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    
    size_type cols = 0;
    is.read((char*) &cols, sizeof(size_type));
    
    if (cols > embedding.cols())
      embedding.conservativeResize(embedding_, cols);
    
    for (size_type i = 0; i != cols; ++ i) {
      size_type word_size = 0;
      is.read((char*) &word_size, sizeof(size_type));
      
      buffer.resize(word_size);
      is.read((char*) &(*buffer.begin()), word_size);
      
      const word_type word(buffer.begin(), buffer.end());

      if (word.id() >= embedding.cols())
	embedding.conservativeResize(embedding_, word.id() + 1);
      
      is.read((char*) embedding.col(word.id()).data(), sizeof(tensor_type::Scalar) * (embedding_));
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
  size_type window_;
  size_type alignment_;

  word_unique_type words_source_;
  word_unique_type words_target_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;

  tensor_type Wc_;
  tensor_type bc_;
  
  //tensor_type Wt_;
  //tensor_type bt_;
  
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
  typedef int    error_type;
  typedef State  state_type;
  typedef cicada::Symbol word_type;
  
  typedef size_t     size_type;
  typedef ptrdiff_t  difference_type;

  static const size_type offset_prev   = 0;
  static const size_type offset_index  = offset_prev  + sizeof(pointer);
  static const size_type offset_score  = offset_index + sizeof(index_type);
  static const size_type offset_error  = offset_score + sizeof(parameter_type);
  static const size_type offset_target = offset_error + sizeof(error_type);
  static const size_type offset_matrix = (offset_target + sizeof(word_type) + 15) & (~15);
  
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

  inline       parameter_type& loss()       { return score(); }
  inline const parameter_type& loss() const { return score(); }

  inline       error_type& error()       { return *reinterpret_cast<error_type*>(base_ + offset_error); }
  inline const error_type& error() const { return *reinterpret_cast<const error_type*>(base_ + offset_error); }

  inline       word_type& target()       { return *reinterpret_cast<word_type*>(base_ + offset_target); }
  inline const word_type& target() const { return *reinterpret_cast<const word_type*>(base_ + offset_target); }
  
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
      return x.score() + x.error() < y.score() + y.error();
    }
  };

  typedef utils::compact_set<word_type,
			     utils::unassigned<word_type>, utils::unassigned<word_type>,
			     boost::hash<word_type>, std::equal_to<word_type>,
			     std::allocator<word_type> > word_set_type;

  typedef Average loss_type;
  
  HMM(const dictionary_type& dict,
      const size_type& sample,
      const size_type& beam)
    : dict_(dict), sample_(sample), beam_(beam) {}
  
  const dictionary_type& dict_;
  
  size_type sample_;
  size_type beam_;
  
  tensor_type delta_alpha_;
  
  state_allocator_type state_allocator_;
  heap_set_type  heaps_;
  heap_set_type  heaps_viterbi_;
  state_map_type states_;

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

  template <typename Embedding>
  void copy_embedding(const sentence_type& source,
		      const sentence_type& target,
		      const model_type& theta,
		      const difference_type source_pos,
		      const difference_type target_pos,
		      Eigen::MatrixBase<Embedding>& embedding)
  {
    const difference_type source_size = source.size();
    const difference_type target_size = target.size();
    
    const difference_type window_size = theta.window_;
    const difference_type embedding_size = theta.embedding_;
    const difference_type embedding_window_size = embedding_size * (theta.window_ * 2 + 1);
    
    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_window_size;
    
    if (source_pos >= source_size + 2) {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	embedding.block(offset_source + embedding_size * i, 0, embedding_size, 1).noalias()
	  = theta.source_.col(vocab_type::EPSILON.id()) * theta.scale_;
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (source_pos + shift <= 0
				 ? vocab_type::BOS
				 : (source_pos + shift > source_size
				    ? vocab_type::EOS
				    : source[source_pos + shift - 1]));
	
	embedding.block(offset_source + embedding_size * i, 0, embedding_size, 1).noalias()
	  = theta.source_.col(word.id()) * theta.scale_;
      }
    }
    
    
    for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
      const difference_type shift = i - window_size;
      const word_type& word = (target_pos + shift <= 0
			       ? vocab_type::BOS
			       : (target_pos + shift > target_size
				  ? vocab_type::EOS
				  : target[target_pos + shift - 1]));
      
      embedding.block(offset_target + embedding_size * i, 0, embedding_size, 1).noalias()
	= theta.target_.col(word.id()) * theta.scale_;
    }    
  }

  template <typename Embedding>
  void propagate_embedding(const sentence_type& source,
			   const sentence_type& target,
			   const model_type& theta,
			   gradient_type& gradient,
			   const difference_type source_pos,
			   const difference_type target_pos,
			   const word_type& target_next,
			   const Eigen::MatrixBase<Embedding>& embedding)
  {
    const difference_type source_size = source.size();
    const difference_type target_size = target.size();
    
    const difference_type window_size = theta.window_;
    const difference_type embedding_size = theta.embedding_;
    const difference_type embedding_window_size = embedding_size * (theta.window_ * 2 + 1);
    
    const difference_type offset_source = 0;
    const difference_type offset_target = embedding_window_size;
    
    if (source_pos >= source_size + 2) {
      tensor_type& dembedding = gradient.source(vocab_type::EPSILON);
      
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i)
	dembedding += embedding.block(offset_source + embedding_size * i, 0, embedding_size, 1);
    } else {
      for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
	const difference_type shift = i - window_size;
	const word_type& word = (source_pos + shift <= 0
				 ? vocab_type::BOS
				 : (source_pos + shift > source_size
				    ? vocab_type::EOS
				    : source[source_pos + shift - 1]));
	
	gradient.source(word) += embedding.block(offset_source + embedding_size * i, 0, embedding_size, 1);
      }
    }
    
    for (difference_type i = 0; i != window_size * 2 + 1; ++ i) {
      if (i == window_size)
	gradient.target(target_next) += embedding.block(offset_target + embedding_size * i, 0, embedding_size, 1);
      else {
	const difference_type shift = i - window_size;
	const word_type& word = (target_pos + shift <= 0
				 ? vocab_type::BOS
				 : (target_pos + shift > target_size
				    ? vocab_type::EOS
				    : target[target_pos + shift - 1]));
	
	gradient.target(word) += embedding.block(offset_target + embedding_size * i, 0, embedding_size, 1);
      }
    }        
  }

  template <typename Gen>
  double forward(const sentence_type& source,
		 const sentence_type& target,
		 const model_type& theta,
		 alignment_type& alignment,
		 Gen& gen)
  {
    typedef Eigen::Map<tensor_type> matrix_type;
   
#if 0
    std::cerr << "forward source: " << source << std::endl
	      << "forward target: " << target << std::endl;
#endif
 
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type embedding_window_size = theta.embedding_ * (theta.window_ * 2 + 1);
    const size_type state_size = embedding_window_size * 2 + theta.hidden_;
    const size_type embedding_size = embedding_window_size * 2;
    
    const size_type offset_embedding = 0;
    const size_type offset_source = 0;
    const size_type offset_target = embedding_window_size;
    const size_type offset_matrix = embedding_window_size * 2;
    
    state_allocator_.assign(state_type::size(state_size));
    
    alignment.clear();
    
    heaps_.clear();
    heaps_.resize(target_size + 2);
    
    heaps_viterbi_.clear();
    heaps_viterbi_.resize(target_size + 2);
    
    sources_.clear();
    targets_.clear();
    
    sources_.insert(source.begin(), source.end());
    targets_.insert(target.begin(), target.end());

    states_.clear();
    states_.resize(target_size + 2);

    state_type state_bos = state_allocator_.allocate();
    
    state_bos.prev() = state_type();
    state_bos.index() = 0;
    state_bos.score() = 0.0;
    state_bos.error() = 0;
    state_bos.target() = vocab_type::BOS;
    
    matrix_type alpha_bos(state_bos.matrix(), state_size, 1);
    
    alpha_bos.block(offset_matrix, 0, theta.hidden_, 1) = theta.bi_;
    
    heaps_[0].push_back(state_bos);
    heaps_viterbi_[0].push_back(state_bos);

    boost::random::uniform_int_distribution<> uniform_source(0, source_size - 1);
    
    for (size_type trg = 1; trg != target_size + 2; ++ trg) {
      const word_type target_next = (trg == 0
				     ? vocab_type::BOS
				     : (trg == target_size + 1
					? vocab_type::EOS
					: target[trg - 1]));
      
      const size_type next_first = utils::bithack::branch(trg == target_size + 1, source_size + 1, size_type(1));
      const size_type next_last  = utils::bithack::branch(trg == target_size + 1, source_size + 2, source_size + 1);

      {
	heap_type& heap      = heaps_[trg - 1];
	heap_type& heap_next = heaps_[trg];
	
	// perform pruning
	heap_type::iterator hiter_begin = heap.begin();
	heap_type::iterator hiter       = heap.end();
	heap_type::iterator hiter_end   = heap.end();
	
	if (heap.size() > beam_) {
	  bool has_error = false;
	  for (/**/; hiter_begin != hiter && hiter_end - hiter != beam_; -- hiter) {
	    std::pop_heap(hiter_begin, hiter, heap_compare());
	    
	    has_error |= ((hiter - 1)->error() > 0);
	  }

	  if (! has_error)
	    for (/**/; hiter_begin != hiter && ! has_error; -- hiter) {
	      std::pop_heap(hiter_begin, hiter, heap_compare());
	      
	      has_error |= ((hiter - 1)->error() > 0);
	    }
	  
	  // deallocate unused states
	  for (/**/; hiter_begin != hiter; ++ hiter_begin)
	    state_allocator_.deallocate(*hiter_begin);
	  
	} else
	  hiter = hiter_begin;
	
	for (/**/; hiter != hiter_end; ++ hiter) {
	  const state_type& state = *hiter;
	  
	  const matrix_type alpha_prev(state.matrix(), state_size, 1);
	  
	  const size_type prev = state.index();
	  
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
	    state_next.error() = state.error();
	    state_next.target() = target_next;
	    
	    matrix_type alpha_next(state_next.matrix(), state_size, 1);
	    
	    copy_embedding(source, target, theta, next, trg, alpha_next);
	    
	    // compute alpha-next
	    alpha_next.block(offset_matrix, 0, theta.hidden_, 1)
	      = (theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
		 + (theta.Wa_.block(theta.hidden_ * shift, offset_embedding, theta.hidden_, embedding_size)
		    * alpha_next.block(offset_embedding, 0, embedding_size, 1))
		 + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
		    * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	    
	    const double score = (theta.Wc_ * alpha_next.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);

	    state_next.score() = score + state.score();
	    
	    heap_next.push_back(state_next);
	    std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	    
	    if (dict_.size(source_next) > 1) {
	      sampled_.clear();

	      for (size_type k = 0; k != sample_; ++ k) {
		word_type target_sampled = dict_.draw(source_next, gen);
		
		while (target_sampled == target_next)
		  target_sampled = dict_.draw(source_next, gen);
		
		if (! sampled_.empty() && sampled_.find(target_sampled) != sampled_.end()) continue;
		
		sampled_.insert(target_sampled);
		
		state_type state_sampled = state_allocator_.allocate();
		state_sampled.prev() = state;
		state_sampled.index() = next;
		state_sampled.target() = target_sampled;
		state_sampled.error() = state.error() + 1;
	      
		matrix_type alpha_sampled(state_sampled.matrix(), state_size, 1);
	      
		// compute alpha-sampled
		alpha_sampled = alpha_next;
	      
		alpha_sampled.block(offset_target + theta.embedding_ * theta.window_, 0, theta.embedding_, 1)
		  = theta.target_.col(target_sampled.id()) * theta.scale_;
		
		alpha_sampled.block(offset_matrix, 0, theta.hidden_, 1)
		  = (theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
		     + (theta.Wa_.block(theta.hidden_ * shift, offset_embedding, theta.hidden_, embedding_size)
			* alpha_sampled.block(offset_embedding, 0, embedding_size, 1))
		     + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
			* alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
		
		const double score = (theta.Wc_ * alpha_sampled.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	      
		state_sampled.score() = score + state.score();
		
		heap_next.push_back(state_sampled);
		std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	      }
	    }
	  }
	  
	  // none...
	  if (trg != target_size + 1) {
	    const size_type next = utils::bithack::branch(prev >= source_size + 2, prev, prev + source_size + 2);
	    
	    state_type state_next = state_allocator_.allocate();
	    state_next.prev() = state;
	    state_next.index() = next;
	    state_next.error() = state.error();
	    state_next.target() = target_next;
	    
	    matrix_type alpha_next(state_next.matrix(), state_size, 1);
	    
	    copy_embedding(source, target, theta, next, trg, alpha_next);
	    
	    // compute alpha-next
	    alpha_next.block(offset_matrix, 0, theta.hidden_, 1)
	      = (theta.bn_
		 + (theta.Wn_.block(0, offset_embedding, theta.hidden_, embedding_size)
		    * alpha_next.block(offset_embedding, 0, embedding_size, 1))
		 + (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
		    * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	    
	    const double score = (theta.Wc_ * alpha_next.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	    
	    state_next.score() = score + state.score();
	    
	    heap_next.push_back(state_next);
	    std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());

	    sampled_.clear();
	    
	    for (size_type k = 0; k != sample_; ++ k) {
	      word_type target_sampled = dict_.draw(source[uniform_source(gen)], gen);
	      
	      //while (targets_.find(target_sampled) != targets_.end())
	      while (target_sampled == target_next)
		target_sampled = dict_.draw(source[uniform_source(gen)], gen);
	      
	      if (! sampled_.empty() && sampled_.find(target_sampled) != sampled_.end()) continue;
		
	      sampled_.insert(target_sampled);
	      
	      state_type state_sampled = state_allocator_.allocate();
	      state_sampled.prev() = state;
	      state_sampled.index() = next;
	      state_sampled.target() = target_sampled;
	      state_sampled.error() = state.error() + 1;
	      
	      matrix_type alpha_sampled(state_sampled.matrix(), state_size, 1);
	      
	      // compute alpha-sampled
	      alpha_sampled = alpha_next;
	      
	      alpha_sampled.block(offset_target + theta.embedding_ * theta.window_, 0, theta.embedding_, 1)
		= theta.target_.col(target_sampled.id()) * theta.scale_;
	      
	      alpha_sampled.block(offset_matrix, 0, theta.hidden_, 1)
		= (theta.bn_
		   + (theta.Wn_.block(0, offset_embedding, theta.hidden_, embedding_size)
		      * alpha_sampled.block(offset_embedding, 0, embedding_size, 1))
		   + (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
		      * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	      
	      const double score = (theta.Wc_ * alpha_sampled.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	      
	      state_sampled.score() = score + state.score();
	      
	      heap_next.push_back(state_sampled);
	      std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	    }
	  }
	}
      }

      // Viterbi iteration...
      {
	heap_type& heap      = heaps_viterbi_[trg - 1];
	heap_type& heap_next = heaps_viterbi_[trg];
	
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
	  
	  const matrix_type alpha_prev(state.matrix(), state_size, 1);
	  
	  const size_type prev = state.index();
	  
	  for (size_type next = next_first; next != next_last; ++ next) {
	    const size_type shift = theta.shift(source_size, target_size, prev, next);
	    
	    state_type state_next = state_allocator_.allocate();
	    state_next.prev() = state;
	    state_next.index() = next;
	    state_next.error() = state.error();
	    state_next.target() = target_next;
	    
	    matrix_type alpha_next(state_next.matrix(), state_size, 1);
	    
	    copy_embedding(source, target, theta, next, trg, alpha_next);
	    
	    // compute alpha-next
	    alpha_next.block(offset_matrix, 0, theta.hidden_, 1)
	      = (theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
		 + (theta.Wa_.block(theta.hidden_ * shift, offset_embedding, theta.hidden_, embedding_size)
		    * alpha_next.block(offset_embedding, 0, embedding_size, 1))
		 + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
		    * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	    
	    const double score = (theta.Wc_ * alpha_next.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	    
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
	    state_next.error() = state.error();
	    state_next.target() = target_next;
	    
	    matrix_type alpha_next(state_next.matrix(), state_size, 1);

	    copy_embedding(source, target, theta, next, trg, alpha_next);
	    
	    // compute alpha-next
	    alpha_next.block(offset_matrix, 0, theta.hidden_, 1)
	      = (theta.bn_
		 + (theta.Wn_.block(0, offset_embedding, theta.hidden_, embedding_size)
		    * alpha_next.block(offset_embedding, 0, embedding_size, 1))
		 + (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
		    * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	    
	    const double score = (theta.Wc_ * alpha_next.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	    
	    state_next.score() = score + state.score();
	    
	    heap_next.push_back(state_next);
	    std::push_heap(heap_next.begin(), heap_next.end(), heap_compare());
	  }
	}
      }
    }
    
    // perform pruning!

    heap_type& heap         = heaps_[target_size + 1];
    heap_type& heap_viterbi = heaps_viterbi_[target_size + 1];
    
    heap_type::iterator hiter_begin = heap.begin();
    heap_type::iterator hiter       = heap.end();
    heap_type::iterator hiter_end   = heap.end();

    heap_type::iterator viter_begin = heap_viterbi.begin();
    heap_type::iterator viter       = heap_viterbi.end();
    heap_type::iterator viter_end   = heap_viterbi.end();
    
    {
      for (/**/; hiter_begin != hiter && hiter_end - hiter != beam_; -- hiter)
	std::pop_heap(hiter_begin, hiter, heap_compare());
      
      // deallocate unused states
      for (/**/; hiter_begin != hiter; ++ hiter_begin)
	state_allocator_.deallocate(*hiter_begin);
    }
    
    {
      for (/**/; viter_begin != viter && viter_end - viter != beam_; -- viter)
	std::pop_heap(viter_begin, viter, heap_compare());
      
      // deallocate unused states
      for (/**/; viter_begin != viter; ++ viter_begin)
	state_allocator_.deallocate(*viter_begin);
    }
      
    // iterate hiter to hiter_end as an instance for error
    state_set_type& states = states_[target_size + 1];
    
    double loss = 0.0;

    size_type num_loss   = 0;
    size_type num_errors = 0;
    for (heap_type::iterator miter = hiter; miter != hiter_end; ++ miter) 
      if (miter->error() > 0) {
	for (heap_type::iterator citer = viter; citer != viter_end; ++ citer)
	  num_loss += double(miter->error()) - (citer->score() - miter->score()) > 0.0;
	
	++ num_errors;
      }
    
    if (num_loss) {
      const double error_factor = 1.0 / (num_errors * (viter_end - viter));
      
      for (heap_type::iterator miter = hiter; miter != hiter_end; ++ miter) 
	if (miter->error() > 0)
	  for (heap_type::iterator citer = viter; citer != viter_end; ++ citer) {
	    const double error = std::max(double(miter->error()) - (citer->score() - miter->score()), 0.0) * error_factor;
	    
	    if (error == 0.0) continue;
	    
	    state_set_type::iterator siter_c = states.find(*citer);
	    if (siter_c != states.end())
	      siter_c->second.loss() += - error_factor;
	    else {
	      state_type buffer = state_allocator_.allocate();
	      
	      buffer.loss() = - error_factor;
	      matrix_type(buffer.matrix(), state_size, 1).setZero();
	      
	      states[*citer] = buffer;
	    }
	    
	    state_set_type::iterator siter_m = states.find(*miter);
	    if (siter_m != states.end())
	      siter_m->second.loss() += error_factor;
	    else {
	      state_type buffer = state_allocator_.allocate();
	      
	      buffer.loss() = error_factor;
	      matrix_type(buffer.matrix(), state_size, 1).setZero();
	      
	      states[*miter] = buffer;
	    }
	    
	    loss += error;
	  }
    }
    
    //std::cerr << "# of pairs: " << pairs << " loss: " << loss << std::endl;
    
    state_type state = heaps_viterbi_[target_size + 1].back();
    
    for (size_type trg = target_size + 1; trg; -- trg) {
      const size_type src = state.index();
      
      if (trg < target_size + 1 && src < source_size + 1)
	alignment.push_back(std::make_pair(src - 1, trg - 1));
      
      state = state.prev();
    }

    return loss;
  }
  
  template <typename Gen>
  double backward(const sentence_type& source,
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

    const size_type embedding_window_size = theta.embedding_ * (theta.window_ * 2 + 1);
    const size_type state_size = embedding_window_size * 2 + theta.hidden_;
    const size_type embedding_size = embedding_window_size * 2;
    
    const size_type offset_embedding = 0;
    const size_type offset_source = 0;
    const size_type offset_target = embedding_window_size;
    const size_type offset_matrix = embedding_window_size * 2;

    // do not accumulate!
    if (states_[target_size + 1].empty()) return 0.0;
    
    double loss = 0.0;
        
    ++ gradient.count_;
    
    for (size_type trg = target_size + 1; trg > 0; -- trg) {
      const state_set_type& states_next = states_[trg];
      state_set_type& states_prev = states_[trg - 1];
      
      state_set_type::const_iterator siter_end = states_next.end();
      for (state_set_type::const_iterator siter = states_next.begin(); siter != siter_end; ++ siter) {
	const state_type state_next(siter->first);
	const state_type state_prev(state_next.prev());
	
	const word_type target_next = state_next.target();
	const word_type target_prev = state_prev.target();
	
	const matrix_type alpha_next(state_next.matrix(), state_size, 1);
	const matrix_type alpha_prev(state_prev.matrix(), state_size, 1);
	
	state_set_type::iterator piter = states_prev.find(state_prev);
	if (piter == states_prev.end()) {
	  state_type buffer = state_allocator_.allocate();
	  
	  buffer.loss() = 0.0;
	  matrix_type(buffer.matrix(), state_size, 1).setZero();
	  
	  piter = states_prev.insert(std::make_pair(state_prev, buffer)).first;
	}
	
	const matrix_type beta_next(siter->second.matrix(), state_size, 1);
	/**/  matrix_type beta_prev(piter->second.matrix(), state_size, 1);

	const parameter_type& loss_next = siter->second.loss();
	/**/  parameter_type& loss_prev = piter->second.loss();
	
	loss_prev += loss_next;
	
	const size_type next = state_next.index();
	const size_type prev = state_prev.index();

	gradient.Wc_         += loss_next * alpha_next.block(offset_matrix, 0, theta.hidden_, 1).transpose();
	gradient.bc_.array() += loss_next;
	
	delta_alpha_ = (alpha_next.block(offset_matrix, 0, theta.hidden_, 1).array().unaryExpr(dhinge())
			* (theta.Wc_.transpose() * loss_next
			   + beta_next.block(offset_matrix, 0, theta.hidden_, 1)).array());
	
	// update Wa or Wn and propagate back to beta_prev...
	if (next >= source_size + 2) {
	  gradient.Wn_.block(0, offset_embedding, theta.hidden_, embedding_size)
	    += delta_alpha_ * alpha_next.block(offset_embedding, 0, embedding_size, 1).transpose();
	  gradient.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
	    += delta_alpha_ * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1).transpose();
	  gradient.bn_
	    += delta_alpha_;
	  
	  propagate_embedding(source, target, theta, gradient, next, trg, target_next,
			      theta.Wn_.block(0, offset_embedding, theta.hidden_, embedding_size).transpose()
			      * delta_alpha_);
	  
	  beta_prev.block(offset_matrix, 0, theta.hidden_, 1)
	    += (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_).transpose()
		* delta_alpha_);
	} else {
	  const size_type shift = theta.shift(source_size, target_size, prev, next); 
	  
	  gradient.Wa_.block(theta.hidden_ * shift, offset_embedding, theta.hidden_, embedding_size)
	    += delta_alpha_ * alpha_next.block(offset_embedding, 0, embedding_size, 1).transpose();
	  gradient.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
	    += delta_alpha_ * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1).transpose();
	  gradient.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
	    += delta_alpha_;
	  
	  propagate_embedding(source, target, theta, gradient, next, trg, target_next,
			      theta.Wa_.block(theta.hidden_ * shift, offset_embedding, theta.hidden_, embedding_size).transpose()
			      * delta_alpha_);
	  
	  beta_prev.block(offset_matrix, 0, theta.hidden_, 1)
	    += (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_).transpose()
		* delta_alpha_);
	}
      }
    }
    
    // final propagation...
    if (states_[0].empty()) return 0.0;
    
    if (states_[0].size() != 1)
      throw std::runtime_error("multiple initial states?");
    
    const matrix_type beta_bos(states_[0].begin()->second.matrix(), state_size, 1);
    
    gradient.bi_ += beta_bos.block(offset_matrix, 0, theta.hidden_, 1);
    
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
    
    const size_type embedding_window_size = theta.embedding_ * (theta.window_ * 2 + 1);
    const size_type state_size = embedding_window_size * 2 + theta.hidden_;
    const size_type embedding_size = embedding_window_size * 2;

    const size_type offset_embedding = 0;
    const size_type offset_source = 0;
    const size_type offset_target = embedding_window_size;
    const size_type offset_matrix = embedding_window_size * 2;
    
    state_allocator_.assign(state_type::size(state_size));
    
    alignment.clear();
    
    heaps_.clear();
    heaps_.resize(target_size + 2);
    
    state_type state_bos = state_allocator_.allocate();
    
    state_bos.prev() = state_type();
    state_bos.index() = 0;
    state_bos.score() = 0.0;
    
    matrix_type alpha_bos(state_bos.matrix(), state_size, 1);
    
    alpha_bos.block(offset_matrix, 0, theta.hidden_, 1) = theta.bi_;
    
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

	const matrix_type alpha_prev(state.matrix(), state_size, 1);
	
	const size_type prev = state.index();

	// alignment into aligned...
	for (size_type next = next_first; next != next_last; ++ next) {
	  const size_type shift = theta.shift(source_size, target_size, prev, next);

	  state_type state_next = state_allocator_.allocate();
	  state_next.prev() = state;
	  state_next.index() = next;
	  
	  matrix_type alpha_next(state_next.matrix(), state_size, 1);
	  
	  copy_embedding(source, target, theta, next, trg, alpha_next);
	  
	  // compute alpha-next
	  alpha_next.block(offset_matrix, 0, theta.hidden_, 1)
	    = (theta.ba_.block(theta.hidden_ * shift, 0, theta.hidden_, 1)
	       + (theta.Wa_.block(theta.hidden_ * shift, offset_embedding, theta.hidden_, embedding_size)
		  * alpha_next.block(offset_embedding, 0, embedding_size, 1))
	       + (theta.Wa_.block(theta.hidden_ * shift, offset_matrix, theta.hidden_, theta.hidden_)
		  * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	  
	  const double score = (theta.Wc_ * alpha_next.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	  
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
	  
	  matrix_type alpha_next(state_next.matrix(), state_size, 1);
	  
	  copy_embedding(source, target, theta, next, trg, alpha_next);
	  
	  // compute alpha-next
	  alpha_next.block(offset_matrix, 0, theta.hidden_, 1)
	    = (theta.bn_
	       + (theta.Wn_.block(0, offset_embedding, theta.hidden_, embedding_size)
		  * alpha_next.block(offset_embedding, 0, embedding_size, 1))
	       + (theta.Wn_.block(0, offset_matrix, theta.hidden_, theta.hidden_)
		  * alpha_prev.block(offset_matrix, 0, theta.hidden_, 1))).array().unaryExpr(hinge());
	  
	  const double score = (theta.Wc_ * alpha_next.block(offset_matrix, 0, theta.hidden_, 1) + theta.bc_)(0, 0);
	  
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
    
    for (size_type trg = target_size + 1; trg; -- trg) {
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
	       const size_type& window,
	       const size_type& alignment,
	       const double& lambda,
	       const double& lambda2,
	       const double& eta0)
    : embedding_(embedding),
      hidden_(hidden),
      window_(window),
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
    
    source_ = tensor_type::Zero(embedding_, vocabulary_size);
    target_ = tensor_type::Zero(embedding_, vocabulary_size);

    const size_type state_size = embedding_ * (window * 2 + 1) * 2 + hidden_;
    
    Wc_ = tensor_type::Zero(1, hidden_);
    bc_ = tensor_type::Zero(1, 1);
    
    //Wt_ = tensor_type::Zero(hidden_, state_size);
    //bt_ = tensor_type::Zero(hidden_, 1);

    Wa_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), state_size);
    ba_ = tensor_type::Zero(hidden_ * (alignment * 2 + 1), 1);
    
    Wn_ = tensor_type::Zero(hidden_, state_size);
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
	     scale);
    
    gradient_embedding_type::const_iterator titer_end = gradient.target_.end();
    for (gradient_embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer)
      update(titer->first,
	     theta.target_,
	     const_cast<tensor_type&>(target_),
	     titer->second,
	     embedding.source_.col(titer->first.id()),
	     scale);

    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, scale, false);
    
    //update(theta.Wt_, const_cast<tensor_type&>(Wt_), gradient.Wt_, scale, lambda_ != 0.0);
    //update(theta.bt_, const_cast<tensor_type&>(bt_), gradient.bt_, scale, false);
    
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
  size_type alignment_;
  
  double lambda_;
  double lambda2_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;

  tensor_type Wc_;
  tensor_type bc_;
  
  //tensor_type Wt_;
  //tensor_type bt_;
  
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

    update(theta.Wc_, gradient.Wc_, scale, lambda_ != 0.0);
    update(theta.bc_, gradient.bc_, scale, false);
    
    //update(theta.Wt_, gradient.Wt_, scale, lambda_ != 0.0);
    //update(theta.bt_, gradient.bt_, scale, false);
    
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

#endif
