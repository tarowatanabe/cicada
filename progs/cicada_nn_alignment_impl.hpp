//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_NN_ALIGNMENT_IMPL__HPP__
#define __CICADA_NN_ALIGNMENT_IMPL__HPP__ 1

#include <cstdlib>
#include <cmath>
#include <climits>

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/fusion/tuple.hpp>
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
#include "utils/mathop.hpp"
#include "utils/unordered_map.hpp"
#include "utils/repository.hpp"
#include "utils/random_seed.hpp"
#include "utils/compress_stream.hpp"

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
    os.write((char*) &count_,     sizeof(size_type));
    
    write(os, Wt_);
    write(os, bt_);

    write(os, Wc_);
    write(os, bc_);
    
    write(os, source_);
    write(os, target_);
  }
  
  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &window_,    sizeof(size_type));
    is.read((char*) &count_,     sizeof(size_type));
    
    read(is, Wt_);
    read(is, bt_);

    read(is, Wc_);
    read(is, bc_);
    
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

  Model& operator+=(const Model& x)
  {
    if (scale_ != x.scale_)
      throw std::runtime_error("different scaling");

    source_ += x.source_;
    target_ += x.target_;
    
    Wt_ += x.Wt_;
    bt_ += x.bt_;
    
    Wc_ += x.Wc_;
    bc_ += x.bc_;
    
    return *this;
  }

  Model& operator*=(const double& x)
  {
    source_ *= x;
    target_ *= x;
    
    Wt_ *= x;
    bt_ *= x;

    Wc_ *= x;
    bc_ *= x;
    
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
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);
    rep["window"]    = utils::lexical_cast<std::string>(window_);
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("source.gz"), rep.path("source.bin"), rep.path("vocab-source"), source_, words_source_);
    write_embedding(rep.path("target.gz"), rep.path("target.bin"), rep.path("vocab-target"), target_, words_target_);
    
    write(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write(rep.path("bt.txt.gz"), rep.path("bt.bin"), bt_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
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
    os.write((char*) &scale_,     sizeof(double));
    
    write(os, Wt_);
    write(os, bt_);

    write(os, Wc_);
    write(os, bc_);
    
    write_embedding(os, source_);
    write_embedding(os, target_);
  }

  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &embedding_, sizeof(size_type));
    is.read((char*) &hidden_,    sizeof(size_type));
    is.read((char*) &window_,    sizeof(size_type));
    is.read((char*) &scale_,     sizeof(double));

    read(is, Wt_);
    read(is, bt_);

    read(is, Wc_);
    read(is, bc_);
    
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
      
      // we treat the loss uniformly, which is not correct (we need to weight by the expectation of alignment scores)
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

    if (! gradient.count_) return;
    
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

    if (! gradient.count_) return;

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

#endif
