//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_RNN_NGRAM_IMPL__HPP__
#define __CICADA_RNN_NGRAM_IMPL__HPP__ 1

#include <cstdlib>
#include <cmath>
#include <climits>

#include <vector>
#include <string>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/mathop.hpp"
#include "utils/unordered_map.hpp"
#include "utils/repository.hpp"
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

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef utils::unordered_map<word_type, tensor_type,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, tensor_type> > >::type embedding_type;
  
  Gradient() : dimension_(0), order_(0), count_(0), shared_(0) {}
  Gradient(const size_type& dimension,
	   const int order) 
    : dimension_(dimension),
      order_(order),
      count_(0),
      shared_(0)
  { initialize(dimension, order); }
  
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
    
    if (! bi_.rows())
      bi_ = tensor_type::Zero(x.bi_.rows(), x.bi_.cols());
    
    Wc_ -= x.Wc_;
    bc_ -= x.bc_;
    
    bi_ -= x.bi_;

    count_ -= x.count_;
    
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

    if (! bi_.rows())
      bi_ = tensor_type::Zero(x.bi_.rows(), x.bi_.cols());
    
    Wc_ += x.Wc_;
    bc_ += x.bc_;
    
    bi_ += x.bi_;

    count_ += x.count_;

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
    
    bi_.setZero();
    
    count_  = 0;
    shared_ = 0;
  }
  
  tensor_type& embedding_input(const word_type& word)
  {
    tensor_type& embedding = embedding_input_[word];
    if (! embedding.cols())
      embedding = tensor_type::Zero(dimension_, 1);
    
    return embedding;
  }

  tensor_type& embedding_output(const word_type& word)
  {
    tensor_type& embedding = embedding_output_[word];
    if (! embedding.cols())
      embedding = tensor_type::Zero(dimension_ + 1, 1);
    
    return embedding;
  }
  
  void initialize(const size_type dimension, const int order)
  {
    if (dimension <= 0)
      throw std::runtime_error("invalid dimension");
    if (order <= 0)
      throw std::runtime_error("invalid order");
    
    dimension_ = dimension;
    order_     = order;
    
    clear();
    
    // initialize...
    Wc_ = tensor_type::Zero(dimension_, dimension_ * 2 * (order - 1));
    bc_ = tensor_type::Zero(dimension_, order - 1).array();
    
    bi_ = tensor_type::Zero(dimension_, 1);

    count_  = 0;
    shared_ = 0;
  }

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
    os.write((char*) &dimension_, sizeof(size_type));
    os.write((char*) &order_, sizeof(size_type));
    os.write((char*) &count_, sizeof(size_type));
    
    write(os, Wc_);
    write(os, bc_);
    write(os, bi_);
    
    write(os, embedding_input_,  false);
    write(os, embedding_output_, true);
  }

  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &dimension_, sizeof(size_type));
    is.read((char*) &order_, sizeof(size_type));
    is.read((char*) &count_, sizeof(size_type));
    
    read(is, Wc_);
    read(is, bc_);
    read(is, bi_);
    
    read(is, embedding_input_,  false);
    read(is, embedding_output_, true);

    // checking...
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
      
      matrix.resize(dimension_ + bias_last, 1);
      
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
  size_type dimension_;
  size_type order_;
  
  // embedding
  embedding_type embedding_input_;
  embedding_type embedding_output_;
  
  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;
  
  // bi for initial context
  tensor_type bi_;
  
  // other variables
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

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  typedef boost::filesystem::path path_type;

  typedef std::vector<bool, std::allocator<bool> > unique_set_type;
  typedef std::vector<word_type, std::allocator<word_type> > word_set_type;
  
  Model() : dimension_(0), order_(0), scale_(1) {}
  template <typename Unigram, typename Gen>
  Model(const size_type& dimension,
	const int order,
	const Unigram& unigram,
	Gen& gen) 
    : dimension_(dimension),
      order_(order),
      scale_(1)
  { initialize(dimension, order, unigram, gen); }
  
  
  void clear()
  {
    // embedding
    embedding_input_.setZero();
    embedding_output_.setZero();
    
    // matrix for context
    Wc_.setZero();
    bc_.setZero();
    
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
  
  template <typename Unigram, typename Gen>
  void initialize(const size_type dimension,
		  const int order,
		  const Unigram& unigram,
		  Gen& gen)
  {
    if (dimension <= 0)
      throw std::runtime_error("invalid dimension");
    if (order <= 0)
      throw std::runtime_error("invalid order");

    dimension_ = dimension;
    order_     = order;
    
    clear();
    
    const size_type vocabulary_size = word_type::allocated();

    const double range_e = std::sqrt(6.0 / (dimension_ + 1));
    const double range_c = std::sqrt(6.0 / (dimension_ + dimension_ + dimension_));
    const double range_i = std::sqrt(6.0 / (dimension_ + 1));
        
    embedding_input_  = tensor_type::Zero(dimension_,     vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    embedding_output_ = tensor_type::Zero(dimension_ + 1, vocabulary_size).array().unaryExpr(randomize<Gen>(gen, range_e));
    embedding_output_.row(dimension_).setZero();
    
    Wc_ = tensor_type::Zero(dimension_, dimension_ * 2 * (order - 1)).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Zero(dimension_, order - 1);
    
    bi_ = tensor_type::Zero(dimension_, 1).array().unaryExpr(randomize<Gen>(gen, range_i));

    uniques_ = unique_set_type(vocabulary_size, false);
    
    // assign id... which is used to compute the final word embeddings...
    words_.clear();
    words_.push_back(vocab_type::BOS);
    words_.push_back(vocab_type::EOS);
    words_.push_back(vocab_type::EPSILON);
    words_.push_back(vocab_type::UNK);
    
    for (size_type pos = 0; pos != unigram.words_.size(); ++ pos)
      if (unigram.words_[pos] != vocab_type::EOS) {
	uniques_[unigram.words_[pos].id()] = true;

	if (unigram.words_[pos] != vocab_type::BOS
	    && unigram.words_[pos] != vocab_type::EPSILON
	    && unigram.words_[pos] != vocab_type::UNK)
	  words_.push_back(unigram.words_[pos]);
      }
    
    word_set_type(words_).swap(words_);
    
    uniques_[vocab_type::BOS.id()] = false;
    uniques_[vocab_type::EOS.id()] = false;

    scale_ = 1.0;
  }

  void finalize()
  {
    // clear unused entries
    embedding_input_.col(vocab_type::EOS.id())  = tensor_type::Zero(dimension_, 1);
    embedding_output_.col(vocab_type::BOS.id()) = tensor_type::Zero(dimension_ + 1, 1);
    embedding_output_.col(vocab_type::EPSILON.id()) = tensor_type::Zero(dimension_ + 1, 1);
    
    const double factor = 1.0 / std::accumulate(uniques_.begin(), uniques_.end(), size_type(0));
    
    tensor_type average_input  = tensor_type::Zero(dimension_, 1);
    tensor_type average_output = tensor_type::Zero(dimension_ + 1, 1);
    
    bool has_unk = false;
    
    for (size_type pos = 0; pos != uniques_.size(); ++ pos) 
      if (uniques_[pos]) {
	average_input  += embedding_input_.col(pos) * factor;
	average_output += embedding_output_.col(pos) * factor;
	
	has_unk |= (word_type(pos) == vocab_type::UNK);
      }
    
    if (! has_unk) {
      embedding_input_.col(vocab_type::UNK.id())  = average_input;
      embedding_output_.col(vocab_type::UNK.id()) = average_output;
    }
    
    if (scale_ != 1.0) {
      embedding_input_.array() *= scale_;
      embedding_output_.block(0, 0, dimension_, embedding_output_.cols()).array() *= scale_;
      
      scale_ = 1.0;
    }
  }

  void rescale()
  {
    if (scale_ == 1.0) return;
    
    embedding_input_.array() *= scale_;
    embedding_output_.block(0, 0, dimension_, embedding_output_.cols()).array() *= scale_;
    
    scale_ = 1.0;
  }
  
  struct real_policy : boost::spirit::karma::real_policies<parameter_type>
  {
    static unsigned int precision(parameter_type)
    {
      return 10;
    }
  };
  
  
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["embedding"] = utils::lexical_cast<std::string>(dimension_);
    rep["order"]     = utils::lexical_cast<std::string>(order_);
    rep["size"]      = utils::lexical_cast<std::string>(words_.size());
    rep["scale"]     = utils::lexical_cast<std::string>(scale_);
    
    write_embedding(rep.path("input.gz"),  rep.path("input.bin"), embedding_input_);
    write_embedding(rep.path("output.gz"), rep.path("output.bin"), embedding_output_);
    
    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
    
    write(rep.path("bi.txt.gz"), rep.path("bi.bin"), bi_);
    
    // vocabulary...
    vocab_type vocab;
    
    vocab.open(rep.path("vocab"), words_.size() >> 1);

    word_set_type::const_iterator witer_end = words_.end();
    for (word_set_type::const_iterator witer = words_.begin(); witer != witer_end; ++ witer)
      vocab.insert(*witer);
    
    vocab.close();
  }
  
  void write_embedding(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix) const
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    const tensor_type::Index rows = matrix.rows();
    const tensor_type::Index cols = matrix.cols();

    karma::real_generator<double, real_policy> float10;
    
    utils::compress_ostream os_txt(path_text, 1024 * 1024);
    utils::compress_ostream os_bin(path_binary, 1024 * 1024);
    std::ostream_iterator<char> iter(os_txt);
    
    word_set_type::const_iterator witer_end = words_.end();
    for (word_set_type::const_iterator witer = words_.begin(); witer != witer_end; ++ witer) {
      karma::generate(iter, standard::string, *witer);
      
      for (difference_type j = 0; j != rows; ++ j)
	karma::generate(iter, karma::lit(' ') << float10, matrix(j, witer->id()));
      
      karma::generate(iter, karma::lit('\n'));

      os_bin.write((char*) matrix.col(witer->id()).data(), sizeof(tensor_type::Scalar) * rows);
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

  Model& operator+=(const Model& x)
  {
    if (scale_ != x.scale_)
      throw std::runtime_error("different scaling");

    embedding_input_  += x.embedding_input_;
    embedding_output_ += x.embedding_output_;
    
    Wc_ += x.Wc_;
    bc_ += x.bc_;
    
    bi_ += x.bi_;

    return *this;
  }

  Model& operator*=(const double& x)
  {
    embedding_input_  *= x;
    embedding_output_ *= x;
    
    Wc_ *= x;
    bc_ *= x;
    
    bi_ *= x;
    
    return *this;
  }

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
    os.write((char*) &dimension_, sizeof(size_type));
    os.write((char*) &order_, sizeof(size_type));
    os.write((char*) &scale_, sizeof(double));
    
    write(os, Wc_);
    write(os, bc_);
    write(os, bi_);
    
    write(os, embedding_input_,  false);
    write(os, embedding_output_, true);
  }

  void read(std::istream& is)
  {
    clear();
    
    is.read((char*) &dimension_, sizeof(size_type));
    is.read((char*) &order_, sizeof(size_type));
    is.read((char*) &scale_, sizeof(double));
    
    read(is, Wc_);
    read(is, bc_);
    read(is, bi_);
    
    read(is, embedding_input_,  false);
    read(is, embedding_output_, true);

    // checking...
  }

  void write(std::ostream& os, const tensor_type& embedding, const bool bias_last) const
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

  void read(std::istream& is, tensor_type& embedding, const bool bias_last)
  {
    buffer_type& buffer = const_cast<buffer_type&>(buffer_);
    
    size_type cols = 0;
    is.read((char*) &cols, sizeof(size_type));
    
    if (cols > embedding.cols())
      embedding.conservativeResize(dimension_ + bias_last, cols);
    
    for (size_type i = 0; i != cols; ++ i) {
      size_type word_size = 0;
      is.read((char*) &word_size, sizeof(size_type));
      
      buffer.resize(word_size);
      is.read((char*) &(*buffer.begin()), word_size);
      
      const word_type word(buffer.begin(), buffer.end());
      
      if (word.id() >= embedding.cols())
	embedding.conservativeResize(dimension_ + bias_last, word.id() + 1);
      
      is.read((char*) embedding.col(word.id()).data(), sizeof(tensor_type::Scalar) * (dimension_ + bias_last));
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
  size_type dimension_;
  size_type order_;
  
  // embedding
  tensor_type embedding_input_;
  tensor_type embedding_output_;

  unique_set_type uniques_;
  word_set_type   words_;
  
  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;
  
  // bi for initial context
  tensor_type bi_;

  // scale for embedding...
  double scale_;
};

struct Unigram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef uint64_t count_type;

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

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
    logprobs_.resize(word_type::allocated(), - std::numeric_limits<double>::infinity());
    
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

struct NGram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;
  typedef Unigram  unigram_type;

  typedef Average log_likelihood_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type    tensor_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
  NGram(const unigram_type& unigram,
	const size_type samples)
    : unigram_(unigram), samples_(samples), log_samples_(std::log(double(samples))) {}

  const unigram_type& unigram_;
  size_type           samples_;
  double              log_samples_;
  
  tensor_type lattice_;
  tensor_type delta_;

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

  template <typename Iterator, typename Gen>
  double learn(Iterator first, Iterator last,
	       const model_type& theta,
	       gradient_type& gradient,
	       Gen& gen)
  {
    // forward path to construct lattice
    forward(first, last, theta);
    
    // backward path to accumlate gradient
    return backward(first, last, theta, gradient, gen);
  }
  
  template <typename Iterator>
  void forward(Iterator first, Iterator last,
	       const model_type& theta)
  {
    const size_type dimension = theta.dimension_;
    const size_type order     = theta.order_;
    
    const size_type offset_embedding = 0;
    const size_type offset_context   = dimension;
    
    if (! lattice_.rows())
      lattice_.resize(dimension, order);
    
    lattice_.col(0) = theta.bi_.array().unaryExpr(hinge());
    
    for (size_type i = 1; i != order; ++ i, ++ first) {
      const size_type shift = (i - 1) * 2 * dimension;
      
      lattice_.col(i) = ((theta.Wc_.block(0, shift + offset_embedding, dimension, dimension)
			  * theta.embedding_input_.col(first->id())
			  * theta.scale_)
			 + (theta.Wc_.block(0, shift + offset_context, dimension, dimension)
			    * lattice_.col(i - 1))
			 + theta.bc_.block(0, i - 1, dimension, 1)).array().unaryExpr(hinge());
    }
  }

  template <typename Iterator, typename Gen>
  double backward(Iterator first, Iterator last,
		  const model_type& theta,
		  gradient_type& gradient,
		  Gen& gen)
  {
    const size_type dimension = theta.dimension_;
    const size_type order     = theta.order_;
    
    const size_type offset_embedding = 0;
    const size_type offset_context   = dimension;
    
    const word_type word(*(last - 1));
    
    const double score = ((theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1).transpose()
			   * lattice_.col(order - 1)
			   * theta.scale_)
			  + theta.embedding_output_.col(word.id()).block(dimension, 0, 1, 1))(0, 0);
    const double score_noise = log_samples_ + unigram_.logprob(word);
    const double z = utils::mathop::logsum(score, score_noise);
    const double logprob = score - z;
    const double logprob_noise = score_noise - z;
    
    const double loss = - 1.0 + std::exp(logprob);
    
    double log_likelihood = logprob;
    
    tensor_type& dembedding = gradient.embedding_output(word);
    
    dembedding.block(0, 0, dimension, 1).noalias() += loss * lattice_.col(order - 1);
    dembedding.block(dimension, 0, 1, 1).array()   += loss;
    
    delta_ = (lattice_.col(order - 1).array().unaryExpr(dhinge())
	      * (theta.embedding_output_.col(word.id()).block(0, 0, dimension, 1) * loss * theta.scale_).array());
    
    // perform sampling here...
    for (size_type k = 0; k != samples_; ++ k) {
      word_type sampled = unigram_.draw(gen);
      
      while (sampled == word)
	sampled = unigram_.draw(gen);
      
      const double score = ((theta.embedding_output_.col(sampled.id()).block(0, 0, dimension, 1).transpose()
			   * lattice_.col(order - 1)
			     * theta.scale_)
			    + theta.embedding_output_.col(sampled.id()).block(dimension, 0, 1, 1))(0, 0);
      const double score_noise = log_samples_ + unigram_.logprob(sampled);
      const double z = utils::mathop::logsum(score, score_noise);
      const double logprob = score - z;
      const double logprob_noise = score_noise - z;
      
      log_likelihood += logprob_noise;
      
      const double loss = std::exp(logprob);
      
      tensor_type& dembedding = gradient.embedding_output(sampled);
      
      dembedding.block(0, 0, dimension, 1).noalias() += loss * lattice_.col(order - 1);
      dembedding.block(dimension, 0, 1, 1).array()   += loss;
      
      delta_.array() += (lattice_.col(order - 1).array().unaryExpr(dhinge())
			 * (theta.embedding_output_.col(sampled.id()).block(0, 0, dimension, 1) * loss * theta.scale_).array());
    }
    
    for (size_type i = order - 1; i != 0; -- i, -- last) {
      const size_type shift = (i - 1) * 2 * dimension;
      const word_type prev(*(last - 2));
      
      gradient.Wc_.block(0, shift + offset_embedding, dimension, dimension)
	+= delta_ * theta.embedding_input_.col(prev.id()).transpose() * theta.scale_;
      gradient.Wc_.block(0, shift + offset_context, dimension, dimension)
	+= delta_ * lattice_.col(i - 1).transpose();
      gradient.bc_.block(0, i - 1, dimension, 1)
	+= delta_;
      
      gradient.embedding_input(prev)
	+= theta.Wc_.block(0, shift + offset_embedding, dimension, dimension).transpose() * delta_;
      
      // propagate this...
      delta_ = (lattice_.col(i - 1).array().unaryExpr(dhinge())
		* (theta.Wc_.block(0, shift + offset_context, dimension, dimension).transpose() * delta_).array());
    }
    
    ++ gradient.count_;
    
    gradient.bi_ += delta_;
    
    return log_likelihood;
  }
};

struct Learn
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  typedef Gradient gradient_type;

  typedef cicada::Symbol   word_type;
  
  typedef model_type::tensor_type tensor_type;
};

struct LearnAdaGrad : public Learn
{
  LearnAdaGrad(const size_type& dimension,
	       const int order,
	       const double& lambda,
	       const double& eta0)
    : dimension_(dimension),
      order_(order),
      lambda_(lambda),
      eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    
    embedding_input_  = tensor_type::Zero(dimension_,     vocabulary_size);
    embedding_output_ = tensor_type::Zero(dimension_ + 1, vocabulary_size);
    
    // initialize...
    Wc_ = tensor_type::Zero(dimension_, dimension_ * 2 * (order - 1));
    bc_ = tensor_type::Zero(dimension_, order - 1).array();
    
    bi_ = tensor_type::Zero(dimension_, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    if (! gradient.count_) return;
    
    embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
    for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
      update(iiter->first,
	     theta.embedding_input_,
	     const_cast<tensor_type&>(embedding_input_),
	     iiter->second,
	     1.0 / gradient.count_,
	     lambda_ != 0.0,
	     false);

    embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
    for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
      update(oiter->first,
	     theta.embedding_output_,
	     const_cast<tensor_type&>(embedding_output_),
	     oiter->second,
	     1.0 / gradient.count_,
	     lambda_ != 0.0,
	     true);
    
    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, 1.0 / gradient.count_, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, 1.0 / gradient.count_, false);

    update(theta.bi_, const_cast<tensor_type&>(bi_), gradient.bi_, 1.0 / gradient.count_, false);
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
	      const bool regularize=true,
	      const bool bias_last=false) const
  {
    if (regularize) {
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
    } else {
      G.col(word.id()).array() += g.array().square() * scale * scale;
      theta.col(word.id()).array() -= eta0_ * scale * g.array() * G.col(word.id()).array().unaryExpr(learning_rate());
    }
  }
  
  size_type dimension_;
  size_type order_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type embedding_input_;
  tensor_type embedding_output_;

  // Wc and bc for context layer
  tensor_type Wc_;
  tensor_type bc_;

  // bi for initial context
  tensor_type bi_;  
};

struct LearnSGD : public Learn
{  
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
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef gradient_type::embedding_type embedding_type;

    if (! gradient.count_) return;
    
    //++ const_cast<size_type&>(epoch_);
    
    const double eta = eta0_ / (epoch_ + 1);
    
    if (lambda_ != 0.0)
      theta.scale_ *= 1.0 - eta * lambda_;
    
    embedding_type::const_iterator iiter_end = gradient.embedding_input_.end();
    for (embedding_type::const_iterator iiter = gradient.embedding_input_.begin(); iiter != iiter_end; ++ iiter)
      update(iiter->first,
	     theta.embedding_input_,
	     iiter->second,
	     1.0 / gradient.count_,
	     theta.scale_,
	     false);    
    
    embedding_type::const_iterator oiter_end = gradient.embedding_output_.end();
    for (embedding_type::const_iterator oiter = gradient.embedding_output_.begin(); oiter != oiter_end; ++ oiter)
      update(oiter->first,
	     theta.embedding_output_,
	     oiter->second,
	     1.0 / gradient.count_,
	     theta.scale_,
	     true);
    
    update(theta.Wc_, gradient.Wc_, 1.0 / gradient.count_, lambda_ != 0.0);
    update(theta.bc_, gradient.bc_, 1.0 / gradient.count_, false);
    
    update(theta.bi_, gradient.bi_, 1.0 / gradient.count_, false);
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
    
    theta -= eta * scale * g;
  }
  
  template <typename Theta, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const double theta_scale,
	      const bool bias_last=false) const
  {
    const double eta = eta0_ / (epoch_ + 1);

    if (bias_last) {
      const size_type rows = g.rows();
      
      theta.col(word.id()).block(0, 0, rows - 1, 1) -= (eta * scale /  theta_scale) * g.block(0, 0, rows - 1, 1);
      theta.col(word.id()).block(rows - 1, 0, 1, 1) -= eta * scale * g.block(rows - 1, 0, 1, 1);
    } else
      theta.col(word.id()) -= (eta * scale /  theta_scale) * g;
  }

  double lambda_;
  double eta0_;
  
  size_type epoch_;
};

struct Data
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef std::vector<word_type, std::allocator<word_type> > data_type;
  
  typedef data_type::const_iterator const_iterator;
  typedef data_type::iterator       iterator;

  Data(const int order) : order_(order) {}
  
  Data& operator+=(const Data& x)
  {
    if (order_ != x.order_)
      throw std::runtime_error("different order?");

    data_.insert(data_.end(), x.data_.begin(), x.data_.end());
    
    return *this;
  };

  void insert(const sentence_type& sent)
  {
    buffer_.resize(order_);
    std::fill(buffer_.begin(), buffer_.end() - 1, vocab_type::EPSILON);
    buffer_.back() = vocab_type::BOS;
    
    sentence_type::const_iterator siter_end = sent.end();
    for (sentence_type::const_iterator siter = sent.begin(); siter != siter_end; ++ siter) {
      buffer_.push_back(*siter);
      data_.insert(data_.end(), buffer_.end() - order_, buffer_.end());
    }
    
    buffer_.push_back(vocab_type::EOS);
    data_.insert(data_.end(), buffer_.end() - order_, buffer_.end());
  }
  
  void reserve(size_type x) { data_.reserve(x * order_); }
  
  void clear()
  {
    buffer_.clear();
    data_.clear();
  };

  void swap(Data& x)
  {
    data_.swap(x.data_);
    std::swap(order_, x.order_);
  }

  void shrink()
  {
    data_type(data_).swap(data_);
  }
  
  size_type size() const { return data_.size() / order_; }
  bool empty() const { return data_.empty(); }
  bool verify() const { return data_.size() % order_ == 0; }
  
  inline const_iterator begin(size_type pos) const { return data_.begin() + pos * order_; }
  inline       iterator begin(size_type pos) { return data_.begin() + pos * order_; }
  inline const_iterator end(size_type pos) const { return data_.begin() + (pos + 1) * order_; }
  inline       iterator end(size_type pos) { return data_.begin() + (pos + 1) * order_; }

  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }

  void erase(iterator first, iterator last)
  {
    data_.erase(first, last);
  }
  
  data_type buffer_;
  data_type data_;
  int       order_;
};

#endif
