// -*- mode: c++ -*-
//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include "tree_rnn.hpp"

#include "utils/map_file.hpp"
#include "utils/unordered_map.hpp"
#include "utils/indexed_set.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/repository.hpp"
#include "utils/spinlock.hpp"

namespace cicada
{
  
  template <typename Value>
  inline
  Value repository_value(const utils::repository& rep, const std::string& key)
  {
    utils::repository::const_iterator iter = rep.find(key);
    if (iter == rep.end())
      throw std::runtime_error("no " + key + "?");
    return utils::lexical_cast<Value>(iter->second);
  }

  class EmbeddingMapped : public TreeRNN::Embedding
  {
  private:
    typedef TreeRNN tree_rnn_type;
    
  public:
    typedef tree_rnn_type::size_type       size_type;
    typedef tree_rnn_type::difference_type difference_type;
    
    typedef tree_rnn_type::word_type  word_type;
    typedef tree_rnn_type::vocab_type vocab_type;

    typedef tree_rnn_type::path_type path_type;
    
    typedef tree_rnn_type::parameter_type parameter_type;
    typedef tree_rnn_type::tensor_type    tensor_type;
    typedef tree_rnn_type::matrix_type    matrix_type;
    typedef tree_rnn_type::embedding_type embedding_type;

  private:
    typedef utils::map_file<parameter_type, std::allocator<parameter_type> > mapped_type;
    
  public:
    EmbeddingMapped() : rows_(0) {}
    EmbeddingMapped(const path_type& path, const size_type rows)
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::read);
      
      rows_ = repository_value<size_type>(rep, "embedding");
      
      mapped_.open(rep.path("input.bin"));
      vocab_.open(rep.path("vocab"));
    }
    
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      
      rep["embedding"] = utils::lexical_cast<std::string>(rows_);
      
      mapped_.write(rep.path("input.bin"));
      vocab_.write(rep.path("vocab"));
    }

    matrix_type operator()(const word_type& word) const
    {
      return matrix_type(const_cast<parameter_type*>(mapped_.begin() + rows_ * vocab_[word]), rows_, 1);
    }

  public:
    static
    boost::shared_ptr<embedding_type> create(const path_type& path, const size_type& rows)
    {
      typedef boost::shared_ptr<embedding_type> shared_type;
      typedef std::string string_type;
      typedef utils::unordered_map<string_type, shared_type,
				   boost::hash<string_type>, std::equal_to<string_type>,
				   std::allocator<std::pair<const string_type, shared_type> > >::type mapped_type;
      
      typedef utils::spinlock         mutex_type;
      typedef mutex_type::scoped_lock lock_type;
      
      static mutex_type  mutex;
      static mapped_type mapped;
      
      lock_type lock(mutex);
      
      const string_type rep = path.string();
      
      mapped_type::iterator miter = mapped.find(rep);
      if (miter == mapped.end())
	miter = mapped.insert(std::make_pair(rep, shared_type(new EmbeddingMapped(path, rows)))).first;
      
      return miter->second;
    }
    
  private:
    vocab_type  vocab_;
    mapped_type mapped_;
    size_type   rows_;
  };
  
  class EmbeddingMemory : public TreeRNN::Embedding
  {
  private:
    typedef TreeRNN tree_rnn_type;
    
  public:
    typedef tree_rnn_type::size_type       size_type;
    typedef tree_rnn_type::difference_type difference_type;
    
    typedef tree_rnn_type::word_type  word_type;
    typedef tree_rnn_type::vocab_type vocab_type;

    typedef tree_rnn_type::path_type path_type;
    
    typedef tree_rnn_type::parameter_type parameter_type;
    typedef tree_rnn_type::tensor_type    tensor_type;
    typedef tree_rnn_type::matrix_type    matrix_type;
    typedef tree_rnn_type::embedding_type embedding_type;

  private:
    typedef utils::indexed_set<word_type, boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<word_type> > word_set_type;
  public:
    EmbeddingMemory() : matrix_(), words_(), oov_(0) {}
    EmbeddingMemory(const path_type& path, const size_type& rows)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;
      typedef boost::fusion::tuple<std::string, parameter_set_type > embedding_parsed_type;
      typedef boost::spirit::istream_iterator iterator_type;
      
      qi::rule<iterator_type, std::string(), standard::blank_type>           word;
      qi::rule<iterator_type, embedding_parsed_type(), standard::blank_type> parser; 
      
      word   %= qi::lexeme[+(standard::char_ - standard::space)];
      parser %= word >> *qi::double_ >> (qi::eol | qi::eoi);
      
      utils::compress_istream is(path, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iterator_type iter(is);
      iterator_type iter_end;
      
      embedding_parsed_type parsed;
      
      oov_ = word_type::id_type(-1);

      while (iter != iter_end) {
	boost::fusion::get<0>(parsed).clear();
	boost::fusion::get<1>(parsed).clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, parsed))
	  if (iter != iter_end)
	    throw std::runtime_error("embedding parsing failed");
	
	if (boost::fusion::get<1>(parsed).size() != rows)
	  throw std::runtime_error("invalid embedding size");
	
	const word_type word = boost::fusion::get<0>(parsed);
	
	word_set_type::const_iterator witer = words_.insert(word).first;
	const word_type::id_type id = witer - words_.begin();

	if (word == vocab_type::UNK)
	  oov_ = id;
	
	if (id >= matrix_.cols())
	  matrix_.conservativeResize(rows, id + 1);
	
	matrix_.col(id) = matrix_type(&(*boost::fusion::get<1>(parsed).begin()), rows, 1);
      }
      
      // estimate OOV by performing averaging...
      if (oov_ == word_type::id_type(-1)) {
	word_set_type::const_iterator witer = words_.insert(vocab_type::UNK).first;
	const word_type::id_type oov_ = witer - words_.begin();
	
	if (oov_ >= matrix_.cols())
	  matrix_.conservativeResize(rows, oov_ + 1);
	
	matrix_.col(oov_).setZero();
	
	const word_type::id_type vocabulary_size = matrix_.cols();
	const double factor = 1.0 / (vocabulary_size - 1);
	for (word_type::id_type id = 0; id != vocabulary_size; ++ id)
	  if (id != oov_)
	    matrix_.col(oov_) += matrix_.col(id) * factor;
      }
    }
    
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      
      rep["embedding"] = utils::lexical_cast<std::string>(matrix_.rows());
      
      utils::compress_ostream os(rep.path("input.bin"), 1024 * 1024);
      
      os.write((char*) matrix_.data(), sizeof(parameter_type) * matrix_.rows() * matrix_.cols());
      
      vocab_type vocab;
      
      vocab.open(rep.path("vocab"), words_.size());
      
      word_set_type::const_iterator witer_end = words_.end();
      for (word_set_type::const_iterator witer = words_.begin(); witer != witer_end; ++ witer)
	vocab.insert(*witer);
    }
    
    matrix_type operator()(const word_type& word) const
    {
      word_set_type::const_iterator witer = words_.find(word);
      
      word_type::id_type id = oov_;
      if (witer != words_.end())
	id = witer - words_.begin();
      
      const size_type rows = matrix_.rows();
      
      return matrix_type(const_cast<parameter_type*>(matrix_.data() + rows * id), rows, 1);
    }

  public:
    static
    boost::shared_ptr<embedding_type> create(const path_type& path, const size_type& rows)
    {
      typedef boost::shared_ptr<embedding_type> shared_type;
      typedef std::string string_type;
      typedef utils::unordered_map<string_type, shared_type,
				   boost::hash<string_type>, std::equal_to<string_type>,
				   std::allocator<std::pair<const string_type, shared_type> > >::type mapped_type;
      
      typedef utils::spinlock         mutex_type;
      typedef mutex_type::scoped_lock lock_type;
      
      static mutex_type  mutex;
      static mapped_type mapped;
      
      lock_type lock(mutex);
      
      const string_type rep = path.string();
      
      mapped_type::iterator miter = mapped.find(rep);
      if (miter == mapped.end())
	miter = mapped.insert(std::make_pair(rep, shared_type(new EmbeddingMemory(path, rows)))).first;
      
      return miter->second;
    }
    
  private:
    tensor_type        matrix_;
    word_set_type      words_;
    word_type::id_type oov_;
  };

  template <typename Path, typename Tensor>
  inline
  void write_matrix(const Path& path_txt,
		    const Path& path_bin,
		    const Tensor& matrix)
  {
    {
      utils::compress_ostream os(path_txt, 1024 * 1024);
      os.precision(10);
      os << matrix;
    }
    
    {
      utils::compress_ostream os(path_bin, 1024 * 1024);
      os.write((char*) matrix.data(), sizeof(typename Tensor::Scalar) * matrix.rows() * matrix.cols());
    }
  }
  
  template <typename Path, typename Tensor>
  inline
  void read_matrix(const Path& path,
		   Tensor& matrix)
  {
    const size_t file_size = boost::filesystem::file_size(path);
    
    if (file_size != sizeof(typename Tensor::Scalar) * matrix.rows() * matrix.cols())
      throw std::runtime_error("file size does not match: " + path.string());
    
    utils::compress_istream is(path, 1024 * 1024);
    
    is.read((char*) matrix.data(), file_size);
  }
  
  void TreeRNN::write(const path_type& path) const
  {
    typedef utils::repository repository_type;
    
    repository_type rep(path, repository_type::write);

    rep["hidden"]    = utils::lexical_cast<std::string>(hidden_);
    rep["embedding"] = utils::lexical_cast<std::string>(embedding_);
    
    write_matrix(rep.path("Wp.txt.gz"), rep.path("Wp.bin"), Wp_);
    write_matrix(rep.path("Bp.txt.gz"), rep.path("Bp.bin"), Bp_);

    write_matrix(rep.path("Wu.txt.gz"), rep.path("Wu.bin"), Wu_);
    write_matrix(rep.path("Bu.txt.gz"), rep.path("Bu.bin"), Bu_);

    write_matrix(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write_matrix(rep.path("Bt.txt.gz"), rep.path("Bt.bin"), Bt_);

    write_matrix(rep.path("Wn.txt.gz"), rep.path("Wn.bin"), Wn_);
    write_matrix(rep.path("Bn.txt.gz"), rep.path("Bn.bin"), Bn_);
    
    if (input_)
      input_->write(rep.path("input"));
  }
  
  void TreeRNN::open(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    if (path.empty() || ! boost::filesystem::exists(path))
      throw std::runtime_error("no file? " + path.string());
    
    repository_type rep(path, repository_type::read);
    
    hidden_    = repository_value<size_type>(rep, "hidden");
    embedding_ = repository_value<size_type>(rep, "embedding");
    
    if (hidden_ == 0)
      throw std::runtime_error("invalid dimension");
    if (embedding_ == 0)
      throw std::runtime_error("invalid dimension");

    Wp_ = tensor_type::Zero(hidden_, embedding_);
    Bp_ = tensor_type::Zero(hidden_, 1);

    Wu_ = tensor_type::Zero(hidden_, hidden_);
    Bu_ = tensor_type::Zero(hidden_, 1);
    
    Wt_ = tensor_type::Zero(hidden_, hidden_ + embedding_);
    Bt_ = tensor_type::Zero(hidden_, 1);
    
    Wn_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    Bn_ = tensor_type::Zero(hidden_, 1);
    
    read_matrix(rep.path("Wp.bin"), Wp_);
    read_matrix(rep.path("Bp.bin"), Bp_);
    
    read_matrix(rep.path("Wu.bin"), Wu_);
    read_matrix(rep.path("Bu.bin"), Bu_);
    
    read_matrix(rep.path("Wt.bin"), Wt_);
    read_matrix(rep.path("Bt.bin"), Bt_);
    
    read_matrix(rep.path("Wn.bin"), Wn_);
    read_matrix(rep.path("Bn.bin"), Bn_);
    
    input_ = EmbeddingMapped::create(rep.path("input"), embedding_);
  }
  
  void TreeRNN::open(const size_type& hidden, const size_type& embedding, const path_type& path)
  {
    if (path.empty() || ! boost::filesystem::exists(path))
      throw std::runtime_error("no embedding? " + path.string());

    hidden_    = hidden;
    embedding_ = embedding;

    if (hidden_ == 0)
      throw std::runtime_error("invalid dimension");
    if (embedding_ == 0)
      throw std::runtime_error("invalid dimension");
    
    Wp_ = tensor_type::Zero(hidden_, embedding_);
    Bp_ = tensor_type::Zero(hidden_, 1);

    Wu_ = tensor_type::Zero(hidden_, hidden_);
    Bu_ = tensor_type::Zero(hidden_, 1);
    
    Wt_ = tensor_type::Zero(hidden_, hidden_ + embedding_);
    Bt_ = tensor_type::Zero(hidden_, 1);
    
    Wn_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    Bn_ = tensor_type::Zero(hidden_, 1);

    input_ = (boost::filesystem::is_directory(path)
	      ? EmbeddingMapped::create(path, embedding_)
	      : EmbeddingMemory::create(path, embedding_));
  }
};
