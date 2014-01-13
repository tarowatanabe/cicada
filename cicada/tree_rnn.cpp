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
    
    write_matrix(rep.path("Wt.txt.gz"), rep.path("Wt.bin"), Wt_);
    write_matrix(rep.path("Bt.txt.gz"), rep.path("Bt.bin"), Bt_);

    write_matrix(rep.path("Wn.txt.gz"), rep.path("Wn.bin"), Wn_);
    write_matrix(rep.path("Bn.txt.gz"), rep.path("Bn.bin"), Bn_);
    
    write_matrix(rep.path("Bi.txt.gz"), rep.path("Bi.bin"), Bi_);
    
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

    Wt_ = tensor_type::Zero(hidden_, hidden_ + embedding_);
    Bt_ = tensor_type::Zero(hidden_, 1);
    
    Wn_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    Bn_ = tensor_type::Zero(hidden_, 1);

    Bi_ = tensor_type::Zero(hidden_, 1);

    read_matrix(rep.path("Wt.bin"), Wt_);
    read_matrix(rep.path("Bt.bin"), Bt_);
    
    read_matrix(rep.path("Wn.bin"), Wn_);
    read_matrix(rep.path("Bn.bin"), Bn_);
    
    read_matrix(rep.path("Bi.bin"), Bi_);
    
    input_ = &embedding_type::create(rep.path("input"));

    if (input_->dimension() != embedding_)
      throw std::runtime_error("invalid dimension for word embedding");
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
    
    Wt_ = tensor_type::Zero(hidden_, hidden_ + embedding_);
    Bt_ = tensor_type::Zero(hidden_, 1);
    
    Wn_ = tensor_type::Zero(hidden_, hidden_ + hidden_);
    Bn_ = tensor_type::Zero(hidden_, 1);

    Bi_ = tensor_type::Zero(hidden_, 1);

    input_ = &embedding_type::create(path);

    if (input_->dimension() != embedding_)
      throw std::runtime_error("invalid dimension for word embedding");
  }
};
