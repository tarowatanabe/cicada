// -*- mode: c++ -*-
//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_RNN__HPP__
#define __CICADA__TREE_RNN__HPP__ 1

// binary tree RNN model
// this is simply a place-holder for actual feature in feature directory...
// + we will share word embedding to save memory space...

#include <stdint.h>

#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/word_embedding.hpp>

#include <utils/bithack.hpp>
#include <utils/mathop.hpp>

#include <Eigen/Core>

namespace cicada
{
  class TreeRNN
  {
  public:
    typedef Symbol                  word_type;
    typedef Vocab                   vocab_type;
    
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    typedef float              logprob_type;
    typedef float              parameter_type;
    typedef double             prob_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic>    tensor_type;
    typedef Eigen::Map<tensor_type>                                          matrix_type;

    typedef WordEmbedding embedding_type;
    
    struct shtanh
    {
      template <typename Tp>
      Tp operator()(const Tp& x) const
      {
	return int(std::min(std::max(x, Tp(- 1)), Tp(1)) * 128) / Tp(128);
      }
    };
    
    struct dshtanh
    {
      template <typename Tp>
      Tp operator()(const Tp& x) const
      {
	return Tp(- 1) < x && x < Tp(1);
      }
    };
    
    
  public:
    TreeRNN()
      : hidden_(0), embedding_(0), input_(0) {}
    TreeRNN(const path_type& path)
    { open(path); }
    TreeRNN(const size_type& hidden, const size_type& embedding)
    { initialize(hidden, embedding); }
    TreeRNN(const size_type& hidden, const size_type& embedding, const path_type& path)
    { open(hidden, embedding, path); }
    
    void write(const path_type& path) const;
    
    void open(const path_type& path);
    void open(const size_type& hidden, const size_type& embedding, const path_type& path);

    void initialize(const size_type& hidden, const size_type& embedding);

  private:
    template <typename Gen>
    struct __randomize
    {
      __randomize(Gen& gen, const double range=0.01) : gen_(gen), range_(range) {}
      
      template <typename Tp>
      Tp operator()(const Tp& x) const
      {
	return boost::random::uniform_real_distribution<Tp>(-range_, range_)(const_cast<Gen&>(gen_));
      }
      
      Gen& gen_;
      double range_;
    };
    
  public:
    template <typename Gen>
    void random(Gen& gen)
    {
      const double range_t = std::sqrt(6.0 / (hidden_ + hidden_ + embedding_));
      const double range_n = std::sqrt(6.0 / (hidden_ + hidden_ + hidden_));
      
      Wt_ = Wt_.array().unaryExpr(__randomize<Gen>(gen, range_t));
      Wn_ = Wn_.array().unaryExpr(__randomize<Gen>(gen, range_n));
    }

    void clear()
    {
      Wt_.setZero();
      Bt_.setZero();

      Wn_.setZero();
      Bn_.setZero();
      
      Bi_.setZero();
    }

    void swap(TreeRNN& x)
    {
      std::swap(hidden_,    x.hidden_);
      std::swap(embedding_, x.embedding_);

      Wt_.swap(x.Wt_);
      Bt_.swap(x.Bt_);

      Wn_.swap(x.Wn_);
      Bn_.swap(x.Bn_);

      Bi_.swap(x.Bi_);

      std::swap(input_, x.input_);
    }

  public:
    friend
    std::ostream& operator<<(std::ostream& os, const TreeRNN& rnn);
    friend
    std::istream& operator>>(std::istream& is, TreeRNN& rnn);

  public:
    template <typename Tp>
    TreeRNN& operator*=(const Tp& x)
    {
      Wt_ *= x;
      Bt_ *= x;
      
      Wn_ *= x;
      Bn_ *= x;
      
      Bi_ *= x;

      return *this;
    }

    TreeRNN& operator+=(const TreeRNN& x)
    {
      Wt_ += x.Wt_;
      Bt_ += x.Bt_;
      
      Wn_ += x.Wn_;
      Bn_ += x.Bn_;
      
      Bi_ += x.Bi_;

      return *this;
    }

    TreeRNN& operator-=(const TreeRNN& x)
    {
      Wt_ -= x.Wt_;
      Bt_ -= x.Bt_;
      
      Wn_ -= x.Wn_;
      Bn_ -= x.Bn_;
      
      Bi_ -= x.Bi_;

      return *this;
    }
    
  public:
    size_type hidden_;
    size_type embedding_;
    
    // binary rule for terminal
    tensor_type Wt_;
    tensor_type Bt_;
    
    // binary rule for non-terminal
    tensor_type Wn_;
    tensor_type Bn_;
    
    // bias for initial state
    tensor_type Bi_;
    
    const embedding_type* input_;
  };
};

namespace std
{
  inline
  void swap(cicada::TreeRNN& x, cicada::TreeRNN& y)
  {
    x.swap(y);
  }
  
};

#endif
