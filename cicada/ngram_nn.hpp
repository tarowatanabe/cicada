// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_NN__HPP__
#define __CICADA__NGRAM_NN__HPP__ 1

// neural network ngram langauge model!

#include <stdint.h>

#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/ngram_cache.hpp>

#include <utils/hashmurmur3.hpp>
#include <utils/array_power2.hpp>
#include <utils/vector2.hpp>
#include <utils/spinlock.hpp>
#include <utils/bithack.hpp>
#include <utils/mathop.hpp>

#include <Eigen/Core>

namespace cicada
{
  class NGramNN : public utils::hashmurmur3<size_t>
  {
  public:
    typedef Symbol                  word_type;
    typedef Vocab                   vocab_type;

    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    typedef word_type::id_type id_type;
    typedef float              logprob_type;
    typedef float              parameter_type;
    typedef double             prob_type;
    
    typedef boost::filesystem::path   path_type;


  private:
    typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic>    tensor_type;
    typedef Eigen::Map<tensor_type>                                          matrix_type;
    
    class MappedMatrix
    {
    public:
      
      typedef utils::map_file<parameter_type, std::allocator<parameter_type> > mapped_type;
      
    public:
      MappedMatrix() : mapped_(), rows_(0), cols_(0) {}
      MappedMatrix(const path_type& path, const size_type rows, const size_type cols)
	: mapped_(path), rows_(rows), cols_(cols)
      {
	if (mapped_.size() != rows * cols)
	  throw std::runtime_error("mapped region does not match");
      }

      void open(const path_type& path, const size_type rows, const size_type cols)
      {
	mapped_.open(path);
	rows_ = rows;
	cols_ = cols;
	
	if (mapped_.size() != rows * cols)
	  throw std::runtime_error("mapped region does not match");
      }
      
      matrix_type operator()() const { return matrix_type(const_cast<parameter_type*>(mapped_.begin()), rows_, cols_); }
      
    public:
      void populate() { mapped_.populate(); }
      void clear()
      {
	mapped_.close();
	rows_ = 0;
	cols_ = 0;
      }
      
    public:      
      mapped_type mapped_;
      size_type rows_;
      size_type cols_;
    };

    typedef MappedMatrix mapped_matrix_type;
    
    struct spinlock_type
    {
      typedef utils::spinlock             mutex_type;
      typedef mutex_type::scoped_lock     lock_type;
      typedef mutex_type::scoped_try_lock trylock_type;
      
      spinlock_type() : mutex_() {}
      spinlock_type(const spinlock_type&) {}
      spinlock_type& operator=(const spinlock_type&) { return *this; }
      
      mutex_type mutex_;
    };
    typedef utils::array_power2<spinlock_type, 16, std::allocator<spinlock_type> > spinlock_set_type;
    
    typedef utils::vector2<parameter_type, std::allocator<parameter_type> > buffer_type;
    typedef cicada::NGramCache<id_type, logprob_type>                       cache_type;
    
  public:
    NGramNN() { clear(); }
    NGramNN(const path_type& path) { open(path); }
    
  public:
   static NGramNN& create(const path_type& path);

  public:
    const vocab_type& vocab() const { return vocab_; }
    
    size_type dimension_embedding() const { return dimension_embedding_; }
    size_type dimension_hidden() const { return dimension_hidden_; }
    const int& order() const { return order_; }
    
    path_type path() const { return path_; }
    bool empty() const { return ! path_.empty(); }
    
    void open(const path_type& path);
    void close() { clear(); }

    void populate()
    {
      embedding_input_.populate();
      embedding_output_.populate();
      
      Wc_.populate();
      bc_.populate();

      Wh_.populate();
      bh_.populate();
      
      vocab_.populate();
    }

    void clear()
    {
      embedding_input_.clear();
      embedding_output_.clear();
      
      Wc_.clear();
      bc_.clear();

      Wh_.clear();
      bh_.clear();
      
      vocab_.clear();

      id_bos_ = id_type(-1);
      id_eos_ = id_type(-1);
      id_eps_ = id_type(-1);
      id_unk_ = id_type(-1);

      embedding_size_      = 0;
      dimension_embedding_ = 0;
      dimension_hidden_    = 0;
      order_               = 0;
      
      path_ = path_type();

      buffer_.clear();
      cache_.clear();
    }

  public:
    template <typename Iterator>
    logprob_type operator()(Iterator first, Iterator last) const
    {
      typedef typename std::iterator_traits<Iterator>::value_type value_type;

      first = std::max(first, last - order_);
      
      if (first == last) return 0;
      
      return logprob_dispatch(first, last, value_type());
    }
    
  private:
    template <typename Iterator, typename __Word>
    logprob_type logprob_dispatch(Iterator first, Iterator last, __Word) const
    {
      typedef std::vector<id_type, std::allocator<id_type> > buffer_type;
      
      buffer_type buffer(last - first);
      
      buffer_type::iterator biter = buffer.begin();
      for (/**/; first != last; ++ first, ++ biter)
	*biter = vocab_[*first];
      
      return logprob_dispatch(buffer.begin(), buffer.end(), id_type());
    }
    
    template <typename Iterator>
    logprob_type logprob_dispatch(Iterator first, Iterator last, id_type) const
    {
      const size_type pos = cache_(first, last);
      
      spinlock_type::lock_type lock(const_cast<spinlock_type&>(locks_[pos % (locks_.size() - 1)]).mutex_);
      
      if (! cache_.equal_to(pos, first, last)) {
	const_cast<cache_type&>(cache_).assign(pos, first, last);
	const_cast<cache_type&>(cache_)[pos] = logprob_buffer(first, last, const_cast<float*>(&(*buffer_.begin(pos % (buffer_.size1() - 1)))));
      }
      
      return cache_[pos];
    }

    struct hinge
    {
      // 50 for numerical stability...
      template <typename Tp>
      Tp operator()(const Tp& x) const
      {
	return std::min(std::max(x, Tp(0)), Tp(50));
      }
    };
    
    template <typename Iterator>
    logprob_type logprob_buffer(Iterator first, Iterator last, void* buffer) const
    {
      matrix_type input(reinterpret_cast<parameter_type*>(buffer), dimension_embedding_ * (order_ - 1), 1);
      
      if (last - first < order_) {
	const id_type id_fill = (last - first == 1 || *first != id_bos_ ? id_eps_ : id_bos_);
	
	size_type i = 0;
	for (/**/; i < order_ - (last - first); ++ i)
	  input.block(dimension_embedding_ * i, 0, dimension_embedding_, 1) = embedding_input_().col(id_fill);
	
	for (/**/; first != last - 1; ++ first, ++ i)
	  input.block(dimension_embedding_ * i, 0, dimension_embedding_, 1) = embedding_input_().col(*first);
      } else {
	size_type i = 0;
	for (/**/; first != last - 1; ++ first, ++ i)
	  input.block(dimension_embedding_ * i, 0, dimension_embedding_, 1) = embedding_input_().col(*first);
      }
      
      return (embedding_output_().col(*(last - 1)).block(0, 0, dimension_embedding_, 1).transpose()
	      * (Wh_() * (Wc_() * input
			  + bc_()).array().unaryExpr(hinge()).matrix()
		 + bh_()).array().unaryExpr(hinge()).matrix()
	      + embedding_output_().col(*(last - 1)).block(dimension_embedding_, 0, 1, 1))(0, 0);
    }

  private:
    // word embedding
    mapped_matrix_type embedding_input_;
    mapped_matrix_type embedding_output_;
    
    // Wc and bc for context layer
    mapped_matrix_type Wc_;
    mapped_matrix_type bc_;
    
    // Wh and bh for hidden layer
    mapped_matrix_type Wh_;
    mapped_matrix_type bh_;
    
    vocab_type vocab_;

    id_type id_bos_;
    id_type id_eos_;
    id_type id_eps_;
    id_type id_unk_;

    size_type embedding_size_;
    size_type dimension_embedding_;
    size_type dimension_hidden_;
    int       order_;
    
    // path to the directory...
    path_type path_;
    
    buffer_type       buffer_;
    cache_type        cache_;
    spinlock_set_type locks_;
  };
};


#endif
