//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/ngram_rnn.hpp"

#include "utils/repository.hpp"

#include "utils/spinlock.hpp"
#include "utils/unordered_map.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/thread_specific_ptr.hpp"

#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>

namespace cicada
{

  template <typename Value>
  inline
  Value repository_value(const utils::repository& rep, const std::string& key)
  {
    utils::repository::const_iterator iter = rep.find(key);
    if (iter == rep.end())
      throw std::runtime_error("no " + key + "?");
    return boost::lexical_cast<Value>(iter->second);
  }
  
  void NGramRNN::open(const path_type& path)
  {
    typedef utils::repository repository_type;

    clear();
    
    repository_type rep(path, repository_type::read);

    path_ = path;
    
    embedding_size_ = repository_value<size_type>(rep, "size");
    dimension_      = repository_value<size_type>(rep, "dimension");
    order_          = repository_value<int>(rep, "order");
    
    embedding_input_.open(rep.path("input.bin"), dimension_, embedding_size_);
    embedding_output_.open(rep.path("output.bin"), dimension_ + 1, embedding_size_);
    
    Wc_.open(rep.path("Wc.bin"), dimension_, dimension_ * 2 * (order_ - 1));
    bc_.open(rep.path("bc.bin"), dimension_, order_ - 1);
    
    bi_.open(rep.path("bi.bin"), dimension_, 1);
    
    vocab_.open(rep.path("vocab"));
    
    id_bos_ = vocab_[vocab_type::BOS];
    id_eos_ = vocab_[vocab_type::EOS];
    id_eps_ = vocab_[vocab_type::EPSILON];
    id_unk_ = vocab_[vocab_type::UNK];
    
    buffer_ = buffer_type(locks_.size(), dimension_);
    
    for (size_type i = 0; i != cache_.size(); ++ i)
      cache_[i] = cache_type(order_);

    // initialize init_ buffers
    
    const size_type offset_embedding = 0;
    const size_type offset_context   = dimension_;
    
    init_ = buffer_type(order_, dimension_);
    
    for (size_type n = 0; n != order_; ++ n) {
      matrix_type context(&(*init_.begin(n)), dimension_, 1);
      
      context = bi_().array().unaryExpr(hinge());
      
      for (size_type i = 0; i != n; ++ i) {
	const size_type shift = i * 2 * dimension_;
	
	context = (Wc_().block(0, shift + offset_embedding, dimension_, dimension_) * embedding_input_().col(id_eps_)
		   + Wc_().block(0, shift + offset_context, dimension_, dimension_) * context
		   + bc_().block(0, i, dimension_, 1)).array().unaryExpr(hinge());
      }
    }
  }
  
  typedef utils::unordered_map<std::string, NGramRNN, boost::hash<utils::piece>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, NGramRNN> > >::type ngram_rnn_map_type;
  
  namespace impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type         __ngram_rnn_mutex;
    static ngram_rnn_map_type __ngram_rnn_map;
  };
  
#ifdef HAVE_TLS
  static __thread ngram_rnn_map_type* __ngram_rnns_tls = 0;
  static utils::thread_specific_ptr<ngram_rnn_map_type> __ngram_rnns;
#else
  static utils::thread_specific_ptr<ngram_rnn_map_type> __ngram_rnns;
#endif
  
  NGramRNN& NGramRNN::create(const path_type& path)
  {
#ifdef HAVE_TLS
    if (! __ngram_rnns_tls) {
      __ngram_rnns.reset(new ngram_rnn_map_type());
      __ngram_rnns_tls = __ngram_rnns.get();
    }
    ngram_rnn_map_type& ngram_rnns_map = *__ngram_rnns_tls;
#else
    if (! __ngram_rnns.get())
      __ngram_rnns.reset(new ngram_rnn_map_type());
    
    ngram_rnn_map_type& ngram_rnns_map = *__ngram_rnns;
#endif
    
    const std::string parameter = path.string();
    
    ngram_rnn_map_type::iterator iter = ngram_rnns_map.find(parameter);
    if (iter == ngram_rnns_map.end()) {
      impl::lock_type lock(impl::__ngram_rnn_mutex);
      
      ngram_rnn_map_type::iterator iter_global = impl::__ngram_rnn_map.find(parameter);
      if (iter_global == impl::__ngram_rnn_map.end())
	iter_global = impl::__ngram_rnn_map.insert(std::make_pair(parameter, NGramRNN(parameter))).first;
      
      iter = ngram_rnns_map.insert(*iter_global).first;
    }
    
    return iter->second;
  }
  
};
