//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/ngram_nn.hpp"

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
  
  void NGramNN::open(const path_type& path)
  {
    typedef utils::repository repository_type;

    clear();
    
    repository_type rep(path, repository_type::read);
    
    embedding_size_      = repository_value<size_type>(rep, "size");
    dimension_embedding_ = repository_value<size_type>(rep, "embedding");
    dimension_hidden_    = repository_value<size_type>(rep, "hidden");
    order_               = repository_value<int>(rep, "order");
    
    embedding_input_.open(rep.path("input.bin"), dimension_embedding_, embedding_size_);
    embedding_output_.open(rep.path("output.bin"), dimension_embedding_ + 1, embedding_size_);
    
    Wc_.open(rep.path("Wc.bin"), dimension_hidden_, dimension_embedding_ * (order_ - 1));
    bc_.open(rep.path("bc.bin"), dimension_hidden_, 1);
    
    Wh_.open(rep.path("Wh.bin"), dimension_embedding_, dimension_hidden_);
    bh_.open(rep.path("bh.bin"), dimension_embedding_, 1);
    
    vocab_.open(rep.path("vocab"));

    id_bos_ = vocab_[vocab_type::BOS];
    id_eos_ = vocab_[vocab_type::EOS];
    id_eps_ = vocab_[vocab_type::EPSILON];
    id_unk_ = vocab_[vocab_type::UNK];
    
    buffer_ = buffer_type(16, dimension_embedding_ * (order_ - 1));
    cache_ = cache_type(order_);
  }
  
  typedef utils::unordered_map<std::string, NGramNN, boost::hash<utils::piece>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, NGramNN> > >::type ngram_nn_map_type;
  
  namespace impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type         __ngram_nn_mutex;
    static ngram_nn_map_type __ngram_nn_map;
  };
  
#ifdef HAVE_TLS
  static __thread ngram_nn_map_type* __ngram_nns_tls = 0;
  static utils::thread_specific_ptr<ngram_nn_map_type> __ngram_nns;
#else
  static utils::thread_specific_ptr<ngram_nn_map_type> __ngram_nns;
#endif
  
  NGramNN& NGramNN::create(const path_type& path)
  {
#ifdef HAVE_TLS
    if (! __ngram_nns_tls) {
      __ngram_nns.reset(new ngram_nn_map_type());
      __ngram_nns_tls = __ngram_nns.get();
    }
    ngram_nn_map_type& ngram_nns_map = *__ngram_nns_tls;
#else
    if (! __ngram_nns.get())
      __ngram_nns.reset(new ngram_nn_map_type());
    
    ngram_nn_map_type& ngram_nns_map = *__ngram_nns;
#endif
    
    const std::string parameter = path.string();
    
    ngram_nn_map_type::iterator iter = ngram_nns_map.find(parameter);
    if (iter == ngram_nns_map.end()) {
      impl::lock_type lock(impl::__ngram_nn_mutex);
      
      ngram_nn_map_type::iterator iter_global = impl::__ngram_nn_map.find(parameter);
      if (iter_global == impl::__ngram_nn_map.end())
	iter_global = impl::__ngram_nn_map.insert(std::make_pair(parameter, NGramNN(parameter))).first;
      
      iter = ngram_nns_map.insert(*iter_global).first;
    }
    
    return iter->second;
  }
  
};
