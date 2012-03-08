//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "ngram_pyp.hpp"

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
  
  void NGramPYP::open(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    clear();
    
    repository_type rep(path, repository_type::read);

    // first, read parameters...
    
    order_ = repository_value<int>(rep, "order");

    if (order_ < 0)
      throw std::runtime_error("negative order?");
    
    p0_      = repository_value<double>(rep, "p0");
    counts0_ = repository_value<count_type>(rep, "counts0");
    
    discount_alpha_ = repository_value<double>(rep, "discount-alpha");
    discount_beta_  = repository_value<double>(rep, "discount-beta");
    strength_shape_ = repository_value<double>(rep, "strength-shape");
    strength_rate_  = repository_value<double>(rep, "strength-rate");
    
    discount_.resize(order_);
    strength_.resize(order_);
    
    for (size_type n = 0; n != discount_.size(); ++ n) {
      discount_[n] = repository_value<double>(rep, "discount" + boost::lexical_cast<std::string>(n));
      strength_[n] = repository_value<double>(rep, "strength" + boost::lexical_cast<std::string>(n));
    }
    
    parameter_set_type(discount_).swap(discount_);
    parameter_set_type(strength_).swap(strength_);
    
    // then, models...
    index_.open(rep.path("index"));
    count_.open(rep.path("count"));
    total_.open(rep.path("total"));
    
    // positions
    positions_.open(rep.path("position"));
    
    // vocab
    vocab_.open(rep.path("vocab"));
    
    // offsets
    offsets_.push_back(0);
    for (int n = 1; n <= order_; ++ n)
      offsets_.push_back(repository_value<size_type>(rep, boost::lexical_cast<std::string>(n) + "-gram-offset"));
    
    offset_set_type(offsets_).swap(offsets_);
  }

  template <typename Tp>
  struct hash_string : public utils::hashmurmur<size_t>
  {
    size_t operator()(const Tp& x) const
    {
      return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
    }
  };

  typedef utils::unordered_map<std::string, NGramPYP, hash_string<std::string>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, NGramPYP> > >::type ngram_pyp_map_type;

  namespace impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type         __ngram_pyp_mutex;
    static ngram_pyp_map_type __ngram_pyp_map;
  };

#ifdef HAVE_TLS
  static __thread ngram_pyp_map_type* __ngram_pyps_tls = 0;
  static boost::thread_specific_ptr<ngram_pyp_map_type> __ngram_pyps;
#else
  static utils::thread_specific_ptr<ngram_pyp_map_type> __ngram_pyps;
#endif
  
  NGramPYP& NGramPYP::create(const path_type& path)
  {
#ifdef HAVE_TLS
    if (! __ngram_pyps_tls) {
      __ngram_pyps.reset(new ngram_pyp_map_type());
      __ngram_pyps_tls = __ngram_pyps.get();
    }
    ngram_pyp_map_type& ngram_pyps_map = *__ngram_pyps_tls;
#else
    if (! __ngram_pyps.get())
      __ngram_pyps.reset(new ngram_pyp_map_type());
    
    ngram_pyp_map_type& ngram_pyps_map = *__ngram_pyps;
#endif

    const std::string parameter = path.string();
    
    ngram_pyp_map_type::iterator iter = ngram_pyps_map.find(parameter);
    if (iter == ngram_pyps_map.end()) {
      impl::lock_type lock(impl::__ngram_pyp_mutex);
      
      ngram_pyp_map_type::iterator iter_global = impl::__ngram_pyp_map.find(parameter);
      if (iter_global == impl::__ngram_pyp_map.end())
	iter_global = impl::__ngram_pyp_map.insert(std::make_pair(parameter, NGramPYP(parameter))).first;
      
      iter = ngram_pyps_map.insert(*iter_global).first;
    }
    
    return iter->second;
  }
  
};
