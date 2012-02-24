//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "ngram.hpp"

#include "utils/spinlock.hpp"
#include "utils/unordered_map.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/thread_specific_ptr.hpp"

#include <boost/thread.hpp>

namespace cicada
{
  void NGram::ShardData::open(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    clear();
    
    repository_type rep(path, repository_type::read);
    
    repository_type::const_iterator oiter = rep.find("offset");
    if (oiter == rep.end())
      throw std::runtime_error("no offset?");
    offset = utils::lexical_cast<size_type>(oiter->second);
    
    if (boost::filesystem::exists(rep.path("quantized"))) {
      quantized.open(rep.path("quantized"));
      
      logprob_map_type logprob_map;
      
      maps.clear();
      maps.push_back(logprob_map);
      for (int n = 1; /**/; ++ n) {
	std::ostringstream stream_map_file;
	stream_map_file << n << "-logprob-map";
	
	if (! boost::filesystem::exists(rep.path(stream_map_file.str()))) break;
	
	std::ifstream is(rep.path(stream_map_file.str()).string().c_str());
	is.read((char*) &(*logprob_map.begin()), sizeof(logprob_type) * logprob_map.size());
	maps.push_back(logprob_map);
      }
    } else
      logprobs.open(rep.path("logprob"));
  }
    

  template <typename Path, typename Shards>
  inline
  void open_shards(const Path& path, Shards& shards, const int debug)
  {
    typedef utils::repository repository_type;
    
    repository_type rep(path, repository_type::read);
    
    repository_type::const_iterator siter = rep.find("shard");
    if (siter == rep.end())
      throw std::runtime_error("no shard size...");

    shards.clear();
    shards.reserve(utils::lexical_cast<size_t>(siter->second));
    shards.resize(utils::lexical_cast<size_t>(siter->second));
    
    for (size_t shard = 0; shard != shards.size(); ++ shard) {
      std::ostringstream stream_shard;
      stream_shard << "ngram-" << std::setfill('0') << std::setw(6) << shard;
      
      shards[shard].open(rep.path(stream_shard.str()));

      if (debug >= 2)
	std::cerr << "ngram data: " << rep.path(stream_shard.str()) << std::endl;
    }
  }

  void NGram::open(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    clear();
    
    repository_type rep(path, repository_type::read);
    
    index.open(rep.path("index"));
    if (boost::filesystem::exists(rep.path("logprob")))
      open_shards(rep.path("logprob"), logprobs, debug);
    if (boost::filesystem::exists(rep.path("backoff")))
      open_shards(rep.path("backoff"), backoffs, debug);
    if (boost::filesystem::exists(rep.path("logbound")))
      open_shards(rep.path("logbound"), logbounds, debug);
    
    repository_type::const_iterator siter = rep.find("smooth");
    if (siter == rep.end())
      throw std::runtime_error("no smoothing parameter...?");
    smooth = atof(siter->second.c_str());
    
    if (debug)
      std::cerr << "ngram: " << path
		<< " # of shards: " << index.size()
		<< " smooth: " << smooth
		<< std::endl;
  }


  template <typename Tp>
  struct hash_string : public utils::hashmurmur<size_t>
  {
    size_t operator()(const Tp& x) const
    {
      return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
    }
  };

  typedef utils::unordered_map<std::string, NGram, hash_string<std::string>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, NGram> > >::type ngram_map_type;

  namespace impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type     __ngram_mutex;
    static ngram_map_type __ngram_map;
  };

#ifdef HAVE_TLS
  static __thread ngram_map_type* __ngrams_tls = 0;
  static boost::thread_specific_ptr<ngram_map_type> __ngrams;
#else
  static utils::thread_specific_ptr<ngram_map_type> __ngrams;
#endif

  NGram& NGram::create(const path_type& path)
  {
#ifdef HAVE_TLS
    if (! __ngrams_tls) {
      __ngrams.reset(new ngram_map_type());
      __ngrams_tls = __ngrams.get();
    }
    ngram_map_type& ngrams_map = *__ngrams_tls;
#else
    if (! __ngrams.get())
      __ngrams.reset(new ngram_map_type());
    
    ngram_map_type& ngrams_map = *__ngrams;
#endif

    const std::string parameter = path.string();
    
    ngram_map_type::iterator iter = ngrams_map.find(parameter);
    if (iter == ngrams_map.end()) {
      impl::lock_type lock(impl::__ngram_mutex);
      
      ngram_map_type::iterator iter_global = impl::__ngram_map.find(parameter);
      if (iter_global == impl::__ngram_map.end())
	iter_global = impl::__ngram_map.insert(std::make_pair(parameter, NGram(parameter))).first;
      
      iter = ngrams_map.insert(*iter_global).first;
    }
    
    return iter->second;
  }
  
};
