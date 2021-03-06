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
  void open_shards(const Path& path, Shards& shards)
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
    }
  }

  void NGram::open(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    clear();

    if (path.empty())
      throw std::runtime_error("no ngram?");
    else if (! boost::filesystem::exists(path))
      throw std::runtime_error("no ngram? " + path.string());
    
    repository_type rep(path, repository_type::read);
    
    index.open(rep.path("index"));
    
    if (boost::filesystem::exists(rep.path("logprob")))
      open_shards(rep.path("logprob"), logprobs);
    
    if (boost::filesystem::exists(rep.path("backoff")))
      open_shards(rep.path("backoff"), backoffs);
    
    if (boost::filesystem::exists(rep.path("logbound")))
      open_shards(rep.path("logbound"), logbounds);
    
    repository_type::const_iterator siter = rep.find("smooth");
    if (siter == rep.end())
      throw std::runtime_error("no smoothing parameter...?");
    smooth = utils::lexical_cast<logprob_type>(siter->second);
    
    if (debug)
      std::cerr << "ngram: " << path
		<< " # of shards: " << index.size()
		<< " smooth: " << smooth
		<< std::endl;
  }

  typedef utils::unordered_map<std::string, NGram, boost::hash<utils::piece>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, NGram> > >::type ngram_map_type;

  namespace impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type     __ngram_mutex;
    static ngram_map_type __ngram_map;
  };


  NGram& NGram::create(const path_type& path)
  {
    const std::string parameter = path.string();
    
    impl::lock_type lock(impl::__ngram_mutex);
    
    ngram_map_type::iterator iter = impl::__ngram_map.find(parameter);
    if (iter == impl::__ngram_map.end())
      iter = impl::__ngram_map.insert(std::make_pair(parameter, NGram(parameter))).first;
    
    return iter->second;
  }
  
};
