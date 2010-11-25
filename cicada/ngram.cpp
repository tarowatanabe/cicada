//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "ngram.hpp"

#include "utils/lexical_cast.hpp"

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
    offset = atoll(oiter->second.c_str());
    
    if (boost::filesystem::exists(rep.path("quantized"))) {
      quantized.open(rep.path("quantized"));
      
      logprob_map_type logprob_map;
      
      maps.clear();
      maps.push_back(logprob_map);
      for (int n = 1; /**/; ++ n) {
	std::ostringstream stream_map_file;
	stream_map_file << n << "-logprob-map";
	
	if (! boost::filesystem::exists(rep.path(stream_map_file.str()))) break;
	
	std::ifstream is(rep.path(stream_map_file.str()).file_string().c_str());
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
    shards.reserve(atoi(siter->second.c_str()));
    shards.resize(atoi(siter->second.c_str()));
    
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
    
    repository_type::const_iterator biter = rep.find("bound-exact");
    if (biter != rep.end())
      bound_exact = utils::lexical_cast<bool>(biter->second);
    else
      bound_exact = false;
    
    if (debug)
      std::cerr << "ngram: " << path
		<< " # of shards: " << index.size()
		<< " smooth: " << smooth
		<< std::endl;
  }
};
