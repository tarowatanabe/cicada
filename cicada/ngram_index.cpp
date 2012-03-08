//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <sstream>
#include <iomanip>

#include "ngram_index.hpp"

#include "utils/lexical_cast.hpp"

namespace cicada
{

  void NGramIndex::Shard::open(const path_type& path)
  {
    typedef utils::repository repository_type;
      
    clear();
      
    repository_type rep(path, repository_type::read);
      
    ids.open(rep.path("index"));
    positions.open(rep.path("position"));
      
    repository_type::const_iterator oiter = rep.find("order");
    if (oiter == rep.end())
      throw std::runtime_error("no order");
    const int order = utils::lexical_cast<int>(oiter->second);
      
    offsets.push_back(0);
    for (int n = 1; n <= order; ++ n) {
      std::ostringstream stream_ngram;
      stream_ngram << n << "-gram-offset";
	
      repository_type::const_iterator iter = rep.find(stream_ngram.str());
      if (iter == rep.end())
	throw std::runtime_error(std::string("no ngram offset? ") + stream_ngram.str());
	
      offsets.push_back(utils::lexical_cast<size_type>(iter->second));
    }

    off_set_type(offsets).swap(offsets);
  }
  
  void NGramIndex::open(const path_type& path)
  {
    typedef utils::repository repository_type;

    close();
    
    repository_type rep(path, repository_type::read);
    
    // shard size
    repository_type::const_iterator siter = rep.find("shard");
    if (siter == rep.end())
      throw std::runtime_error("no shard size...");
    __shards.resize(utils::lexical_cast<size_type>(siter->second));
    
    // order
    repository_type::const_iterator oiter = rep.find("order");
    if (oiter == rep.end())
      throw std::runtime_error("no order");
    __order = utils::lexical_cast<int>(oiter->second);
    
    // vocabulary...
    __vocab.open(rep.path("vocab"));
    
    for (size_t shard = 0; shard != __shards.size(); ++ shard) {
      std::ostringstream stream_shard;
      stream_shard << "ngram-" << std::setfill('0') << std::setw(6) << shard;
      
      __shards[shard].open(rep.path(stream_shard.str()));
    }
    
    __path = path;
  }
  

};
