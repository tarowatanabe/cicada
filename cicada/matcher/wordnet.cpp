#include <cstdlib>
#include <stdexcept>

#include "cicada/matcher/wordnet.hpp"

#include "utils/spinlock.hpp"
#include "utils/simple_vector.hpp"
#include "utils/array_power2.hpp"

#include "wn/wordnet.hpp"

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

namespace cicada
{
  namespace matcher
  {
    WordNet::WordNet() { wn::WordNet __wn; }
    WordNet::WordNet(const std::string& path) { wn::WordNet __wn(path); }
    
    // do we cache...?
    
    typedef cicada::Symbol symbol_type;
    typedef utils::simple_vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    struct __wordnet_cache
    {
      symbol_type     word;
      symbol_set_type synsets;

      __wordnet_cache() : word(), synsets() {}
    };

    typedef __wordnet_cache cache_type;
    
    typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;
    
    
    bool WordNet::operator()(const symbol_type& x, const symbol_type& y) const
    {
#ifdef HAVE_TLS
      static __thread cache_set_type* __caches_tls = 0;
      static boost::thread_specific_ptr<cache_set_type> __caches;
      
      if (! __caches_tls) {
	__caches.reset(new cache_set_type());
	__caches_tls = __caches.get();
      }
      cache_set_type& caches = *__caches_tls;    
#else
      static boost::thread_specific_ptr<cache_set_type> __caches;
      
      if (! __caches.get())
	__caches.reset(new cache_set_type());
      
      cache_set_type& caches = *__caches;
#endif
      
      const size_type pos_x = hash_value(x) & (caches.size() - 1);
      const size_type pos_y = hash_value(y) & (caches.size() - 1);

      if (pos_x != pos_y) {
	cache_type& cache_x = caches[pos_x];
	cache_type& cache_y = caches[pos_y];
	
	if (cache_x.word != x) {
	  wordnet_type::synset_set_type synsets;
	  wordnet(x, synsets);

	  cache_x.word = x;
	  cache_x.synsets.resize(synsets.size());
	  for (size_t i = 0; i != synsets.size(); ++ i)
	    cache_x.synsets[i] = synsets[i].word;
	  
	}

	if (cache_y.word != y) {
	  wordnet_type::synset_set_type synsets;
	  wordnet(y, synsets);
	  
	  cache_y.synsets.resize(synsets.size());
	  for (size_t i = 0; i != synsets.size(); ++ i)
	    cache_y.synsets[i] = synsets[i].word;
	}
	
	// check if we have matching between cache_x and cache_y
	
	
	
      } else {
	// check if we have an entry in cache... otherwise...
	
	
      }
    }
  };
};
