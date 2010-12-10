#include <cstdlib>
#include <stdexcept>

#include "cicada/matcher/wordnet.hpp"

#include "utils/spinlock.hpp"

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include "wn/wn.h"

namespace cicada
{
  namespace matcher
  {
    // this is the global guard!
    static utils::spinlock __wordnet_mutex;
    
    
    
    
    
    static boost::once_flag __wordnet_init_once = BOOST_ONCE_INIT;
    static std::string      __wordnet_path;
    
    static void __wordnet_init()
    {
      if (wninit()) {
	if (__wordnet_path.empty() || ! boost::filesystem::exists(__wordnet_path))
	  throw std::runtime_error("no wordnet database? check WNHOME, WNSEARCHDIR or supply path with file=\"path to db\"");
	
	setenv("WNHOME", __wordnet_path.c_str(), 1);
	if (wninit()) {
	  setenv("WNSEARCHDIR", __wordnet_path.c_str(), 1);
	  
	  if (wninit())
	    throw std::runtime_error("no wordnet database? " + __wordnet_path);
	}
      }
    }
    
    void WordNet::initialize(const std::string& path)
    {
      // we will assure locking + once to make sure this will be called only once...
      utils::spinlock::scoped_lock lock(__wordnet_mutex);
      
      __wordnet_path = path;
      
      boost::call_once(__wordnet_init_once, __wordnet_init);
    }
  };
};
