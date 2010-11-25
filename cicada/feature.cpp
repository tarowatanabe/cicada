//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "feature.hpp"

namespace cicada
{
  Feature::feature_map_type& Feature::__feature_maps()
  {
#ifdef HAVE_TLS
    static __thread feature_map_type* __maps_tls = 0;
    static boost::thread_specific_ptr<feature_map_type> __maps;
    
    if (! __maps_tls) {
      __maps.reset(new feature_map_type());
      __maps->reserve(allocated());
      
      __maps_tls = __maps.get();
    }
      
    return *__maps_tls;
#else
    static boost::thread_specific_ptr<feature_map_type> __maps;
      
    if (! __maps.get()) {
      __maps.reset(new feature_map_type());
      __maps->reserve(allocated());
    }
    
    return *__maps;
#endif
  }
  
  Feature::mutex_type    Feature::__mutex;
  
};
