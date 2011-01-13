//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "feature.hpp"

namespace cicada
{
  Feature::mutex_type    Feature::__mutex;
  
  class __feature_set_instance
  {
  public:
    __feature_set_instance() : instance(0) {}
    ~__feature_set_instance() { if (instance) delete instance; }
    
    Feature::feature_set_type* instance;
  };
  
  static boost::once_flag      __features_once = BOOST_ONCE_INIT;
  static __feature_set_instance __features_instance;
  
  static void __features_init()
  {
    __features_instance.instance = new Feature::feature_set_type();
  }
  
  Feature::feature_set_type& Feature::__features()
  {
    boost::call_once(__features_once, __features_init);

    return *__features_instance.instance;
  }
  
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
    
};
