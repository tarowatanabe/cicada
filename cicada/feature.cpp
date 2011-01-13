//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "feature.hpp"

namespace cicada
{

  struct FeatureImpl
  {
    typedef Feature::feature_map_type feature_map_type;
  };
  
  Feature::mutex_type    Feature::__mutex;
  
  Feature::feature_set_type& Feature::__features()
  {
    static feature_set_type features;
    
    return features;
  }

#ifdef HAVE_TLS
  static __thread FeatureImpl::feature_map_type* feature_maps_tls = 0;
#endif
  static boost::thread_specific_ptr<FeatureImpl::feature_map_type> feature_maps;

  Feature::feature_map_type& Feature::__feature_maps()
  {
#ifdef HAVE_TLS
    if (! feature_maps_tls) {
      feature_maps.reset(new feature_map_type());
      feature_maps->reserve(allocated());
      
      feature_maps_tls = feature_maps.get();
    }
    
    return *feature_maps_tls;
#else
    if (! feature_maps.get()) {
      feature_maps.reset(new feature_map_type());
      feature_maps->reserve(allocated());
    }
    
    return *feature_maps;
#endif
  }
};
