//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "feature.hpp"

namespace cicada
{

  struct FeatureImpl
  {
    typedef Feature::feature_set_type feature_set_type;

    static boost::once_flag once;

    static feature_set_type* features;
    
    static void initialize()
    {
      features = new feature_set_type();
    }
  };

  boost::once_flag FeatureImpl::once = BOOST_ONCE_INIT;
  
  FeatureImpl::feature_set_type* FeatureImpl::features = 0;
  
  Feature::mutex_type    Feature::__mutex;
  
  Feature::feature_set_type& Feature::__features()
  {
    static feature_set_type features;
    
    return features;

#if 0
    boost::call_once(FeatureImpl::once, FeatureImpl::initialize);
    
    return *FeatureImpl::features;
#endif
  }

  Feature::feature_map_type& Feature::__feature_maps()
  {
#ifdef HAVE_TLS
    static __thread feature_map_type* feature_maps_tls = 0;
#endif
    static boost::thread_specific_ptr<feature_map_type> feature_maps;

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
