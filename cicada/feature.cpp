//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "feature.hpp"

namespace cicada
{

  Feature::mutex_type    Feature::__mutex;

  struct FeatureImpl
  {
    typedef Feature::feature_set_type feature_set_type;
    typedef Feature::feature_map_type feature_map_type;

    static boost::once_flag once;

    static feature_set_type* features;
    
#ifdef HAVE_TLS
    static __thread feature_map_type* feature_maps_tls;
#endif
    static boost::thread_specific_ptr<feature_map_type> feature_maps;

    static void initialize()
    {
      features = new feature_set_type();
    }
  };

  boost::once_flag FeatureImpl::once = BOOST_ONCE_INIT;
  
  FeatureImpl::feature_set_type* FeatureImpl::features = 0;
  
#ifdef HAVE_TLS
  __thread FeatureImpl::feature_map_type* FeatureImpl::feature_maps_tls = 0;
#endif
  
  Feature::feature_set_type& Feature::__features()
  {
    boost::call_once(FeatureImpl::once, FeatureImpl::initialize);
    
    return *FeatureImpl::features;
  }

  Feature::feature_map_type& Feature::__feature_maps()
  {
#ifdef HAVE_TLS
    if (! FeatureImpl::feature_maps_tls) {
      FeatureImpl::feature_maps.reset(new feature_map_type());
      FeatureImpl::feature_maps->reserve(allocated());
      
      FeatureImpl::feature_maps_tls = FeatureImpl::feature_maps.get();
    }
    
    return *FeatureImpl::feature_maps_tls;
#else
    if (! FeatureImpl::feature_maps.get()) {
      FeatureImpl::feature_maps.reset(new feature_map_type());
      FeatureImpl::feature_maps->reserve(allocated());
    }
    
    return *FeatureImpl::feature_maps;
#endif
  }
};
