//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "attribute.hpp"

namespace cicada
{
  Attribute::attribute_map_type& Attribute::__attribute_maps()
  {
#ifdef HAVE_TLS
    static __thread attribute_map_type* __maps_tls = 0;
    static boost::thread_specific_ptr<attribute_map_type> __maps;
    
    if (! __maps_tls) {
      __maps.reset(new attribute_map_type());
      __maps->reserve(allocated());
      
      __maps_tls = __maps.get();
    }
      
    return *__maps_tls;
#else
    static boost::thread_specific_ptr<attribute_map_type> __maps;
      
    if (! __maps.get()) {
      __maps.reset(new attribute_map_type());
      __maps->reserve(allocated());
    }
    
    return *__maps;
#endif
  }
  
  Attribute::mutex_type    Attribute::__mutex;
  
};
