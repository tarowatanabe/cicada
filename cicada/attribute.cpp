//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "attribute.hpp"

namespace cicada
{
  struct AttributeImpl
  {
    typedef Attribute::attribute_map_type attribute_map_type;
  };
  
  Attribute::mutex_type    Attribute::__mutex;
  
  Attribute::attribute_set_type& Attribute::__attributes()
  {
    static attribute_set_type attributes;
    
    return attributes;
  }
  
#ifdef HAVE_TLS
  static __thread AttributeImpl::attribute_map_type* attribute_maps_tls = 0;
#endif
  static boost::thread_specific_ptr<AttributeImpl::attribute_map_type> attribute_maps;
  
  Attribute::attribute_map_type& Attribute::__attribute_maps()
  {

#ifdef HAVE_TLS
    if (! attribute_maps_tls) {
      attribute_maps.reset(new attribute_map_type());
      attribute_maps->reserve(allocated());
      
      attribute_maps_tls = attribute_maps.get();
    }
    
    return *attribute_maps_tls;
#else
    if (! attribute_maps.get()) {
      attribute_maps.reset(new attribute_map_type());
      attribute_maps->reserve(allocated());
    }
    
    return *attribute_maps;
#endif
  }
};
