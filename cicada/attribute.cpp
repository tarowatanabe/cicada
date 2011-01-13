//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "attribute.hpp"

namespace cicada
{

  Attribute::mutex_type    Attribute::__mutex;

  struct AttributeImpl
  {
    typedef Attribute::attribute_set_type attribute_set_type;
    typedef Attribute::attribute_map_type attribute_map_type;

    static boost::once_flag once;

    static attribute_set_type* attributes;
    
#ifdef HAVE_TLS
    static __thread attribute_map_type* attribute_maps_tls;
#endif
    static boost::thread_specific_ptr<attribute_map_type> attribute_maps;

    static void initialize()
    {
      attributes = new attribute_set_type();
    }
  };

  boost::once_flag AttributeImpl::once = BOOST_ONCE_INIT;
  
  AttributeImpl::attribute_set_type* AttributeImpl::attributes = 0;
  
#ifdef HAVE_TLS
  __thread AttributeImpl::attribute_map_type* AttributeImpl::attribute_maps_tls = 0;
#endif
  boost::thread_specific_ptr<AttributeImpl::attribute_map_type> AttributeImpl::attribute_maps;
  
  Attribute::attribute_set_type& Attribute::__attributes()
  {
    boost::call_once(AttributeImpl::once, AttributeImpl::initialize);
    
    return *AttributeImpl::attributes;
  }

  Attribute::attribute_map_type& Attribute::__attribute_maps()
  {
#ifdef HAVE_TLS
    if (! AttributeImpl::attribute_maps_tls) {
      AttributeImpl::attribute_maps.reset(new attribute_map_type());
      AttributeImpl::attribute_maps->reserve(allocated());
      
      AttributeImpl::attribute_maps_tls = AttributeImpl::attribute_maps.get();
    }
    
    return *AttributeImpl::attribute_maps_tls;
#else
    if (! AttributeImpl::attribute_maps.get()) {
      AttributeImpl::attribute_maps.reset(new attribute_map_type());
      AttributeImpl::attribute_maps->reserve(allocated());
    }
    
    return *AttributeImpl::attribute_maps;
#endif
  }
};
