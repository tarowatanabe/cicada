// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_FEATURE_VECTOR_LINEAR__HPP__
#define __CICADA__MSGPACK_FEATURE_VECTOR_LINEAR__HPP__ 1

#include <utils/config.hpp>

#include <cicada/feature_vector_linear.hpp>
#include <cicada/msgpack/feature.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>
#include <msgpack/type/float.hpp>

namespace msgpack
{
  template <typename T, typename A>
  inline
  cicada::FeatureVectorLinear<T,A>& operator>>(msgpack::object o, cicada::FeatureVectorLinear<T,A>& v)
  {
    if (o.type != msgpack::type::MAP)
      throw msgpack::type_error();
    
    msgpack::object_kv* p(o.via.map.ptr);
    msgpack::object_kv* const pend(o.via.map.ptr + o.via.map.size);
    
    v.clear();
    
    cicada::Feature feature;
    for (/**/; p != pend; ++ p) {
      p->key.convert(&feature);
      p->val.convert(&v[feature]);
    }
      
    return v;
  }
    
  template <typename Stream, typename T, typename A>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::FeatureVectorLinear<T,A>& v)
  {
    o.pack_map(v.size());
      
    typename cicada::FeatureVectorLinear<T,A>::const_iterator it_end = v.end();
    for (typename cicada::FeatureVectorLinear<T,A>::const_iterator it(v.begin()); it != it_end; ++ it) {
      o.pack(it->first);
      o.pack(it->second);
    }
      
    return o;      
  }
    
  template <typename T, typename A>
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::FeatureVectorLinear<T,A>& v)
  {
    o.type = msgpack::type::MAP;
	
    if (v.empty()) {
      o.via.map.ptr  = NULL;
      o.via.map.size = 0;
    } else {
      msgpack::object_kv* p = (msgpack::object_kv*) o.zone->malloc(sizeof(msgpack::object_kv) * v.size());
      msgpack::object_kv* const pend = p + v.size();
      
      o.via.map.ptr  = p;
      o.via.map.size = v.size();
      
      for (typename cicada::FeatureVectorLinear<T,A>::const_iterator it(v.begin()); p != pend; ++ p, ++ it) {
	p->key = msgpack::object(it->first, o.zone);
	p->val = msgpack::object(it->second, o.zone);
      }
    }
  }
};

#endif
#endif
