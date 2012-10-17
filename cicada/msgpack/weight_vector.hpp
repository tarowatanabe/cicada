// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_WEIGHT_VECTOR__HPP__
#define __CICADA__MSGPACK_WEIGHT_VECTOR__HPP__ 1

#include <utils/config.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <cicada/weight_vector.hpp>

#include <cicada/msgpack/feature.hpp>

#include <msgpack/object.hpp>
#include <msgpack/type/float.hpp>

namespace cicada
{
  namespace msgpack
  {
    template <typename T, typename A>
    inline
    cicada::WeightVector<T,A>& operator>>(::msgpack::object o, cicada::WeightVector<T,A>& v)
    {
      if (o.type != ::msgpack::type::MAP)
	throw ::msgpack::type_error();
            
      ::msgpack::object_kv* p(o.via.map.ptr);
      ::msgpack::object_kv* const pend(o.via.map.ptr + o.via.map.size);
      
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
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::WeightVector<T,A>& v)
    {
      o.pack_map(v.size());
      
      for (cicada::Feature::id_type id = 0; id != v.size(); ++ id) {
	const cicada::Feature feature(id);
	
	o.pack(feature);
	o.pack(v[feature]);
      }
      
      return o;      
    }
    
    template <typename T, typename A>
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::WeightVector<T,A>& v)
    {
      o.type = ::msgpack::type::MAP;
      
      if (v.empty()) {
	o.via.map.ptr  = NULL;
	o.via.map.size = 0;
      } else {
	::msgpack::object_kv* p = (::msgpack::object_kv*) o.zone->malloc(sizeof(::msgpack::object_kv) * v.size());
	::msgpack::object_kv* const pend = p + v.size();
	
	o.via.map.ptr  = p;
	o.via.map.size = v.size();
	
	for (cicada::Feature::id_type id = 0; id != v.size(); ++ id, ++ p) {
	  const cicada::Feature feature(id);
	  
	  p->key = ::msgpack::object(feature, o.zone);
	  p->val = ::msgpack::object(v[feature], o.zone);
	}
      }
    }
  };
};

#endif
#endif
