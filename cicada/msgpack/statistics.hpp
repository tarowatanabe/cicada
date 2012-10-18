// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_STATISTICS__HPP__
#define __CICADA__MSGPACK_STATISTICS__HPP__ 1

#include <utils/config.hpp>

#include <cicada/statistics.hpp>
#include <cicada/msgpack/attribute.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>
#include <msgpack/type/float.hpp>

namespace msgpack
{
  inline
  cicada::Statistics::stat_type& operator>>(msgpack::object o, cicada::Statistics::stat_type& v)
  {
    if (o.type != msgpack::type::ARRAY)
      throw msgpack::type_error();
    if (o.via.array.size != 6)
      throw msgpack::type_error();
      
    o.via.array.ptr[0].convert(&v.count);
    o.via.array.ptr[1].convert(&v.node);
    o.via.array.ptr[2].convert(&v.edge);
    o.via.array.ptr[3].convert(&v.user_time);
    o.via.array.ptr[4].convert(&v.cpu_time);
    o.via.array.ptr[5].convert(&v.thread_time);
      
    return v;
  }
    
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::Statistics::stat_type& v)
  {
    o.pack_array(6);
    o.pack(v.count);
    o.pack(v.node);
    o.pack(v.edge);
    o.pack(v.user_time);
    o.pack(v.cpu_time);
    o.pack(v.thread_time);
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::Statistics::stat_type& v)
  {
    o.type = msgpack::type::ARRAY;
      
    msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * 6);
      
    o.via.array.ptr = p;
    o.via.array.size = 6;
      
    p[0] = msgpack::object(v.count,       o.zone);
    p[1] = msgpack::object(v.node,        o.zone);
    p[2] = msgpack::object(v.edge,        o.zone);
    p[3] = msgpack::object(v.user_time,   o.zone);
    p[4] = msgpack::object(v.cpu_time,    o.zone);
    p[5] = msgpack::object(v.thread_time, o.zone);
  }


  inline
  cicada::Statistics& operator>>(msgpack::object o, cicada::Statistics& v)
  {
    if (o.type != msgpack::type::MAP)
      throw msgpack::type_error();
    
    msgpack::object_kv* p(o.via.map.ptr);
    msgpack::object_kv* const pend(o.via.map.ptr + o.via.map.size);
    
    v.clear();
    
    cicada::Statistics::attribute_type attr
    for (/**/; p != pend; ++ p) {
      p->key.convert(&attr);
      p->val.convert(&v[attr]);
    }
    
    return v;
  }
  
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::Statistics& v)
  {
    o.pack_map(v.size());
    
    typename cicada::Statistics::const_iterator it_end = v.end();
    for (typename cicada::Statistics::const_iterator it(v.begin()); it != it_end; ++ it) {
      o.pack(it->first);
      o.pack(it->second);
    }
    
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::Statistics& v)
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
	
      for (typename cicada::Statistics::const_iterator it(v.begin()); p != pend; ++ p, ++ it) {
	p->key = msgpack::object(it->first, o.zone);
	p->val = msgpack::object(it->second, o.zone);
      }
    }
  }
};

#endif
#endif
