// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_NGRAM_COUNT_SET__HPP__
#define __CICADA__MSGPACK_NGRAM_COUNT_SET__HPP__ 1

#include <utils/config.hpp>

#include <cicada/ngram_count_set.hpp>
#include <cicada/msgpack/symbol_vector.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>

namespace msgpack
{
  inline
  cicada::NGramCountSet& operator>>(msgpack::object o, cicada::NGramCountSet& v)
  {
    if (o.type != msgpack::type::MAP)
      throw msgpack::type_error();
    
    msgpack::object_kv* p(o.via.map.ptr);
    msgpack::object_kv* const pend(o.via.map.ptr + o.via.map.size);
    
    v.clear();
    
    cicada::NGramCountSet::ngram_type ngram;
    for (/**/; p != pend; ++ p) {
      p->key.convert(&ngram);
      p->val.convert(&v[ngram]);
    }
      
    return v;
  }
    
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::NGramCountSet& v)
  {
    o.pack_map(v.size());
    
    cicada::NGramCountSet::const_iterator it_end = v.end();
    for (cicada::NGramCountSet::const_iterator it(v.begin()); it != it_end; ++ it) {
      o.pack(it->first);
      o.pack(it->second);
    }
      
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::NGramCountSet& v)
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
	
      for (cicada::NGramCountSet::const_iterator it(v.begin()); p != pend; ++ p, ++ it) {
	p->key = msgpack::object(it->first, o.zone);
	p->val = msgpack::object(it->second, o.zone);
      }
    }
  }
};

#endif
#endif
