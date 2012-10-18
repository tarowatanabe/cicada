// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_SENTENCE_VECTOR__HPP__
#define __CICADA__MSGPACK_SENTENCE_VECTOR__HPP__ 1

#include <utils/config.hpp>

#include <cicada/sentence_vector.hpp>
#include <cicada/msgpack/sentence.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <msgpack/object.hpp>

namespace msgpack
{
  inline
  cicada::SentenceVector& operator>>(msgpack::object o, cicada::SentenceVector& v)
  {
    if (o.type != msgpack::type::ARRAY)
      throw msgpack::type_error();
    
    v.resize(o.via.array.size);
    
    if (o.via.array.size) {
      msgpack::object* p = o.via.array.ptr;
      msgpack::object* const pend = o.via.array.ptr + o.via.array.size;
      
      for (cicada::SentenceVector::iterator it = v.begin(); p != pend; ++ p, ++ it)
	p->convert(&(*it));
    }
    
    return v;
  }
  
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::SentenceVector& v)
  {
    o.pack_array(v.size());
    
    cicada::SentenceVector::const_iterator it_end = v.end();
    for (cicada::SentenceVector::const_iterator it(v.begin()); it != it_end; ++ it) 
      o.pack(*it);
    
    return o;
  }
  
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::SentenceVector& v)
  {
    o.type = msgpack::type::ARRAY;
    
    if (v.empty()) {
      o.via.array.ptr = NULL;
      o.via.array.size = 0;
    } else {
      msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * v.size());
      msgpack::object* const pend = p + v.size();
      
      o.via.array.ptr = p;
      o.via.array.size = v.size();
      
      for (cicada::SentenceVector::const_iterator it(v.begin()); p != pend; ++ p, ++ it)
	*p = msgpack::object(*it, o.zone);
    }
  }
};

#endif
#endif
