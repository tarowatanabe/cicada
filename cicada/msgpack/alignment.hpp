// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_ALIGNMENT__HPP__
#define __CICADA__MSGPACK_ALIGNMENT__HPP__ 1

#include <cicada/alignment.hpp>

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>

namespace cicada
{
  namespace msgpack
  {
    inline
    cicada::Alignment::point_type& operator>>(::msgpack::object o, cicada::Alignment::point_type& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      if (o.via.array.size != 2)
	throw ::msgpack::type_error();
      
      o.via.array.ptr[0].convert(&v.source);
      o.via.array.ptr[1].convert(&v.target);
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::Alignment::point_type& v)
    {
      o.pack_array(2);
      o.pack(v.source);
      o.pack(v.target);
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::Alignment::point_type& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      ::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * 2);
      
      o.via.array.ptr = p;
      o.via.array.size = 2;
      
      p[0] = ::msgpack::object(v.source, o.zone);
      p[1] = ::msgpack::object(v.target, o.zone);
    }


    inline
    cicada::Alignment& operator>>(::msgpack::object o, cicada::Alignment& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      
      v.resize(o.via.array.size);
      
      if (o.via.array.size) {
	::msgpack::object* p = o.via.array.ptr;
	::msgpack::object* const pend = o.via.array.ptr + o.via.array.size;
	
	for (cicada::Alignment::iterator it = v.begin(); p != pend; ++ p, ++ it)
	  p->convert(&(*it));
      }
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::Alignment& v)
    {
      o.pack_array(v.size());
      
      cicada::Alignment::const_iterator it_end = v.end();
      for (cicada::Alignment::const_iterator it(v.begin()); it != it_end; ++ it) 
	o.pack(*it);
      
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::Alignment& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      if (v.empty()) {
	o.via.array.ptr = NULL;
	o.via.array.size = 0;
      } else {
	::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * v.size());
	::msgpack::object* const pend = p + v.size();
	
	o.via.array.ptr = p;
	o.via.array.size = v.size();
	
	for (cicada::Alignment::const_iterator it(v.begin()); p != pend; ++ p, ++ it)
	  *p = ::msgpack::object(*it, o.zone);
      }
    }
  };
};

#endif
