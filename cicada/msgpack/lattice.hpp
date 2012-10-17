// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_LATTICE__HPP__
#define __CICADA__MSGPACK_LATTICE__HPP__ 1

#include <cicada/lattice.hpp>

#include <cicada/msgpack/symbol.hpp>
#include <cicada/msgpack/feature_vector.hpp>

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>

namespace cicada
{
  namespace msgpack
  {
    inline
    cicada::Lattice& operator>>(::msgpack::object o, cicada::Lattice::arc_type& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      if (o.via.array.size != 3)
	throw ::msgpack::type_error();
      
      o.via.array.ptr[0].convert(&v.feature);
      o.via.array.ptr[1].convert(&v.label);
      o.via.array.ptr[2].convert(&v.distance);
      
      return v;
    }

    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::Lattice::arc_type& v)
    {
      o.pack_array(3);
      o.pack(v.features);
      o.pack(v.label);
      o.pack(v.distance);
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::Lattice::arc_type& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      ::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * 3);
      
      o.via.array.ptr = p;
      o.via.array.size = 3;
      
      p[0] = ::msgpack::object(v.features, o.zone);
      p[1] = ::msgpack::object(v.label,    o.zone);
      p[2] = ::msgpack::object(v.distance, o.zone);
    }

    inline
    cicada::Lattice& operator>>(::msgpack::object o, cicada::Lattice::arc_set_type& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      
      v.resize(o.via.array.size);
      
      if (o.via.array.size) {
	::msgpack::object* p = o.via.array.ptr;
	::msgpack::object* const pend = o.via.array.ptr + o.via.array.size;
	
	for (cicada::Lattice::iterator it = v.begin(); p != pend; ++ p, ++ it)
	  p->convert(&(*it));
      }
      
      return v;
    }

    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::Lattice::arc_set_type& v)
    {
      o.pack_array(v.size());
      
      cicada::Lattice::const_iterator it_end = v.end();
      for (cicada::Lattice::const_iterator it(v.begin()); it != it_end; ++ it) 
	o.pack(*it);
      
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::Lattice::arc_set_type& v)
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
	
	for (cicada::Lattice::const_iterator it(v.begin()); p != pend; ++ p, ++ it)
	  *p = ::msgpack::object(*it, o.zone);
      }
    }
    
    inline
    cicada::Lattice& operator>>(::msgpack::object o, cicada::Lattice& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      
      v.resize(o.via.array.size);
      
      if (o.via.array.size) {
	::msgpack::object* p = o.via.array.ptr;
	::msgpack::object* const pend = o.via.array.ptr + o.via.array.size;
	
	for (cicada::Lattice::iterator it = v.begin(); p != pend; ++ p, ++ it)
	  p->convert(&(*it));
      }
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::Lattice& v)
    {
      o.pack_array(v.size());
      
      cicada::Lattice::const_iterator it_end = v.end();
      for (cicada::Lattice::const_iterator it(v.begin()); it != it_end; ++ it) 
	o.pack(*it);
      
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::Lattice& v)
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
	
	for (cicada::Lattice::const_iterator it(v.begin()); p != pend; ++ p, ++ it)
	  *p = ::msgpack::object(*it, o.zone);
      }
    }
  };
};

#endif
