// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_RULE__HPP__
#define __CICADA__MSGPACK_RULE__HPP__ 1

#include <utils/config.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <cicada/rule.hpp>

#include <cicada/msgpack/symbol.hpp>
#include <cicada/msgpack/symbol_vector.hpp>

#include <msgpack/object.hpp>

namespace cicada
{
  namespace msgpack
  {
    inline
    cicada::Rule& operator>>(::msgpack::object o, cicada::Rule& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      if (o.via.array.size != 2)
	throw ::msgpack::type_error();
      
      o.via.array.ptr[0].convert(&v.lhs);
      o.via.array.ptr[1].convert(&v.rhs);
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::Rule& v)
    {
      o.pack_array(2);
      o.pack(v.lhs);
      o.pack(v.rhs);
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::Rule& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      ::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * 2);
      
      o.via.array.ptr = p;
      o.via.array.size = 2;
      
      p[0] = ::msgpack::object(v.lhs, o.zone);
      p[1] = ::msgpack::object(v.rhs, o.zone);
    }
  };
};

#endif
#endif
