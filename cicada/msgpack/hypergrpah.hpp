// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_HYPERGRAPH__HPP__
#define __CICADA__MSGPACK_HYPERGRAPH__HPP__ 1

#include <cicada/hypergraph.hpp>

#include <cicada/msgpack/rule.hpp>
#include <cicada/msgpack/feature_vector.hpp>
#include <cicada/msgpack/attribute_vector.hpp>

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>
#include <msgpack/type/vector.hpp>

namespace cicada
{
  namespace msgpack
  {
    inline
    cicada::HyperGraph::node_type& operator>>(::msgpack::object o, cicada::HyperGraph::node_type& v)
    {
       if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      if (o.via.array.size != 2)
	throw ::msgpack::type_error();
      
      o.via.array.ptr[0].convert(&v.edges);
      o.via.array.ptr[1].convert(&v.id);
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::HyperGraph::node_type& v)
    {
      o.pack_array(2);
      o.pack(v.edges);
      o.pack(v.id);
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::HyperGraph::node_type& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      ::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * 2);
      
      o.via.array.ptr = p;
      o.via.array.size = 2;
      
      p[0] = ::msgpack::object(v.edges, o.zone);
      p[1] = ::msgpack::object(v.id,    o.zone);
    }

    inline
    cicada::HyperGraph::edge_type& operator>>(::msgpack::object o, cicada::HyperGraph::edge_type& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      if (o.via.array.size != 6)
	throw ::msgpack::type_error();
      
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::HyperGraph::edge_type& v)
    {
      o.pack_array(6);
      
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::HyperGraph::edge_type& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      ::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * 6);
      
      o.via.array.ptr = p;
      o.via.array.size = 6;
      
      
    }
    
    inline
    cicada::HyperGraph& operator>>(::msgpack::object o, cicada::HyperGraph& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();
      if (o.via.array.size != 3)
	throw ::msgpack::type_error();
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::HyperGraph& v)
    {
      o.pack_array(3);
      
      
      
      
      return o;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::HyperGraph& v)
    {
      o.type = ::msgpack::type::ARRAY;
      
      ::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * 3);
      
      o.via.array.ptr = p;
      o.via.array.size = 3;
      
      
    }
  };
};

#endif
