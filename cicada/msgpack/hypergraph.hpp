// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_HYPERGRAPH__HPP__
#define __CICADA__MSGPACK_HYPERGRAPH__HPP__ 1

#include <utils/config.hpp>

#include <cicada/hypergraph.hpp>

#include <cicada/msgpack/rule.hpp>
#include <cicada/msgpack/feature_vector.hpp>
#include <cicada/msgpack/attribute_vector.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>
#include <msgpack/type/vector.hpp>

namespace msgpack
{
  namespace detail
  {
    template <typename Storage>
    inline
    void decode_array(msgpack::object& o, Storage& storage)
    {
      if (o.type != msgpack::type::ARRAY)
	throw msgpack::type_error();
	
      storage.resize(o.via.array.size);
	
      if (o.via.array.size) {
	msgpack::object* p = o.via.array.ptr;
	msgpack::object* const pend = o.via.array.ptr + o.via.array.size;
	  
	for (typename Storage::iterator it = storage.begin(); p != pend; ++ p, ++ it)
	  p->convert(&(*it));
      }
    }

    template <typename Stream, typename Iterator>
    inline
    void encode_array(msgpack::packer<Stream>& o, Iterator first, Iterator last)
    {
      o.pack_array(std::distance(first, last));
      for (/**/; first != last; ++ first)
	o.pack(*first);
    }
      
    template <typename Iterator>
    inline
    void encode_array(msgpack::object& t, msgpack::object::with_zone& o, Iterator first, Iterator last)
    {
      const size_t size = std::distance(first, last);
	
      t.type = msgpack::type::ARRAY;
      if (! size) {
	t.via.array.ptr = NULL;
	t.via.array.size = 0;
      } else {
	msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * size);
	msgpack::object* const pend = p + size;
	  
	t.via.array.ptr  = p;
	t.via.array.size = size;
	  
	for (/**/; p != pend; ++ p, ++ first)
	  *p = msgpack::object(*first, o.zone);
      }
    }
  };

  inline
  cicada::HyperGraph::node_type& operator>>(msgpack::object o, cicada::HyperGraph::node_type& v)
  {
    if (o.type != msgpack::type::ARRAY)
      throw msgpack::type_error();
    if (o.via.array.size != 2)
      throw msgpack::type_error();
    
    detail::decode_array(o.via.array.ptr[0], v.edges);
    o.via.array.ptr[1].convert(&v.id);
      
    return v;
  }
    
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::HyperGraph::node_type& v)
  {
    o.pack_array(2);
    detail::encode_array(o, v.edges.begin(), v.edges.end());
    o.pack(v.id);
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::HyperGraph::node_type& v)
  {
    o.type = msgpack::type::ARRAY;
      
    msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * 2);
      
    o.via.array.ptr = p;
    o.via.array.size = 2;
    
    detail::encode_array(p[0], o, v.edges.begin(), v.edges.end());
    p[1] = msgpack::object(v.id,    o.zone);
  }

  inline
  cicada::HyperGraph::edge_type& operator>>(msgpack::object o, cicada::HyperGraph::edge_type& v)
  {
    if (o.type != msgpack::type::ARRAY)
      throw msgpack::type_error();
    if (o.via.array.size != 6)
      throw msgpack::type_error();
      
    // head
    o.via.array.ptr[0].convert(&v.head);
    // tails
    detail::decode_array(o.via.array.ptr[1], v.tails);
      
    // features/attributes
    o.via.array.ptr[2].convert(&v.features);
    o.via.array.ptr[3].convert(&v.attributes);
      
    // rule/id
    cicada::Rule rule;
    o.via.array.ptr[4].convert(&rule);
    o.via.array.ptr[5].convert(&v.id);
      
    v.rule = cicada::Rule::create(rule);
      
    return v;
  }
    
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::HyperGraph::edge_type& v)
  {
    o.pack_array(6);
    o.pack(v.head);
      
    // tails...
    detail::encode_array(o, v.tails.begin(), v.tails.end());
      
    o.pack(v.features);
    o.pack(v.attributes);
    o.pack(*v.rule);
    o.pack(v.id);
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::HyperGraph::edge_type& v)
  {
    o.type = msgpack::type::ARRAY;
      
    msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * 6);
      
    o.via.array.ptr = p;
    o.via.array.size = 6;
      
    p[0] = msgpack::object(v.head,       o.zone);
      
    detail::encode_array(p[1], o, v.tails.begin(), v.tails.end());
      
    p[2] = msgpack::object(v.features,   o.zone);
    p[3] = msgpack::object(v.attributes, o.zone);
    p[4] = msgpack::object(*v.rule,      o.zone);
    p[5] = msgpack::object(v.id,         o.zone);
  }
    
  inline
  cicada::HyperGraph& operator>>(msgpack::object o, cicada::HyperGraph& v)
  {
    if (o.type != msgpack::type::ARRAY)
      throw msgpack::type_error();
    if (o.via.array.size != 3)
      throw msgpack::type_error();
      
    detail::decode_array(o.via.array.ptr[0], v.nodes);
    detail::decode_array(o.via.array.ptr[1], v.edges);
    o.via.array.ptr[2].convert(&v.goal);
      
    return v;
  }
    
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::HyperGraph& v)
  {
    o.pack_array(3);
      
    detail::encode_array(o, v.nodes.begin(), v.nodes.end());
    detail::encode_array(o, v.edges.begin(), v.edges.end());
    o.pack(v.goal);
      
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::HyperGraph& v)
  {
    o.type = msgpack::type::ARRAY;
      
    msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * 3);
      
    o.via.array.ptr = p;
    o.via.array.size = 3;
      
    detail::encode_array(p[0], o, v.nodes.begin(), v.nodes.end());
    detail::encode_array(p[1], o, v.edges.begin(), v.edges.end());
    p[2] = msgpack::object(v.goal,   o.zone);
  }
};

#endif
#endif
