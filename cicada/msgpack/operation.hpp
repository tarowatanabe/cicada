// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_OPERATION__HPP__
#define __CICADA__MSGPACK_OPERATION__HPP__ 1

#include <utils/config.hpp>

#include <cicada/operation.hpp>

#include <cicada/msgpack/hypergraph.hpp>
#include <cicada/msgpack/lattice.hpp>
#include <cicada/msgpack/span_vector.hpp>
#include <cicada/msgpack/alignment.hpp>
#include <cicada/msgpack/dependency.hpp>
#include <cicada/msgpack/sentence_vector.hpp>
#include <cicada/msgpack/ngram_count_set.hpp>
#include <cicada/msgpack/statistics.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>

namespace msgpack
{

  inline
  cicada::Operation::data_type& operator>>(msgpack::object o, cicada::Operation::data_type& v)
  {
    if (o.type != msgpack::type::ARRAY)
      throw msgpack::type_error();
    if (o.via.array.size != 9)
      throw msgpack::type_error();
    
    
    o.via.array.ptr[0].convert(&v.id);
    o.via.array.ptr[1].convert(&v.hypergraph);
    o.via.array.ptr[2].convert(&v.lattice);
    o.via.array.ptr[3].convert(&v.spans);
    o.via.array.ptr[4].convert(&v.alignment);
    o.via.array.ptr[5].convert(&v.dependency);
    o.via.array.ptr[6].convert(&v.targets);
    o.via.array.ptr[7].convert(&v.ngrma_counts);
    o.via.array.ptr[8].convert(&v.statistics);
    
    return v;
  }
    
  template <typename Stream>
  inline
  msgpack::packer<Stream>& operator<<(msgpack::packer<Stream>& o, const cicada::Operation::data_type& v)
  {
    o.pack_array(9);
    o.pack(v.id);
    o.pack(v.hypergraph);
    o.pack(v.lattice);
    o.pack(v.spans);
    o.pack(v.alignment);
    o.pack(v.dependency);
    o.pack(v.targets);
    o.pack(v.ngram_counts);
    o.pack(v.statistics);
    return o;
  }
    
  inline
  void operator<<(msgpack::object::with_zone& o, const cicada::Operation::data_type& v)
  {
    o.type = msgpack::type::ARRAY;
      
    msgpack::object* p = (msgpack::object*) o.zone->malloc(sizeof(msgpack::object) * 9);
    
    o.via.array.ptr = p;
    o.via.array.size = 9;
    
    p[0] = msgpack::object(v.id,           o.zone);
    p[1] = msgpack::object(v.hypergraph,   o.zone);
    p[2] = msgpack::object(v.lattice,      o.zone);
    p[3] = msgpack::object(v.spans,        o.zone);
    p[4] = msgpack::object(v.alignment,    o.zone);
    p[5] = msgpack::object(v.dependency,   o.zone);
    p[6] = msgpack::object(v.targets,      o.zone);
    p[7] = msgpack::object(v.ngram_counts, o.zone);
    p[8] = msgpack::object(v.statistics,   o.zone);
  }
};

#endif
#endif
