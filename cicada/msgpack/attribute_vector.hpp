// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_ATTRIBUTE_VECTOR__HPP__
#define __CICADA__MSGPACK_ATTRIBUTE_VECTOR__HPP__ 1

#include <utils/config.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <cicada/attribute_vector.hpp>

#include <cicada/msgpack/attribute.hpp>

#include <msgpack/object.hpp>
#include <msgpack/type/string.hpp>
#include <msgpack/type/int.hpp>
#include <msgpack/type/float.hpp>

namespace cicada
{
  namespace msgpack
  {
    namespace detail
    {
      template <typename Object>
      struct object_visotor : public boost::static_visitor<bool>
      {
	typedef cicada::AttributeVector::int_type int_type;
	typedef cicada::AttributeVector::float_type float_type;
	typedef cicada::AttributeVector::string_type string_type;

	object_visotor(Object& __o) : o(__o) {}
	
	bool operator()(const int_type& x) const { o << x; return true; }
	bool operator()(const float_type& x) const { o << x; return true; }
	bool operator()(const string_type& x) const { o << x; return true; }
	
	Object& o;
      };
      
    };

    inline
    cicada::AttributeVector& operator>>(::msgpack::object o, cicada::AttributeVector::data_type& v)
    {
      switch (o.type) {
      case ::msgpack::type::POSITIVE_INTEGER:
      case ::msgpack::type::NEGATIVE_INTEGER:
	cicada::AttributeVector::int_type vi;
	o >> vi;
	v = vi;
	break;
      case ::msgpack::type::DOUBLE:
	cicada::AttributeVector::float_type vf;
	o >> vf;
	v = vf;
	break;
      case ::msgpack::type::RAW:
	cicada::AttributeVector::string_type vs;
	o >> vs;
	v = vs;
	break;
      default:
	throw ::msgpack::type_error();
      }
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::AttributeVector::data_type& v)
    {
      boost::apply_visitor(detail::object_visotor<::msgpack::packer<Stream> >(o), v);
      
      return o;      
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::AttributeVector::data_type& v)
    {
      boost::apply_visitor(detail::object_visotor<::msgpack::object::with_zone>(o), v);
    }
    
    inline
    void operator<< (::msgpack::object& o, const cicada::AttributeVector::data_type& v)
    {
      boost::apply_visitor(detail::object_visotor<::msgpack::object>(o), v);
    }
    
    inline
    cicada::AttributeVector& operator>>(::msgpack::object o, cicada::AttributeVector& v)
    {
      if (o.type != ::msgpack::type::MAP)
	throw ::msgpack::type_error();
            
      ::msgpack::object_kv* p(o.via.map.ptr);
      ::msgpack::object_kv* const pend(o.via.map.ptr + o.via.map.size);
      
      v.clear();
      v.rehash(o.via.map.size);
      
      cicada::Attribute attribute;
      for (/**/; p != pend; ++ p) {
	p->key.convert(&attribute);
	p->val.convert(&v[attribute]);
      }
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::AttributeVector& v)
    {
      o.pack_map(v.size());
      
      cicada::AttributeVector::const_iterator it_end = v.end();
      for (cicada::AttributeVector::const_iterator it(v.begin()); it != it_end; ++ it) {
	o.pack(it->first);
	o.pack(it->second);
      }
      
      return o;      
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::AttributeVector& v)
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
	
	for (cicada::AttributeVector::const_iterator it(v.begin()); p != pend; ++ p, ++ it) {
	  p->key = ::msgpack::object(it->first, o.zone);
	  p->val = ::msgpack::object(it->second, o.zone);
	}
      }
    }
  };
};

#endif
#endif
