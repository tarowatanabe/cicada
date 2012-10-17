// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MSGPACK_TREE_RULE__HPP__
#define __CICADA__MSGPACK_TREE_RULE__HPP__ 1

#include <utils/config.hpp>

#ifdef HAVE_MSGPACK_HPP

#include <vector>

#include <cicada/tree_rule.hpp>

#include <cicada/msgpack/symbol.hpp>

#include <msgpack/object.hpp>
#include <msgpack/type/int.hpp>

namespace cicada
{
  namespace msgpack
  {
    namespace detail
    {
      template <typename BufferLabel, typename BufferSize>
      inline
      void encode(const cicada::TreeRule& rule, BufferLabel& labels, BufferSize& sizes)
      {
	labels.push_back(rule.label);
	sizes.push_back(rule.antecedents.size());
	
	for (size_t i = 0; i != rule.antecedents.size(); ++ i)
	  encode(rule.antecedents[i], labels, sizes);
      }
      
      typedef ::msgpack::object* object_ptr_type;
      
      inline
      void decode(const cicada::TreeRule& rule, object_ptr_type& p)
      {
	uint32_t size = 0;
	
	p->convert(&rule.label);
	++ p;
	p->convert(&size);
	++ p;
	
	rule.antecedents.resize(size);
	
	for (size_t i = 0; i != rule.antecedents.size(); ++ i)
	  decode(rule.antecedents[i], p);
      }
      
    };
    
    inline
    cicada::TreeRule& operator>>(::msgpack::object o, cicada::TreeRule& v)
    {
      if (o.type != ::msgpack::type::ARRAY)
	throw ::msgpack::type_error();

      v.clear();
      
      if (o.via.array.size) {
	::msgpack::object* p = o.via.array.ptr;
	
	detail::decode(v, p);
      }
      
      return v;
    }
    
    template <typename Stream>
    inline
    ::msgpack::packer<Stream>& operator<<(::msgpack::packer<Stream>& o, const cicada::TreeRule& v)
    {
      typedef uint32_t       size_type;
      typedef cicada::Symbol symbol_type;
      typedef std::vector<size_type, std::allocator<size_type> > size_set_type;
      typedef std::vector<symbol_type, std::allocator<symbol_type> > label_set_type;
      
      label_set_type labels;
      size_set_type  sizes;
      
      detail::encode(v, labels, sizes);
      
      o.pack_array(sizes.size() * 2);
      
      for (size_t i = 0; i != sizes.size(); ++ i) {
	o.pack(labels[i]);
	o.pack(sizes[i]);
      }
      
      return 0;
    }
    
    inline
    void operator<<(::msgpack::object::with_zone& o, const cicada::TreeRule& v)
    {
      typedef uint32_t       size_type;
      typedef cicada::Symbol symbol_type;
      typedef std::vector<size_type, std::allocator<size_type> > size_set_type;
      typedef std::vector<symbol_type, std::allocator<symbol_type> > label_set_type;
      
      label_set_type labels;
      size_set_type  sizes;
      
      detail::encode(v, labels, sizes);
      
      o.type = ::msgpack::type::ARRAY;
      
      if (sizes.empty()) {
	o.via.array.ptr = NULL;
	o.via.array.size = 0;
      } else {
	::msgpack::object* p = (::msgpack::object*) o.zone->malloc(sizeof(::msgpack::object) * sizes.size() * 2);
	
	o.via.array.ptr = p;
	o.via.array.size = sizes.size() * 2;
	
	for (size_t i = 0; i != sizes.size(); ++ i) {
	  *p = ::msgpack::object(labels[i], o.zone);
	  ++ p;
	  *p = ::msgpack::object(sizes[i], o.zone);
	  ++ p;
	}
      }
    }
  };
};

#endif
#endif
