// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__ATTRIBUTE_VECTOR_CODEC__HPP__
#define __CICADA__ATTRIBUTE_VECTOR_CODEC__HPP__ 1

#include <memory>
#include <utility>
#include <algorithm>
#include <iterator>

#include <cicada/attribute.hpp>
#include <cicada/attribute_vector.hpp>

#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>
#include <utils/simple_vector.hpp>

namespace cicada
{
  // a compact attribute vector representation which do not allow any modification, and
  // uses input-iterator, not bidirectional/random-access iterator
  // we use double as our underlying stroage..
  
  struct AttributeVectorCODEC
  {
    typedef cicada::Attribute       attribute_type;
    typedef cicada::AttributeVector attribute_set_type;
    
    typedef char byte_type;
    typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;
    
    struct __attribute_vector_attribute_codec
    {
      typedef cicada::Attribute attribute_type;
      typedef uint8_t byte_type;

      static size_t encode(const attribute_type::id_type& value)
      {
	return utils::__byte_aligned_code<attribute_type::id_type, sizeof(attribute_type::id_type)>::byte_size(value);
      }

      static size_t encode(const attribute_type& value)
      {
	return encode(value.id());
      }

      static size_t encode(byte_type* buffer, const attribute_type::id_type& value)
      {
	return utils::byte_aligned_encode(value, reinterpret_cast<char*>(buffer));
      }
      
      static size_t encode(byte_type* buffer, const attribute_type& value)
      {
	return encode(buffer, value.id());
      }
      
      static size_t decode(const byte_type* buffer, attribute_type::id_type& value)
      {
	return utils::byte_aligned_decode(value, reinterpret_cast<const char*>(buffer));
      }
      static size_t decode(const byte_type* buffer, attribute_type& value)
      {
	attribute_type::id_type value_id = 0;
	const size_t ret = utils::byte_aligned_decode(value_id, reinterpret_cast<const char*>(buffer));
	value = attribute_type(value_id);
	return ret;
      }
    };
  
    struct __attribute_vector_data_codec
    {    
      typedef uint8_t byte_type;
      typedef AttributeVector attribute_set_type;
      
      static const uint8_t mask_float    = 1 << (4 + 0);
      static const uint8_t mask_unsigned = 1 << (4 + 1);
      static const uint8_t mask_signed   = 1 << (4 + 2);
      static const uint8_t mask_string   = 1 << (4 + 3);
      static const uint8_t mask_size     = 0x0f;

      template <typename __Tp>
      static byte_type* cast(__Tp& x)
      {
	return (byte_type*) &x;
      }
    
      template <typename __Tp>
      static const byte_type* cast(const __Tp& x)
      {
	return (const byte_type*) &x;
      }
    
      static 
      size_t byte_size(const uint64_t& x)
      {
	return (1 
		+ bool(x & 0xffffffffffffff00ull)
		+ bool(x & 0xffffffffffff0000ull)
		+ bool(x & 0xffffffffff000000ull)
		+ bool(x & 0xffffffff00000000ull)
		+ bool(x & 0xffffff0000000000ull)
		+ bool(x & 0xffff000000000000ull)
		+ bool(x & 0xff00000000000000ull));
      }

      struct __encoder_size : public boost::static_visitor<size_t>
      {
	size_t operator()(const attribute_set_type::int_type& x) const
	{
	  return byte_size(utils::bithack::branch(x < 0, - x, x)) + 1;
	}
	
	size_t operator()(const attribute_set_type::float_type& x) const
	{
	  return sizeof(attribute_set_type::float_type) + 1;
	}
	
	size_t operator()(const attribute_set_type::string_type& x) const
	{
	  return x.size() + (x.size() / mask_size) + 1;
	}
	
      };
      
      struct __encoder : public boost::static_visitor<size_t>
      {
	__encoder(byte_type* __buffer) : buffer(__buffer) {}
	
	size_t operator()(const attribute_set_type::int_type& x) const
	{
	  const uint64_t value_encode = utils::bithack::branch(x < 0, - x, x);
	  const size_t   value_size = byte_size(value_encode);
	  
	  *buffer = utils::bithack::branch(x < 0, mask_signed, mask_unsigned) | (value_size & mask_size);
	  
	  byte_type* data = buffer + 1;
	  
	  switch (value_size) {
	  case 8: data[value_size - 8] = (value_encode >> 56);
	  case 7: data[value_size - 7] = (value_encode >> 48);
	  case 6: data[value_size - 6] = (value_encode >> 40);
	  case 5: data[value_size - 5] = (value_encode >> 32);
	  case 4: data[value_size - 4] = (value_encode >> 24);
	  case 3: data[value_size - 3] = (value_encode >> 16);
	  case 2: data[value_size - 2] = (value_encode >> 8);
	  case 1: data[value_size - 1] = (value_encode);
	  }
	  
	  return value_size + 1;
	}
	
	size_t operator()(const attribute_set_type::float_type& x) const
	{
	  *buffer = (mask_float | (sizeof(attribute_set_type::float_type) & mask_size));
	  
	  std::copy(cast(x), cast(x) + sizeof(attribute_set_type::float_type), buffer + 1);
	  
	  return sizeof(attribute_set_type::float_type) + 1;
	}
	
	size_t operator()(const attribute_set_type::string_type& x) const
	{
	  byte_type* biter = buffer;
	  size_t copied = 0;
	  attribute_set_type::string_type::const_iterator iter_end = x.end();
	  for (attribute_set_type::string_type::const_iterator iter = x.begin(); iter <= iter_end; iter += mask_size) {
	    const size_t block_size = std::distance(iter, std::min(iter + mask_size, iter_end));
	    
	    *biter = (mask_string | (block_size & mask_size));
	    ++ biter;
	    
	    std::copy(iter, iter + block_size, biter);
	    copied += block_size + 1;
	  }
	  return copied;
	}
	
	byte_type* buffer;
      };
      
      static size_t encode(const attribute_set_type::data_type& value)
      {
	return boost::apply_visitor(__encoder_size(), value);
      }

      static size_t encode(byte_type* buffer, const attribute_set_type::data_type& value)
      {
	return boost::apply_visitor(__encoder(buffer), value);
      }
      
      static
      size_t decode(const byte_type* buffer, attribute_set_type::data_type& value) 
      {
	if (*buffer & mask_float) {
	  attribute_set_type::float_type value_decode;
	  ++ buffer;
	  std::copy(buffer, buffer + sizeof(attribute_set_type::float_type), cast(value_decode));
	  value = value_decode;
	  return sizeof(attribute_set_type::float_type) + 1;
	} else if ((*buffer & mask_signed) || (*buffer & mask_unsigned)) {
	  const bool value_signed = (*buffer & mask_signed);
	  const size_t value_size = (*buffer & mask_size);

	  ++ buffer;
	
	  const uint64_t mask = 0xff;
	  attribute_set_type::int_type value_decode = 0;
	  switch (value_size) {
	  case 8: value_decode |= ((uint64_t(buffer[value_size - 8]) & mask) << 56);
	  case 7: value_decode |= ((uint64_t(buffer[value_size - 7]) & mask) << 48);
	  case 6: value_decode |= ((uint64_t(buffer[value_size - 6]) & mask) << 40);
	  case 5: value_decode |= ((uint64_t(buffer[value_size - 5]) & mask) << 32);
	  case 4: value_decode |= ((uint64_t(buffer[value_size - 4]) & mask) << 24);
	  case 3: value_decode |= ((uint64_t(buffer[value_size - 3]) & mask) << 16);
	  case 2: value_decode |= ((uint64_t(buffer[value_size - 2]) & mask) << 8);
	  case 1: value_decode |= ((uint64_t(buffer[value_size - 1]) & mask));
	  }
	  
	  value = utils::bithack::branch(value_signed, - int64_t(value_decode), int64_t(value_decode));
	  
	  return value_size + 1;
	} else {
	  attribute_set_type::string_type value_decode;
	  size_t ret = 0;
	  for (;;) {
	    const size_t value_size = (*buffer & mask_size);
	    
	    ++ buffer;
	    
	    std::copy(buffer, buffer + value_size, std::back_inserter(value_decode));
	    ret += value_size + 1;
	    
	    if (value_size != mask_size) break;
	  }

	  value = value_decode;
	  
	  return ret;
	}
      }
    };

    typedef __attribute_vector_attribute_codec codec_attribute_type;
    typedef __attribute_vector_data_codec      codec_data_type;

    struct __attribute_set_encoder
    {
      template <typename Iterator>
      size_t operator()(Iterator first, Iterator last) const
      {
	attribute_type::id_type id_prev = 0;
	size_t coded_size = 0;
	for (/**/; first != last; ++ first) {
	  const attribute_type::id_type id = attribute_type(first->first).id();
	  
	  coded_size += codec_attribute_type::encode(id - id_prev);
	  coded_size += codec_data_type::encode(first->second);
	  
	  id_prev = id;
	}
	return coded_size;
      }
      
      template <typename Iterator, typename Output>
      Output operator()(Iterator first, Iterator last, Output iter) const
      {
	attribute_type::id_type id_prev = 0;
	for (/**/; first != last; ++ first) {
	  const attribute_type::id_type id = attribute_type(first->first).id();
	  
	  std::advance(iter, codec_attribute_type::encode(reinterpret_cast<codec_attribute_type::byte_type*>(&(*iter)), id - id_prev));
	  std::advance(iter, codec_data_type::encode(reinterpret_cast<codec_data_type::byte_type*>(&(*iter)), first->second));
	  
	  id_prev = id;
	}
	return iter;
      }
    };
    
    template <typename Iterator>
    void encode(const attribute_set_type& x, Iterator result)
    {
      buffer.clear();
      buffer.resize(__attribute_set_encoder()(x.begin(), x.end()));
      
      std::copy(buffer.begin(), __attribute_set_encoder()(x.begin(), x.end(), buffer.begin()), result);
    }
    
    template <typename Iterator>
    void decode(Iterator first, Iterator last, attribute_set_type& x)
    {
      x.clear();
      
      attribute_type::id_type id = 0;
      attribute_type::id_type inc = 0;
      attribute_set_type::data_type data;
      while (first != last) {
	std::advance(first, codec_attribute_type::decode(reinterpret_cast<const codec_attribute_type::byte_type*>(&(*first)), inc));
	std::advance(first, codec_data_type::decode(reinterpret_cast<const codec_data_type::byte_type*>(&(*first)), data));
	
	id += inc;
	x[id] = data;
      }
    }
    
    buffer_type buffer;
  };
  
  template <typename Iterator>
  inline
  void attribute_vector_encode(const AttributeVector& x, Iterator result)
  {
    AttributeVectorCODEC codec;
    codec.encode(x, result);
  }
  
  template <typename Iterator>
  void attribute_vector_decode(Iterator first, Iterator last, AttributeVector& x)
  {
    AttributeVectorCODEC codec;
    codec.decode(first, last, x);
  }
  
};

#endif
