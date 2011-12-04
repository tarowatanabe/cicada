// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE_VECTOR_CODEC__HPP__
#define __CICADA__FEATURE_VECTOR_CODEC__HPP__ 1

#include <memory>
#include <utility>
#include <algorithm>
#include <iterator>

#include <cicada/feature.hpp>

#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>
#include <utils/simple_vector.hpp>

namespace cicada
{
  template <typename __Tp, typename __Alloc >
  class FeatureVector;

  // a compact feature vector representation which do not allow any modification, and
  // uses input-iterator, not bidirectional/random-access iterator
  // we use double as our underlying stroage..
  
  struct FeatureVectorCODEC
  {
    struct __feature_vector_feature_codec
    {
      typedef cicada::Feature feature_type;
      typedef uint8_t byte_type;

      static size_t encode(byte_type* buffer, const feature_type::id_type& value)
      {
	return utils::byte_aligned_encode(value, reinterpret_cast<char*>(buffer));
      }
      
      static size_t encode(byte_type* buffer, const feature_type& value)
      {
	return encode(buffer, value.id());
      }
      
      static size_t decode(const byte_type* buffer, feature_type::id_type& value)
      {
	return utils::byte_aligned_decode(value, reinterpret_cast<const char*>(buffer));
      }
      static size_t decode(const byte_type* buffer, feature_type& value)
      {
	feature_type::id_type value_id = 0;
	const size_t ret = utils::byte_aligned_decode(value_id, reinterpret_cast<const char*>(buffer));
	value = feature_type(value_id);
	return ret;
      }
    };
  
    struct __feature_vector_data_codec
    {    
      typedef uint8_t byte_type;

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

      static size_t encode(byte_type* buffer, const double& value)
      {
	static const uint8_t mask_float    = 1 << (4 + 0);
	static const uint8_t mask_unsigned = 1 << (4 + 1);
	static const uint8_t mask_signed   = 1 << (4 + 2);
	static const uint8_t mask_size     = 0x0f;
      
	if (::fmod(value, 1.0) == 0.0) {
	  const int64_t  val = value;
	  const uint64_t value_encode = utils::bithack::branch(val < 0, - val, val);
	  const size_t   value_size = byte_size(value_encode);
	
	  *buffer = utils::bithack::branch(val < 0, mask_signed, mask_unsigned) | (value_size & mask_size);
	  ++ buffer;
	
	  switch (value_size) {
	  case 8: buffer[value_size - 8] = (value_encode >> 56);
	  case 7: buffer[value_size - 7] = (value_encode >> 48);
	  case 6: buffer[value_size - 6] = (value_encode >> 40);
	  case 5: buffer[value_size - 5] = (value_encode >> 32);
	  case 4: buffer[value_size - 4] = (value_encode >> 24);
	  case 3: buffer[value_size - 3] = (value_encode >> 16);
	  case 2: buffer[value_size - 2] = (value_encode >> 8);
	  case 1: buffer[value_size - 1] = (value_encode);
	  }
	
	  return value_size + 1;
	} else {
	  *buffer = (mask_float | (sizeof(double) & mask_size));
	  ++ buffer;
	
	  std::copy(cast(value), cast(value) + sizeof(double), buffer);
	
	  return sizeof(double) + 1;
	}
      }
    
      static
      size_t decode(const byte_type* buffer, double& value) 
      {
	static const uint8_t mask_float    = 1 << (4 + 0);
	static const uint8_t mask_unsigned = 1 << (4 + 1);
	static const uint8_t mask_signed   = 1 << (4 + 2);
	static const uint8_t mask_size     = 0x0f;

	if (*buffer & mask_float) {
	  ++ buffer;
	  std::copy(buffer, buffer + sizeof(double), cast(value));
	  return sizeof(double) + 1;
	} else if ((*buffer & mask_signed) || (*buffer & mask_unsigned)) {
	  const bool value_signed = (*buffer & mask_signed);
	  const size_t value_size = (*buffer & mask_size);

	  ++ buffer;
	
	  const uint64_t mask = 0xff;
	  uint64_t value_decode = 0;
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
	} else
	  throw std::runtime_error("invalid type for decoding");
      }
    };

    typedef __feature_vector_feature_codec codec_feature_type;
    typedef __feature_vector_data_codec    codec_data_type;

    typedef cicada::Feature feature_type;
    typedef char byte_type;
    typedef std::vector<byte_type, std::allocator<byte_type> > buffer_type;

    struct __encoder
    {
      template <typename Iterator, typename Output>
      Output operator()(Iterator first, Iterator last, Output iter) const
      {
	feature_type::id_type id_prev = 0;
	for (/**/; first != last; ++ first) {
	  const feature_type::id_type id = feature_type(first->first).id();
	  
	  std::advance(iter, codec_feature_type::encode(reinterpret_cast<codec_feature_type::byte_type*>(&(*iter)), id - id_prev));
	  std::advance(iter, codec_data_type::encode(reinterpret_cast<codec_data_type::byte_type*>(&(*iter)), first->second));
	  
	  id_prev = id;
	}
	
	return iter;
      }
    };
    
    template <typename Tp, typename Alloc, typename Iterator>
    void encode(const FeatureVector<Tp, Alloc>& x, Iterator result)
    {
      __encoder encoder;

      buffer.clear();
      buffer.reserve(x.size() * 16);
      buffer.resize(x.size() * 16);
      
      std::copy(buffer.begin(), encoder(x.begin(), x.end(), buffer.begin()), result);
    }
    
    template <typename Iterator, typename Tp, typename Alloc>
    void decode(Iterator first, Iterator last, FeatureVector<Tp, Alloc>& x)
    {
      x.clear();
      
      feature_type::id_type id = 0;
      feature_type::id_type inc = 0;
      double data;
      while (first != last) {
	std::advance(first, codec_feature_type::decode(reinterpret_cast<const codec_feature_type::byte_type*>(&(*first)), inc));
	std::advance(first, codec_data_type::decode(reinterpret_cast<const codec_data_type::byte_type*>(&(*first)), data));
	
	id += inc;
	x[id] = data;
      }
    }

    template <typename Iterator, typename _Vocab, typename Tp, typename Alloc>
    void decode(Iterator first, Iterator last, const _Vocab& vocab, FeatureVector<Tp, Alloc>& x)
    {
      x.clear();
      
      feature_type::id_type id = 0;
      feature_type::id_type inc = 0;
      double data;
      while (first != last) {
	std::advance(first, codec_feature_type::decode(reinterpret_cast<const codec_feature_type::byte_type*>(&(*first)), inc));
	std::advance(first, codec_data_type::decode(reinterpret_cast<const codec_data_type::byte_type*>(&(*first)), data));
	
	id += inc;
	x[vocab[id]] = data;
      }
    }
    
    buffer_type buffer;
  };
  
  
  template <typename Tp, typename Alloc, typename Iterator>
  inline
  void feature_vector_encode(const FeatureVector<Tp, Alloc>& x, Iterator result)
  {
    FeatureVectorCODEC codec;
    codec.encode(x, result);
  }
  
  template <typename Iterator, typename Tp, typename Alloc>
  inline
  void feature_vector_decode(Iterator first, Iterator last, FeatureVector<Tp, Alloc>& x)
  {
    FeatureVectorCODEC codec;
    codec.decode(first, last, x);
  }

  template <typename Iterator, typename _Vocab, typename Tp, typename Alloc>
  inline
  void feature_vector_decode(Iterator first, Iterator last, const _Vocab& vocab, FeatureVector<Tp, Alloc>& x)
  {
    FeatureVectorCODEC codec;
    codec.decode(first, last, vocab, x);
  }
  
};


#include <cicada/feature_vector.hpp>

#endif
