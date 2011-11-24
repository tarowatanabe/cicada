// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE_VECTOR_COMPACT__HPP__
#define __CICADA__FEATURE_VECTOR_COMPACT__HPP__ 1

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

      const int64_t val = value;
      
      if (double(val) == value) {
	const uint64_t value_encode = utils::bithack::branch(val < 0, - val, val);
	const size_t value_size = byte_size(value_encode);
	
	*buffer = utils::bithack::branch(val < 0, mask_signed, mask_unsigned) | (value_size & mask_size);
	++ buffer;
	
	switch (value_size) {
	case 8: buffer[7] = (value_encode >> 56);
	case 7: buffer[6] = (value_encode >> 48);
	case 6: buffer[5] = (value_encode >> 40);
	case 5: buffer[4] = (value_encode >> 32);
	case 4: buffer[3] = (value_encode >> 24);
	case 3: buffer[2] = (value_encode >> 16);
	case 2: buffer[1] = (value_encode >> 8);
	case 1: buffer[0] = (value_encode);
	}
	
	return value_size + 1;
      } else {
	const float val = value;
	
	if (double(val) == value) {
	  *buffer = (mask_float | (sizeof(float) & mask_size));
	  ++ buffer;
	  
	  std::copy(cast(val), cast(val) + sizeof(float), buffer);
	  
	  return sizeof(float) + 1;
	} else {
	  *buffer = (mask_float | (sizeof(double) & mask_size));
	  ++ buffer;
	  
	  std::copy(cast(value), cast(value) + sizeof(double), buffer);
	  
	  return sizeof(double) + 1;
	}
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
	if ((*buffer & mask_size) == sizeof(float)) {
	  ++ buffer;
	  
	  float value_decode;
	  std::copy(buffer, buffer + sizeof(float), cast(value_decode));
	  value = value_decode;
	  
	  return sizeof(float) + 1;
	} else if ((*buffer & mask_size) == sizeof(double)) {
	  ++ buffer;
	  
	  std::copy(buffer, buffer + sizeof(double), cast(value));
	  
	  return sizeof(double) + 1;
	} else
	  throw std::runtime_error("invalid size for float");
      } else if ((*buffer & mask_signed) || (*buffer & mask_unsigned)) {
	const bool value_signed = (*buffer & mask_signed);
	const size_t value_size = (*buffer & mask_size);

	++ buffer;
	
	const uint64_t mask = 0xff;
	uint64_t value_decode = 0;
	switch (value_size) {
	case 8: value_decode |= ((uint64_t(buffer[7]) & mask) << 56);
	case 7: value_decode |= ((uint64_t(buffer[6]) & mask) << 48);
	case 6: value_decode |= ((uint64_t(buffer[5]) & mask) << 40);
	case 5: value_decode |= ((uint64_t(buffer[4]) & mask) << 32);
	case 4: value_decode |= ((uint64_t(buffer[3]) & mask) << 24);
	case 3: value_decode |= ((uint64_t(buffer[2]) & mask) << 16);
	case 2: value_decode |= ((uint64_t(buffer[1]) & mask) << 8);
	case 1: value_decode |= ((uint64_t(buffer[0]) & mask));
	}
	
	value = utils::bithack::branch(value_signed, - int64_t(value_decode), int64_t(value_decode));
	
	return value_size + 1;
      } else
	throw std::runtime_error("invalid type for decoding");
    }
  };

  
  class FeatureVectorCompact
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Feature feature_type;
    typedef cicada::Feature key_type;
    typedef double mapped_type;
    typedef double data_type;
    
    typedef std::pair<const feature_type, data_type> value_type;

    typedef uint8_t byte_type;
        
  private:
    typedef utils::simple_vector<byte_type, std::allocator<byte_type> > storage_type;
    
  public:
    typedef __feature_vector_feature_codec codec_key_type;
    typedef __feature_vector_feature_codec codec_feature_type;
    typedef __feature_vector_data_codec    codec_data_type;
    typedef __feature_vector_data_codec    codec_mapped_type;
    
  public:
    struct iterator
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      typedef std::input_iterator_tag   iterator_category;
      typedef std::pair<const feature_type, data_type> value_type;
      typedef const value_type* pointer;
      typedef const value_type& reference;
      
      typedef const byte_type* ptr_type;
      
    public:
      iterator(ptr_type iter, ptr_type last) : __iter(iter), __last(last), __impl()
      {
	if (__iter != __last) {
	  const size_type feature_size = codec_feature_type::decode(&(*__iter), const_cast<feature_type&>(__impl.first));
	  __iter += feature_size;
	  
	  const size_type data_size = codec_data_type::decode(&(*__iter), const_cast<data_type&>(__impl.second));
	  __iter += data_size;
	} else {
	  __iter = 0;
	  __last = 0;
	}
      }
      
      iterator() : __iter(0), __last(0), __impl() {}
      iterator(const iterator& x) : __iter(x.__iter), __last(x.__last), __impl(x.__impl) {}
      iterator& operator=(const iterator& x)
      {
	__iter = x.__iter;
	__last = x.__last;
	const_cast<feature_type&>(__impl.first) = x.__impl.first;
	const_cast<data_type&>(__impl.second) = x.__impl.second;
	
	return *this;
      }
      
      const value_type& operator*() const { return __impl; }
      const value_type* operator->() const { return &__impl; }
      
      iterator& operator++()
      {
	increment();
	return *this;
      }
      
      iterator operator++(int)
      {
	iterator tmp = *this;
	increment();
	return tmp;
      }
      
      friend
      bool operator==(const iterator& x, const iterator& y)
      {
	return (x.__iter == y.__iter) && (x.__last == y.__last);
      }
      
      friend
      bool operator!=(const iterator& x, const iterator& y)
      {
	return (x.__iter != y.__iter) || (x.__last != y.__last);
      }
      
    private:
      void increment()
      {
	if (__iter == __last) {
	  __iter = 0;
	  __last = 0;
	} else {
	  feature_type::id_type feature_inc = 0;
	  const size_type feature_size = codec_feature_type::decode(&(*__iter), feature_inc);
	  const_cast<feature_type&>(__impl.first) = feature_type(__impl.first.id() + feature_inc);
	  __iter += feature_size;
	  
	  const size_type data_size = codec_data_type::decode(&(*__iter), const_cast<data_type&>(__impl.second));
	  __iter += data_size;
	}
      }

    private:
      ptr_type   __iter;
      ptr_type   __last;
      value_type __impl;
    };
    
    typedef iterator const_iterator;
    
    typedef const value_type& reference;
    typedef const value_type& const_reference;
    
    struct encoder_type
    {
      template <typename Iterator, typename Output>
      Output operator()(Iterator first, Iterator last, Output iter) const
      {
	feature_type::id_type id_prev = 0;
	for (/**/; first != last; ++ first) {
	  const feature_type::id_type id = feature_type(first->first).id();
	  
	  std::advance(iter, codec_feature_type::encode(&(*iter), id - id_prev));
	  std::advance(iter, codec_data_type::encode(&(*iter), first->second));
	  
	  id_prev = id;
	}
	
	return iter;
      }
    };

  private:
    template <typename Tp>
    struct less_first
    {
      bool operator()(const Tp& x, const Tp& y) const
      {
	return x.first < y.first;
      }
    };
    
  public:
    FeatureVectorCompact() {}
    
    template <typename Iterator>
    FeatureVectorCompact(Iterator first, Iterator last, const bool sorted=false) { assign(first, last, sorted); }
    
    template <typename T, typename A>
    FeatureVectorCompact(const FeatureVector<T, A>& x) { assign(x); }
    
    FeatureVectorCompact(const FeatureVectorCompact& x) : storage(x.storage) {}

    FeatureVectorCompact& operator=(const FeatureVectorCompact& x)
    {
      storage = x.storage;
      return *this;
    }
    
    template <typename T, typename A>
    FeatureVectorCompact& operator=(const FeatureVector<T, A>& x)
    {
      assign(x);
      return *this;
    }

    void assign(const FeatureVectorCompact& x)
    {
      storage.assign(x.storage);
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last, const bool sorted=false)
    {
      typedef std::pair<feature_type, double> pair_type;
      typedef std::vector<pair_type, std::allocator<pair_type> > raw_type;
      typedef std::vector<byte_type, std::allocator<byte_type> > compressed_type;

      encoder_type encoder;

      if (sorted) {
	compressed_type compressed(std::distance(first, last) * 16);
	
	storage.assign(compressed.begin(), encoder(first, last, compressed.begin()));
      } else {
	raw_type raw(first, last);
	std::sort(raw.begin(), raw.end(), less_first<pair_type>());
	
	compressed_type compressed(raw.size() * 16);
	
	storage.assign(compressed.begin(), encoder(raw.begin(), raw.end(), compressed.begin()));
      }
    }
    
    template <typename T, typename A>
    void assign(const FeatureVector<T, A>& x)
    {
      typedef std::vector<byte_type, std::allocator<byte_type> > compressed_type;

      encoder_type encoder;
      compressed_type compressed(x.size() * 16);
      
      storage.assign(compressed.begin(), encoder(x.begin(), x.end(), compressed.begin()));
    }

    
  public:
    const_iterator begin() const { return const_iterator(&(*storage.begin()), &(*storage.end())); }
    const_iterator end() const { return const_iterator(); }
    
    bool empty() const { return storage.empty(); }
    size_type size_compressed() const  { return storage.size(); }

    void swap(FeatureVectorCompact& x)
    {
      storage.swap(x.storage);
    }
    
  public:
    friend bool operator==(const FeatureVectorCompact& x, const FeatureVectorCompact& y) { return x.storage == y.storage; }
    friend bool operator!=(const FeatureVectorCompact& x, const FeatureVectorCompact& y) { return x.storage != y.storage; }
    friend bool operator<(const FeatureVectorCompact& x, const FeatureVectorCompact& y) { return x.storage < y.storage; }
    friend bool operator>(const FeatureVectorCompact& x, const FeatureVectorCompact& y) { return x.storage > y.storage; }
    friend bool operator<=(const FeatureVectorCompact& x, const FeatureVectorCompact& y) { return x.storage <= y.storage; }
    friend bool operator>=(const FeatureVectorCompact& x, const FeatureVectorCompact& y) { return x.storage >= y.storage; }
    
  private:
    storage_type storage;
  }; 
};

namespace std
{
  inline
  void swap(cicada::FeatureVectorCompact& x, cicada::FeatureVectorCompact& y)
  {
    x.swap(y);
  }
};


#include <cicada/feature_vector.hpp>

#endif
