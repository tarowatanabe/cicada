// -*- mode: c++ -*-

#ifndef __CODEC__CODEC_IMPL_HPP__
#define __CODEC__CODEC_IMPL_HPP__ 1

#include <boost/detail/endian.hpp>

namespace codec
{
  namespace impl
  {
#ifdef BOOST_LITTLE_ENDIAN
    inline
    void write_size(const uint32_t size, char* buffer)
    {
      *((uint32_t*) buffer) = size;
    }

    inline
    uint32_t read_size(const char* buffer)
    {
      return *((const uint32_t*) buffer);
    }
    
#else
    inline
    void write_size(const uint32_t size, char* buffer)
    {
      buffer[0] = size & 0xff;
      buffer[1] = (size >> 8) & 0xff;
      buffer[2] = (size >> 16) & 0xff;
      buffer[3] = (size >> 24) & 0xff;
    }

    inline
    uint32_t read_size(const char* buffer)
    {
      return (uint32_t((uint8_t) buffer[0]) 
	      | (uint32_t((uint8_t) buffer[1]) << 8)
	      | (uint32_t((uint8_t) buffer[2]) << 16)
	      | (uint32_t((uint8_t) buffer[3]) << 24));
    }
    
#endif
  };
};
#endif
