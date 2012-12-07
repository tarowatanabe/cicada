//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "codec/fastlz.hpp"
#include "codec/fastlz.h"
#include "codec/codec_impl.hpp"

#include "utils/bithack.hpp"

namespace codec
{
  namespace detail
  {
    const fastlz_param::size_type fastlz_param::chunk_size = 1024 * 1024;
    const fastlz_param::size_type fastlz_param::bound_size = (1024 * 1024 + ((1024 * 1024) >> 4)) + 4;
  };

  fastlz_compressor_impl::fastlz_compressor_impl()
    : buffer(new byte_type[chunk_size]),
      pos(0),
      buffer_compressed(new byte_type[bound_size]),
      size_compressed(0),
      pos_compressed(0)
  { }
  
  fastlz_compressor_impl::fastlz_compressor_impl(const fastlz_compressor_impl&)
    : buffer(new byte_type[chunk_size]),
      pos(0),
      buffer_compressed(new byte_type[bound_size]),
      size_compressed(0),
      pos_compressed(0)
  { }
  
  bool fastlz_compressor_impl::filter(const char*& src_begin,
				   const char* src_end,
				   char*& dest_begin,
				   char* dest_end,
				   bool flush)
  {
    const size_type src_copied = utils::bithack::min(size_type(src_end - src_begin), chunk_size - pos);
      
    // copy into buffer
    if (src_copied) {
      std::copy(src_begin, src_begin + src_copied, buffer + pos);
      src_begin += src_copied;
      pos += src_copied;
    }
      
    // perform compression, if all the compressed buffer is dumped, pos == chunk or flush
    if (pos_compressed == size_compressed && (pos == chunk_size || (pos && flush))) {
      size_compressed = fastlz_compress(buffer, pos, buffer_compressed + 4);
      
      impl::write_size(size_compressed, buffer_compressed);
      
      size_compressed += 4;
      pos_compressed = 0;
      pos = 0;
    }
    
    // perform copy into dest
    const size_type dest_copied = utils::bithack::min(size_type(dest_end - dest_begin), size_compressed - pos_compressed);
    
    if (dest_copied) {
      std::copy(buffer_compressed + pos_compressed, buffer_compressed + pos_compressed + dest_copied, dest_begin);
      dest_begin += dest_copied;
      pos_compressed += dest_copied;
    }
    
    return (size_compressed != pos_compressed) || pos;
  }


  fastlz_decompressor_impl::fastlz_decompressor_impl()
    : buffer_compressed(new byte_type[bound_size]),
      size_compressed(0),
      pos_compressed(0),
      buffer(new byte_type[chunk_size]),
      size(0),
      pos(0)
  { }
  
  fastlz_decompressor_impl::fastlz_decompressor_impl(const fastlz_decompressor_impl&)
    : buffer_compressed(new byte_type[bound_size]),
      size_compressed(0),
      pos_compressed(0),
      buffer(new byte_type[chunk_size]),
      size(0),
      pos(0)
  { }
  
  bool fastlz_decompressor_impl::filter(const char*& src_begin,
				     const char* src_end,
				     char*& dest_begin,
				     char* dest_end,
				     bool flush)
  {
    // copy into buffer as much as possible
    
    const size_type src_copied = utils::bithack::min(size_type(src_end - src_begin),
						     utils::bithack::branch(size_compressed,
									    size_compressed + 4,
									    bound_size)
						     - pos_compressed);
    
    if (src_copied) {
      std::copy(src_begin, src_begin + src_copied, buffer_compressed + pos_compressed);
      src_begin += src_copied;
      pos_compressed += src_copied;
    }

    // assig size-compressed if possible
    if ((! size_compressed) && pos_compressed >= 4)
      size_compressed = impl::read_size(buffer_compressed);

    // perform actual uncompression...
    if (pos == size && size_compressed && pos_compressed >= size_compressed + 4) {
      size = fastlz_decompress(buffer_compressed + 4, size_compressed, buffer, chunk_size);
      pos = 0;
      
      // copy with potential overlap...
      std::copy(buffer_compressed + size_compressed + 4, buffer_compressed + pos_compressed, buffer_compressed);
      pos_compressed -= size_compressed + 4;
      size_compressed = 0;
      
      // assig size-compressed if possible
      if ((! size_compressed) && pos_compressed >= 4)
	size_compressed = impl::read_size(buffer_compressed);
    }
    
    // dump into dest
    const size_type dest_copied = utils::bithack::min(size_type(dest_end - dest_begin), size - pos);
    
    if (dest_copied) {
      std::copy(buffer + pos, buffer + pos + dest_copied, dest_begin);
      dest_begin += dest_copied;
      pos += dest_copied;
    }
    
    return (pos != size) || pos_compressed;
  }

}
