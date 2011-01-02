// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__ICU_FILTER__HPP__
#define __UTILS__ICU_FILTER__HPP__ 1

#include <stdint.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>

#include <unicode/uchar.h>
#include <unicode/ucnv.h>

#include <boost/iostreams/filter/symmetric.hpp>
#include <boost/iostreams/constants.hpp>
#include <boost/iostreams/pipeline.hpp>
#include <boost/iostreams/detail/ios.hpp>

#include <string>

namespace utils
{

  struct icu_filter_param
  {
    typedef enum {
      substitute = 0,
      skip,
      stop,
      escape,
      escape_icu,
      escape_java,
      escape_c,
      escape_xml,
      escape_xml_hex,
      escape_xml_dec,
      escape_unicode,
      __callback_size,
    } callback_type;
    
  };
  
  class __icu_filter_impl : public icu_filter_param
  {
  public:
    typedef char      char_type;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint64_t  off_type;
    
    typedef icu_filter_param::callback_type callback_type;

    __icu_filter_impl(const std::string& codepage_from,
		      const std::string& codepage_to,
		      const callback_type callback=stop)
      : ucnv_from(0), ucnv_to(0), pivot_start(0) { __initialize(codepage_from, codepage_to, callback); }
    ~__icu_filter_impl() { close(); __clear(); }
    
  public:
    bool filter(const char_type*& src_begin, const char_type* src_end,
		char_type*& dest_begin, char_type* dest_end,
		bool flush)
    {
      if (! ucnv_from && ! ucnv_to) {
	// no converter... simply copy!
	
	const size_t copy_size = std::min(src_end - src_begin, dest_end - dest_begin);
	
	std::copy(src_begin, src_begin + copy_size, dest_begin);
	
	src_begin += copy_size;
	dest_begin += copy_size;
	
	return false;
      }

      UChar* pivot_end = pivot_start + boost::iostreams::default_device_buffer_size;
      
      UErrorCode status = U_ZERO_ERROR;
      
      if (pivot_target != pivot_end) {
	const char_type* src_begin_prev = src_begin;
	UChar* pivot_target_prev = pivot_target;
	
	status = U_ZERO_ERROR;
	ucnv_toUnicode(ucnv_from, &pivot_target, pivot_end, &src_begin, src_end, 0, flush, &status);
	
	offset_from += src_begin - src_begin_prev;
	offset_pivot_target += pivot_target - pivot_target_prev;
	
	if (status != U_BUFFER_OVERFLOW_ERROR && U_FAILURE(status)) {
	  UErrorCode status_getname = U_ZERO_ERROR;
	  const char* encoding = ucnv_getName(ucnv_from, &status_getname);
	  
	  std::ostringstream offset_stream;
	  offset_stream << offset_from;

	  std::ostringstream offset_stream_unicode;
	  offset_stream_unicode << offset_pivot_target;
	  
	  message_from = (std::string("ucnv_toUnicode(): ") + u_errorName(status)
			  + " from " + encoding
			  + " offset: " + offset_stream.str()
			  + " unicode offset: " + offset_stream_unicode.str());
	  throw BOOST_IOSTREAMS_FAILURE(message_from);
	}
      }
      
      char_type*   dest_begin_prev = dest_begin;
      const UChar* pivot_source_prev = pivot_source;
      
      status = U_ZERO_ERROR;
      ucnv_fromUnicode(ucnv_to, &dest_begin, dest_end, &pivot_source, pivot_target, 0, flush, &status);
      
      offset_to += dest_begin - dest_begin_prev;
      offset_pivot_source += pivot_source - pivot_source_prev;
      
      if (status != U_BUFFER_OVERFLOW_ERROR && U_FAILURE(status)) {
	UErrorCode status_getname = U_ZERO_ERROR;
	const char* encoding = ucnv_getName(ucnv_to, &status_getname);
	
	std::ostringstream offset_stream;
	offset_stream << offset_to;
	
	std::ostringstream offset_stream_unicode;
	offset_stream_unicode << offset_pivot_source;
	
	message_to = (std::string("ucnv_fromUnicode(): ") + u_errorName(status) 
		      + " to " + encoding
		      + " offset: " + offset_stream.str()
		      + " unicode offset: " + offset_stream_unicode.str());
	throw BOOST_IOSTREAMS_FAILURE(message_to);
      }
      
      if (pivot_source == pivot_target) {
	pivot_source = pivot_start;
	pivot_target = pivot_start;
      }
      
      return status == U_BUFFER_OVERFLOW_ERROR;
    }
    
    void close()
    {
      __close();
    }
    
  private:
    void __initialize(const std::string& codepage_from,
		      const std::string& codepage_to,
		      const callback_type callback=stop);
    void __clear();    
    void __close();
    
  private:
    UConverter* ucnv_from;
    UConverter* ucnv_to;
    
    UChar*       pivot_start;
    const UChar* pivot_source;
    UChar*       pivot_target;
    
    off_type    offset_from;
    off_type    offset_pivot_source;
    off_type    offset_pivot_target;
    off_type    offset_to;

    std::string message_from;
    std::string message_to;
  };
  
  template <typename Alloc = std::allocator<char> >
  struct basic_icu_filter
    : boost::iostreams::symmetric_filter<__icu_filter_impl, Alloc>,
    icu_filter_param
  {
  private:
    typedef __icu_filter_impl impl_type;
    typedef boost::iostreams::symmetric_filter<__icu_filter_impl, Alloc> base_type;
    
  public:
    typedef typename base_type::char_type char_type;
    typedef typename base_type::category  category;
    
    typedef typename impl_type::callback_type callback_type;
    

    basic_icu_filter(const std::string codepage_from = "",
		     const std::string codepage_to = "",
		     const callback_type callback = impl_type::stop,
		     int buffer_size = boost::iostreams::default_device_buffer_size)
      : base_type(buffer_size, codepage_from, codepage_to, callback) {}
  };
  
  typedef basic_icu_filter<> icu_filter;

};

#endif
