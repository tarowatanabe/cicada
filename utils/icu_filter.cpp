//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <stdexcept>

#include <unicode/unistr.h>
#include <unicode/uclean.h>

#include "icu_filter.hpp"

namespace utils
{
  struct callback_struct
  {
    // from-u-callback
    UConverterFromUCallback callback_from;
    const void* context_from;
    
    // to-u-callback
    UConverterToUCallback   callback_to;
    const void* context_to;
  };
  
  static callback_struct callbacks[__icu_filter_impl::__callback_size] = {
    { UCNV_FROM_U_CALLBACK_SUBSTITUTE, 0,
      UCNV_TO_U_CALLBACK_SUBSTITUTE, 0,},
    { UCNV_FROM_U_CALLBACK_SKIP, 0,
      UCNV_TO_U_CALLBACK_SKIP, 0 },
    { UCNV_FROM_U_CALLBACK_STOP, 0,
      UCNV_TO_U_CALLBACK_STOP, 0 },
    { UCNV_FROM_U_CALLBACK_ESCAPE, 0,
      UCNV_TO_U_CALLBACK_ESCAPE, 0},
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_ICU,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_ICU },
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_JAVA,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_JAVA },
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_C,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_C },
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_XML_HEX,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_XML_HEX },
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_XML_HEX,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_XML_HEX },
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_XML_DEC,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_XML_DEC },
    { UCNV_FROM_U_CALLBACK_ESCAPE, UCNV_ESCAPE_UNICODE,
      UCNV_TO_U_CALLBACK_ESCAPE, UCNV_ESCAPE_UNICODE }
  };
  
  void __icu_filter_impl::__initialize(const std::string& codepage_from,
				       const std::string& codepage_to,
				       const callback_type callback)
  {
    __clear();
    
    if (codepage_from == codepage_to)
      return;
    
    // conversin from
    UErrorCode status = U_ZERO_ERROR;
    ucnv_from = ucnv_open(codepage_from.empty() ? 0 : codepage_from.c_str(), &status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucnv_open(): ") + codepage_from + " " + u_errorName(status));
    // setup callback
    status = U_ZERO_ERROR;
    ucnv_setToUCallBack(ucnv_from,
			callbacks[callback].callback_to, callbacks[callback].context_to,
			0, 0, &status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucnv_setToUCallBack(): ") + u_errorName(status));
    
    // convertion to
    status = U_ZERO_ERROR;
    ucnv_to = ucnv_open(codepage_to.empty() ? 0 : codepage_to.c_str(), &status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucnv_open(): ") + codepage_to + " " + u_errorName(status));
    // setup callback
    status = U_ZERO_ERROR;
    ucnv_setFromUCallBack(ucnv_to,
			  callbacks[callback].callback_from, callbacks[callback].context_from,
			  0, 0, &status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucnv_setFromUCallBack(): ") + u_errorName(status));
    
    pivot_start = new UChar[boost::iostreams::default_device_buffer_size];
    pivot_source = pivot_start;
    pivot_target = pivot_start;

    status = U_ZERO_ERROR;
    const char* encoding_from = ucnv_getName(ucnv_from, &status);
    status = U_ZERO_ERROR;
    const char* encoding_to   = ucnv_getName(ucnv_to, &status);
    
    // if the same encoding, clear!
    if (strcasecmp(encoding_from, encoding_to) == 0)
      __clear();

    message_from.clear();
    message_to.clear();
  }

  void __icu_filter_impl::__close()
  {
    if (ucnv_from)
      ucnv_reset(ucnv_from);
    if (ucnv_to)
      ucnv_reset(ucnv_to);
    
    // dump error message...
    if (! message_from.empty())
      std::cerr << message_from << std::endl;
    if (! message_to.empty())
      std::cerr << message_to << std::endl;
    
    offset_from = 0;
    offset_pivot_source = 0;
    offset_pivot_target = 0;
    offset_to = 0;
    
    message_from.clear();
    message_to.clear();
  }
  
  void __icu_filter_impl::__clear()
  {
    if (ucnv_from)
      ucnv_close(ucnv_from);
    if (ucnv_to)
      ucnv_close(ucnv_to);
    if (pivot_start)
      delete [] pivot_start;
 
    ucnv_from = 0;
    ucnv_to = 0;
    
    pivot_start = 0;
    pivot_source = 0;
    pivot_target = 0;

    offset_from = 0;
    offset_pivot_source = 0;
    offset_pivot_target = 0;
    offset_to = 0;
    
    message_from.clear();
    message_to.clear();
  }
  
}
