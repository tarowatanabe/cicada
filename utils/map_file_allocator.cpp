//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "map_file_allocator.hpp"

namespace utils
{
  
  boost::mutex __map_alloc_impl::global_alloc::mutex;
  
  __map_alloc_impl::global_alloc __map_alloc_impl::__global_alloc = __map_alloc_impl::global_alloc();
  
#ifdef HAVE_TLS
  __thread __map_file_allocator_base::map_alloc_type* __map_file_allocator_base::local_alloc_thread = 0;
#endif
  boost::thread_specific_ptr<__map_file_allocator_base::map_alloc_type> __map_file_allocator_base::local_alloc;
  
};
