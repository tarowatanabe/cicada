//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "operation.hpp"

#include "utils/hashmurmur.hpp"
#include "utils/unordered_map.hpp"
#include "utils/compress_stream.hpp"
#include "utils/thread_specific_ptr.hpp"

#include <boost/filesystem.hpp>

namespace cicada
{

  namespace operation_detail {
    
    typedef Operation::weight_set_type   weight_set_type;
    typedef Operation::weights_path_type weights_path_type;
    
    struct hash_string : public utils::hashmurmur<size_t>
    {
      size_t operator()(const std::string& x) const
      {
	return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
      }
    };
    
    typedef utils::unordered_map<std::string, weights_path_type, hash_string, std::equal_to<std::string>,
				 std::allocator<std::pair<const std::string, weights_path_type> > >::type weight_map_type;

#ifdef HAVE_TLS
    static __thread weight_map_type* __weights_tls = 0;
    static boost::thread_specific_ptr<weight_map_type> __weights;
#else
    static utils::thread_specific_ptr<weight_map_type> __weights;
#endif

  };
  
  const Operation::weights_path_type& Operation::weights()
  {
    return weights(path_type());
  }

  const Operation::weights_path_type& Operation::weights(const path_type& path)
  {
    typedef operation_detail::weight_map_type weight_map_type;
    
#ifdef HAVE_TLS
    if (! operation_detail::__weights_tls) {
      operation_detail::__weights.reset(new weight_map_type());
      operation_detail::__weights_tls = operation_detail::__weights.get();
    }
    weight_map_type& weights_map = *operation_detail::__weights_tls;
#else
    if (! operation_detail::__weights.get())
      operation_detail::__weights.reset(new weight_map_type());
    
    weight_map_type& weights_map = *operation_detail::__weights;
#endif
    
    weight_map_type::iterator iter = weights_map.find(path.string());
    if (iter == weights_map.end()) {
      iter = weights_map.insert(std::make_pair(path.string(), weights_path_type(path))).first;
      
      if (! path.empty()) {
	if (path != "-" && ! boost::filesystem::exists(path))
	  throw std::runtime_error("no feture weights? " + path.string());
	
	utils::compress_istream is(path, 1024 * 1024);
	is >> iter->second.weights;
      }
    }
    return iter->second;
  }
};
