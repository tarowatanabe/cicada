//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "signature.hpp"
#include "signature/chinese.hpp"
#include "signature/english.hpp"

#include "parameter.hpp"

#include <utils/sgi_hash_map.hpp>
#include <utils/thread_specific_ptr.hpp>
#include <utils/piece.hpp>
#include <utils/lexical_cast.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

namespace cicada
{
  const char* Signature::lists()
  {
    static const char* desc = "\
english: English signature\n\
";
    
    return desc;
  }
  
  template <typename Tp>
  struct hash_string : public utils::hashmurmur<size_t>
  {
    size_t operator()(const Tp& x) const
    {
      return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
    }
  };

  typedef boost::shared_ptr<Signature> signature_ptr_type;
  
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, signature_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, signature_ptr_type> > > signature_map_type;
#else
  typedef sgi::hash_map<std::string, signature_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, signature_ptr_type> > > signature_map_type;
#endif


#ifdef HAVE_TLS
  static __thread signature_map_type* __signatures_tls = 0;
  static boost::thread_specific_ptr<signature_map_type> __signatures;
#else
  static utils::thread_specific_ptr<signature_map_type> __signatures;
#endif


  Signature& Signature::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    if (! __signatures_tls) {
      __signatures.reset(new signature_map_type());
      __signatures_tls = __signatures.get();
    }
    signature_map_type& signatures_map = *__signatures_tls;    
#else
    if (! __signatures.get())
      __signatures.reset(new signature_map_type());
    
    signature_map_type& signatures_map = *__signatures;
#endif
    
    const parameter_type param(parameter);
    
    if (utils::ipiece(param.name()) == "chinese") {
      const std::string name("chinese");
      
      signature_map_type::iterator iter = signatures_map.find(name);
      if (iter == signatures_map.end()) {
	iter = signatures_map.insert(std::make_pair(name, signature_ptr_type(new signature::Chinese()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "english") {
      const std::string name("english");
      
      signature_map_type::iterator iter = signatures_map.find(name);
      if (iter == signatures_map.end()) {
	iter = signatures_map.insert(std::make_pair(name, signature_ptr_type(new signature::English()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
    
  }
  
};
