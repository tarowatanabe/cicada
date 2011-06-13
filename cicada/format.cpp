//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "format.hpp"
#include "parameter.hpp"

#include "format/number.hpp"

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/thread_specific_ptr.hpp>

namespace cicada
{
  const char* Format::lists()
  {
    static const char* desc = "\
number: number format\n\
\tsource=[locale] parser locale\n\
\ttarget=[locale] generator locale\n\
\tparser=[locale] parser locale\n\
\tgenerator=[locale] generator locale\n\
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
  
  typedef boost::shared_ptr<Format> format_ptr_type;
  
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, format_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, format_ptr_type> > > format_map_type;
#else
  typedef sgi::hash_map<std::string, format_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, format_ptr_type> > > format_map_type;
#endif
  
#ifdef HAVE_TLS
  static __thread format_map_type* __formats_tls = 0;
  static boost::thread_specific_ptr<format_map_type> __formats;
#else
  static utils::thread_specific_ptr<format_map_type> __formats;
#endif

  Format& Format::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    if (! __formats_tls) {
      __formats.reset(new format_map_type());
      __formats_tls = __formats.get();
    }
    format_map_type& formats_map = *__formats_tls;    
#else
    if (! __formats.get())
      __formats.reset(new format_map_type());
    
    format_map_type& formats_map = *__formats;
#endif
    
    const parameter_type param(parameter);

    if (utils::ipiece(param.name()) == "number") {
      
      format_map_type::iterator iter = formats_map.find(parameter);
      if (iter == formats_map.end()) {
	std::string locale_parser;
	std::string locale_generator;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "parser" || utils::ipiece(piter->first) == "source")
	    locale_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "generator" || utils::ipiece(piter->first) == "target")
	    locale_generator = piter->second;
	  else
	    throw std::runtime_error("unsupported parameter: " + parameter);
	}
	
	iter = formats_map.insert(std::make_pair(parameter, format_ptr_type(new format::Number(locale_parser, locale_generator)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
};
