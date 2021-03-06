//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "format.hpp"
#include "parameter.hpp"

#include "format/region.hpp"
#include "format/language.hpp"
#include "format/date.hpp"
#include "format/number.hpp"

#include <utils/unordered_map.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/thread_specific_ptr.hpp>

namespace cicada
{
  const char* Format::lists()
  {
    static const char* desc = "\
date: date/time format\n\
\tsource=[file] parser file\n\
\ttarget=[file] generator file\n\
\tparser=[file] parser file\n\
\tgenerator=[file] generator file\n\
\tlocale-source=[locale] parser locale\n\
\tlocale-target=[locale] generator locale\n\
\tlocale-parser=[locale] parser locale\n\
\tlocale-generator=[locale] generator locale\n\
language: language format\n\
\tlocale-source=[locale] parser locale\n\
\tlocale-target=[locale] generator locale\n\
\tlocale-parser=[locale] parser locale\n\
\tlocale-generator=[locale] generator locale\n\
number: number format\n\
\tsource=[file] parser file\n\
\ttarget=[file] generator file\n\
\tparser=[file] parser file\n\
\tgenerator=[file] generator file\n\
\tlocale-source=[locale] parser locale\n\
\tlocale-target=[locale] generator locale\n\
\tlocale-parser=[locale] parser locale\n\
\tlocale-generator=[locale] generator locale\n\
region: region format\n\
\tlocale-source=[locale] parser locale\n\
\tlocale-target=[locale] generator locale\n\
\tlocale-parser=[locale] parser locale\n\
\tlocale-generator=[locale] generator locale\n\
";
    
    return desc;
  }

  typedef boost::shared_ptr<Format> format_ptr_type;
  
  typedef utils::unordered_map<std::string, format_ptr_type, boost::hash<utils::piece>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, format_ptr_type> > >::type format_map_type;
  
#ifdef HAVE_TLS
  static __thread format_map_type* __formats_tls = 0;
  static utils::thread_specific_ptr<format_map_type> __formats;
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
	std::string file_parser;
	std::string file_generator;
	std::string locale_parser;
	std::string locale_generator;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "locale-parser" || utils::ipiece(piter->first) == "locale-source")
	    locale_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "locale-generator" || utils::ipiece(piter->first) == "locale-target")
	    locale_generator = piter->second;
	  else if (utils::ipiece(piter->first) == "parser" || utils::ipiece(piter->first) == "source")
	    file_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "generator" || utils::ipiece(piter->first) == "target")
	    file_generator = piter->second;
	  else
	    throw std::runtime_error("unsupported parameter: " + parameter);
	}
	
	iter = formats_map.insert(std::make_pair(parameter, format_ptr_type(new format::Number(file_parser, file_generator, locale_parser, locale_generator)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "date") {
      
      format_map_type::iterator iter = formats_map.find(parameter);
      if (iter == formats_map.end()) {
	std::string file_parser;
	std::string file_generator;
	std::string locale_parser;
	std::string locale_generator;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "locale-parser" || utils::ipiece(piter->first) == "locale-source")
	    locale_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "locale-generator" || utils::ipiece(piter->first) == "locale-target")
	    locale_generator = piter->second;
	  else if (utils::ipiece(piter->first) == "parser" || utils::ipiece(piter->first) == "source")
	    file_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "generator" || utils::ipiece(piter->first) == "target")
	    file_generator = piter->second;
	  else
	    throw std::runtime_error("unsupported parameter: " + parameter);
	}
	
	iter = formats_map.insert(std::make_pair(parameter, format_ptr_type(new format::Date(file_parser, file_generator, locale_parser, locale_generator)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "region") {
      
      format_map_type::iterator iter = formats_map.find(parameter);
      if (iter == formats_map.end()) {
	std::string locale_parser;
	std::string locale_generator;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "locale-parser" || utils::ipiece(piter->first) == "locale-source")
	    locale_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "locale-generator" || utils::ipiece(piter->first) == "locale-target")
	    locale_generator = piter->second;
	  else
	    throw std::runtime_error("unsupported parameter: " + parameter);
	}
	
	iter = formats_map.insert(std::make_pair(parameter, format_ptr_type(new format::Region(locale_parser, locale_generator)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "language") {
      
      format_map_type::iterator iter = formats_map.find(parameter);
      if (iter == formats_map.end()) {
	std::string locale_parser;
	std::string locale_generator;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "locale-parser" || utils::ipiece(piter->first) == "locale-source")
	    locale_parser = piter->second;
	  else if (utils::ipiece(piter->first) == "locale-generator" || utils::ipiece(piter->first) == "locale-target")
	    locale_generator = piter->second;
	  else
	    throw std::runtime_error("unsupported parameter: " + parameter);
	}
	
	iter = formats_map.insert(std::make_pair(parameter, format_ptr_type(new format::Language(locale_parser, locale_generator)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
};
