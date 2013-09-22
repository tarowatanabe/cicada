//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

#include "language.hpp"

#include <unicode/locid.h>
#include <unicode/locdspnm.h>

#include "utils/unordered_map.hpp"
#include "utils/piece.hpp"

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace format
  {
    class LanguageImpl
    {
    public:
      typedef Format::phrase_type     phrase_type;
      typedef Format::phrase_tag_type phrase_tag_type;
      typedef Format::phrase_set_type phrase_set_type;

private:
      typedef utils::unordered_map<std::string, std::string, boost::hash<utils::piece>, std::equal_to<std::string>,
				   std::allocator<std::pair<const std::string, std::string> > >::type language_map_type;

      
    public:
      LanguageImpl(const std::string& locale_str_source,
		   const std::string& locale_str_target)
      {
	const icu::Locale locale_source(locale_str_source.c_str());
	const icu::Locale locale_target(locale_str_target.c_str());
	
	if (locale_source.isBogus())
	  throw std::runtime_error("invalid locale: " + locale_str_source);
	if (locale_target.isBogus())
	  throw std::runtime_error("invalid locale: " + locale_str_target);

	icu::UnicodeString usource;
	icu::UnicodeString utarget;
	
	std::string source;
	std::string target;

	std::auto_ptr<LocaleDisplayNames> lsource(LocaleDisplayNames::createInstance(locale_source));
	std::auto_ptr<LocaleDisplayNames> ltarget(LocaleDisplayNames::createInstance(locale_target));
	
	for (const char* const* iter =  icu::Locale::getISOLanguages(); *iter; ++ iter) {
	  lsource->languageDisplayName(*iter, usource);
	  ltarget->languageDisplayName(*iter, utarget);
	  
	  source.clear();
	  target.clear();

	  usource.toUTF8String(source);
	  utarget.toUTF8String(target);

	  // skip untrnlatated langauge
	  if (source == *iter || target == *iter) continue;

	  // std::cout << "lang: " << *iter << " source: " << source << " target: " << target << std::endl;
	  
	  language[source] = target;
	}
      }
      
      
      void generate(const phrase_type& phrase, phrase_set_type& generated) const
      {
	generated.clear();
	
	language_map_type::const_iterator iter = language.find(phrase);
	if (iter != language.end())
	  generated.push_back(phrase_tag_type(iter->second, "language"));
      }

      language_map_type language;
    };
    
    
    void Language::operator()(const phrase_type& phrase, phrase_set_type& generated) const
    {
      pimpl->generate(phrase, generated);
    }
    
    
    Language::Language(const std::string& locale_str_source,
		       const std::string& locale_str_target)
      : pimpl(new LanguageImpl(locale_str_source, locale_str_target)) {}
    
    Language::~Language()
    {
      delete pimpl;
    }
    
  };
};
