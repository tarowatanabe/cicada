//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

#include "region.hpp"

#include <unicode/locid.h>
#include <unicode/locdspnm.h>

#include "utils/unordered_map.hpp"
#include "utils/piece.hpp"

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace format
  {
    class RegionImpl
    {
    public:
      typedef Format::phrase_type     phrase_type;
      typedef Format::phrase_tag_type phrase_tag_type;
      typedef Format::phrase_set_type phrase_set_type;

    private:
      typedef utils::unordered_map<std::string, std::string, boost::hash<utils::piece>, std::equal_to<std::string>,
				   std::allocator<std::pair<const std::string, std::string> > >::type region_map_type;

      
    public:
      RegionImpl(const std::string& locale_str_source,
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

	for (const char* const* iter =  icu::Locale::getISOCountries(); *iter; ++ iter) {
	  lsource->regionDisplayName(*iter, usource);
	  ltarget->regionDisplayName(*iter, utarget);
	  
	  source.clear();
	  target.clear();

	  usource.toUTF8String(source);
	  utarget.toUTF8String(target);
	  
	  // skip untrnlatated region
	  if (source == *iter || target == *iter) continue;
	  
          // std::cout << "region: " << *iter << " source: " << source << " target: " << target << std::endl;

	  region[source] = target;
	}
      }
      
      
      void generate(const phrase_type& phrase, phrase_set_type& generated) const
      {
	generated.clear();
	
	region_map_type::const_iterator iter = region.find(phrase);
	if (iter != region.end())
	  generated.push_back(phrase_tag_type(iter->second, "region"));
      }

      region_map_type region;
    };
    
    
    void Region::operator()(const phrase_type& phrase, phrase_set_type& generated) const
    {
      pimpl->generate(phrase, generated);
    }
    
    
    Region::Region(const std::string& locale_str_source,
		   const std::string& locale_str_target)
      : pimpl(new RegionImpl(locale_str_source, locale_str_target)) {}
    
    Region::~Region()
    {
      delete pimpl;
    }
    
  };
};
