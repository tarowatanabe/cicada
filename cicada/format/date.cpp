//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <vector>
#include <string>
#include <set>
#include <memory>
#include <stdexcept>

#include <unicode/format.h>
#include <unicode/datefmt.h>
#include <unicode/unistr.h>
#include <unicode/locid.h>
#include <unicode/bytestream.h>

#include "date.hpp"

#include "utils/unordered_map.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/array_power2.hpp"

#include <boost/functional/hash.hpp>


namespace cicada
{
  namespace format
  {
    struct DateImpl : public utils::hashmurmur<size_t>
    {
    public:
      typedef icu::DateFormat parser_type;
      typedef icu::DateFormat generator_type;
      
      typedef std::vector<parser_type*, std::allocator<parser_type*> >       parser_set_type;
      typedef std::vector<generator_type*, std::allocator<generator_type*> > generator_set_type;
      
      typedef cicada::Format::phrase_type     phrase_type;
      typedef cicada::Format::phrase_set_type phrase_set_type;
      
      typedef utils::hashmurmur<size_t> hasher_type;

      struct cache_type
      {
	phrase_type     key;
	phrase_set_type value;
	
	cache_type() : key(), value() {}
      };
      typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;

    public:
      DateImpl() {}
      DateImpl(const DateImpl& x) { assign(x); }
      DateImpl& operator=(const DateImpl& x)
      {
	assign(x);
	return *this;
      }
      ~DateImpl() { clear(); }
      
    public:
      const phrase_set_type& operator()(const phrase_type& input) const
      {
	typedef std::set<std::string, std::less<std::string>, std::allocator<std::string> > phrase_unique_type;

	cache_type& cache = const_cast<cache_type&>(caches[hasher_type::operator()(input.begin(), input.end(), 0) & (caches.size() - 1)]);
	if (cache.key != input) {
	  icu::UnicodeString uphrase = icu::UnicodeString::fromUTF8(input);
	  icu::UnicodeString ugenerated;
	  std::string   generated;
	  
	  phrase_unique_type uniques;
	
	  parser_set_type::const_iterator piter_end = parsers.end();
	  for (parser_set_type::const_iterator piter = parsers.begin(); piter != piter_end; ++ piter) {
	    const parser_type* parser = *piter;
	    icu::Formattable   formattable;
	    icu::ParsePosition pos(0);
	  
	    parser->parseObject(uphrase, formattable, pos);
	  
	    if (pos.getErrorIndex() >= 0 || pos.getIndex() != uphrase.length()) continue;
	  
	    generator_set_type::const_iterator giter_end = generators.end();
	    for (generator_set_type::const_iterator giter = generators.begin(); giter != giter_end; ++ giter) {
	      const generator_type* generator = *giter;
	      UErrorCode status(U_ZERO_ERROR);
	    
	      ugenerated.remove();
	      generator->format(formattable, ugenerated, status);
	    
	      if (U_FAILURE(status)) continue;
	    
	      generated.clear();
	      ugenerated.toUTF8String(generated);
	      uniques.insert(generated);
	    }
	  }
	  
	  cache.key = input;
	  cache.value.clear();
	  cache.value.insert(cache.value.end(), uniques.begin(), uniques.end());
	}
	
	return cache.value;
      }
      
    public:
      void assign(const DateImpl& x)
      {
	clear();
	
	parser_set_type::const_iterator piter_end = x.parsers.end();
	for (parser_set_type::const_iterator piter = x.parsers.begin(); piter != piter_end; ++ piter)
	  parsers.push_back(dynamic_cast<parser_type*>((*piter)->clone()));
	
	generator_set_type::const_iterator giter_end = x.generators.end();
	for (generator_set_type::const_iterator giter = x.generators.begin(); giter != giter_end; ++ giter)
	  generators.push_back(dynamic_cast<generator_type*>((*giter)->clone()));
      }
      
      void clear()
      {
	parser_set_type::iterator piter_end = parsers.end();
	for (parser_set_type::iterator piter = parsers.begin(); piter != piter_end; ++ piter)
	  delete *piter;
	
	generator_set_type::iterator giter_end = generators.end();
	for (generator_set_type::iterator giter = generators.begin(); giter != giter_end; ++ giter)
	  delete *giter;
	
	parsers.clear();
	generators.clear();

	caches.clear();
      }
      
    public:
      cache_set_type     caches;
      parser_set_type    parsers;
      generator_set_type generators;
      std::string name;
    };
    
    void Date::operator()(const phrase_type& phrase, phrase_set_type& generated) const
    {
      typedef std::set<std::string, std::less<std::string>, std::allocator<std::string> > phrase_unique_type;
      
      phrase_unique_type uniques;
      
      for (pimpl_set_type::const_iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter) {
	const phrase_set_type& results = (*iter)->operator()(phrase);
	
	uniques.insert(results.begin(), results.end());
      }
      
      generated.clear();
      generated.insert(generated.end(), uniques.begin(), uniques.end());
    }

    
    Date::Date(const std::string& locale_str_source,
	       const std::string& locale_str_target)
    {
      // pre-defined rule-set
      typedef utils::unordered_map<std::string, impl_type, boost::hash<std::string>, std::equal_to<std::string>,
				   std::allocator<std::pair<const std::string, impl_type> > >::type impl_map_type;
      
      const icu::Locale locale_source(locale_str_source.c_str());
      const icu::Locale locale_target(locale_str_target.c_str());
      
      if (locale_source.isBogus())
	throw std::runtime_error("invalid locale: " + locale_str_source);
      if (locale_target.isBogus())
	throw std::runtime_error("invalid locale: " + locale_str_target);
      
      impl_map_type sources;
      impl_map_type targets;
      
      const icu::DateFormat::EStyle date_styles[4] = { icu::DateFormat::SHORT, icu::DateFormat::MEDIUM, icu::DateFormat::LONG, icu::DateFormat::FULL};
      const char* names[4] = {"short", "medium", "long", "full"};
      
      for (int i = 0; i != 4; ++ i) {
	std::auto_ptr<icu::DateFormat> date(icu::DateFormat::createDateInstance(date_styles[i], locale_source));
	std::auto_ptr<icu::DateFormat> time(icu::DateFormat::createTimeInstance(date_styles[i], locale_source));

	date->setLenient(true);
	time->setLenient(true);
	
	sources[std::string("date-") + names[i]].parsers.push_back(date.release());
	sources[std::string("time-") + names[i]].parsers.push_back(time.release());
      }

      for (int i = 0; i != 4; ++ i) {
	targets[std::string("date-") + names[i]].generators.push_back(icu::DateFormat::createDateInstance(date_styles[i], locale_target));
	targets[std::string("time-") + names[i]].generators.push_back(icu::DateFormat::createTimeInstance(date_styles[i], locale_target));
      }
      
      // try match!
      impl_map_type::iterator siter_end = sources.end();
      for (impl_map_type::iterator siter = sources.begin(); siter != siter_end; ++ siter) {
	impl_map_type::iterator titer = targets.find(siter->first);
	
	if (titer == targets.end()) continue;
	
	pimpls.push_back(new impl_type());
	
	pimpls.back()->parsers.swap(siter->second.parsers);
	pimpls.back()->generators.swap(titer->second.generators);
	pimpls.back()->name = siter->first;
      }
    }
    
    Date::Date(const Date& x)
      : pimpls()
    {
      for (pimpl_set_type::const_iterator iter = x.pimpls.begin(); iter != x.pimpls.end(); ++ iter)
	pimpls.push_back(new impl_type(*(*iter)));
    }
    
    Date& Date::operator=(const Date& x)
    {
      for (pimpl_set_type::iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter)
	delete *iter;
      pimpls.clear();

      for (pimpl_set_type::const_iterator iter = x.pimpls.begin(); iter != x.pimpls.end(); ++ iter)
	pimpls.push_back(new impl_type(*(*iter)));

      return *this;
    }

    Date::~Date()
    {
      for (pimpl_set_type::iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter)
	delete *iter;
      pimpls.clear();
    }
    
  };
};
