//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

#include <unicode/calendar.h>
#include <unicode/format.h>
#include <unicode/datefmt.h>
#include <unicode/unistr.h>
#include <unicode/locid.h>
#include <unicode/bytestream.h>
#include <unicode/smpdtfmt.h>

#include "date.hpp"

#include "utils/unordered_map.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/array_power2.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/getline.hpp"

#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

namespace cicada
{
  namespace format
  {
    struct DateImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef icu::DateFormat parser_type;
      typedef icu::DateFormat generator_type;
      
      typedef std::vector<parser_type*, std::allocator<parser_type*> >       parser_set_type;
      typedef std::vector<generator_type*, std::allocator<generator_type*> > generator_set_type;
      
      typedef cicada::Format::phrase_type     phrase_type;
      typedef cicada::Format::phrase_set_type phrase_set_type;
      
      typedef utils::hashmurmur3<size_t> hasher_type;

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
	cache_type& cache = const_cast<cache_type&>(caches[hasher_type::operator()(input.begin(), input.end(), 0) & (caches.size() - 1)]);
	if (cache.key != input) {
	  icu::UnicodeString uphrase = icu::UnicodeString::fromUTF8(input);
	  icu::UnicodeString ugenerated;
	  std::string   generated;
	  
	  cache.key = input;
	  cache.value.clear();
	  
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
	      
	      icu::Formattable   reparsed;
	      icu::ParsePosition pos(0);
	      
	      generator->parseObject(ugenerated, reparsed, pos);
	      
	      if (pos.getErrorIndex() >= 0 || pos.getIndex() != ugenerated.length() || reparsed != formattable) continue;
	      
	      generated.clear();
	      ugenerated.toUTF8String(generated);
	      
	      cache.value.push_back(std::make_pair(generated, name));
	    }
	  }
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
      generated.clear();
      
      for (pimpl_set_type::const_iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter) {
	const phrase_set_type& results = (*iter)->operator()(phrase);
	
	generated.insert(generated.end(), results.begin(), results.end());
      }
    }

    
    void Date::initialize(const path_type& path_source,
			  const path_type& path_target,
			  const std::string& locale_str_source,
			  const std::string& locale_str_target)
    {
      // pre-defined rule-set
      typedef utils::unordered_map<std::string, impl_type, boost::hash<utils::piece>, std::equal_to<std::string>,
				   std::allocator<std::pair<const std::string, impl_type> > >::type impl_map_type;
      
      const icu::Locale locale_source(locale_str_source.c_str());
      const icu::Locale locale_target(locale_str_target.c_str());
      
      if (locale_source.isBogus())
	throw std::runtime_error("invalid locale: " + locale_str_source);
      if (locale_target.isBogus())
	throw std::runtime_error("invalid locale: " + locale_str_target);
      
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<icu::Calendar> calendar(icu::Calendar::createInstance(locale_source, status));
      if (U_FAILURE(status))
        throw std::runtime_error(std::string("Calendar::createInstance(): ") + u_errorName(status));
      
      status = U_ZERO_ERROR;
      calendar->setTime(icu::Calendar::getNow(), status);
      if (U_FAILURE(status))
        throw std::runtime_error(std::string("Calendar::setTime(): ") + u_errorName(status));
      
      impl_map_type sources;
      impl_map_type targets;
      
      const icu::DateFormat::EStyle date_styles[8] = { icu::DateFormat::kShort,
						       icu::DateFormat::kMedium,
						       icu::DateFormat::kLong,
						       icu::DateFormat::kFull,
						       icu::DateFormat::kShortRelative,
						       icu::DateFormat::kMediumRelative,
						       icu::DateFormat::kLongRelative,
						       icu::DateFormat::kFullRelative,};
      const char* names[8] = {"short",
			      "medium",
			      "long",
			      "full",
			      "short-relative",
			      "medium-relative",
			      "long-relative",
			      "full-relative"};
      
      for (int i = 0; i != sizeof(date_styles) / sizeof(icu::DateFormat::EStyle); ++ i) {
	std::auto_ptr<icu::DateFormat> date(icu::DateFormat::createDateInstance(date_styles[i], locale_source));
	std::auto_ptr<icu::DateFormat> time(icu::DateFormat::createTimeInstance(date_styles[i], locale_source));
	
	date->setLenient(true);
	time->setLenient(true);
	
	date->setCalendar(*calendar);
	time->setCalendar(*calendar);
	
	sources[std::string("date-") + names[i]].parsers.push_back(date.release());
	sources[std::string("time-") + names[i]].parsers.push_back(time.release());
      }

      for (int i = 0; i != sizeof(date_styles) / sizeof(icu::DateFormat::EStyle); ++ i) {
	std::auto_ptr<icu::DateFormat> date(icu::DateFormat::createDateInstance(date_styles[i], locale_target));
	std::auto_ptr<icu::DateFormat> time(icu::DateFormat::createTimeInstance(date_styles[i], locale_target));
	
	date->setLenient(true);
	time->setLenient(true);
	
	date->setCalendar(*calendar);
	time->setCalendar(*calendar);
	
	targets[std::string("date-") + names[i]].generators.push_back(date.release());
	targets[std::string("time-") + names[i]].generators.push_back(time.release());
      }
      
      // user-defined..
      if (! path_source.empty() && ! path_target.empty()) {
	if (! boost::filesystem::exists(path_source))
	  throw std::runtime_error("no file? " + path_source.string());
	
	if (! boost::filesystem::exists(path_target))
	  throw std::runtime_error("no file? " + path_target.string());
	
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	typedef std::string::const_iterator iter_type;
	
	qi::rule<iter_type, std::string()> rule_tag(qi::lexeme[+(standard::char_ - standard::space)]);
	qi::rule<iter_type, std::string()> pattern_tag(qi::lexeme[+(standard::char_)]);
	
	{
	  utils::compress_istream is(path_source);
	  std::string line;
	  
	  while (utils::getline(is, line)) {
	    std::string::const_iterator iter(line.begin());
	    std::string::const_iterator end(line.end());

	    std::string tag;
	    std::string pattern;
	    
	    if (! qi::parse(iter, end, rule_tag >> qi::omit[+standard::space] >> pattern_tag, tag, pattern))
	      continue;
	    
	    boost::algorithm::trim(pattern);
	    
	    UErrorCode status = U_ZERO_ERROR;
	    std::auto_ptr<icu::DateFormat> date(new icu::SimpleDateFormat(icu::UnicodeString::fromUTF8(pattern), locale_source, status));
	    if (U_FAILURE(status))
	      throw std::runtime_error(std::string("SimpleDateFormat(): ") + u_errorName(status));
	    
	    
	    date->setLenient(true);
	    date->setCalendar(*calendar);
	    
	    sources[tag].parsers.push_back(date.release());
	  }
	}
	
	{
	  utils::compress_istream is(path_target);
	  std::string line;
	  
	  while (utils::getline(is, line)) {
	    std::string::const_iterator iter(line.begin());
	    std::string::const_iterator end(line.end());

	    std::string tag;
	    std::string pattern;
	    
	    if (! qi::parse(iter, end, rule_tag >> qi::omit[+standard::space] >> pattern_tag, tag, pattern))
	      continue;
	    
	    boost::algorithm::trim(pattern);
	    
	    UErrorCode status = U_ZERO_ERROR;
	    std::auto_ptr<icu::DateFormat> date(new icu::SimpleDateFormat(icu::UnicodeString::fromUTF8(pattern), locale_target, status));
	    if (U_FAILURE(status))
	      throw std::runtime_error(std::string("SimpleDateFormat(): ") + u_errorName(status));
	    
	    
	    date->setLenient(true);
	    date->setCalendar(*calendar);
	    
	    targets[tag].generators.push_back(date.release());
	  }
	}
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

    Date::~Date()
    {
      for (pimpl_set_type::iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter)
	delete *iter;
      pimpls.clear();
    }
    
  };
};
