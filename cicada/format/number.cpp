//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <vector>
#include <string>
#include <set>
#include <memory>
#include <stdexcept>

#include <unicode/curramt.h>
#include <unicode/format.h>
#include <unicode/numfmt.h>
#include <unicode/rbnf.h>
#include <unicode/unistr.h>
#include <unicode/locid.h>
#include <unicode/bytestream.h>
#include <unicode/resbund.h>
#include <unicode/ures.h>
#include <unicode/udata.h>

#define U_ICUDATA_RBNF U_ICUDATA_NAME U_TREE_SEPARATOR_STRING "rbnf"

#include "number.hpp"

#include "utils/sgi_hash_map.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/array_power2.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>

namespace cicada
{
  namespace format
  {
    struct NumberImpl : public utils::hashmurmur<size_t>
    {
    public:
      typedef icu::NumberFormat parser_type;
      typedef icu::NumberFormat generator_type;
      
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
      NumberImpl() :currency(false){}
      NumberImpl(const NumberImpl& x) { assign(x); }
      NumberImpl& operator=(const NumberImpl& x)
      {
	assign(x);
	return *this;
      }
      ~NumberImpl() { clear(); }
      
    public:
      const phrase_set_type& operator()(const phrase_type& input) const
      {
	typedef std::set<std::string, std::less<std::string>, std::allocator<std::string> > phrase_unique_type;
	
	cache_type& cache = const_cast<cache_type&>(caches[hasher_type::operator()(input.begin(), input.end(), 0) & (caches.size() - 1)]);
	if (cache.key != input) {
	  UnicodeString uphrase = UnicodeString::fromUTF8(input);
	  UnicodeString ugenerated;
	  std::string   generated;
	
	  phrase_unique_type uniques;
	
	  parser_set_type::const_iterator piter_end = parsers.end();
	  for (parser_set_type::const_iterator piter = parsers.begin(); piter != piter_end; ++ piter) {
	    const parser_type* parser = *piter;
	    icu::Formattable   formattable;
	    icu::ParsePosition pos(0);
	  
	    if (currency)
	      parser->parseCurrency(uphrase, formattable, pos);
	    else
	      parser->parse(uphrase, formattable, pos);
	  
	    if (pos.getErrorIndex() >= 0 || pos.getIndex() != uphrase.length()) continue;
	  
	    generator_set_type::const_iterator giter_end = generators.end();
	    for (generator_set_type::const_iterator giter = generators.begin(); giter != giter_end; ++ giter) {
	      const generator_type* generator = *giter;
	      UErrorCode    status(U_ZERO_ERROR);
	    
	      ugenerated.remove();
	      generator->format(formattable, ugenerated, status);
	    
	      if (U_FAILURE(status)) continue;
	    
	      generated.clear();
	      ugenerated.toUTF8String(generated);
	      uniques.insert(generated);

	      //std::cerr << name << ": " << generated << std::endl;
	    }
	  }
	  
	  cache.key = input;
	  cache.value.clear();
	  cache.value.insert(cache.value.end(), uniques.begin(), uniques.end());
	}
	
	return cache.value;
      }
      
    public:
      void assign(const NumberImpl& x)
      {
	clear();
	
	parser_set_type::const_iterator piter_end = x.parsers.end();
	for (parser_set_type::const_iterator piter = x.parsers.begin(); piter != piter_end; ++ piter)
	  parsers.push_back(dynamic_cast<parser_type*>((*piter)->clone()));
	
	generator_set_type::const_iterator giter_end = x.generators.end();
	for (generator_set_type::const_iterator giter = x.generators.begin(); giter != giter_end; ++ giter)
	  generators.push_back(dynamic_cast<generator_type*>((*giter)->clone()));
	
	name = x.name;
	currency = x.currency;
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
      bool currency;
    };
    
    void Number::operator()(const phrase_type& phrase, phrase_set_type& generated) const
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

    template <typename Format>
    inline
    bool has_rbnf_rule_set(const Format& format, const char* name)
    {
      const int num_rules = format.getNumberOfRuleSetNames();
      for (int i = 0; i < num_rules; ++ i)
	if (format.getRuleSetName(i).indexOf(name) >= 0)
	  return true;
      return false;
    }

    template <typename Locale, typename Name>
    inline
    icu::NumberFormat* create_rbnf_instance(const Locale& locale, const Name& name)
    {
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(URBNF_SPELLOUT, locale, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RuleBasedNumberFormat: ") + u_errorName(status));
      
      status = U_ZERO_ERROR;
      rbnf->setDefaultRuleSet(name, status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RuleBasedNumberFormat::setDefaultRuleSet: ") + u_errorName(status));
      
      rbnf->setLenient(true);
      
      return rbnf.release();
    }
    
    
    
    void Number::initialize(const path_type& path_source,
			    const path_type& path_target,
			    const std::string& locale_str_source,
			    const std::string& locale_str_target)
    {
      // pre-defined rule-set
#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<std::string, impl_type, boost::hash<std::string>, std::equal_to<std::string>,
				      std::allocator<std::pair<const std::string, impl_type> > > impl_map_type;
#else
      typedef sgi::hash_map<std::string, impl_type, boost::hash<std::string>, std::equal_to<std::string>,
			    std::allocator<std::pair<const std::string, impl_type> > > impl_map_type;
#endif
      
      const icu::Locale locale_source(locale_str_source.c_str());
      const icu::Locale locale_target(locale_str_target.c_str());
      
      if (locale_source.isBogus())
	throw std::runtime_error("invalid locale: " + locale_str_source);
      if (locale_target.isBogus())
	throw std::runtime_error("invalid locale: " + locale_str_target);
      
      impl_map_type sources;
      impl_map_type targets;
      
      UErrorCode status = U_ZERO_ERROR;
      std::auto_ptr<icu::RuleBasedNumberFormat> rbnf_source(new icu::RuleBasedNumberFormat(URBNF_SPELLOUT, locale_source, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
      
      UnicodeString rules_rbnf_source;
      const char* rules_tag = "RBNFRules";
      const char* fmt_tag   = "SpelloutRules";
      
      status = U_ZERO_ERROR;
      icu::ResourceBundle bundle(U_ICUDATA_RBNF, locale_source, status);
      
      if (U_SUCCESS(status)) {
	icu::ResourceBundle bundle_rbnf = bundle.getWithFallback(rules_tag, status);
        if (U_FAILURE(status))
	  throw std::runtime_error("no rbnf?");
	
	icu::ResourceBundle bundle_rules = bundle_rbnf.getWithFallback(fmt_tag, status);
        if (U_FAILURE(status))
	  throw std::runtime_error("no rules?");
	
	while (bundle_rules.hasNext() && U_SUCCESS(status))
	  rules_rbnf_source += bundle_rules.getNextString(status);
	
	if (U_FAILURE(status))
	  throw std::runtime_error("no string?");
      } else 
	throw std::runtime_error("no RBNF rules?");
      
      if (has_rbnf_rule_set(*rbnf_source, "cardinal")) {
	UnicodeString rules(rules_rbnf_source);
	
	rules.findAndReplace("%spellout-numbering", "%%spellout-numbering");
	rules.findAndReplace("%spellout-ordinal",   "%%spellout-ordinal");
	
	UErrorCode status = U_ZERO_ERROR;
	UParseError perror;
	std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(rules, locale_source, perror, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
	
	rbnf->setLenient(true);
	sources["cardinal"].parsers.push_back(rbnf.release());
      }
      
      if (has_rbnf_rule_set(*rbnf_source, "ordinal")) {
	UnicodeString rules(rules_rbnf_source);
	
	rules.findAndReplace("%spellout-numbering", "%%spellout-numbering");
	rules.findAndReplace("%spellout-cardinal",  "%%spellout-cardinal");
	
	UErrorCode status = U_ZERO_ERROR;
	UParseError perror;
	std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(rules, locale_source, perror, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
	
	rbnf->setLenient(true);
	sources["ordinal"].parsers.push_back(rbnf.release());
      }
      
      if (has_rbnf_rule_set(*rbnf_source, "numbering")) {
	UnicodeString rules(rules_rbnf_source);
	
	rules.findAndReplace("%spellout-numbering-year",  "%%spellout-numbering-year");
	rules.findAndReplace("%spellout-ordinal",  "%%spellout-ordinal");
	rules.findAndReplace("%spellout-cardinal", "%%spellout-cardinal");

	UErrorCode status = U_ZERO_ERROR;
	UParseError perror;
	std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(rules, locale_source, perror, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
	
	rbnf->setLenient(true);
	sources["numbering"].parsers.push_back(rbnf.release());
      }
      
      // user defined rules for source side
      if (! path_source.empty()) {
	if (! boost::filesystem::exists(path_source))
	  throw std::runtime_error("no file? " + path_source.string());
	
	UnicodeString rules;
	
	utils::compress_istream is(path_target);
	std::string line;
	while (std::getline(is, line)) {
	  rules += icu::UnicodeString::fromUTF8(line);
	  rules += '\n';
	}
	
	status = U_ZERO_ERROR;
	UParseError perror;
	std::auto_ptr<icu::RuleBasedNumberFormat> nf_rule(new icu::RuleBasedNumberFormat(rules, locale_source, perror, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("RuleBasedNumberFormat: ") + u_errorName(status)
				   + std::string(" offset: ") + utils::lexical_cast<std::string>(perror.offset));
	
	if (has_rbnf_rule_set(*nf_rule, "cardinal")) {
	  UnicodeString rules_local(rules);
	  
	  rules_local.findAndReplace("%spellout-numbering", "%%spellout-numbering");
	  rules_local.findAndReplace("%spellout-ordinal",   "%%spellout-ordinal");
	  
	  UErrorCode status = U_ZERO_ERROR;
	  UParseError perror;
	  std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(rules_local, locale_source, perror, status));
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
	  
	  rbnf->setLenient(true);
	  sources["cardinal"].parsers.push_back(rbnf.release());
	}
	
	if (has_rbnf_rule_set(*nf_rule, "ordinal")) {
	  UnicodeString rules_local(rules);
	  
	  rules_local.findAndReplace("%spellout-numbering", "%%spellout-numbering");
	  rules_local.findAndReplace("%spellout-cardinal",  "%%spellout-cardinal");
	  
	  UErrorCode status = U_ZERO_ERROR;
	  UParseError perror;
	  std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(rules_local, locale_source, perror, status));
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
	  
	  rbnf->setLenient(true);
	  sources["ordinal"].parsers.push_back(rbnf.release());
	}
	
	if (has_rbnf_rule_set(*nf_rule, "numbering")) {
	  UnicodeString rules_local(rules);
	
	  rules_local.findAndReplace("%spellout-numbering-year",  "%%spellout-numbering-year");
	  rules_local.findAndReplace("%spellout-ordinal",  "%%spellout-ordinal");
	  rules_local.findAndReplace("%spellout-cardinal", "%%spellout-cardinal");

	  UErrorCode status = U_ZERO_ERROR;
	  UParseError perror;
	  std::auto_ptr<icu::RuleBasedNumberFormat> rbnf(new icu::RuleBasedNumberFormat(rules_local, locale_source, perror, status));
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
	  
	  rbnf->setLenient(true);
	  sources["numbering"].parsers.push_back(rbnf.release());
	}
      }

    
      status = U_ZERO_ERROR;
      std::auto_ptr<NumberFormat> nf_source(icu::NumberFormat::createInstance(locale_source, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("NumberFormat::createInstance: ") + u_errorName(status));

      sources["numbering"].parsers.push_back(dynamic_cast<impl_type::parser_type*>(nf_source->clone()));
      sources["cardinal"].parsers.push_back(dynamic_cast<impl_type::parser_type*>(nf_source->clone()));
      sources["any"].parsers.push_back(dynamic_cast<impl_type::parser_type*>(nf_source->clone()));
      
      status = U_ZERO_ERROR;
      std::auto_ptr<icu::RuleBasedNumberFormat> rbnf_target(new icu::RuleBasedNumberFormat(URBNF_SPELLOUT, locale_target, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
      
      const int num_rules = rbnf_target->getNumberOfRuleSetNames();
      for (int i = 0; i < num_rules; ++ i) {
	const icu::UnicodeString uname = rbnf_target->getRuleSetName(i);
	
	if (uname.indexOf("numbering") >= 0 && uname.indexOf("year") < 0)
	  targets["numbering"].generators.push_back(create_rbnf_instance(locale_target, uname));
	else if (uname.indexOf("ordinal") >= 0)
	  targets["ordinal"].generators.push_back(create_rbnf_instance(locale_target, uname));
	else if (uname.indexOf("cardinal") >= 0)
	  targets["cardinal"].generators.push_back(create_rbnf_instance(locale_target, uname));
	else
	  targets["any"].generators.push_back(create_rbnf_instance(locale_target, uname));
      }

      
      // user defined rules for target side
      if (! path_target.empty()) {
	if (! boost::filesystem::exists(path_target))
	  throw std::runtime_error("no file? " + path_target.string());
	
	UnicodeString rules;
	
	utils::compress_istream is(path_target);
	std::string line;
	while (std::getline(is, line)) {
	  rules += icu::UnicodeString::fromUTF8(line);
	  rules += '\n';
	}
	
	status = U_ZERO_ERROR;
	UParseError perror;
	std::auto_ptr<icu::RuleBasedNumberFormat> nf_rule(new icu::RuleBasedNumberFormat(rules, locale_target, perror, status));
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("RuleBasedNumberFormat: ") + u_errorName(status)
				   + std::string(" offset: ") + utils::lexical_cast<std::string>(perror.offset));
	
	int32_t num_rule_set = nf_rule->getNumberOfRuleSetNames();
	for (int32_t i = 0; i < num_rule_set; ++ i) {
	  const icu::UnicodeString uname = nf_rule->getRuleSetName(i);

	  status = U_ZERO_ERROR;
	  UParseError perror;
	  std::auto_ptr<icu::RuleBasedNumberFormat> formatter(new icu::RuleBasedNumberFormat(rules, locale_target, perror, status));
	  
	  status = U_ZERO_ERROR;
	  formatter->setDefaultRuleSet(uname, status);
	  
	  if (uname.indexOf("numbering") >= 0 && uname.indexOf("year") < 0)
	    targets["numbering"].generators.push_back(formatter.release());
	  else if (uname.indexOf("ordinal") >= 0)
	    targets["ordinal"].generators.push_back(formatter.release());
	  else if (uname.indexOf("cardinal") >= 0)
	    targets["cardinal"].generators.push_back(formatter.release());
	  else
	    targets["any"].generators.push_back(formatter.release());
	}
      }
      
      status = U_ZERO_ERROR;
      std::auto_ptr<NumberFormat> nf_target(icu::NumberFormat::createInstance(locale_target, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("NumberFormat::createInstance: ") + u_errorName(status));
      
      targets["numbering"].generators.push_back(dynamic_cast<impl_type::generator_type*>(nf_target->clone()));
      targets["cardinal"].generators.push_back(dynamic_cast<impl_type::generator_type*>(nf_target->clone()));
      targets["any"].generators.push_back(dynamic_cast<impl_type::generator_type*>(nf_target->clone()));
      
#if 0
      // currently, we will disable currency parsing/generation
      status = U_ZERO_ERROR;
      std::auto_ptr<NumberFormat> curr_source(icu::NumberFormat::createCurrencyInstance(locale_source, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("NumberFormat::createCurrencyInstance: ") + u_errorName(status));
      
      status = U_ZERO_ERROR;
      std::auto_ptr<NumberFormat> curr_target(icu::NumberFormat::createCurrencyInstance(locale_target, status));
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("NumberFormat::createCurrencyInstance: ") + u_errorName(status));
      
      sources["currency"].parsers.push_back(curr_source.release());
      targets["currency"].generators.push_back(curr_target.release());
#endif

      // try match!
      impl_map_type::iterator siter_end = sources.end();
      for (impl_map_type::iterator siter = sources.begin(); siter != siter_end; ++ siter) {
	impl_map_type::iterator titer = targets.find(siter->first);
	
	if (titer == targets.end()) {
	  if (siter->first == "ordinal") {
	    titer = targets.find("cardinal");
	    if (titer != targets.end()) {
	      pimpls.push_back(new impl_type());
	  
	      pimpls.back()->parsers.swap(siter->second.parsers);
	      pimpls.back()->name = siter->first;
	      
	      impl_type::generator_set_type::const_iterator giter_end = titer->second.generators.end();
	      for (impl_type::generator_set_type::const_iterator giter = titer->second.generators.begin(); giter != giter_end; ++ giter)
		pimpls.back()->generators.push_back(dynamic_cast<impl_type::generator_type*>((*giter)->clone()));
	    }
	  }
	} else {
	  pimpls.push_back(new impl_type());
	  
	  pimpls.back()->parsers.swap(siter->second.parsers);
	  pimpls.back()->generators.swap(titer->second.generators);
	  pimpls.back()->name = siter->first;
	  pimpls.back()->currency = (siter->first == "currency");
	}
      }
    }
    
    Number::Number(const Number& x)
      : pimpls()
    {
      for (pimpl_set_type::const_iterator iter = x.pimpls.begin(); iter != x.pimpls.end(); ++ iter)
	pimpls.push_back(new impl_type(*(*iter)));
    }
    
    Number& Number::operator=(const Number& x)
    {
      for (pimpl_set_type::iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter)
	delete *iter;
      pimpls.clear();

      for (pimpl_set_type::const_iterator iter = x.pimpls.begin(); iter != x.pimpls.end(); ++ iter)
	pimpls.push_back(new impl_type(*(*iter)));

      return *this;
    }

    Number::~Number()
    {
      for (pimpl_set_type::iterator iter = pimpls.begin(); iter != pimpls.end(); ++ iter)
	delete *iter;
      pimpls.clear();
    }
    
  };
};
