//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "stemmer.hpp"
#include "stemmer/arabic.hpp"
#include "stemmer/prefix.hpp"
#include "stemmer/suffix.hpp"
#include "stemmer/katakana.hpp"
#include "stemmer/latin.hpp"
#include "stemmer/digit.hpp"
#include "stemmer/lower.hpp"
#include "stemmer/upper.hpp"
#include "stemmer/halfwidth.hpp"
#include "stemmer/fullwidth.hpp"
#include "stemmer/simplified.hpp"
#include "stemmer/traditional.hpp"
#include "stemmer/nfkc.hpp"
#include "stemmer/nfkd.hpp"
#include "stemmer/snowball.hpp"

#include "parameter.hpp"

#include <utils/unordered_map.hpp>
#include <utils/thread_specific_ptr.hpp>
#include <utils/piece.hpp>
#include "utils/lexical_cast.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

namespace cicada
{

  const char* Stemmer::lists()
  {
    static const char* desc = "\
snowball: snowball stemming\n\
\tlanguage=[language] stemming algorithm (en, de etc.)\n\
arabic: Arabic stemming\n\
prefix: taking prefix of letters\n\
\tsize=[int] prefix size\n\
suffix: taking suffix of letters\n\
\tsize=[int] suffix size\n\
digit: digits normalized to @\n\
katakana: katakana\n\
latin: romanization\n\
lower: lower casing\n\
upper: upper casing\n\
halfwidth: Fullwidth-Halfwidth\n\
fullwidth: Halfwidth-Fullwidth\n\
traditional: Simplified-Traditional\n\
simplified: Traditioanl-Simplified\n\
nfkc: NFKC\n\
nfkd: NFKD\n\
";
    return desc;
  }
  
  
  typedef boost::shared_ptr<Stemmer> stemmer_ptr_type;

  typedef utils::unordered_map<std::string, stemmer_ptr_type, boost::hash<utils::piece>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, stemmer_ptr_type> > >::type stemmer_map_type;

#ifdef HAVE_TLS
  static __thread stemmer_map_type* __stemmers_tls = 0;
  static utils::thread_specific_ptr<stemmer_map_type> __stemmers;
#else
  static utils::thread_specific_ptr<stemmer_map_type> __stemmers;
#endif
  
  Stemmer& Stemmer::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    if (! __stemmers_tls) {
      __stemmers.reset(new stemmer_map_type());
      __stemmers_tls = __stemmers.get();
    }
    stemmer_map_type& stemmers_map = *__stemmers_tls;    
#else
    if (! __stemmers.get())
      __stemmers.reset(new stemmer_map_type());
    
    stemmer_map_type& stemmers_map = *__stemmers;
#endif
    
    const parameter_type param(parameter);
    
    if (utils::ipiece(param.name()) == "prefix") {
      int size = 0;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "size")
	  size = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "unsupported parameter for prefix stemmer: " << piter->first << "=" << piter->second << std::endl;
      }

      if (size <= 0)
	throw std::runtime_error("invalid prefix size: " + utils::lexical_cast<std::string>(size));
      
      const std::string name = "prefix:" + utils::lexical_cast<std::string>(size);
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Prefix(size)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "suffix") {
      int size = 0;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "size")
	  size = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "unsupported parameter for suffix stemmer: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (size <= 0)
	throw std::runtime_error("invalid suffix size: " + utils::lexical_cast<std::string>(size));
      
      const std::string name = "suffix:" + utils::lexical_cast<std::string>(size);
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Suffix(size)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "digit" || utils::ipiece(param.name()) == "digits") {
      const std::string name("digit");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Digit()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "katakana") {
      const std::string name("katakana");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Katakana()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "latin") {
      const std::string name("latin");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Latin()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "lower") {
      const std::string name("lower");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Lower()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "upper") {
      const std::string name("upper");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Upper()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "simplified") {
      const std::string name("simplified");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Simplified()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "traditional") {
      const std::string name("traditional");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Traditional()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "halfwidth" || utils::ipiece(param.name()) == "half") {
      const std::string name("halfwidth");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Halfwidth()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "fullwidth" || utils::ipiece(param.name()) == "full") {
      const std::string name("fullwidth");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Fullwidth()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "nfkc") {
      const std::string name("nfkc");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::NFKC()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "nfkd") {
      const std::string name("nfkd");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::NFKD()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "arabic") {
      const std::string name("arabic");
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Arabic()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "snowball") {
      
      std::string language;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "lang"
	    || utils::ipiece(piter->first) == "language"
	    || utils::ipiece(piter->first) == "algorithm")
	  language = piter->second;
	else
	  std::cerr << "unsupported parameter for snowball stemmer: " << piter->first << "=" << piter->second << std::endl;
      }

      if (language.empty())
	throw std::runtime_error("no stemming algorithm?");
      
      const std::string name = "snowball:" + language;
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end()) {
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new stemmer::Snowball(language)))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
};
