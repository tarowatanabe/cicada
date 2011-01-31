//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "tokenizer.hpp"

#include "parameter.hpp"

#include "tokenizer/lower.hpp"
#include "tokenizer/nonascii.hpp"
#include "tokenizer/tokenize.hpp"
#include "tokenizer/nist.hpp"
#include "tokenizer/penntreebank.hpp"
#include "tokenizer/stemmer.hpp"

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/thread_specific_ptr.hpp>
#include <utils/piece.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

namespace cicada
{

  const char* Tokenizer::lists()
  {
    static const char* desc = "\
lower: lower casing\n\
nonascii: split non ascii characters\n\
nist: NIST mteval style tokenization\n\
penn: penn-treebank stye tokenization\n\
stemmer: tokenize by stemming algorithm\n\
\talgorithm=[stemmer spec]\n\
tokenize: use the chain of tokenization\n\
\tlower=[true|false] perform lower casing\n\
\tnonascii=[true|false] perform non ascii character splitting\n\
\tnist=[true|false] perform NIST tokenization\n\
\tpenn=[true|false] perform Penn-treebank tokenization\n\
\tstemmer=[stemmer spec] perform by stemmring tokenization\n\
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
  
  typedef boost::shared_ptr<Tokenizer> tokenizer_ptr_type;
  
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, tokenizer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, tokenizer_ptr_type> > > tokenizer_map_type;
#else
  typedef sgi::hash_map<std::string, tokenizer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, tokenizer_ptr_type> > > tokenizer_map_type;
#endif

#ifdef HAVE_TLS
  static __thread tokenizer_map_type* __tokenizers_tls = 0;
  static boost::thread_specific_ptr<tokenizer_map_type> __tokenizers;
#else
  static utils::thread_specific_ptr<tokenizer_map_type> __tokenizers;
#endif
  
  Tokenizer& Tokenizer::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    if (! __tokenizers_tls) {
      __tokenizers.reset(new tokenizer_map_type());
      __tokenizers_tls = __tokenizers.get();
    }
    tokenizer_map_type& tokenizers_map = *__tokenizers_tls;    
#else
    if (! __tokenizers.get())
      __tokenizers.reset(new tokenizer_map_type());
    
    tokenizer_map_type& tokenizers_map = *__tokenizers;
#endif
    
    const parameter_type param(parameter);
    
    if (utils::ipiece(param.name()) == "lower") {
      const std::string name("lower");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Lower()))).first;
	iter->second->__algorithm = parameter;
      }

      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "nist") {
      const std::string name("nist");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Nist()))).first;
	iter->second->__algorithm = parameter;
      }

      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "penn" || utils::ipiece(param.name()) == "penntreebank") {
      const std::string name("penn");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Penntreebank()))).first;
	iter->second->__algorithm = parameter;
      }

      return *(iter->second);
    } else if (utils::ipiece(param.name() == "nonascii")) {
      const std::string name("nonascii");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::NonAscii()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "stemmer") {
      std::string algorithm;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "algorithm")
	  algorithm = piter->second;
	else
	  std::cerr << "unsupported parameter for stemming tokenizer: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (algorithm.empty())
	throw std::runtime_error("no stemming algorithm?");

      const std::string name = "stemmer:" + algorithm;
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Stemmer(&cicada::Stemmer::create(algorithm))))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (utils::ipiece(param.name()) == "tokenize") {
      tokenizer_map_type::iterator iter = tokenizers_map.find(parameter);
      if (iter == tokenizers_map.end()) {
	std::auto_ptr<tokenizer::Tokenize> tokenize(new tokenizer::Tokenize());
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "lower") {
	    if (utils::lexical_cast<bool>(piter->second))
	      tokenize->insert(create("lower"));
	  } else if (utils::ipiece(piter->first) == "nist") {
	    if (utils::lexical_cast<bool>(piter->second))
	      tokenize->insert(create("nist"));
	  } else if (utils::ipiece(piter->first) == "penn") {
	    if (utils::lexical_cast<bool>(piter->second))
	      tokenize->insert(create("penn"));
	  } else if (utils::ipiece(piter->first) == "nonascii") {
	    if (utils::lexical_cast<bool>(piter->second))
	      tokenize->insert(create("nonascii"));
	  } else if (utils::ipiece(piter->first) == "stemmer")
	    tokenize->insert(create("stemmer:algorithm=" + piter->second));
	  else
	    std::cerr << "unsupported parameter for combined tokenizer: " << piter->first << "=" << piter->second << std::endl;
	}

	if (tokenize->empty())
	  throw std::runtime_error("no tokenization for combined tokenizer?");
	
	iter = tokenizers_map.insert(std::make_pair(parameter, tokenizer_ptr_type(tokenize.release()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
  
};
