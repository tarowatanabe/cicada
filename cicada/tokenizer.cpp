#include "tokenizer.hpp"

#include "parameter.hpp"

#include "tokenizer/lower.hpp"
#include "tokenizer/nonascii.hpp"
#include "tokenizer/tokenize.hpp"
#include "tokenizer/nist.hpp"
#include "tokenizer/penntreebank.hpp"

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/lexical_cast.hpp>

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
tokenize: use the chain of tokenization\n\
\tlower=[true|false] perform lower casing\n\
\tnonascii=[true|false] perform non ascii character splitting\n\
\tnist=[true|false] perform NIST tokenization\n\
\tpenn=[true|false] perform Penn-treebank tokenization\n\
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
  
  
  Tokenizer& Tokenizer::create(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    static __thread tokenizer_map_type* __tokenizers_tls = 0;
    static boost::thread_specific_ptr<tokenizer_map_type> __tokenizers;
    
    if (! __tokenizers_tls) {
      __tokenizers.reset(new tokenizer_map_type());
      __tokenizers_tls = __tokenizers.get();
    }
    tokenizer_map_type& tokenizers_map = *__tokenizers_tls;    
#else
    static boost::thread_specific_ptr<tokenizer_map_type> __tokenizers;
    
    if (! __tokenizers.get())
      __tokenizers.reset(new tokenizer_map_type());
    
    tokenizer_map_type& tokenizers_map = *__tokenizers;
#endif
    
    const parameter_type param(parameter);
    
    if (param.name() == "lower") {
      const std::string name("lower");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Lower()))).first;
	iter->second->__algorithm = parameter;
      }

      return *(iter->second);
    } else if (param.name() == "nist") {
      const std::string name("nist");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Nist()))).first;
	iter->second->__algorithm = parameter;
      }

      return *(iter->second);
    } else if (param.name() == "penn") {
      const std::string name("penn");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Penntreebank()))).first;
	iter->second->__algorithm = parameter;
      }

      return *(iter->second);
    } else if (param.name() == "nonascii") {
      const std::string name("nonascii");
      
      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
      if (iter == tokenizers_map.end()) {
	iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::NonAscii()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else if (param.name() == "tokenize") {
      tokenizer_map_type::iterator iter = tokenizers_map.find(parameter);
      if (iter == tokenizers_map.end()) {
	std::auto_ptr<tokenizer::Tokenize> tokenize(new tokenizer::Tokenize());
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "lower") == 0) {
	    if (utils::lexical_cast<bool>(piter->second)) {
	      const std::string name("lower");
      
	      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
	      if (iter == tokenizers_map.end()) {
		iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Lower()))).first;
		iter->second->__algorithm = "lower";
	      }
	      
	      tokenize->insert(*(iter->second));
	    }
	  } else if (strcasecmp(piter->first.c_str(), "nist") == 0) {
	    if (utils::lexical_cast<bool>(piter->second)) {
	      const std::string name("nist");
      
	      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
	      if (iter == tokenizers_map.end()) {
		iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Nist()))).first;
		iter->second->__algorithm = "nist";
	      }
	      
	      tokenize->insert(*(iter->second));
	    }
	  } else if (strcasecmp(piter->first.c_str(), "penn") == 0) {
	    if (utils::lexical_cast<bool>(piter->second)) {
	      const std::string name("penn");
      
	      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
	      if (iter == tokenizers_map.end()) {
		iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::Penntreebank()))).first;
		iter->second->__algorithm = "penn";
	      }
	      
	      tokenize->insert(*(iter->second));
	    }
	  } else if (strcasecmp(piter->first.c_str(), "nonascii") == 0) {
	    if (utils::lexical_cast<bool>(piter->second)) {
	      const std::string name("nonascii");
	      
	      tokenizer_map_type::iterator iter = tokenizers_map.find(name);
	      if (iter == tokenizers_map.end()) {
		iter = tokenizers_map.insert(std::make_pair(name, tokenizer_ptr_type(new tokenizer::NonAscii()))).first;
		iter->second->__algorithm = "nonascii";
	      }
	      
	      tokenize->insert(*(iter->second));
	    }
	  } else
	    std::cerr << "unsupported parameter for combined tokenizer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	iter = tokenizers_map.insert(std::make_pair(parameter, tokenizer_ptr_type(tokenize.release()))).first;
	iter->second->__algorithm = parameter;
      }
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
  
};
