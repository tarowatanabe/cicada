#include "tokenizer.hpp"

#include "parameter.hpp"

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

namespace cicada
{

  const char* Tokenizer::lists()
  {
    static const char* desc = "\
lower: lower casing\n\
nonascii: split non ascii characters\n\
tokenize: use the chain of tokenization\n\
\tlower=[true|false] perform lower casing\n\
\tnonascii=[true|false] perform non ascii character splitting\n\
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
      
      
    } else if (param.name() == "nonascii") {
      
    } else if (param.name() == "tokenize") {
      
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
  
};
