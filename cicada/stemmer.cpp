
#include "stemmer.hpp"
#include "stemmer_simple.hpp"
#include "parameter.hpp"

#include <utils/sgi_hash_map.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

namespace cicada
{

  template <typename Tp>
  struct hash_string : public utils::hashmurmur<size_t>
  {
    size_t operator()(const Tp& x) const
    {
      return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
    }
  };

  typedef boost::shared_ptr<Stemmer> stemmer_ptr_type;

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, stemmer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, stemmer_ptr_type> > > stemmer_map_type;
#else
  typedef sgi::hash_map<std::string, stemmer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, stemmer_ptr_type> > > stemmer_map_type;
#endif
  


  Stemmer& Stemmer::create(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    static __thread stemmer_map_type* __stemmers_tls = 0;
    static boost::thread_specific_ptr<stemmer_map_type> __stemmers;
    
    if (! __stemmers_tls) {
      __stemmers.reset(new stemmer_map_type());
      __stemmers_tls = __stemmers.get();
    }
    stemmer_map_type& stemmers_map = *__stemmers_tls;    
#else
    static boost::thread_specific_ptr<stemmer_map_type> __stemmers;
    
    if (! __stemmers.get())
      __stemmers.reset(new stemmer_map_type());
    
    stemmer_map_type& stemmers_map = *__stemmers;
#endif
    
    const parameter_type param(parameter);
    
    if (param.name() == "prefix") {
      int size = 0;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "size") == 0)
	  size = boost::lexical_cast<int>(piter->second);
	else
	  std::cerr << "unsupported parameter for prefix stemmer: " << piter->first << "=" << piter->second << std::endl;
      }

      if (size <= 0)
	throw std::runtime_error("invalid prefix size: " + boost::lexical_cast<std::string>(size));
      
      const std::string name = "prefix:" + boost::lexical_cast<std::string>(size);
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end())
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new StemmerPrefix(size)))).first;
      
      return *(iter->second);
    } else if (param.name() == "suffix") {
      int size = 0;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "size") == 0)
	  size = boost::lexical_cast<int>(piter->second);
	else
	  std::cerr << "unsupported parameter for suffix stemmer: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (size <= 0)
	throw std::runtime_error("invalid suffix size: " + boost::lexical_cast<std::string>(size));
      
      const std::string name = "suffix:" + boost::lexical_cast<std::string>(size);
      
      stemmer_map_type::iterator iter = stemmers_map.find(name);
      if (iter == stemmers_map.end())
	iter = stemmers_map.insert(std::make_pair(name, stemmer_ptr_type(new StemmerSuffix(size)))).first;
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
};
