//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "parameter.hpp"
#include "head_finder.hpp"

#include "head/chinese.hpp"
#include "head/collins.hpp"
#include "head/collins_modified.hpp"

#include <utils/sgi_hash_map.hpp>
#include <utils/hashmurmur.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

namespace cicada
{
  
  const char* HeadFinder::lists()
  {
    static const char* desc = "\
collins: Collins head finder\n\
collins-modified: modified Collins head finder (found in Stanford parser)\n\
chinese: Chinese head finder\n\
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

  typedef boost::shared_ptr<HeadFinder> finder_ptr_type;
  
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, finder_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, finder_ptr_type> > > finder_map_type;
#else
  typedef sgi::hash_map<std::string, finder_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, finder_ptr_type> > > finder_map_type;
#endif
  
  HeadFinder& HeadFinder::create(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;

#ifdef HAVE_TLS
    static __thread finder_map_type* __finders_tls = 0;
    static boost::thread_specific_ptr<finder_map_type> __finders;
    
    if (! __finders_tls) {
      __finders.reset(new finder_map_type());
      __finders_tls = __finders.get();
    }
    finder_map_type& finders_map = *__finders_tls;    
#else
    static boost::thread_specific_ptr<finder_map_type> __finders;
    
    if (! __finders.get())
      __finders.reset(new finder_map_type());
    
    finder_map_type& finders_map = *__finders;
#endif
    
    const parameter_type param(parameter);
    if (param.name() == "collins") {
      const std::string name("collins");
      
      finder_map_type::iterator iter = finders_map.find(name);
      if (iter == finders_map.end())
	iter = finders_map.insert(std::make_pair(name, finder_ptr_type(new head::Collins()))).first;
      
      return *(iter->second);
    } else if (param.name() == "collins-modified" || param.name() == "modified-collins") {
      const std::string name("collins");
      
      finder_map_type::iterator iter = finders_map.find(name);
      if (iter == finders_map.end())
	iter = finders_map.insert(std::make_pair(name, finder_ptr_type(new head::CollinsModified()))).first;
      
      return *(iter->second);
    } else if (param.name() == "chinese") {
      const std::string name("chinese");
      
      finder_map_type::iterator iter = finders_map.find(name);
      if (iter == finders_map.end())
	iter = finders_map.insert(std::make_pair(name, finder_ptr_type(new head::Chinese()))).first;
      
      return *(iter->second);
    } else
      throw std::runtime_error("invalid parameter: " + parameter);
  }
  
};
