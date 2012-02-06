//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cicada/tree_transducer.hpp>
#include <cicada/tree_grammar_mutable.hpp>
#include <cicada/tree_grammar_static.hpp>
#include <cicada/tree_grammar_shared.hpp>
#include <cicada/tree_grammar_simple.hpp>

#include <cicada/parameter.hpp>

#include <boost/filesystem.hpp>

#include "utils/hashmurmur.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/compress_stream.hpp"
#include "utils/spinlock.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/thread_specific_ptr.hpp"

namespace cicada
{
  
  const char* TreeTransducer::lists()
  {
    static const char* desc ="\
file-name: indexed tree grammar or plain text tree grammar\n\
\tcky|cyk=[true|false] indexing for CKY|CYK parsing/composition\n\
\tkey-value=[true|false] store key-value format of features/attributes\n\
\tfeature0=[feature-name]\n\
\tfeature1=[feature-name]\n\
\t...\n\
\tattribute0=[attribute-name]\n\
\tattribute1=[attribute-name]\n\
\t...\n\
fallback: fallback source-to-target transfer rule\n\
\tnon-terminal=[defaut non-terminal] target side non-terminal\n\
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

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, TreeTransducer::transducer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, TreeTransducer::transducer_ptr_type> > > tree_transducer_map_type;
#else
  typedef sgi::hash_map<std::string, TreeTransducer::transducer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, TreeTransducer::transducer_ptr_type> > > tree_transducer_map_type;
#endif
  
  
  namespace impl
  {
    typedef utils::spinlock             mutex_type;
    typedef mutex_type::scoped_lock     lock_type;
    
    static mutex_type               __tree_transducer_mutex;
    static tree_transducer_map_type __tree_transducer_map;
  };
  
#ifdef HAVE_TLS
  static __thread tree_transducer_map_type* __tree_transducers_tls = 0;
  static boost::thread_specific_ptr<tree_transducer_map_type> __tree_transducers;
#else
  static utils::thread_specific_ptr<tree_transducer_map_type> __tree_transducers;
#endif

  TreeTransducer::transducer_ptr_type TreeTransducer::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;
    typedef boost::filesystem::path path_type;
    
    const parameter_type param(parameter);
    
    if (utils::ipiece(param.name()) == "fallback") {
      symbol_type non_terminal;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else
	  throw std::runtime_error("unsupported parameter for fallback grammar: " + piter->first + "=" + piter->second);
      }
      
      if (! non_terminal.empty() && ! non_terminal.is_non_terminal())
	throw std::runtime_error("invalid non-terminal for fallback grammar: " + static_cast<const std::string&>(non_terminal));
      
      return transducer_ptr_type(new TreeGrammarFallback(non_terminal));
    } else {
#ifdef HAVE_TLS
      if (! __tree_transducers_tls) {
	__tree_transducers.reset(new tree_transducer_map_type());
	__tree_transducers_tls = __tree_transducers.get();
      }
      tree_transducer_map_type& tree_transducers_map = *__tree_transducers_tls;
#else
      if (! __tree_transducers.get())
	__tree_transducers.reset(new tree_transducer_map_type());
      
      tree_transducer_map_type& tree_transducers_map = *__tree_transducers;
#endif
      
      tree_transducer_map_type::iterator iter = tree_transducers_map.find(parameter);
      if (iter == tree_transducers_map.end()) {	
	impl::lock_type lock(impl::__tree_transducer_mutex);
	
	tree_transducer_map_type::iterator iter_global = impl::__tree_transducer_map.find(parameter);
	if (iter_global == impl::__tree_transducer_map.end()) {
	  const path_type path = param.name();
	  if (path != "-" && ! boost::filesystem::exists(path))
	    throw std::runtime_error("invalid parameter: " + parameter);
	  
	  transducer_ptr_type ptr;
	  if (path != "-" && boost::filesystem::is_directory(path))
	    ptr.reset(new TreeGrammarStatic(parameter));
	  else
	    ptr.reset(new TreeGrammarShared(parameter));
	  
	  iter_global = impl::__tree_transducer_map.insert(std::make_pair(parameter, ptr)).first;
	}
	
	iter = tree_transducers_map.insert(*iter_global).first;
      }
      
      return iter->second;
    }
  }
  
};
