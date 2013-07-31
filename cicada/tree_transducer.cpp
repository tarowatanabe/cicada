//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cicada/tree_transducer.hpp>
#include <cicada/tree_grammar_mutable.hpp>
#include <cicada/tree_grammar_static.hpp>
#include <cicada/tree_grammar_shared.hpp>
#include <cicada/tree_grammar_simple.hpp>

#include <cicada/parameter.hpp>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include "utils/lexical_cast.hpp"
#include "utils/compress_stream.hpp"
#include "utils/spinlock.hpp"
#include "utils/unordered_map.hpp"
#include "utils/thread_specific_ptr.hpp"

namespace cicada
{
  
  const char* TreeTransducer::lists()
  {
    static const char* desc ="\
file-name: indexed tree grammar or plain text tree grammar\n\
\tmax-span=[int] maximum span (<=0 for no-constraint)\n\
\tcky|cyk=[true|false] indexing for CKY|CYK parsing/composition\n\
\tkey-value=[true|false] store key-value format of features/attributes\n\
\tpopulate=[true|false] \"populate\" by pre-fetching\n\
\tfeature-prefix=[prefix for feature name] add prefix to the default feature name: tree-rule-table\n\
\tattribute-prefix=[prefix for attribute name] add prefix to the default attribute name: tree-rule-table\n\
\tfeature0=[feature-name]\n\
\tfeature1=[feature-name]\n\
\t...\n\
\tattribute0=[attribute-name]\n\
\tattribute1=[attribute-name]\n\
\t...\n\
fallback: fallback source-to-target, tree-to-{string,tree} transfer rule\n\
\tgoal=[default goal label] target side goal\n\
\tnon-terminal=[defaut non-terminal] target side non-terminal\n\
";
    return desc;
  }

  typedef utils::unordered_map<std::string, TreeTransducer::transducer_ptr_type, boost::hash<utils::piece>, std::equal_to<std::string>,
			       std::allocator<std::pair<const std::string, TreeTransducer::transducer_ptr_type> > >::type tree_transducer_map_type;
  
  namespace impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type               __tree_transducer_mutex;
    static tree_transducer_map_type __tree_transducer_map;
  };
  
#ifdef HAVE_TLS
  static __thread tree_transducer_map_type* __tree_transducers_tls = 0;
  static utils::thread_specific_ptr<tree_transducer_map_type> __tree_transducers;
#else
  static utils::thread_specific_ptr<tree_transducer_map_type> __tree_transducers;
#endif

  TreeTransducer::transducer_ptr_type TreeTransducer::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;
    typedef boost::filesystem::path path_type;
    
    const parameter_type param(parameter);
    
    if (utils::ipiece(param.name()) == "fallback") {
      symbol_type goal;
      symbol_type non_terminal;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else
	  throw std::runtime_error("unsupported parameter for fallback grammar: " + piter->first + "=" + piter->second);
      }
      
      if (! goal.empty() && ! goal.is_non_terminal())
	throw std::runtime_error("invalid goal for fallback grammar: " + static_cast<const std::string&>(goal));
      
      if (! non_terminal.empty() && ! non_terminal.is_non_terminal())
	throw std::runtime_error("invalid non-terminal for fallback grammar: " + static_cast<const std::string&>(non_terminal));

      if (! goal.empty() || ! non_terminal.empty())
	if (goal.empty() || non_terminal.empty())
	  throw std::runtime_error("fallback grammar should specify both of goal and non-terminal (or nothing)");
	
      return transducer_ptr_type(new TreeGrammarFallback(goal, non_terminal));
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
