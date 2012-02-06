//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <iterator>
#include <iostream>

#include <cicada/transducer.hpp>
#include <cicada/grammar_format.hpp>
#include <cicada/grammar_mutable.hpp>
#include <cicada/grammar_shared.hpp>
#include <cicada/grammar_static.hpp>
#include <cicada/grammar_simple.hpp>
#include <cicada/grammar_unknown.hpp>

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
  const char* Transducer::lists()
  {
    static const char* desc ="\
file-name: indexed grammar or plain text grammar\n\
\tmax-span=[int] maximum span (<=0 for no-constraint)\n\
\tkey-value=[true|false] store key-value format of features/attributes\n\
\tfeature0=[feature-name]\n\
\tfeature1=[feature-name]\n\
\t...\n\
\tattribute0=[attribute-name]\n\
\tattribute1=[attribute-name]\n\
\t...\n\
glue: glue rules\n\
\tgoal=[goal non-terminal]\n\
\tnon-terminal=[default non-terminal]\n\
\tfallback=[file] the list of fallback non-terminal\n\
\tstraight=[true|false] straight glue-rule\n\
\tinverted=[true|false] inverted glue-rule\n\
insertion: terminal insertion rule\n\
\tnon-terminal=[defaut non-terminal]\n\
deletion: terminal deletion rule\n\
\tnon-terminal=[defaut non-terminal]\n\
pair: terminal pair rule (for alignment composition)\n\
\tnon-terminal=[defaut non-terminal]\n\
pos: terminal pos rule (for POS annotated input) \n\
unknown: pos assignment by signature\n\
\tsignature=[signature for OOV]\n\
\tfile=[file-name] lexical grammar\n\
\tcharacter=[file-name] character model\n\
format: ICU's number/date format rules\n\
\tnon-terminal=[default non-terminal]\n\
\tformat=[formtter spec]\n\
\tremove-space=[true|false] remove space (like Chinese/Japanese)\n\
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
  typedef std::tr1::unordered_map<std::string, Transducer::transducer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, Transducer::transducer_ptr_type> > > transducer_map_type;
#else
  typedef sgi::hash_map<std::string, Transducer::transducer_ptr_type, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, Transducer::transducer_ptr_type> > > transducer_map_type;
#endif
  
  
  namespace impl
  {
    typedef utils::spinlock             mutex_type;
    typedef mutex_type::scoped_lock     lock_type;
    
    static mutex_type          __transducer_mutex;
    static transducer_map_type __transducer_map;
  };
  
#ifdef HAVE_TLS
  static __thread transducer_map_type* __transducers_tls = 0;
  static boost::thread_specific_ptr<transducer_map_type> __transducers;
#else
  static utils::thread_specific_ptr<transducer_map_type> __transducers;
#endif

  Transducer::transducer_ptr_type Transducer::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;
    typedef boost::filesystem::path path_type;
    
    const parameter_type param(parameter);
    
    if (utils::ipiece(param.name()) == "glue") {
      symbol_type goal;
      symbol_type non_terminal;
      path_type fallback_file;
      bool straight = false;
      bool inverted = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "goal")
	  goal = piter->second;
	else if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else if (utils::ipiece(piter->first) == "fallback")
	  fallback_file = piter->second;
	else if (utils::ipiece(piter->first) == "straight")
	  straight = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "inverted" || utils::ipiece(piter->first) == "invert")
	  inverted = utils::lexical_cast<bool>(piter->second);
	else
	  throw std::runtime_error("unsupported parameter for glue grammar: " + piter->first + "=" + piter->second);
      }
      
      if (int(inverted) + straight == 0)
	throw std::runtime_error("no insetion or straight glue rules?");
      if (goal.empty() || ! goal.is_non_terminal())
	throw std::runtime_error("invalid goal for glue rules? " + static_cast<const std::string&>(goal));
      
      if (! fallback_file.empty()) {
	if (fallback_file != "-" && boost::filesystem::exists(fallback_file))
	  throw std::runtime_error("invalid fallback for glue rules?" + fallback_file.string());
	
	if (! non_terminal.empty() && ! non_terminal.is_non_terminal())
	  throw std::runtime_error("invalid non_terminal for glue rules? " + static_cast<const std::string&>(non_terminal));
	
	utils::compress_istream is(fallback_file, 1024 * 1024);
	
	return transducer_ptr_type(new GrammarGlue(goal,
						   non_terminal,
						   std::istream_iterator<std::string>(is),
						   std::istream_iterator<std::string>(),
						   straight,
						   inverted));
      } else {
	if (non_terminal.empty() || ! non_terminal.is_non_terminal())
	  throw std::runtime_error("invalid non_terminal for glue rules? " + static_cast<const std::string&>(non_terminal));
	
	return transducer_ptr_type(new GrammarGlue(goal, non_terminal, straight, inverted));
      }
      
    } else if (utils::ipiece(param.name()) == "insertion") {
      symbol_type non_terminal;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else
	  throw std::runtime_error("unsupported parameter for insertion grammar: " + piter->first + "=" + piter->second);
      }
      
      if (non_terminal.empty() || ! non_terminal.is_non_terminal())
	throw std::runtime_error("invalid non-terminal for insertion grammar: " + static_cast<const std::string&>(non_terminal));
      
      return transducer_ptr_type(new GrammarInsertion(non_terminal));
    } else if (utils::ipiece(param.name()) == "deletion") {
      symbol_type non_terminal;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else
	  throw std::runtime_error("unsupported parameter for deletion grammar: " + piter->first + "=" + piter->second);
      }
      
      if (non_terminal.empty() || ! non_terminal.is_non_terminal())
	throw std::runtime_error("invalid non-terminal for deletion grammar: " + static_cast<const std::string&>(non_terminal));
      
      return transducer_ptr_type(new GrammarDeletion(non_terminal));
    } else if (utils::ipiece(param.name()) == "pair") {
      symbol_type non_terminal;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else
	  throw std::runtime_error("unsupported parameter for pair grammar: " + piter->first + "=" + piter->second);
      }
      
      if (non_terminal.empty() || ! non_terminal.is_non_terminal())
	throw std::runtime_error("invalid non-terminal for pair grammar: " + static_cast<const std::string&>(non_terminal));
      
      return transducer_ptr_type(new GrammarPair(non_terminal));
    } else if (utils::ipiece(param.name()) == "pos") {
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	throw std::runtime_error("unsupported parameter for POS grammar: " + piter->first + "=" + piter->second);
      
      return transducer_ptr_type(new GrammarPOS());
    } else if (utils::ipiece(param.name()) == "format") {
      std::string format;
      symbol_type non_terminal;
      bool remove_space = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "non-terminal")
	  non_terminal = piter->second;
	else if (utils::ipiece(piter->first) == "format")
	  format = piter->second;
	else if (utils::ipiece(piter->first) == "remove-space")
	  remove_space = utils::lexical_cast<bool>(piter->second);
	else
	  throw std::runtime_error("unsupported parameter for format grammar: " + piter->first + "=" + piter->second);
      }
      
      if (non_terminal.empty() || ! non_terminal.is_non_terminal())
	throw std::runtime_error("invalid non-terminal for format grammar: " + static_cast<const std::string&>(non_terminal));
      
      return transducer_ptr_type(new GrammarFormat(non_terminal, format, remove_space));
    } else if (utils::ipiece(param.name()) == "unknown") {
      std::string signature;
      std::string file;
      std::string character;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "signature")
	  signature = piter->second;	
	else if (utils::ipiece(piter->first) == "file")
	  file = piter->second;
	else if (utils::ipiece(piter->first) == "character" || utils::ipiece(piter->first) == "char")
	  character = piter->second;
	else
	  throw std::runtime_error("unsupported parameter for unknown grammar: " + piter->first + "=" + piter->second);
      }
      
      if (file.empty())
	throw std::runtime_error("no file? ");
      
      if (! character.empty())
	if (character != "-" && ! boost::filesystem::exists(character))
	  throw std::runtime_error("no character model file? " + character);
      
      if (signature.empty())
	throw std::runtime_error("no signature?");
      
      if (! character.empty())
	return transducer_ptr_type(new GrammarUnknown(signature, file, character));
      else
	return transducer_ptr_type(new GrammarUnknown(signature, file));
    } else {
      
#ifdef HAVE_TLS
      if (! __transducers_tls) {
	__transducers.reset(new transducer_map_type());
	__transducers_tls = __transducers.get();
      }
      transducer_map_type& transducers_map = *__transducers_tls;
#else
      if (! __transducers.get())
	__transducers.reset(new transducer_map_type());
      
      transducer_map_type& transducers_map = *__transducers;
#endif
      
      transducer_map_type::iterator iter = transducers_map.find(parameter);
      if (iter == transducers_map.end()) {	
	impl::lock_type lock(impl::__transducer_mutex);
	
	transducer_map_type::iterator iter_global = impl::__transducer_map.find(parameter);
	if (iter_global == impl::__transducer_map.end()) {
	  const path_type path = param.name();
	  if (path != "-" && ! boost::filesystem::exists(path))
	    throw std::runtime_error("invalid parameter: " + parameter);
	  
	  transducer_ptr_type ptr;
	  if (path != "-" && boost::filesystem::is_directory(path))
	    ptr.reset(new GrammarStatic(parameter));
	  else
	    ptr.reset(new GrammarShared(parameter));
	  
	  iter_global = impl::__transducer_map.insert(std::make_pair(parameter, ptr)).first;
	}
	
	iter = transducers_map.insert(*iter_global).first;
      }
      
      return iter->second;
    }
  }
};
