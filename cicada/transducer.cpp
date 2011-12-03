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

#include "utils/lexical_cast.hpp"
#include "utils/compress_stream.hpp"

namespace cicada
{
  const char* Transducer::lists()
  {
    static const char* desc ="\
file-name: indexed grammar or plain text grammar\n\
\tmax-span=[int] maximum span (<=0 for no-constraint)\n\
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
      const path_type path = param.name();
      if (path != "-" && ! boost::filesystem::exists(path))
	throw std::runtime_error("invalid parameter: " + parameter);
      
      if (path != "-" && boost::filesystem::is_directory(path))
	return transducer_ptr_type(new GrammarStatic(parameter));
      else
	return transducer_ptr_type(new GrammarShared(parameter));
    }
  }
};
