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

namespace cicada
{
  
  const char* TreeTransducer::lists()
  {
    static const char* desc ="\
file-name: indexed tree grammar or plain text tree grammar\n\
\tcky|cyk=[true|false] indexing for CKY|CYK parsing/composition\n\
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
      const path_type path = param.name();
      if (path != "-" && ! boost::filesystem::exists(path))
	throw std::runtime_error("invalid parameter: " + parameter);
      
      if (path != "-" && boost::filesystem::is_directory(path))
	return transducer_ptr_type(new TreeGrammarStatic(parameter));
      else
	return transducer_ptr_type(new TreeGrammarShared(parameter));
    }
  }
  
};
