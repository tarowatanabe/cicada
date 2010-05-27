
#include "grammar_static.hpp"

#include <boost/tokenizer.hpp>

#include <utils/space_separator.hpp>

int main(int argc, char** argv)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;

  if (argc == 1) {
    std::cout << argv[0] << " [indexed-grammar]" << std::endl;
    return 1;
  }
  
  cicada::GrammarStatic grammar(argv[1]);

  std::string line;
  while (std::getline(std::cin, line)) {
    tokenizer_type tokenizer(line);
    
    cicada::GrammarStatic::id_type node = grammar.root();
    
    for (tokenizer_type::iterator iter = tokenizer.begin(); iter != tokenizer.end(); ++ iter)
      node = grammar.next(node, *iter);
    
    if (! grammar.rules(node).empty()) {
      const cicada::GrammarStatic::rule_set_type& rules = grammar.rules(node);
      
      for (cicada::GrammarStatic::rule_set_type::const_iterator riter = rules.begin(); riter != rules.end(); ++ riter)
	std::cout << *(*riter) << std::endl;
    }
  }
}

