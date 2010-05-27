
#include "grammar_mutable.hpp"

int main(int argc, char** argv)
{
  cicada::GrammarMutable grammar("-,feature0=prob-soruce-target,feature2=lex-target-source,feature2=prob-target-source,feature3=lex-target-source");
  
  grammar.insert("[x] ||| good boy |||  bad boy ||| feature1=0.4");
  
  typedef cicada::GrammarMutable::id_type id_type;
  typedef cicada::GrammarMutable::rule_type rule_type;
  typedef cicada::GrammarMutable::rule_set_type rule_set_type;
  
  id_type node = grammar.root();
  node = grammar.next(node, std::string("good"));
  
  if (! grammar.rules(node).empty()) {
    const rule_set_type& rules = grammar.rules(node);
    
    for (rule_set_type::const_iterator riter = rules.begin(); riter != rules.end(); ++ riter)
      std::cout << *(*riter) << std::endl;
  }

  if (grammar.has_next(node)) {
    node = grammar.next(node, std::string("boy"));
    
    if (! grammar.rules(node).empty()) {
      const rule_set_type& rules = grammar.rules(node);
      
      for (rule_set_type::const_iterator riter = rules.begin(); riter != rules.end(); ++ riter)
	std::cout << *(*riter) << std::endl;
    }
  }

}
