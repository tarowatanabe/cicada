//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "tree_grammar_static.hpp"


int main(int argc, char** argv)
{
  typedef std::vector<cicada::Symbol> node_set_type;
  typedef std::vector<node_set_type> hyperpath_type;

  typedef cicada::TreeGrammarStatic::edge_type edge_type;
  typedef cicada::TreeGrammarStatic::id_type id_type;

  if (argc < 2) {
    std::cout << argv[0] << " grammar-file" << std::endl;
    return 1;
  }
  
  cicada::TreeGrammarStatic grammar(argv[1]);

  
  std::string line;
  cicada::TreeRule rule;

  hyperpath_type hyperpath;
  node_set_type  nodes;
  
  while (std::getline(std::cin, line)) {
    if (line.empty()) continue;
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator iter_end = line.end();
    
    if (! rule.assign(iter, iter_end)) continue;
    
    rule.hyperpath(hyperpath);

    id_type id = grammar.root();
    const edge_type edge_none = grammar.edge(cicada::Vocab::NONE);
    
    hyperpath_type::const_iterator piter_end = hyperpath.end();
    for (hyperpath_type::const_iterator piter = hyperpath.begin(); piter != piter_end; ++ piter) {
      
      nodes.clear();
      
      node_set_type::const_iterator niter_end = piter->end();
      for (node_set_type::const_iterator niter = piter->begin(); niter != niter_end; ++ niter) 
	if (*niter == cicada::Vocab::NONE) {
	  const edge_type edge_id = grammar.edge(&(*nodes.begin()), &(*nodes.end()));
	  
	  if (edge_id == edge_type())
	    id = grammar.root();
	  else
	    id = grammar.next(id, edge_id);
	  
	  if (id == grammar.root()) break;
	  
	  nodes.clear();
	} else
	  nodes.push_back(*niter);
      
      if (id == grammar.root()) break;
      
      id = grammar.next(id, edge_none);
    }
    
    if (id == grammar.root()) {
      std::cout << "no rule: " << rule << std::endl;
    } else {
      const cicada::TreeGrammarStatic::rule_pair_set_type& rules = grammar.rules(id);
      
      std::cout << "rule: " << rule << " # of entries: " << rules.size() << std::endl;
      
      cicada::TreeGrammarStatic::rule_pair_set_type::const_iterator riter_end = rules.end();
      for (cicada::TreeGrammarStatic::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter)
	std::cout << "source: " << (*riter->source) << " ||| " << (*riter->target) << std::endl;
    }
    
  }
}
