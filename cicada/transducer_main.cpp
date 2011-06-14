#include <iostream>
#include <vector>

#include <cicada/transducer.hpp>
#include <cicada/sentence.hpp>

int main(int argc, char** argv)
{
  typedef cicada::Transducer transducer_type;
  typedef transducer_type::transducer_ptr_type transducer_ptr_type;
  typedef cicada::Sentence sentence_type;
  typedef std::vector<transducer_type::id_type> id_set_type;

  transducer_ptr_type grammar(transducer_type::create(argv[1]));
  
  sentence_type sentence;
  
  while (std::cin >> sentence) {
    
    id_set_type nodes;
    id_set_type nodes_next;
    
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
      
      nodes.push_back(grammar->root());
      nodes_next.clear();
      
      id_set_type::const_iterator niter_end = nodes.end();
      for (id_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	transducer_type::id_type node = grammar->next(*niter, *siter);
	
	if (node == grammar->root()) {
	  std::cerr << "invalid node?" << std::endl;
	  continue;
	}
	
	nodes_next.push_back(node);
	
	const transducer_type::rule_pair_set_type& rules = grammar->rules(node);
	if (rules.empty())
	  std::cerr << "no rules?" << std::endl;
	
	transducer_type::rule_pair_set_type::const_iterator riter_end = rules.end();
	for (transducer_type::rule_pair_set_type::const_iterator riter = rules.begin(); riter != riter_end; ++ riter) {
	  std::cout << *(riter->source) << " ||| " << *(riter->target) << " |||";
	  transducer_type::feature_set_type::const_iterator fiter_end = riter->features.end();
	  for (transducer_type::feature_set_type::const_iterator fiter = riter->features.begin(); fiter != fiter_end; ++ fiter)
	    std::cout << ' ' << fiter->first << '=' << fiter->second;
	  std::cout << " ||| " << riter->attributes << std::endl;
	}
      }
      
      nodes.swap(nodes_next);
    }
  }
}
