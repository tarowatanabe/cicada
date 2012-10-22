#include <iostream>

#include <cicada/transducer.hpp>
#include <cicada/sentence.hpp>
#include <cicada/grammar_unknown.hpp>

int main(int argc, char** argv)
{
  try {
    typedef cicada::Transducer transducer_type;
    typedef transducer_type::transducer_ptr_type transducer_ptr_type;
    typedef cicada::Sentence sentence_type;

    transducer_ptr_type grammar(transducer_type::create(argv[1]));
  
    sentence_type sentence;
  
    while (std::cin >> sentence) {
      grammar->assign(sentence);
    
      sentence_type::const_iterator siter_end = sentence.end();
      for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	transducer_type::id_type node = grammar->next(grammar->root(), *siter);
	if (node == grammar->root()) {
	  cicada::GrammarUnknown* unk = dynamic_cast<cicada::GrammarUnknown*>(&(*grammar));

	  if (! unk) {
	    std::cerr << "no unk grammar?" << std::endl;
	    continue;
	  }
	  
	  cicada::GrammarMutable* oov = dynamic_cast<cicada::GrammarMutable*>(&(*unk->grammar_oov()));

	  if (! oov) {
	    std::cerr << "no oov grammar?" << std::endl;
	    continue;
	  }
	  
	  node = oov->next(oov->root(), *siter);
	  
	  const transducer_type::rule_pair_set_type& rules = oov->rules(node);
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
	} else {
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
      }
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
}
