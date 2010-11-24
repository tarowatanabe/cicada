
#include "hypergraph.hpp"
#include "tree_grammar_mutable.hpp"

#include <boost/tokenizer.hpp>

#include <utils/space_separator.hpp>

typedef cicada::HyperGraph hypergraph_type;


int main(int argc, char** argv)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;
  
  if (argc < 2) {
    std::cout << argv[0] << " grammar-file" << std::endl;
    return 1;
  }
  
  cicada::TreeGrammarMutable grammar(argv[1]);
  
  cicada::TreeRule rule;
}
