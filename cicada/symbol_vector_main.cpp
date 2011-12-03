#include <cstdlib>
#include <iostream>
#include <vector>

#include "symbol.hpp"
#include "symbol_vector.hpp"
#include "symbol_vector_compact.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"

typedef cicada::Symbol symbol_type;
typedef cicada::SymbolVector symbol_set_type;
typedef cicada::SymbolVectorCompact symbol_compact_type;

typedef std::vector<symbol_type, std::allocator<symbol_type> > tmp_type;

int main(int argc, char** argv)
{
  srandom(utils::random_seed());

  tmp_type tmptmp;
  
  std::string word;
  while (std::cin >> word) {
    
    tmptmp.push_back(word);
    
    if ((random() & 0x0f) == 0) {
      symbol_set_type     symbols(tmptmp.begin(), tmptmp.end());

      symbol_compact_type compacts1(symbols.begin(), symbols.end());
      symbol_compact_type compacts2(tmptmp.begin(), tmptmp.end());

      symbol_set_type symbols1(compacts1.begin(), compacts1.end());
      symbol_set_type symbols2(compacts2.begin(), compacts2.end());

      if (symbols != symbols1)
	std::cerr << "differ1" << std::endl;
      if (symbols != symbols2)
	std::cerr << "differ2" << std::endl;
      
      std::cerr << "raw: " << (symbols.size() * sizeof(symbol_type)) << " compressed: " << compacts1.size_compressed() << std::endl;
      
      
      tmptmp.clear();
    }
  }
  
}
