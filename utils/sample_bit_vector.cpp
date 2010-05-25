
#include <cstdlib>

#include <iostream>
#include <vector>

#include "bit_vector.hpp"

int main(int argc, char** argv)
{
  srandom(time(0) * getpid());
  
  for (int sample = 0; sample < 8; ++ sample) {
    utils::bit_vector<1024> bvec;
    std::vector<char> stdvec(1024);
    
    for (int i = 0; i < 512; ++ i) {
      const int pos = random() % (1024);
      bvec.set(pos); 
      stdvec[pos] = true;
    }
    
    size_t rank1 = 0;
    size_t rank0 = 0;
    for (int i = 0; i < bvec.size(); ++ i) {
      if (bvec[i] != stdvec[i])
	std::cout << "DIFFERENT" << std::endl;
      
      if (bvec.test(i)) {
	++ rank1;
	
	const size_t select1 = bvec.select(rank1, true);
	if (select1 != i)
	  std::cout << "DIFFER for select1: i = " << i << std::endl;
	
      } else {
	++ rank0;

	size_t select0 = bvec.select(rank0, false);
	if (select0 != i)
	  std::cout << "DIFFER for select0: i = " << i << std::endl;
      }

      const size_t rank0_decoded = bvec.rank(i, false);
      const size_t rank1_decoded = bvec.rank(i, true);
      
      if (rank0 != rank0_decoded)
	std::cout << "DIFFERENT RANK" << std::endl;
      if (rank1 != rank1_decoded)
	std::cout << "DIFFERENT RANK" << std::endl;
    }
  }
  
  
}
