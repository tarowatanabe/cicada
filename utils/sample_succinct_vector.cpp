//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>

#include "utils/succinct_vector.hpp"

int main(int argc, char** argv)
{
  srandom(time(0) * getpid());

  for (int samples = 0; samples < 10; ++ samples) {
    
    std::vector<bool> stdvector(2048 + 7);
    utils::succinct_vector<std::allocator<char> > bvector;
    for (int i = 0; i < 2048 + 7; ++ i) {
      const int pos = random() % (2048 + 7);
      bvector.set(pos); 
      stdvector[pos] = true;
    }
    
    std::cout << "build..." << std::endl;
    bvector.build();

    bvector.write("tmptmp-succinct");
    utils::succinct_vector_mapped<std::allocator<char> > bvector_mapped;
    bvector_mapped.open("tmptmp-succinct");
    
    utils::succinct_vector<std::allocator<char> > bvector_copied(bvector_mapped);

    uint32_t rank1 = 0;
    uint32_t rank0 = 0;
    for (int i = 0; i < bvector.size(); ++ i) {
      std::cout << "i = " << i << " value: " << bvector.test(i) << std::endl;
      
      if (bvector.test(i) != stdvector[i])
	std::cout << "DIFFERENT from STDVECTOR!" << std::endl;
      
      if (bvector.test(i)) {
	++ rank1;
	
	if (! bvector_mapped.test(i))
	  std::cout << "DIFFERENT from mapped" << std::endl;
	if (! bvector_copied.test(i))
	  std::cout << "DIFFERENT from copied" << std::endl;
	
	size_t select1 = bvector.select(rank1, true);
	std::cout << "select1: " << select1 << std::endl;
	if (select1 != i)
	  std::cout << "DIFFER for select1: i = " << i << std::endl;

	if (select1 != bvector_mapped.select(rank1, true))
	  std::cout << "DIFFERENT from mapped" << std::endl;
	if (select1 != bvector_copied.select(rank1, true))
	  std::cout << "DIFFERENT from copied" << std::endl;
		
      } else {
	++ rank0;
	
	if (bvector_mapped.test(i))
	  std::cout << "DIFFERENT from mapped" << std::endl;
	if (bvector_copied.test(i))
	  std::cout << "DIFFERENT from copied" << std::endl;

	size_t select0 = bvector.select(rank0, false);
	std::cout << "select0: " << select0 << std::endl;
	if (select0 != i)
	  std::cout << "DIFFER for select0: i = " << i << std::endl;	

	if (select0 != bvector_mapped.select(rank0, false))
	  std::cout << "DIFFERENT from mapped" << std::endl;
	if (select0 != bvector_copied.select(rank0, false))
	  std::cout << "DIFFERENT from copied" << std::endl;
      }
      
      const size_t rank0_decoded = bvector.rank(i, false);
      const size_t rank1_decoded = bvector.rank(i, true);
      
      const size_t rank0_decoded_mapped = bvector_mapped.rank(i, false);
      const size_t rank1_decoded_mapped = bvector_mapped.rank(i, true);
      
      const size_t rank0_decoded_copied = bvector_copied.rank(i, false);
      const size_t rank1_decoded_copied = bvector_copied.rank(i, true);
      
      std::cout << "\trank1 = " << rank1_decoded << std::endl;
      if (rank1 != rank1_decoded)
	std::cout << "DIFFER for rank1!" << std::endl;      
      if (rank1 != rank1_decoded_mapped)
	std::cout << "DIFFER for rank1! (mapped)" << std::endl;      
      if (rank1 != rank1_decoded_copied)
	std::cout << "DIFFER for rank1! (copied)" << std::endl;      

      std::cout << "\trank0 = " << rank0_decoded << std::endl;
      if (rank0 != rank0_decoded)
	std::cout << "DIFFER for rank0! i = " << i << " " << rank0 << " " << rank0_decoded << std::endl;
      
      if (rank0 != rank0_decoded_mapped)
	std::cout << "DIFFER for rank0! (mapped)" << std::endl;      
      if (rank0 != rank0_decoded_copied)
	std::cout << "DIFFER for rank0! (copied)" << std::endl;      
    }
  }
}
