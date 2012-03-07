//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>

#include "snappy_vector.hpp"
#include "snappy_device.hpp"

int main(int argc, char**argv)
{
  srandom(time(0) * getpid());

  for (int i = 0; i < 16; ++ i) {
    std::cerr << "iteration: " << i << std::endl;
    
    std::vector<double> floats;
    for (int j = 0; j < 1024 * 1024 + 10; ++ j)
      floats.push_back(random() / random());
    
    boost::iostreams::filtering_ostream os;
    os.push(utils::snappy_sink<>("tmptmp.snappy.double"));

    for (int j = 0; j != floats.size(); ++ j)
      os.write((char*) &floats[j], sizeof(double));
    os.pop();
    
    utils::snappy_vector_mapped<double> snapped("tmptmp.snappy.double");

    std::cerr << "raw: " << snapped.size_bytes()
	      << " compressed: " << snapped.size_compressed()
	      << " cached: " << snapped.size_cache()
	      << std::endl;
    
    for (int j = 0; j != floats.size(); ++ j)
      if (floats[j] != snapped[j])
	std::cerr << "differ: " << "j = " << j << " " << floats[j] << " " << snapped[j] << std::endl;
    
    for (int j = 0; j != floats.size(); ++ j) {
      const int pos = random() % floats.size();
      
      if (floats[pos] != snapped[pos])
	std::cerr << "differ: " << "j = " << pos << " " << floats[pos] << " " << snapped[pos] << std::endl;
    }
    
  }
}
