#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <base64.hpp>

int main(int argc, char** argv)
{
  srandom(time(0) * getpid());
  
  for (int i = 0; i != 1024 * 32; ++ i) {
    const double value =  double(random()) / random();
    
    std::string stream;
    utils::encode_base64(value, std::back_inserter(stream));

    std::string encoded = utils::encode_base64(value);
    
    const double decoded = utils::decode_base64<double>(encoded);

    if (stream != encoded)
      std::cerr << "stream based encoding failed" << std::endl;

    if (decoded != value)
      std::cerr << "encoding/decoding failed" << std::endl;
  }
  
}
