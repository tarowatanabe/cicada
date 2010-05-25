
#include "utils/compress_stream.hpp"

int main(int argc, char** argv)
{
  if (argc == 2) {
    utils::compress_istream is(argv[1]);
    
    char buffer[4096];
    
    do {
      is.read(buffer, 4096);
      if (is.gcount() > 0)
	std::cout.write(buffer, is.gcount());
    } while (is);
  } else if (argc == 3) {
    utils::compress_istream is(argv[1]);
    utils::compress_ostream os(argv[2]);
    
    char buffer[4096];
    
    do {
      is.read(buffer, 4096);
      if (is.gcount() > 0)
	os.write(buffer, is.gcount());
    } while (is);
  }
}
