//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>

#include "utils/snappy_device.hpp"
#include "utils/snappy_file.hpp"
#include "utils/resource.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << argv[0] << " [output file]" << std::endl;
    return 1;
  }
  
  boost::iostreams::filtering_ostream os;
  os.push(utils::snappy_sink<>(argv[1]));
  
  utils::resource compress_start;

  char buffer[4096];
  do {
    std::cin.read(buffer, 4096);
    if (std::cin.gcount() > 0)
      os.write(buffer, std::cin.gcount());
  } while(std::cin);

  os.pop();

  utils::resource compress_end;

  utils::snappy_file<> file(argv[1]);

  std::cerr << "size-bytes: " << file.size_bytes()
	    << " size-compressed: " << file.size_compressed()
	    << " size-cache: " << file.size_cache() 
	    << std::endl;

  utils::resource decompress_start;
  
  const size_t last = file.size();
  for (size_t first = 0; first != last; /**/) {
    const size_t read_size = std::min(last - first, size_t(4096));
    
    file.read(buffer, read_size, first);
    
    std::cout.write(buffer, read_size);
    
    first += read_size;
  }

  utils::resource decompress_end;
  
  std::cerr << "compress:"
	    << " user time: " << (compress_end.user_time() - compress_start.user_time())
	    << " cpu time: " << (compress_end.cpu_time() - compress_start.cpu_time())
	    << " decompress:"
	    << " user time: " << (decompress_end.user_time() - decompress_start.user_time())
	    << " cpu time: " << (decompress_end.cpu_time() - decompress_start.cpu_time())
	    << std::endl;

}
