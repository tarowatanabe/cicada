//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

#include "codec/snappy.hpp"

int main(int argc, char** argv)
{
  boost::iostreams::filtering_ostream os;
  
  os.push(codec::snappy_compressor());
  os.push(boost::iostreams::file_descriptor_sink(::dup(STDOUT_FILENO), boost::iostreams::close_handle));
  
  char buffer[4096];
  
  do {
    std::cin.read(buffer, 4096);
    if (std::cin.gcount() > 0)
      os.write(buffer, std::cin.gcount());
  } while (std::cin);
}
