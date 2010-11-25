//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <string>

#include "icu_filter.hpp"

#include <boost/version.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/device/file.hpp>

int main(int argc, char** argv)
{
  try {
    boost::iostreams::filtering_istream is;
    is.push(utils::icu_filter("utf-8", "euc-jp", utils::icu_filter_param::stop));
#if BOOST_VERSION >= 104400
    is.push(boost::iostreams::file_descriptor_source(::dup(STDIN_FILENO), boost::iostreams::close_handle));
#else
    is.push(boost::iostreams::file_descriptor_source(::dup(STDIN_FILENO), true));
#endif
    
    char buffer[4096];
    
    size_t total = 0;
    do {
      is.read(buffer, 4096);
      
      int read_size = is.gcount();
      std::cerr << "read size: " << read_size << std::endl;
      
      total += (is.gcount() > 0 ? is.gcount() : 0);
      
      if (is.gcount() > 0)
        std::cout.write(buffer, is.gcount());
    } while (is);
    
    std::cerr << "total size: " << total << std::endl;
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
