//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <string>

#include "utils/subprocess.hpp"
#include "utils/async_device.hpp"

#include <boost/thread.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>

struct Task
{
  int fd;
  
  Task(const int _fd) : fd(_fd) {}
  
  void operator()()
  {
    boost::iostreams::filtering_ostream os;
    os.push(utils::async_sink(fd, true));

    std::string line;
    while (std::getline(std::cin, line))
      os << line << std::endl;
  }
};

int main(int argc, char** argv)
{
  utils::subprocess run(boost::filesystem::path("cat"));

  boost::thread thread(Task(run.desc_write()));
  
  boost::iostreams::filtering_istream is;
  is.push(utils::async_source(run.desc_read(), true));
  
  std::string line;
  while (std::getline(is, line))
    std::cout << line << std::endl;

  thread.join();
}
