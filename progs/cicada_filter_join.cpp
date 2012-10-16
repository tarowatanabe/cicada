//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// join lines with "|||"

#include <boost/spirit/include/karma.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <vector>
#include <string>
#include <iterator>

#include "utils/compress_stream.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type files;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    typedef std::vector<std::istream*, std::allocator<std::istream*> > istream_set_type;
    typedef std::vector<std::string, std::allocator<std::string> > line_set_type;

    options(argc, argv);

    if (files.empty())
      files.push_back("-");
    
    istream_set_type istreams(files.size());

    for (size_t i = 0; i != files.size(); ++ i) {
      if (files[i] != "-" && ! boost::filesystem::exists(files[i]))
	throw std::runtime_error("no file? " + files[i].string());
      
      istreams[i] = new utils::compress_istream(files[i], 1024 * 1024);
    }
    
    std::ostream_iterator<char> oiter(std::cout);
    
    line_set_type lines(istreams.size());
    for (;;) {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      bool terminated = false;
      
      for (size_t i = 0; i != istreams.size(); ++ i)
	terminated |= (! std::getline(*istreams[i], lines[i]));
      
      if (terminated) break;
      
      if (! karma::generate(oiter, (standard::string % " ||| ") << '\n', lines))
	throw std::runtime_error("generation failed");
    }
    
    bool incompatible = false;
    
    istream_set_type::iterator iiter_end = istreams.end();
    for (istream_set_type::iterator iiter = istreams.begin(); iiter != iiter_end; ++ iiter) {
      if (**iiter)
	incompatible = true;
      
      delete *iiter;
    }
    
    if (incompatible)
      throw std::runtime_error("# of lines do not match");
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("help", "help message");
  
  po::options_description hidden;
  hidden.add_options()
    ("file", po::value<path_set_type>(&files), "files");
  
  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);
  
  po::positional_options_description pos;
  pos.add("file", -1); // all the files
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " file(s)" << '\n' << desc << '\n';
    exit(0);
  }
}
