//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// tee, but support compressed output...
//

#include <boost/filesystem.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

typedef boost::filesystem::path path_type;

path_type tee_file;

int options(int argc, char** argv);

int main(int argc, char** argv)
{
  if (options(argc, argv) != 0)
    return 1;
  
  try {
    if (tee_file.empty())
      throw std::runtime_error("no path for tee output?");
    
    utils::compress_ostream os(tee_file, 1024 * 1024);
    
    char buffer[4096];

    do {
      std::cin.read(buffer, 4096);
      if (std::cin.gcount() > 0) {
	os.write(buffer, std::cin.gcount());
	std::cout.write(buffer, std::cin.gcount());
      }
    } while (std::cin);
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

int options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description hidden;
  hidden.add_options()
    ("tee", po::value<path_type>(&tee_file), "tee file");
  
  po::options_description cmdline_options;
  cmdline_options.add(hidden);

  po::positional_options_description pos;
  pos.add("tee", 1); // only one...

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " output-file" << '\n';
    return 1;
  }
  
  return 0;
}
