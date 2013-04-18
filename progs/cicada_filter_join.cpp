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
#include "utils/getline.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

path_set_type input_files;
path_type     output_file = "-";

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    typedef std::vector<std::istream*, std::allocator<std::istream*> > istream_set_type;
    typedef std::vector<std::string, std::allocator<std::string> > line_set_type;

    options(argc, argv);

    if (input_files.empty())
      input_files.push_back("-");
    
    istream_set_type istreams(input_files.size());

    for (size_t i = 0; i != input_files.size(); ++ i) {
      if (input_files[i] != "-" && ! boost::filesystem::exists(input_files[i]))
	throw std::runtime_error("no file? " + input_files[i].string());
      
      istreams[i] = new utils::compress_istream(input_files[i], 1024 * 1024);
    }

    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    line_set_type lines(istreams.size());
    for (;;) {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      bool terminated = false;
      
      for (size_t i = 0; i != istreams.size(); ++ i)
	terminated |= (! utils::getline(*istreams[i], lines[i]));
      
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
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    ("help", "help message");
  
  po::options_description hidden;
  hidden.add_options()
    ("input", po::value<path_set_type>(&input_files), "input file");

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::positional_options_description pos;
  pos.add("input", -1); // all the files

  po::command_line_parser parser(argc, argv);
  parser.style(po::command_line_style::unix_style & (~po::command_line_style::allow_guessing));
  parser.options(cmdline_options);
  parser.positional(pos);
  
  po::variables_map vm;
  po::store(parser.run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " file(s)" << '\n' << desc << '\n';
    exit(0);
  }
}
