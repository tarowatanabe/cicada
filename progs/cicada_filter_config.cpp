//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//


// a simple cicada.config filter
// we will replace the occurencce of
//
// ${weights} ${weight_file} ${directory} ${file}
//
// where
//   weight_file refers to weight file
//   weights refers to weights operations (weights=weight_file or weights-one=true)
//   directory will be the directory where output will be dumped
//   file will be the file where output will be dumped
//

#include <string>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <boost/xpressive/xpressive.hpp>

#include "utils/compress_stream.hpp"

typedef boost::filesystem::path path_type;

path_type input_file = "-";
path_type output_file = "-";

std::string weights;
std::string weight_file;
std::string kbest;
std::string directory;
std::string file;

bool remove_operation = false;
bool remove_feature_function = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    std::string config;

    utils::compress_istream is(input_file, 1024 * 1024);
    
    char buffer[4096];
    do {
      is.read(buffer, 4096);
      if (is.gcount() > 0)
	config.append(buffer, is.gcount());
    } while (is);
    
    // replace!
    if (! weights.empty())
      boost::algorithm::replace_all(config, "${weights}", weights);
    if (! weight_file.empty())
      boost::algorithm::replace_all(config, "${weight_file}", weight_file);
    if (! kbest.empty())
      boost::algorithm::replace_all(config, "${kbest}", kbest);
    if (! directory.empty())
      boost::algorithm::replace_all(config, "${directory}", directory);
    if (! file.empty())
      boost::algorithm::replace_all(config, "${file}", file);

    if (remove_operation) {
      namespace xpressive = boost::xpressive;
      
      const xpressive::sregex re = ((xpressive::s1= (xpressive::bos | xpressive::_ln))
				    >> (xpressive::s2= -*(xpressive::blank) >> "operation" >> *(xpressive::_s) >> '=' >> *(~xpressive::_ln)));
      
      
      config = xpressive::regex_replace(config, re, "$1# $2");
    }

    if (remove_feature_function) {
      namespace xpressive = boost::xpressive;
      
      const xpressive::sregex re = ((xpressive::s1= (xpressive::bos | xpressive::_ln))
				    >> (xpressive::s2= -*(xpressive::blank) >> "feature-function" >> *(xpressive::_s) >> '=' >> *(~xpressive::_ln)));
					
      
      config = xpressive::regex_replace(config, re, "$1# $2");
    }
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    os << config;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")
    
    ("weights",     po::value<std::string>(&weights),     "substitute ${weights}")
    ("weight-file", po::value<std::string>(&weight_file), "substitute ${weight_file}")
    ("kbest",       po::value<std::string>(&kbest),       "substitute ${kbest}")
    ("directory",   po::value<std::string>(&directory),   "substitute ${directory}")
    ("file",        po::value<std::string>(&file),        "substitute ${file}")
    
    ("remove-operation",        po::bool_switch(&remove_operation),        "remove operation(s)")
    ("remove-feature-function", po::bool_switch(&remove_feature_function), "remove feature function(s)")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

