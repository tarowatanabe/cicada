//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// clear temporary files

#include <cstdlib>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/find.hpp>

#include <utils/tempfile.hpp>
#include <utils/filesystem.hpp>

typedef boost::filesystem::path path_type;
typedef boost::program_options::variables_map variable_set_type;

void options(int argc, char** argv, variable_set_type& variables);

int main(int argc, char** argv)
{
  try {
    variable_set_type variables;
    options(argc, argv, variables);
    
    if (variables.count("temporary") && ! variables["temporary"].as<path_type>().empty())
      ::setenv("TMPDIR_SPEC", variables["temporary"].as<path_type>().string().data(), 1);
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();

    static const char* cicada_suffix[] = {
      "extract",
      "score",
      "attr",
      "edge",
      "cluster",
      "lexicon",
      "rule",
      "source",
      "target",
      "feature",
      "feature-data",
      "feature-vocab",
      "attribute",
      "attribute-data",
      "attribute-vocab",
      "vocab",
    };
    
    static const char* succinct_suffix[] = {
      "size",
      "key-data",
    };
    
    const int cicada_suffix_size   = sizeof(cicada_suffix) / sizeof(const char*);
    const int succinct_suffix_size = sizeof(succinct_suffix) / sizeof(const char*);
    
    boost::filesystem::directory_iterator iter_end;
    for (boost::filesystem::directory_iterator iter(tmp_dir); iter != iter_end; ++ iter) {
      const path_type path = *iter;

#if BOOST_FILESYSTEM_VERSION == 2
      for (int i = 0; i < cicada_suffix_size; ++ i)
	if (path_type(path.filename()).file_string().find(std::string("cicada.") + cicada_suffix[i]) != std::string::npos)
	  utils::filesystem::remove_all(path);
      for (int i = 0; i < succinct_suffix_size; ++ i)
	if (path_type(path.filename()).file_string().find(std::string("succinct-db.") + succinct_suffix[i]) != std::string::npos)
	  utils::filesystem::remove_all(path);
#else
      for (int i = 0; i < cicada_suffix_size; ++ i)
	if (path_type(path.filename()).string().find(std::string("cicada.") + cicada_suffix[i]) != std::string::npos)
	  utils::filesystem::remove_all(path);
      for (int i = 0; i < succinct_suffix_size; ++ i)
	if (path_type(path.filename()).string().find(std::string("succinct-db.") + succinct_suffix[i]) != std::string::npos)
	  utils::filesystem::remove_all(path);
#endif
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}

void options(int argc, char** argv, variable_set_type& variables)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("temporary", po::value<path_type>(), "temporary directory")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
