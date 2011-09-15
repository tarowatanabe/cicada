//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// dependency filter to project source dependency into target dependency using the word alignment
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <memory>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <cicada/alignment.hpp>
#include <cicada/dependency.hpp>

#include "utils/bithack.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lexical_cast.hpp"

typedef cicada::Alignment  alignment_type;
typedef cicada::Dependency dependency_type;

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef std::vector<std::string, std::allocator<std::string> > sentence_type;

template <typename Iterator>
struct sentence_parser : boost::spirit::qi::grammar<Iterator, sentence_type(), boost::spirit::standard::blank_type>
{
  sentence_parser() : sentence_parser::base_type(tokens)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"] >> (qi::eoi | qi::eol);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, sentence_type(), blank_type> tokens;
};

path_set_type source_files;
path_set_type target_files;
path_set_type alignment_files;
path_set_type dependency_files;

path_type list_source_file;
path_type list_target_file;
path_type list_alignment_file;
path_type list_dependency_file;

bool projective_mode = false;

// dependency output
path_type output_file = "-";

void read_list(const path_type& path, path_set_type& files);
void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    read_list(list_source_file, source_files);
    read_list(list_target_file, target_files);
    read_list(list_alignment_file, alignment_files);
    read_list(list_dependency_file, dependency_files);
    
    if (source_files.empty())
      source_files.push_back("-");
    if (target_files.empty())
      target_files.push_back("-");
    if (alignment_files.empty())
      alignment_files.push_back("-");
    if (dependency_files.empty())
      dependency_files.push_back("-");
    
    if (source_files.size() != target_files.size())
      throw std::runtime_error("# of files do not match");
    
    if (source_files.size() != alignment_files.size())
      throw std::runtime_error("# of alignment files do not match");
    if (target_files.size() != alignment_files.size())
      throw std::runtime_error("# of alignment files do not match");
    if (source_files.size() != dependency_files.size())
      throw std::runtime_error("# of dependency files do not match");
    if (target_files.size() != dependency_files.size())
      throw std::runtime_error("# of dependency files do not match");
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    typedef boost::spirit::istream_iterator iiter_type;
    
    sentence_parser<iiter_type>    parser;

    for (size_t i = 0; i != source_files.size(); ++ i) {
      utils::compress_istream is_src(source_files[i], 1024 * 1024);
      utils::compress_istream is_trg(target_files[i], 1024 * 1024);
      utils::compress_istream is_align(alignment_files[i], 1024 * 1024);
      utils::compress_istream is_dep(dependency_files[i], 1024 * 1024);
      
      is_src.unsetf(std::ios::skipws);
      is_trg.unsetf(std::ios::skipws);
      
      iiter_type siter(is_src);
      iiter_type titer(is_trg);
      iiter_type siter_end;
      iiter_type titer_end;
      
      sentence_type source;
      sentence_type target;
      alignment_type alignment;
      dependency_type dependency;
      dependency_type projected;
      
      for (size_t line_no = 0; siter != siter_end && titer != titer_end; ++ line_no) {
	source.clear();
	target.clear();
	alignment.clear();
	dependency.clear();
	
	if (! boost::spirit::qi::phrase_parse(siter, siter_end, parser, boost::spirit::standard::blank, source))
	  throw std::runtime_error("source sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	if (! boost::spirit::qi::phrase_parse(titer, titer_end, parser, boost::spirit::standard::blank, target))
	  throw std::runtime_error("target sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	
	if (! (is_align >> alignment))
	  throw std::runtime_error("no alignment?");
	if (! (is_dep >> dependency))
	  throw std::runtime_error("no dependency?");
	
	const int source_size = source.size();
	const int target_size = target.size();

	if (source_size != static_cast<int>(dependency.size()))
	  throw std::runtime_error("source size do not match with dependency size at # " + utils::lexical_cast<std::string>(line_no));
	
	if (source_size == 0 || target_size == 0) {
	  os << '\n';
	  continue;
	};
	
	projected.clear();
	projected.resize(target_size, -1);
	
	// perform projection....
	
	
	os << projected << '\n';
      }
    }
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}


void read_list(const path_type& path, path_set_type& files)
{
  if (path.empty()) return;
  if (path != "-" && ! boost::filesystem::exists(path))
    throw std::runtime_error("no file? " + path.string());
  
  utils::compress_istream is(path);
  std::string file;
  while (std::getline(is, file)) {
    if (file.empty()) continue;
    if (! boost::filesystem::exists(file))
      throw std::runtime_error("no file? " + file);
    files.push_back(file);
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("source",     po::value<path_set_type>(&source_files)->multitoken(),     "source file(s)")
    ("target",     po::value<path_set_type>(&target_files)->multitoken(),     "target file(s)")
    ("alignment",  po::value<path_set_type>(&alignment_files)->multitoken(),  "alignment file(s)")
    ("dependency", po::value<path_set_type>(&dependency_files)->multitoken(), "dependency file(s)")
    
    ("list-source",     po::value<path_type>(&list_source_file),     "source list file")
    ("list-target",     po::value<path_type>(&list_target_file),     "target list file")
    ("list-alignment",  po::value<path_type>(&list_alignment_file),  "alignment list file")
    ("list-dependency", po::value<path_type>(&list_dependency_file), "dependency list file")
    
    ("projective", po::bool_switch(&projective_mode), "project into projective dependency")
    
    ("output",    po::value<path_type>(&output_file)->default_value(output_file),       "output file")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
