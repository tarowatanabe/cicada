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
#include "utils/lexical_cast.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

class Reader
{
public:
  Reader() : finished(false) {}
  virtual ~Reader() {}
  
  virtual Reader& read(std::string& line) = 0;
  
  virtual operator bool() const { return ! finished; }
  
  bool finished;
};

class ReaderFile : public Reader
{
public:
  ReaderFile(const path_type& path)
    : is(new utils::compress_istream(path, 1024 * 1024)) {}
  
  ~ReaderFile()
  {
    if (is)
      delete is;
  }
  
  Reader& read(std::string& line)
  {
    if (! utils::getline(*is, line)) {
      line.clear();
      finished = true;
    }
    
    return *this;
  }
  
  std::istream* is;
};

class ReaderDirectory : public Reader
{
public:
  ReaderDirectory(const path_type& __path)
    : path(__path), id(0) {}
  
  Reader& read(std::string& line)
  {
    const path_type file = path / (utils::lexical_cast<std::string>(id) + ".gz");
    
    if (! boost::filesystem::exists(file)) {
      line.clear();
      finished = true;
      return *this;
    }
    
    utils::compress_istream is(file, 1024 * 1024);
    if (! utils::getline(is, line)) {
      line.clear();
      finished = true;
      return *this;
    }
    
    ++ id;

    return *this;
  }
  
  path_type path;
  uint64_t id;
};

path_set_type input_files;
path_type     output_file = "-";

bool id_mode = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    typedef Reader reader_type;
    typedef std::vector<reader_type*, std::allocator<reader_type*> > reader_set_type;
    typedef std::vector<std::string, std::allocator<std::string> > line_set_type;
    
    options(argc, argv);

    if (input_files.empty())
      input_files.push_back("-");
    
    reader_set_type readers(input_files.size());
    line_set_type lines(input_files.size());
    
    for (size_t i = 0; i != input_files.size(); ++ i) {
      if (input_files[i] != "-" && ! boost::filesystem::exists(input_files[i]))
	throw std::runtime_error("no file? " + input_files[i].string());

      if (boost::filesystem::is_directory(input_files[i]))
	readers[i] = new ReaderDirectory(input_files[i]);
      else
	readers[i] = new ReaderFile(input_files[i]);
    }
    
    uint64_t id = 0;
    karma::uint_generator<uint64_t> id_gen;

    utils::compress_ostream os(output_file, 1024 * 1024);
    std::ostream_iterator<char> oiter(os);
    
    for (;;) {
      bool terminated = false;
      
      for (size_t i = 0; i != readers.size(); ++ i)
	terminated |= ! readers[i]->read(lines[i]);
      
      if (terminated) break;
      
      if (id_mode) {
	if (! karma::generate(oiter, id_gen << " ||| " << (standard::string % " ||| ") << '\n', id, lines))
	  throw std::runtime_error("generation failed");
      } else {
	if (! karma::generate(oiter, (standard::string % " ||| ") << '\n', lines))
	  throw std::runtime_error("generation failed");
      }

      ++ id;
    }
    
    int incompatible = 0;
    
    reader_set_type::iterator iiter_end = readers.end();
    for (reader_set_type::iterator iiter = readers.begin(); iiter != iiter_end; ++ iiter) {
      if (**iiter)
	++ incompatible;
      
      delete *iiter;
    }
    
    if (incompatible)
      throw std::runtime_error("# of lines do not match: " + utils::lexical_cast<std::string>(incompatible));
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
    ("id", po::bool_switch(&id_mode), "output id")
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
