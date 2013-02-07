//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <utils/compress_stream.hpp>
#include <utils/program_options.hpp>

#include <unicode/ucsdet.h>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef std::vector<char, std::allocator<char> > buffer_type;

struct Detector
{
  struct result_type
  {
    std::string lang;
    std::string name;
    int         confidence;
    
    result_type() : lang("unknown"), name(), confidence(0) {}
    result_type(const std::string& __lang,
		const std::string& __name,
		const int& __confidence)
      : lang(__lang), name(__name), confidence(__confidence) {}
  };

  Detector(bool sgml_filter=false) : csd(0)
  {
    UErrorCode status = U_ZERO_ERROR;
    csd = ucsdet_open(&status);
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucsdet_open(): ") + u_errorName(status));
    
    if (sgml_filter)
      ucsdet_enableInputFilter(csd, true);
  }
  ~Detector()
  {
    if (csd)
      ucsdet_close(csd);
  }

  template <typename Buffer>
  result_type operator()(const Buffer& buffer)
  {
    UErrorCode status = U_ZERO_ERROR;
    ucsdet_setText(csd, &(*buffer.begin()), buffer.size(), &status);
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucsdet_setText(): ") + u_errorName(status));
    
    int32_t match_count = 0;
    status = U_ZERO_ERROR;
    const UCharsetMatch **csm = ucsdet_detectAll(csd, &match_count, &status);
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("ucsdet_detectAll(): ") + u_errorName(status));
    
    if (match_count > 0) {
      status = U_ZERO_ERROR;
      const std::string name = ucsdet_getName(csm[0], &status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("ucsdet_getName(): ") + u_errorName(status));

      status = U_ZERO_ERROR;
      const std::string lang = ucsdet_getLanguage(csm[0], &status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("ucsdet_getLanguage(): ") + u_errorName(status));

      status = U_ZERO_ERROR;
      const int confidence = ucsdet_getConfidence(csm[0], &status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("ucsdet_getConfidence(): ") + u_errorName(status));

      return result_type(lang.empty() ? std::string("unknown") : lang, name, confidence);
    } else
      return result_type();
  }
  
private:
  UCharsetDetector* csd;
};

path_type     list_file;
path_type     output_file = "-";
path_set_type input_files;

int confidence_threshold = 80;
int buffer_size = 1024;
bool sgml_filter = false;
bool simple = false;

int getoptions(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    if (getoptions(argc, argv) != 0) 
      return 1;

    if (confidence_threshold <= 0 || confidence_threshold > 100)
      throw std::runtime_error("confidence must be > 0 and <= 100");
    
    if (! list_file.empty()) {
      if (list_file != "-" && ! boost::filesystem::exists(list_file))
	throw std::runtime_error("no list file? " + list_file.string());

      const path_type parent_path = list_file.parent_path();
      
      utils::compress_istream is(list_file);
      
      std::string file;
      while (std::getline(is, file)) {
	if (boost::filesystem::exists(file))
	  input_files.push_back(file);
	else if (boost::filesystem::exists(parent_path / file))
	  input_files.push_back(parent_path / file);
	else
	  throw std::runtime_error("no file? " + file);
      }
    }
    
    if (input_files.empty())
      input_files.push_back("-");

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    buffer_type buffer;
    buffer_type buffer_read(buffer_size, 0);

    Detector detector(sgml_filter);
    
    for (path_set_type::const_iterator iter = input_files.begin(); iter != input_files.end(); ++ iter) {
      utils::compress_istream is(*iter);
      const bool dump_filename = (*iter != "-");
      
      buffer.clear();
      
      Detector::result_type detected;
      
      while (detected.confidence < confidence_threshold) {
	is.read(&(*buffer_read.begin()), buffer_read.size());
	
	if (is.gcount() == 0) break;
	
	buffer.clear();
	buffer.insert(buffer.end(), buffer_read.begin(), buffer_read.begin() + is.gcount());
	
	const Detector::result_type result = detector(buffer);
	
	if (result.confidence > detected.confidence)
	  detected = result;
      }
      
      if (simple)
	os << detected.name << '\n';
      else {
	if (! dump_filename)
	  os << "lang: " << detected.lang
	     << " name: " << detected.name
	     << " confidence: " << detected.confidence
	     << '\n';
	else
	  os << iter->string()
	     << " lang: " << detected.lang
	     << " name: " << detected.name
	     << " confidence: " << detected.confidence
	     << '\n';
      }
    }
  } 
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

int getoptions(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("list",        po::value<path_type>(&list_file),                               "file list")
    ("output",      po::value<path_type>(&output_file)->default_value(output_file), "output file")
    ("confidence",  po::value<int>(&confidence_threshold)->default_value(confidence_threshold), "confidence threshold (> 0 and <= 100)")
    ("buffer",      po::value<int>(&buffer_size)->default_value(buffer_size),  "buffer size")
    ("sgml",        po::bool_switch(&sgml_filter), "ignore sgml-like tags")
    ("simple",      po::bool_switch(&simple),      "simple output")
    
    ("help", "help message");
  
  po::options_description hidden;
  hidden.add_options()
    ("input-file", po::value<path_set_type>(&input_files), "input file");
  
  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);
  
  po::positional_options_description pos;
  pos.add("input-file", -1); // all the files
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options] file(s)" << '\n' << desc << '\n';
    return 1;
  }
  
  return 0;
}
