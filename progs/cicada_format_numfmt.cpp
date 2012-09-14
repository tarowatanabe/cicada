//
// ICU RBNF
//

#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <memory>

#include <unicode/numfmt.h>
#include <unicode/unistr.h>
#include <unicode/locid.h>
#include <unicode/bytestream.h>
#include <unicode/resbund.h>
#include <unicode/ures.h>
#include <unicode/udata.h>

#define U_ICUDATA_RBNF U_ICUDATA_NAME U_TREE_SEPARATOR_STRING "rbnf"

#include "utils/compress_stream.hpp"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

typedef boost::filesystem::path path_type;

struct ostream_sink : public ByteSink
{
  
  ostream_sink(std::ostream& _os) : os(_os) {}
  
  virtual void Append(const char* data, int32_t n) 
  {
    os.write((char*) data, n);
  }
  
  void write(char c)
  {
    os.write((char*) & c, sizeof(c));
  }
  
  std::ostream& os;
};

path_type input_file = "-";
path_type output_file = "-";

std::string locale;
std::string rule;

bool resource = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (locale.empty())
      throw std::runtime_error("no locale?");
    
    icu::Locale loc(locale.c_str());
    
    if (loc.isBogus())
      throw std::runtime_error("invalid ocale: " + locale);

    typedef std::pair<int, icu::NumberFormat*> formatter_type;
    typedef std::vector<formatter_type, std::allocator<formatter_type> > formatter_set_type;
    
    formatter_set_type formatters;
    
    for (int style = 0; style != UNUM_FORMAT_STYLE_COUNT; ++ style) {
      UErrorCode status = U_ZERO_ERROR;  
      std::auto_ptr<icu::NumberFormat> format(icu::NumberFormat::createInstance(loc, UNumberFormatStyle(style), status));
      
      if (U_FAILURE(status)) {
	std::cerr << "no style: " << style << std::endl;
	continue;
      }
      
      formatters.push_back(std::make_pair(style, format.release()));
    }
    
    const bool flush_output = (output_file == "-"
                               || (boost::filesystem::exists(output_file)
                                   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_istream is(input_file);
    utils::compress_ostream os(output_file, 4096 * (! flush_output));
    ostream_sink sink(os);
    
    int64_t integer;
    UnicodeString formatted;
    while (is >> integer) {
      for (size_t i = 0; i != formatters.size(); ++ i) {
	{
	  UErrorCode status(U_ZERO_ERROR);
	  UChar uCurrency[4];
	  u_charsToUChars("USD", uCurrency, 4);
	  formatters[i].second->setCurrency(uCurrency, status);
	  
	  if (U_FAILURE(status))
	    std::cerr << "not a currency!" << std::endl;
	}

	os << formatters[i].first << ": ";
	
	{
	  FieldPosition pos;
	  formatted.remove();
	  formatters[i].second->format(integer, formatted, pos);
	  formatted.toUTF8(sink);
	}
	
	sink.write('\n');
	icu::ParsePosition pos(0);
	std::auto_ptr<icu::CurrencyAmount> curr(formatters[i].second->parseCurrency(formatted, pos));
	
	if (curr.get() && pos.getIndex() == formatted.length()) {
	  os << "re-format currency: ";

	  UErrorCode status = U_ZERO_ERROR;
	  formatters[i].second->setCurrency(curr->getISOCurrency(), status);
	  
	  formatted.remove();
	  status = U_ZERO_ERROR;
	  formatters[i].second->format(curr->getNumber(), formatted, status);
	  
	  if (U_FAILURE(status))
	    os << "failed";
	  else
	    formatted.toUTF8(sink);
	  
	  sink.write('\n');
	} else {
	  UErrorCode status = U_ZERO_ERROR;
	  icu::Formattable parsed;
	  formatters[i].second->parse(formatted, parsed, status);
	  
	  if (U_FAILURE(status))
	    std::cerr << "parsing failed" << std::endl;
	  else {
	    os << "re-format: ";
	    
	    UErrorCode status = U_ZERO_ERROR;
	    formatted.remove();
	    formatters[i].second->format(parsed, formatted, status);
	    
	    if (U_FAILURE(status))
	      os << "failed";
	    else
	      formatted.toUTF8(sink);
	    
	    sink.write('\n');
	  }
	}
      }
    }
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
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output")
    
    ("locale",   po::value<std::string>(&locale), "locale")
    ("rule",     po::value<std::string>(&rule),   "rule set")
    ("resource", po::bool_switch(&resource),      "dump resource")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}
