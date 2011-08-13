//
// ICU RBNF
//

#include <iostream>
#include <string>

#include <unicode/rbnf.h>
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
    
    Locale loc(locale.c_str());
    
    if (loc.isBogus())
      throw std::runtime_error("invalid ocale: " + locale);

    if (resource) {
      const char* rules_tag = "RBNFRules";
      const char* fmt_tag   = "SpelloutRules";

      utils::compress_ostream os(output_file);
      ostream_sink sink(os);
      
      UErrorCode status = U_ZERO_ERROR;
      icu::ResourceBundle bundle(U_ICUDATA_RBNF, loc, status);
      
      if (U_SUCCESS(status)) {
	icu::ResourceBundle bundle_rbnf = bundle.getWithFallback(rules_tag, status);
        if (U_FAILURE(status))
	  throw std::runtime_error("no rbnf?");
	
	icu::ResourceBundle bundle_rules = bundle_rbnf.getWithFallback(fmt_tag, status);
        if (U_FAILURE(status))
	  throw std::runtime_error("no rules?");
	
	while (bundle_rules.hasNext() && U_SUCCESS(status)) {
	  bundle_rules.getNextString(status).toUTF8(sink);
	  sink.write('\n');
	}
	
	if (U_FAILURE(status))
	  throw std::runtime_error("no string?");
      } else 
	throw std::runtime_error("no RBNF rules?");
      
      return 0;
    }
    
    UErrorCode status = U_ZERO_ERROR;
    std::auto_ptr<RuleBasedNumberFormat> formatter(new RuleBasedNumberFormat(URBNF_SPELLOUT, loc, status));
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
    
    if (! rule.empty()) {
      UErrorCode status = U_ZERO_ERROR;
      formatter->setDefaultRuleSet(UnicodeString(rule.c_str()), status);
      
      if (U_FAILURE(status)) {
	std::cerr << "# of rules: " << formatter->getNumberOfRuleSetNames() << std::endl;
	for (int i = 0; i < formatter->getNumberOfRuleSetNames(); ++ i) {
	  std::string name;
	  formatter->getRuleSetName(i).toUTF8String(name);
	  
	  std::cerr << "rule set: " << name << std::endl;
	}
	
	throw std::runtime_error("unsupported rule set: " + rule);
      }
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
      FieldPosition pos;
      formatted.remove();
      formatter->format(integer, formatted, pos);
      formatted.toUTF8(sink);
      sink.write('\n');
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
