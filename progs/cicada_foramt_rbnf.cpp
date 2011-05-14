//
// ICU RBNF
//

#include <iostream>
#include <string>

#include <unicode/rbnf.h>
#include <unicode/unistr.h>
#include <unicode/locid.h>

#include <boost/program_options.hpp>

typedef boost::filesystem::path path_type;

path_type input_file = "-";
path_type output_file = "-";

std::string locale;
std::string rule;

int main(int argc, char** argv)
{
  options(argc, argv);

  try {
    if (local.empty())
      throw std::runtime_error("no locale?");
    
    Locale loc(locale.c_str());
    
    if (loc.isBogus())
      throw std::runtime_error("invlaid ocale: " + locale);
    
    UErrorCode status = U_ZERO_ERROR;
    std::auto_ptr<RuleBasedNumberFormat> parser_rule_spell(new RuleBasedNumberFormat(URBNF_SPELLOUT, loc, status));
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RuleBasedNumberFormat::spell_out: ") + u_errorName(status));
    
    
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
    
    ("locale", po::value<std::string>(&locale), "locale")
    ("rule",   po::value<std::string>(&rule),   "rule set")
    
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
