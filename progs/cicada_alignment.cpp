// cicada alignment tool:
//
// we can produce intersection/union/grow-{,diag}-{final,final-and}/source/target/itg/max-match from GIZA++ alingment
//        invert alignment
//

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>

#include <cicada/alignment.hpp>

#include "utils/program_options.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;

typedef cicada::Alignment alignment_type;

path_type source_target_file;
path_type target_source_file;
path_type span_source_file;
path_type span_target_file;
path_type input_file;
path_type output_file = "-";

bool source_target_mode = false;
bool target_source_mode = false;

bool itg_mode = false;
bool max_match_mode = false;

bool intersection_mode = false;
bool union_mode = false;
bool grow_mode = false;
bool final_mode = false;
bool diag_mode = false;
bool final_and_mode = false;
bool invert_mode = false;

int threads = 1;
int debug = 0;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    
    
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
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source-target", po::value<path_type>(&source_target_file), "P(target | source) viterbi output")
    ("target-source", po::value<path_type>(&target_source_file), "P(source | target) viterbi output")
    ("span-source",   po::value<path_type>(&span_source_file),   "source span data")
    ("span-target",   po::value<path_type>(&span_target_file),   "target span data")
    ("input",         po::value<path_type>(&input_file),                      "input alignement")
    ("output",        po::value<path_type>(&output_file)->default_value("-"), "output alignment")
    
    ("f2e", po::bool_switch(&source_target_mode), "source target")
    ("e2f", po::bool_switch(&target_source_mode), "target source")
    
    ("itg",          po::bool_switch(&itg_mode),          "itg")
    ("max-match",    po::bool_switch(&max_match_mode),    "max-match")
    ("intersection", po::bool_switch(&intersection_mode), "intersection")
    ("union",        po::bool_switch(&union_mode),        "union")
    ("grow",         po::bool_switch(&grow_mode),         "grow")
    ("diag",         po::bool_switch(&diag_mode),         "diag")
    ("final",        po::bool_switch(&final_mode),        "final")
    ("final-and",    po::bool_switch(&final_and_mode),    "final-and")
    ("invert",       po::bool_switch(&invert_mode),       "invert alignment")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
