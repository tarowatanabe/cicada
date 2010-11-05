
#include "cicada_extract_score_impl.hpp"

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <utility>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>

#include <utils/resource.hpp>

typedef boost::filesystem::path                                    path_type;
typedef std::vector<path_type, std::allocator<path_type> >         path_set_type;
typedef std::vector<path_set_type, std::allocator<path_set_type> > path_map_type;

path_set_type counts_files;

path_type lexicon_source_target_file;
path_type lexicon_target_source_file;

path_type output_file;

bool score_phrase = false;
bool score_scfg   = false;
bool score_tree   = false;

double max_malloc = 8; // 8 GB
int    threads = 1;
double discount_dp = 0.0;

int debug = 0;

int main(int argc, char** argv)
{
  
  
}
