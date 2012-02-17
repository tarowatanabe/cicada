//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// alignment filter

#include <stdexcept>
#include <iostream>

#include <algorithm>
#include <memory>

#include <cicada/alignment.hpp>
#include <cicada/dependency.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/bithack.hpp"
#include "utils/mathop.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Alignment  alignment_type;
typedef cicada::Dependency permutation_type;


path_set_type input_files;
path_set_type permutation_source_files;
path_set_type permutation_target_files;
path_type output_file = "-";

bool inverse_mode = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
        
    if (input_files.empty())
      input_files.push_back("-");
    
    if (! permutation_source_files.empty())
      if (permutation_source_files.size() != input_files.size())
	throw std::runtime_error("# of permutation files does not match");
    
    if (! permutation_target_files.empty())
      if (permutation_target_files.size() != input_files.size())
	throw std::runtime_error("# of permutation files does not match");
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    if (! permutation_source_files.empty() || ! permutation_target_files.empty()) {
      permutation_type permutation_source;
      permutation_type permutation_target;
      alignment_type   alignment;
      
      const bool has_permutation_source = ! permutation_source_files.empty();
      const bool has_permutation_target = ! permutation_target_files.empty();
      
      for (size_t i = 0; i != input_files.size(); ++ i) {
	std::auto_ptr<std::istream> ps_source(has_permutation_source ? new utils::compress_istream(permutation_source_files[i], 1024 * 1024) : 0);
	std::auto_ptr<std::istream> ps_target(has_permutation_target ? new utils::compress_istream(permutation_target_files[i], 1024 * 1024) : 0);
	
	utils::compress_istream is(input_files[i], 1024 * 1024);
	
	for (;;) {
	  is >> alignment;
	  
	  if (has_permutation_source)
	    *ps_source >> permutation_source;
	  if (has_permutation_target)
	    *ps_target >> permutation_target;
	  
	  if (! is || (ps_source.get() && ! *ps_source) || (ps_target.get() && ! *ps_target)) break;
	  
	  alignment_type::iterator aiter_end = alignment.end();
	  for (alignment_type::iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
	    if (has_permutation_source) {
	      if (aiter->source >= static_cast<int>(permutation_source.size()))
		throw std::runtime_error("invalid source permutation");
	      
	      aiter->source = permutation_source[aiter->source];
	    }
	    
	    if (has_permutation_target) {
	      if (aiter->target >= static_cast<int>(permutation_target.size()))
		throw std::runtime_error("invalid target permutation");
		      
	      aiter->target = permutation_target[aiter->target];
	    }
	  }
	  
	  if (inverse_mode) {
	    alignment.inverse();
	    std::sort(alignment.begin(), alignment.end());
	  }
	  
	  os << alignment << '\n';
	}
	
	if (is || (ps_source.get() && *ps_source) || (ps_target.get() && *ps_target))
	  throw std::runtime_error("# of samples do not match");
      }
    } else {
      for (path_set_type::const_iterator fiter = input_files.begin(); fiter != input_files.end(); ++ fiter) {
	utils::compress_istream is(*fiter, 1024 * 1024);
	
	alignment_type alignment;
	
	while (is >> alignment) {
	  
	  if (inverse_mode) {
	    alignment.inverse();
	    std::sort(alignment.begin(), alignment.end());
	  }
	  
	  os << alignment << '\n';
	}
      }
    }
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
    ("input",              po::value<path_set_type>(&input_files)->multitoken(),              "input file(s)")
    ("permutation-source", po::value<path_set_type>(&permutation_source_files)->multitoken(), "source side permutation file(s)")
    ("permutation-target", po::value<path_set_type>(&permutation_target_files)->multitoken(), "target side permutation file(s)")
    ("output",             po::value<path_type>(&output_file)->default_value(output_file),    "output file")
    
    ("inverse",   po::bool_switch(&inverse_mode), "inverse alignment")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}


