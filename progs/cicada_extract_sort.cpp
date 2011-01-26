//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>

#include "cicada_extract_score_impl.hpp"

#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <utils/resource.hpp>



typedef PhrasePair       rule_pair_type;
typedef PhrasePairParser rule_pair_parser_type;

#ifdef HAVE_TR1_UNORDERED_SET
typedef std::tr1::unordered_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
				std::allocator<rule_pair_type> > rule_pair_set_type;
#else
typedef sgi::hash_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
		      std::allocator<rule_pair_type> > rule_pair_set_type;
#endif

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

template <typename Tp>
struct less_ptr
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return *x < *y;
  }
};

path_set_type input_files;
path_type output_file = "-";

int debug;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (input_files.empty())
      input_files.push_back("-");
    
    std::string           line;
    rule_pair_type        rule_pair;
    rule_pair_set_type    rule_pairs;
    rule_pair_parser_type parser;
    
    path_set_type::const_iterator fiter_end = input_files.end();
    for (path_set_type::const_iterator fiter = input_files.begin(); fiter != fiter_end; ++ fiter) {
      utils::compress_istream is (*fiter, 1024 * 1024);
      
      while (std::getline(is, line)) {
	if (! parser(line, rule_pair)) continue;
	
	std::pair<rule_pair_set_type::iterator, bool> result = rule_pairs.insert(rule_pair);
	if (! result.second)
	  const_cast<rule_pair_type&>(*result.first).increment(rule_pair.counts.begin(), rule_pair.counts.end());
      }
    }
    
    typedef std::vector<const rule_pair_type*, std::allocator<const rule_pair_type*> > sorted_type;
   
    sorted_type sorted(rule_pairs.size());
    {
      sorted_type::iterator siter = sorted.begin();
      rule_pair_set_type::const_iterator citer_end = rule_pairs.end();
      for (rule_pair_set_type::const_iterator citer = rule_pairs.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    
    std::sort(sorted.begin(), sorted.end(), less_ptr<rule_pair_type>());
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    os.precision(20);
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << *(*siter) << '\n';
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("options");
  
  opts_config.add_options()
    ("input",  po::value<path_set_type>(&input_files)->multitoken(), "input file(s)")
    ("output", po::value<path_type>(&output_file),                   "output file")
    ("help", "help message");
  
  po::options_description desc;
  desc.add(opts_config);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}


