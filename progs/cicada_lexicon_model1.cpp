//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/alloc_vector.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <google/dense_hash_map>


typedef cicada::Symbol   word_type;
typedef cicada::Sentence sentence_type;
typedef cicada::Vocab    vocab_type;
typedef boost::filesystem::path path_type;

typedef double count_type;
typedef double prob_type;

struct count_map_type
{
  typedef google::dense_hash_map<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type> > counts_type;

  typedef counts_type::value_type      value_type;
  typedef counts_type::size_type       size_type;
  typedef counts_type::difference_type difference_type;
      
  typedef counts_type::mapped_type     mapped_type;
  typedef counts_type::key_type        key_type;
  
  typedef counts_type::const_iterator const_iterator;
  typedef counts_type::iterator       iterator;
  
  typedef counts_type::const_reference const_reference;
  typedef counts_type::reference       reference;
  
  count_map_type() { counts.set_empty_key(word_type()); }

  inline const_iterator begin() const { return counts.begin(); }
  inline       iterator begin()       { return counts.begin(); }
  inline const_iterator end() const { return counts.end(); }
  inline       iterator end()       { return counts.end(); }

  mapped_type& operator[](const key_type& key) { return counts[key]; }
      
  size_type size() const { return counts.size(); }
  bool empty() const { return counts.empty(); }

  void swap(count_map_type& x) { counts.swap(x.counts); }
  void clear() { counts.clear(); }
  
  counts_type counts;
};

typedef utils::alloc_vector<count_map_type, std::allocator<count_map_type> > count_dict_type;

path_type source_file = "-";
path_type target_file = "-";
path_type output_source_target_file = "-";
path_type output_target_source_file = "-";

bool individual_mode = false;
bool symmetric_mode = false;
bool posterior_mode = false;

// parameter...
double p0    = 1e-4;
double prior = 1e-4;

int threads = 2;

int debug = 0;

struct LearnIndividual;
struct LearnSymmetric;
struct LearnPosterior;

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (int(individual_mode) + symmetric_mode + posterior_mode > 1)
      throw std::runtime_error("specify either individual|symmetric|posterior");
    
    if (int(individual_mode) + symmetric_mode + posterior_mode == 0)
      individual_mode = true;
    
    threads = std::max(threads, 1);
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct LearnBase
{
  LearnBase(const count_dict_type& __ttable_source_target)
    : ttable_source_target(__ttable_target_source),
      ttable_target_source(__ttable_source_target)
  {}
  
  const count_dict_type& ttable_source_target;
  const count_dict_type& ttable_target_source;
  count_dict_type counts_source_target;
  count_dict_type counts_target_source;
};

struct LearnIndividual : public LearnBase
{
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    
    
    
  }
};

struct LearnSymmetric : public LearnBase
{

};

struct LearnPosterior : public LearnBase
{
  
    
  
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    
    ("output_source_targert", po::value<path_type>(&output_source_target_file), "output for P(target | source)")
    ("output_targert_source", po::value<path_type>(&output_source_target_file), "output for P(source | target)")
    
    ("individual", po::bool_switch(&individual_mode), "individual model1 training")
    ("symmetric",  po::bool_switch(&symmetric_mode),  "symmetric model1 training")
    ("posterior",  po::bool_switch(&posterior_mode),  "posterior constrained model1 training")
    
    ("p0",    po::value<double>(&p0),    "parameter for NULL alignment")
    ("prior", po::value<double>(&prior), "prior for Dirichlet smoothing")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
