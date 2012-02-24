//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <map>
#include <cmath>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/dense_hash_map.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef cicada::Symbol   word_type;
typedef cicada::Sentence sentence_type;
typedef boost::filesystem::path path_type;

typedef uint64_t count_type;
typedef std::vector<count_type, std::allocator<count_type> > count_unigram_type;

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
path_type output_file = "-";
bool inverse = false;

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    count_unigram_type sources;
    count_unigram_type targets;
    count_dict_type    dict;

    utils::compress_istream is_src(! inverse ? source_file : target_file, 1024 * 1024);
    utils::compress_istream is_trg(! inverse ? target_file : source_file, 1024 * 1024);

    sentence_type source;
    sentence_type target;
    
    while (1) {
      is_src >> source;
      is_trg >> target;

      if (! is_src || ! is_trg) break;

      if (source.empty() || target.empty()) continue;
      
      sentence_type::const_iterator siter_begin = source.begin();
      sentence_type::const_iterator siter_end   = source.end();
      sentence_type::const_iterator titer_begin = target.begin();
      sentence_type::const_iterator titer_end   = target.end();

      for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	if (titer->id() >= targets.size())
	  targets.resize(titer->id() + 1, 0);
	++ targets[titer->id()];
      }
      
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	if (siter->id() >= sources.size())
	  sources.resize(siter->id() + 1, 0);
	++ sources[siter->id()];
	
	for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	  ++ dict[siter->id()][*titer];
      }
    }
    if (is_src || is_trg)
      throw std::runtime_error("# of lines do not match");
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    os.precision(10);
    
    typedef std::multimap<double, word_type, std::greater<double>, std::allocator<std::pair<const double, word_type> > > sorted_type;
    
    sorted_type sorted;
    
    // dump..
    count_dict_type::const_iterator siter_begin = dict.begin();
    count_dict_type::const_iterator siter_end   = dict.end();
    for (count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
      if (*siter) {
	const word_type source(word_type::id_type(siter - siter_begin));
	const count_map_type& target = *(*siter);

	sorted.clear();
	
	count_map_type::const_iterator titer_end = target.end();
	for (count_map_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) {
	  const double score = ((2.0 * titer->second) / (sources[source.id()] + targets[titer->first.id()]));
	  
	  sorted.insert(std::make_pair(score, titer->first));
	}
	
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  os << iter->second << ' ' << source << ' ' << iter->first << '\n';
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
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    ("output", po::value<path_type>(&output_file)->default_value("-"), "dice: for target source")
    ("inverse", po::bool_switch(&inverse), "inverse source/target")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
