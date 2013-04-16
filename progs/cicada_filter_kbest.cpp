//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// kbest filter:
//
// id ||| sentence ||| features etc.
//
// TODO: add subprocess filter so that we can uncover the kbest-format again...
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <deque>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/variant.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <cicada/feature.hpp>
#include <cicada/symbol.hpp>
#include <cicada/lattice.hpp>
#include <cicada/sentence.hpp>
#include <cicada/unite.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/subprocess.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/simple_vector.hpp"
#include "utils/unordered_set.hpp"
#include "utils/lexical_cast.hpp"

#include "cicada_output_impl.hpp"

typedef boost::filesystem::path path_type;

typedef size_t size_type;
typedef cicada::Symbol  word_type;
typedef cicada::Feature feature_type;
typedef std::pair<feature_type, double> feature_value_type;

typedef std::vector<word_type, std::allocator<word_type> > tokens_type;
typedef std::pair<feature_type, double> feature_value_type;
typedef std::vector<feature_value_type, std::allocator<feature_value_type> > features_type;
typedef boost::fusion::tuple<size_type, tokens_type, features_type> kbest_type;

typedef utils::unordered_set<feature_type, boost::hash<feature_type>, std::equal_to<feature_type>,
			     std::allocator<feature_type> >::type feature_unique_type;
typedef std::vector<feature_type, std::allocator<feature_type> > feature_list_type;

struct hypothesis_type
{
  typedef utils::simple_vector<word_type, std::allocator<word_type> >                   sentence_type;
  typedef utils::simple_vector<feature_value_type, std::allocator<feature_value_type> > feature_set_type;
  
  hypothesis_type() : sentence(), features() {}
  hypothesis_type(const kbest_type& x)
    : sentence(boost::fusion::get<1>(x).begin(), boost::fusion::get<1>(x).end()),
      features(boost::fusion::get<2>(x).begin(), boost::fusion::get<2>(x).end())
  {
    std::sort(features.begin(), features.end());
  }
  
  sentence_type    sentence;
  feature_set_type features;
};

inline
size_t hash_value(hypothesis_type const& x)
{
  typedef utils::hashmurmur3<size_t> hasher_type;
  
  return hasher_type()(x.sentence.begin(), x.sentence.end(), hasher_type()(x.features.begin(), x.features.end(), 0));
}

inline
bool operator==(const hypothesis_type& x, const hypothesis_type& y)
{
  return x.sentence == y.sentence && x.features == y.features;
}

template <typename Iterator>
struct kbest_parser : boost::spirit::qi::grammar<Iterator, kbest_type(), boost::spirit::standard::blank_type>
{
  kbest_parser() : kbest_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    word    %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];
    tokens  %= *word;
    
    feature %= qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi)) >> (standard::char_ - standard::space))] >> '=' >> qi::double_;
    features %= -(feature % (+standard::space));
    
    remains %= *qi::lexeme[+(standard::char_ - standard::space)];
    
    kbest %= size >> "|||" >> tokens >> "|||" >> features >> -qi::omit["|||" >> remains] >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type> word;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  
  boost::spirit::qi::rule<Iterator, std::pair<std::string, double>()> feature;
  boost::spirit::qi::rule<Iterator, features_type()>                  features;

  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  
  boost::spirit::qi::rule<Iterator, kbest_type(), blank_type>  kbest;
};


typedef boost::fusion::tuple<size_type, std::string, features_type> kbest_json_type;

template <typename Iterator>
struct kbest_json_parser : boost::spirit::qi::grammar<Iterator, kbest_json_type(), boost::spirit::standard::blank_type>
{
  typedef int64_t     int_type;
  typedef double      float_type;
  typedef std::string string_type;
  
  kbest_json_parser() : kbest_json_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    
    
    feature %= qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi)) >> (standard::char_ - standard::space))] >> '=' >> qi::double_;
    features %= -(feature % (+standard::space));
    
    remains %= *qi::lexeme[+(standard::char_ - standard::space)];
    
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1> size;
  
  boost::spirit::qi::rule<Iterator, std::pair<std::string, double>()> feature;
  boost::spirit::qi::rule<Iterator, features_type()>                  features;
  
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  
  boost::spirit::qi::rule<Iterator, kbest_json_type(), blank_type>  kbest;
};

struct real_precision20 : boost::spirit::karma::real_policies<double>
{
  static unsigned int precision(double) 
  { 
    return 20;
  }
};


struct Task
{
  typedef utils::lockfree_list_queue<kbest_type, std::allocator<kbest_type> > queue_type;
  typedef utils::subprocess subprocess_type;

  Task(queue_type&      __queue,
       subprocess_type& __subprocess,
       std::ostream&    __os)
    : queue(__queue),
      subprocess(__subprocess),
      os(__os) {}
  
  queue_type& queue;
  subprocess_type& subprocess;
  std::ostream& os;

  void operator()()
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    boost::iostreams::filtering_istream is;
#if BOOST_VERSION >= 104400
    is.push(boost::iostreams::file_descriptor_source(subprocess.desc_read(), boost::iostreams::close_handle));
#else
    is.push(boost::iostreams::file_descriptor_source(subprocess.desc_read(), true));
#endif
    
    karma::real_generator<double, real_precision20> double20;
    
    kbest_type kbest;
    std::string line;
    while (queue.pop(kbest)) {
      if (boost::fusion::get<0>(kbest) == size_type(-1)) break;
      
      // read from is...
      if (! std::getline(is, line))
	throw std::runtime_error("# of lines do not match");
      
      os << boost::fusion::get<0>(kbest) << " ||| ";
      os << line;
      os << " ||| ";
      if (! karma::generate(std::ostream_iterator<char>(os),
			    -((standard::string << '=' << double20) % ' '),
			    boost::fusion::get<2>(kbest)))
	throw std::runtime_error("tokens generation failed...?");
      
      os << '\n';
    }
  }
};


path_type input_file = "-";
path_type output_file = "-";

std::string filter;

feature_list_type features_removes;

// alternative mode...
bool merge_mode = false;
bool lattice_mode = false;
bool directory_mode = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (int(merge_mode) + lattice_mode + (! filter.empty()) > 1)
      throw std::runtime_error("you can either --{merge, lattice} or --filter");
    
    feature_unique_type removes(features_removes.begin(), features_removes.end());
    
    const bool directory_input_mode = (boost::filesystem::exists(input_file)
				       && boost::filesystem::is_directory(input_file));
    const bool directory_output_mode = directory_mode;
    
    if (merge_mode) {
      namespace qi = boost::spirit::qi;
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      typedef utils::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
				   std::allocator<hypothesis_type> >::type hypothesis_set_type;
      typedef std::deque<hypothesis_set_type, std::allocator<hypothesis_set_type> > hypothesis_map_type;

      typedef boost::spirit::istream_iterator iter_type;

      hypothesis_map_type hypotheses;
      
      kbest_parser<iter_type> parser;
      karma::real_generator<double, real_precision20> double20;
      
      
      kbest_type kbest;
      features_type features_removed;
      
      if (directory_input_mode) {
	for (size_t i = 0; /**/; ++ i) {
	  const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	  
	  const path_type path_input = input_file / file_name;
	  
	  if (! boost::filesystem::exists(path_input)) break;
	  
	  utils::compress_istream is(path_input, 1024 * 1024);
	  is.unsetf(std::ios::skipws);
	  
	  iter_type iter(is);
	  iter_type iter_end;
	  
	  while (iter != iter_end) {
	    boost::fusion::get<1>(kbest).clear();
	    boost::fusion::get<2>(kbest).clear();
	    
	    if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	      if (iter != iter_end)
		throw std::runtime_error("kbest parsing failed: merge+directory");
	  
	    if (! removes.empty()) {
	      features_removed.clear();
	      
	      features_type::const_iterator fiter_end = boost::fusion::get<2>(kbest).end();
	      for (features_type::const_iterator fiter = boost::fusion::get<2>(kbest).begin(); fiter != fiter_end; ++ fiter)
		if (removes.find(fiter->first) == removes.end())
		  features_removed.push_back(*fiter);
	    
	      boost::fusion::get<2>(kbest).swap(features_removed);
	    }
	  
	    const size_t& id = boost::fusion::get<0>(kbest);

	    if (id != i)
	      throw std::runtime_error("invalid directory input format");
	    
	    if (id >= hypotheses.size())
	      hypotheses.resize(id + 1);
	    
	    hypotheses[id].insert(hypothesis_type(kbest));
	  }	  
	}
      } else {
	utils::compress_istream is(input_file, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed: merge");
	  
	  if (! removes.empty()) {
	    features_removed.clear();
	    
	    features_type::const_iterator fiter_end = boost::fusion::get<2>(kbest).end();
	    for (features_type::const_iterator fiter = boost::fusion::get<2>(kbest).begin(); fiter != fiter_end; ++ fiter)
	      if (removes.find(fiter->first) == removes.end())
		features_removed.push_back(*fiter);
	    
	    boost::fusion::get<2>(kbest).swap(features_removed);
	  }
	  
	  const size_t& id = boost::fusion::get<0>(kbest);
	  
	  if (id >= hypotheses.size())
	    hypotheses.resize(id + 1);
	  
	  hypotheses[id].insert(hypothesis_type(kbest));
	}
      }
      
      if (directory_output_mode) {
	prepare_directory(output_file);
	
	for (size_t id = 0; id != hypotheses.size(); ++ id) {
	  const std::string file_name = utils::lexical_cast<std::string>(id) + ".gz";
	  
	  utils::compress_ostream os(output_file / file_name, 1024 * 1024);
	  
	  hypothesis_set_type::const_iterator hiter_end = hypotheses[id].end();
	  for (hypothesis_set_type::const_iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter) {
	    os << id << " ||| ";
	    
	    if (! karma::generate(std::ostream_iterator<char>(os), -(standard::string % ' '), hiter->sentence))
	      throw std::runtime_error("tokens generation failed...?");
	    os << " ||| ";
	    if (! karma::generate(std::ostream_iterator<char>(os),
				  -((standard::string << '=' << double20) % ' '),
				  hiter->features))
	      throw std::runtime_error("tokens generation failed...?");
	    os << '\n';
	  }
	}
      } else {
	utils::compress_ostream os(output_file, 1024 * 1024);
	
	for (size_t id = 0; id != hypotheses.size(); ++ id) {
	  hypothesis_set_type::const_iterator hiter_end = hypotheses[id].end();
	  for (hypothesis_set_type::const_iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter) {
	    os << id << " ||| ";
	    
	    if (! karma::generate(std::ostream_iterator<char>(os), -(standard::string % ' '), hiter->sentence))
	      throw std::runtime_error("tokens generation failed...?");
	    os << " ||| ";
	    if (! karma::generate(std::ostream_iterator<char>(os),
				  -((standard::string << '=' << double20) % ' '),
				  hiter->features))
	      throw std::runtime_error("tokens generation failed...?");
	    os << '\n';
	  }
	}
      }
      
    } else if (lattice_mode) {
      namespace qi = boost::spirit::qi;
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      typedef std::deque<hypothesis_type, std::allocator<hypothesis_type> > hypothesis_set_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice  lattice_type;
      
      typedef boost::spirit::istream_iterator iter_type;

      const bool flush_output = (output_file == "-"
				 || (boost::filesystem::exists(output_file)
				     && ! boost::filesystem::is_regular_file(output_file)));
      
      utils::compress_istream is(input_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
      
      kbest_parser<iter_type> parser;
      iter_type iter(is);
      iter_type iter_end;
      
      kbest_type kbest;
      features_type features_removed;

      size_t id = size_t(-1);
      hypothesis_set_type hypotheses;
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
      
	if (! qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed: lattice");
	
	if (! removes.empty()) {
	  features_removed.clear();
	  
	  features_type::const_iterator fiter_end = boost::fusion::get<2>(kbest).end();
	  for (features_type::const_iterator fiter = boost::fusion::get<2>(kbest).begin(); fiter != fiter_end; ++ fiter)
	    if (removes.find(fiter->first) == removes.end())
	      features_removed.push_back(*fiter);
	  
	  boost::fusion::get<2>(kbest).swap(features_removed);
	}
	
	if (boost::fusion::get<0>(kbest) != id) {
	  
	  if (! hypotheses.empty()) {
	    if (debug)
	      std::cerr << id << " merge lattice" << std::endl;
	    
	    lattice_type lattice;
	    
	    hypothesis_set_type::const_iterator hiter_end = hypotheses.end();
	    for (hypothesis_set_type::const_iterator hiter = hypotheses.begin(); hiter != hiter_end; ++ hiter) {
	      lattice_type lattice_local(sentence_type(hiter->sentence.begin(), hiter->sentence.end()));
	      
	      lattice_local.front().front().features.assign(hiter->features.begin(), hiter->features.end());
	      
	      cicada::unite(lattice, lattice_local);
	    }
	    
	    os << id << " ||| " << lattice << '\n';
	  }
	  
	  hypotheses.clear();
	}
	
	id = boost::fusion::get<0>(kbest);
	hypotheses.push_back(hypothesis_type(kbest));
      }
      
      if (! hypotheses.empty()) {
	if (debug)
	  std::cerr << id << " merge lattice" << std::endl;
	
	lattice_type lattice;
	
	hypothesis_set_type::const_iterator hiter_end = hypotheses.end();
	for (hypothesis_set_type::const_iterator hiter = hypotheses.begin(); hiter != hiter_end; ++ hiter) {
	  lattice_type lattice_local(sentence_type(hiter->sentence.begin(), hiter->sentence.end()));
	  
	  lattice_local.front().front().features.assign(hiter->features.begin(), hiter->features.end());
	  
	  cicada::unite(lattice, lattice_local);
	}
	
	os << id << " ||| " << lattice << '\n';
	
	hypotheses.clear();
      }
      
    } else if (! filter.empty()) {
      namespace qi = boost::spirit::qi;
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      typedef Task task_type;
      typedef boost::spirit::istream_iterator iter_type;

      const bool flush_output = (output_file == "-"
				 || (boost::filesystem::exists(output_file)
				     && ! boost::filesystem::is_regular_file(output_file)));
      
      utils::compress_istream is(input_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
      
      kbest_parser<iter_type> parser;
      iter_type iter(is);
      iter_type iter_end;
      
      kbest_type kbest;
      features_type features_removed;
      
      task_type::subprocess_type subprocess(filter);
      task_type::queue_type      queue;
      
      boost::thread thread(task_type(queue, subprocess, os));

      {
	boost::iostreams::filtering_ostream os_filter;
#if BOOST_VERSION >= 104400
	os_filter.push(boost::iostreams::file_descriptor_sink(subprocess.desc_write(), boost::iostreams::close_handle));
#else
	os_filter.push(boost::iostreams::file_descriptor_sink(subprocess.desc_write(), true));
#endif
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! qi::phrase_parse(iter, iter_end, parser, standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed: filter");
	  
	  if (! removes.empty()) {
	    features_removed.clear();
	    
	    features_type::const_iterator fiter_end = boost::fusion::get<2>(kbest).end();
	    for (features_type::const_iterator fiter = boost::fusion::get<2>(kbest).begin(); fiter != fiter_end; ++ fiter)
	      if (removes.find(fiter->first) == removes.end())
		features_removed.push_back(*fiter);
	    
	    boost::fusion::get<2>(kbest).swap(features_removed);
	  }
	  
	  if (! karma::generate(std::ostream_iterator<char>(os_filter),
				-(standard::string % ' '),
				boost::fusion::get<1>(kbest)))
	    throw std::runtime_error("tokens generation failed...?");
	  os_filter << '\n';
	  
	  queue.push(kbest);
	}
      }
      
      boost::fusion::get<0>(kbest) = size_type(-1);
      queue.push(kbest);
      
      thread.join();
    } else {
      namespace qi = boost::spirit::qi;
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      typedef boost::spirit::istream_iterator iter_type;

      const bool flush_output = (output_file == "-"
				 || (boost::filesystem::exists(output_file)
				     && ! boost::filesystem::is_regular_file(output_file)));
      
      utils::compress_istream is(input_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));

      karma::real_generator<double, real_precision20> double20;
      
      kbest_parser<iter_type> parser;
      iter_type iter(is);
      iter_type iter_end;
      
      kbest_type kbest;
      features_type features_removed;
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
      
	if (! qi::phrase_parse(iter, iter_end, parser, standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
      
	if (! removes.empty()) {
	  features_removed.clear();
	  
	  features_type::const_iterator fiter_end = boost::fusion::get<2>(kbest).end();
	  for (features_type::const_iterator fiter = boost::fusion::get<2>(kbest).begin(); fiter != fiter_end; ++ fiter)
	    if (removes.find(fiter->first) == removes.end())
	      features_removed.push_back(*fiter);
	  
	  boost::fusion::get<2>(kbest).swap(features_removed);
	}
	
	os << boost::fusion::get<0>(kbest) << " ||| ";
	
	if (! karma::generate(std::ostream_iterator<char>(os),
			      -(standard::string % ' '),
			      boost::fusion::get<1>(kbest)))
	  throw std::runtime_error("tokens generation failed...?");
	
	os << " ||| ";
	if (! karma::generate(std::ostream_iterator<char>(os),
			      -((standard::string << '=' << double20) % ' '),
			      boost::fusion::get<2>(kbest)))
	  throw std::runtime_error("tokens generation failed...?");
	
	os << '\n';
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
    
    ("filter", po::value<std::string>(&filter), "filter for sentences")
    ("feature-remove", po::value<feature_list_type>(&features_removes)->multitoken(), "remove featureso")

    ("merge",     po::bool_switch(&merge_mode),     "merge features")
    ("lattice",   po::bool_switch(&lattice_mode),   "output merged lattice")
    ("directory", po::bool_switch(&directory_mode), "output in directory")
    
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
