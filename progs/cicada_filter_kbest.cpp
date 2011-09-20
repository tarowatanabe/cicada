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
#include "utils/hashmurmur.hpp"
#include "utils/simple_vector.hpp"
#include "utils/sgi_hash_set.hpp"
#include "utils/lexical_cast.hpp"

typedef boost::filesystem::path path_type;

typedef size_t size_type;
typedef std::vector<std::string, std::allocator<std::string> > tokens_type;
typedef boost::fusion::tuple<size_type, tokens_type, tokens_type> kbest_type;

typedef std::pair<std::string, double> feature_type;
typedef std::vector<feature_type, std::allocator<feature_type> > features_type;

typedef boost::fusion::tuple<size_type, tokens_type, features_type> kbest_feature_type;

struct hypothesis_type
{
  typedef cicada::Symbol  word_type;
  typedef cicada::Feature feature_type;
  typedef std::pair<feature_type, double> feature_value_type;
  
  typedef utils::simple_vector<word_type, std::allocator<word_type> >                   sentence_type;
  typedef utils::simple_vector<feature_value_type, std::allocator<feature_value_type> > feature_set_type;
  
  hypothesis_type() : sentence(), features() {}
  hypothesis_type(const kbest_feature_type& x)
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
  typedef utils::hashmurmur<size_t> hasher_type;
  
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
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    remains %= *qi::lexeme[+(standard::char_ - standard::space)];
    
    kbest %= size >> "|||" >> tokens >> -("|||" >> remains) >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  boost::spirit::qi::rule<Iterator, kbest_type(), blank_type>  kbest;
};

template <typename Iterator>
struct kbest_feature_parser : boost::spirit::qi::grammar<Iterator, kbest_feature_type(), boost::spirit::standard::blank_type>
{
  kbest_feature_parser() : kbest_feature_parser::base_type(kbest)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
    remains %= *qi::lexeme[+(standard::char_ - standard::space)];
    
    feature %= qi::lexeme[+(standard::char_ - standard::space - '=')] >> '=' >> qi::double_;
    features %= *feature;
    
    kbest %= size >> "|||" >> tokens >> "|||" >> features >> -("|||" >> remains) >> (qi::eol | qi::eoi);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::uint_parser<size_type, 10, 1, -1>         size;
  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> tokens;
  
  boost::spirit::qi::rule<Iterator, feature_type(), blank_type>  feature;
  boost::spirit::qi::rule<Iterator, features_type(), blank_type> features;

  boost::spirit::qi::rule<Iterator, tokens_type(), blank_type> remains;
  
  boost::spirit::qi::rule<Iterator, kbest_feature_type(), blank_type>  kbest;
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
       std::ostream&    __os,
       const bool __id_mode,
       const bool __features_mode)
    : queue(__queue),
      subprocess(__subprocess),
      os(__os),
      id_mode(__id_mode),
      features_mode(__features_mode) {}
  
  queue_type& queue;
  subprocess_type& subprocess;
  std::ostream& os;
  
  bool id_mode;
  bool features_mode;

  void operator()()
  {
    boost::iostreams::filtering_istream is;
#if BOOST_VERSION >= 104400
    is.push(boost::iostreams::file_descriptor_source(subprocess.desc_read(), boost::iostreams::close_handle));
#else
    is.push(boost::iostreams::file_descriptor_source(subprocess.desc_read(), true));
#endif
    
    kbest_type kbest;
    std::string line;
    while (queue.pop(kbest)) {
      if (boost::fusion::get<0>(kbest) == size_type(-1)) break;
      
      // read from is...
      if (! std::getline(is, line))
	throw std::runtime_error("# of lines do not match");
      
      if (id_mode)
	os << boost::fusion::get<0>(kbest) << " ||| ";
      
      os << line;
      
      if (features_mode) {
	const tokens_type& features = boost::fusion::get<2>(kbest);
	
	os << " ||| ";
	if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -(boost::spirit::standard::string % ' '), features))
	  throw std::runtime_error("tokens generation failed...?");
      }
      
      os << '\n';
    }
  }
};

void moses_to_cicada(tokens_type& features)
{
  std::string feature_name = "";
  int id = 0;
  tokens_type features_new;
  
  tokens_type::const_iterator fiter_end = features.end();
  for (tokens_type::const_iterator fiter = features.begin(); fiter != fiter_end; ++ fiter) {
    const size_t feature_size = fiter->size();
    
    if (fiter->operator[](feature_size - 1) == ':') {
      feature_name = *fiter;
      id = 0;
    } else {
      features_new.push_back(feature_name + utils::lexical_cast<std::string>(id) + "=" + *fiter);
      ++ id;
    }
  }
  features.swap(features_new);
}

path_type input_file = "-";
path_type output_file = "-";

std::string filter;

bool id_mode = false;
bool features_mode = false;
bool moses_mode = false;

// alternative mode...
bool merge_mode = false;
bool lattice_mode = false;

int debug = 0;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (int(merge_mode) + lattice_mode + (! filter.empty()) > 1)
      throw std::runtime_error("you can either --{merge, lattice} or --filter");
    
    if (merge_mode) {
#ifdef HAVE_TR1_UNORDERED_SET
      typedef std::tr1::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
	std::allocator<hypothesis_type> > hypothesis_set_type;
#else
      typedef sgi::hash_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
	std::allocator<hypothesis_type> > hypothesis_set_type;
#endif
      typedef std::deque<hypothesis_set_type, std::allocator<hypothesis_set_type> > hypothesis_map_type;

      typedef boost::spirit::istream_iterator iter_type;
      
      utils::compress_istream is(input_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);

      boost::spirit::karma::real_generator<double, real_precision20> double20;
      
      kbest_feature_parser<iter_type> parser;
      iter_type iter(is);
      iter_type iter_end;
      
      kbest_feature_type kbest;

      hypothesis_map_type hypotheses;
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
      
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
	
	const size_t& id = boost::fusion::get<0>(kbest);
	
	if (id >= hypotheses.size())
	  hypotheses.resize(id + 1);
	
	hypotheses[id].insert(hypothesis_type(kbest));
      }
      
      utils::compress_ostream os(output_file, 1024 * 1024);
      
      for (size_t id = 0; id != hypotheses.size(); ++ id) {
	namespace karma = boost::spirit::karma;
	namespace standard = boost::spirit::standard;
	
	hypothesis_set_type::const_iterator hiter_end = hypotheses[id].end();
	for (hypothesis_set_type::const_iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter) {
	  os << id << " ||| ";
	  
	  if (! karma::generate(std::ostream_iterator<char>(os), -(standard::string % ' '), hiter->sentence))
	    throw std::runtime_error("tokens generation failed...?");
	  os << " ||| ";
	  if (! karma::generate(std::ostream_iterator<char>(os), -((standard::string << '=' << double20) % ' '), hiter->features))
	    throw std::runtime_error("tokens generation failed...?");
	  os << '\n';
	}
      }
      
    } else if (lattice_mode) {
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
      
      kbest_feature_parser<iter_type> parser;
      iter_type iter(is);
      iter_type iter_end;
      
      kbest_feature_type kbest;

      size_t id = size_t(-1);
      hypothesis_set_type hypotheses;
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
      
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
	
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
      
      task_type::subprocess_type subprocess(filter);
      task_type::queue_type      queue;
      
      boost::thread thread(task_type(queue, subprocess, os, id_mode, features_mode));

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
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  const tokens_type& tokens = boost::fusion::get<1>(kbest);
	  
	  if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os_filter), -(boost::spirit::standard::string % ' '), tokens))
	    throw std::runtime_error("tokens generation failed...?");
	  os_filter << '\n';
	  
	  queue.push(kbest);
	}
      }
      
      boost::fusion::get<0>(kbest) = size_type(-1);
      queue.push(kbest);
      
      thread.join();
    } else {
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
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
      
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
      
	if (id_mode)
	  os << boost::fusion::get<0>(kbest) << " ||| ";
      
	const tokens_type& tokens = boost::fusion::get<1>(kbest);
      
	if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -(boost::spirit::standard::string % ' '), tokens))
	  throw std::runtime_error("tokens generation failed...?");
      
	if (features_mode) {
	  const tokens_type& features = boost::fusion::get<2>(kbest);
	  
	  os << " ||| ";
	  if (! boost::spirit::karma::generate(std::ostream_iterator<char>(os), -(boost::spirit::standard::string % ' '), features))
	    throw std::runtime_error("tokens generation failed...?");
	}

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

    ("id",       po::bool_switch(&id_mode),       "output id")
    ("features", po::bool_switch(&features_mode), "output features")
    ("moses",    po::bool_switch(&moses_mode),    "features in moses format")
    ("merge",    po::bool_switch(&merge_mode),    "merge features")
    ("lattice",  po::bool_switch(&lattice_mode),  "output merged lattice")
    
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
