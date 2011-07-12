//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// bitext filter

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <memory>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>


#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <cicada/vocab.hpp>
#include <cicada/stemmer.hpp>

#include "utils/bithack.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/subprocess.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"

typedef cicada::Vocab     vocab_type;

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef std::vector<std::string, std::allocator<std::string> > sentence_type;

template <typename Iterator>
struct sentence_parser : boost::spirit::qi::grammar<Iterator, sentence_type(), boost::spirit::standard::blank_type>
{
  sentence_parser() : sentence_parser::base_type(tokens)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    tokens  %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"] >> (qi::eoi | qi::eol);
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, sentence_type(), blank_type> tokens;
};

template <typename Iterator>
struct non_terminal_parser : boost::spirit::qi::grammar<Iterator, std::string()>
{
  non_terminal_parser() : non_terminal_parser::base_type(string)  
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    string %= '[' >> +(standard::char_ - ',' - ']') >> -(',' >> qi::int_) >> ']';
  }
  
  boost::spirit::qi::rule<Iterator, std::string()> string;
};

template <typename Iterator>
struct sentence_generator : boost::spirit::karma::grammar<Iterator, sentence_type()>
{
  sentence_generator() : sentence_generator::base_type(tokens)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    tokens  %= -(standard::string % ' ') << '\n';
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::karma::rule<Iterator, sentence_type()> tokens;
};


path_set_type input_files;
path_type list_file;
path_type output_file = "-";

std::string stemmer_spec;
bool stemmer_list = false;

int max_length = 0;

bool add_bos_eos = false;

void read_list(const path_type& path, path_set_type& files);
void options(int argc, char** argv);

template <typename Iterator, typename Grammar>
inline
bool verify(Iterator first, Iterator last, const Grammar& grammar)
{
  namespace qi = boost::spirit::qi;
  
  for (/**/; first != last; ++ first) {
    std::string::const_iterator iter = first->begin();
    std::string::const_iterator end = first->end();
    
    if (qi::parse(iter, end, grammar) && iter == end) return false;
  }
  return true;
}

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (stemmer_list) {
      std::cout << cicada::Stemmer::lists();
      return 1;
    }

    cicada::Stemmer* stemmer = (! stemmer_spec.empty() ? &cicada::Stemmer::create(stemmer_spec) : 0);
    
    read_list(list_file, input_files);

    if (input_files.empty())
      input_files.push_back("-");
    
    const std::string bos = static_cast<const std::string&>(vocab_type::BOS);
    const std::string eos = static_cast<const std::string&>(vocab_type::EOS);

    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));

    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    typedef boost::spirit::istream_iterator iiter_type;
    typedef std::ostream_iterator<char>     oiter_type;
    
    sentence_parser<iiter_type>    parser;
    sentence_generator<oiter_type> generator;
    non_terminal_parser<std::string::const_iterator> non_terminal_parser;

    for (path_set_type::const_iterator fiter = input_files.begin(); fiter != input_files.end(); ++ fiter) {
      utils::compress_istream is(*fiter, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iiter_type iter(is);
      iiter_type iter_end;
      
      sentence_type sentence;
      
      for (size_t line_no = 0; iter != iter_end; ++ line_no) {
	sentence.clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, sentence))
	  throw std::runtime_error("sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	
	if (! verify(sentence.begin(), sentence.end(), non_terminal_parser))
	  throw std::runtime_error("sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	
	if (sentence.size() == 0) continue;
	if (max_length > 0 && sentence.size() > max_length) continue;
	
	if (stemmer) {
	  sentence_type::iterator siter_end = sentence.end();
	  for (sentence_type::iterator siter = sentence.begin(); siter != siter_end; ++ siter)
	    *siter = stemmer->operator()(*siter);
	}
	
	if (add_bos_eos) {
	  sentence.insert(sentence.begin(), bos);
	  sentence.push_back(eos);
	}
	
	if (! boost::spirit::karma::generate(oiter_type(os), generator, sentence))
	  throw std::runtime_error("source sentence generation failed at # " + utils::lexical_cast<std::string>(line_no));
      }
    }
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}

void read_list(const path_type& path, path_set_type& files)
{
  if (path.empty()) return;
  if (path != "-" && ! boost::filesystem::exists(path))
    throw std::runtime_error("no file? " + path.string());
  
  utils::compress_istream is(path);
  std::string file;
  while (std::getline(is, file)) {
    if (file.empty()) continue;
    if (! boost::filesystem::exists(file))
      throw std::runtime_error("no file? " + file);
    files.push_back(file);
  }
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",  po::value<path_set_type>(&input_files)->multitoken(),           "input file(s)")
    ("list",   po::value<path_type>(&list_file),                               "list file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("stemmer",      po::value<std::string>(&stemmer_spec), "stemmer")
    ("stemmer-list", po::bool_switch(&stemmer_list),        "list of stemmers")
    
    ("max-length",    po::value<int>(&max_length)->default_value(max_length),          "maximum length")
    ("add-bos-eos",   po::bool_switch(&add_bos_eos), "add BOS/EOS for each sentence")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

