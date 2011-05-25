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

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <cicada/vocab.hpp>
#include <cicada/alignment.hpp>

#include "utils/bithack.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/subprocess.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"

typedef cicada::Vocab     vocab_type;
typedef cicada::Alignment alignment_type;

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


path_set_type source_files;
path_set_type target_files;
path_set_type alignment_files;

path_type list_source_file;
path_type list_target_file;
path_type list_alignment_file;

path_type output_source_file = "-";
path_type output_target_file = "-";
path_type output_alignment_file = "-";

int max_length = 100;
double max_fertility = 10;

bool add_bos_eos = false;

void read_list(const path_type& path, path_set_type& files);
void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    read_list(list_source_file, source_files);
    read_list(list_target_file, target_files);
    read_list(list_alignment_file, alignment_files);

    if (source_files.empty())
      source_files.push_back("-");
    if (target_files.empty())
      targetfiles.push_back("-");
    
    const bool alignment_mode = ! alignment_files.empty();

    if (source_files.size() != target_files.size())
      throw std::runtime_error("# of files do not match");
    if (alignment_mode) {
      if (source_files.size() != alignment_files.size())
	throw std::runtime_error("# of alignemnt files do not match");
      if (target_files.size() != alignment_files.size())
	throw std::runtime_error("# of alignemnt files do not match");
    }
    
    const std::string bos = static_cast<const std::string&>(vocab_type::BOS);
    const std::string eos = static_cast<const std::string&>(vocab_type::EOS);

    utils::compress_ostream os_src(output_source_file, 1024 * 1024);
    utils::compress_ostream os_trg(output_target_file, 1024 * 1024);
    std::auto_ptr<utils::compress_ostream> os_align(alignment_mode ? new utils::compress_ostream(output_alignment_file, 1024 * 1024) : 0);

    typedef boost::spirit::istream_iterator iiter_type;
    typedef std::ostream_iterator<char>     oiter_type;
    
    sentence_parser<iiter_type>    parser;
    sentence_generator<oiter_type> generator;
    
    for (size_t i = 0; i != source_files.size(); ++ i) {
      utils::compress_istream is_src(source_files[i], 1024 * 1024);
      utils::compress_istream is_trg(target_files[i], 1024 * 1024);
      std::auto_ptr<utils::compress_istream> is_align(alignment_mode ? new utils::compress_istream(alignment_files[i], 1024 * 1024) : 0);

      is_src.unsetf(std::ios::skipws);
      is_trg.unsetf(std::ios::skipws);
      
      iiter_type siter(is_src);
      iiter_type titer(is_trg);
      iiter_type siter_end;
      iiter_type titer_end;
      
      sentence_type source;
      sentence_type target;
      alignment_type alignment;
      alignment_type alignment_new;
    
      for (size_t line_no = 0; siter != siter_end && titer != titer_end; ++ line_no) {
	source.clear();
	target.clear();
	alignment.clear();
      
	if (! boost::spirit::qi::phrase_parse(siter, siter_end, parser, boost::spirit::standard::blank, source))
	  throw std::runtime_error("source sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	if (! boost::spirit::qi::phrase_parse(titer, titer_end, parser, boost::spirit::standard::blank, target))
	  throw std::runtime_error("target sentence parsing failed at # " + utils::lexical_cast<std::string>(line_no));
	
	if (alignment_mode) {
	  *is_align >> alignment;
	  if (! *is_align)
	    throw std::runtime_error("no alignment?");
	}
      
	const int source_size = source.size();
	const int target_size = target.size();
      
	if (source_size == 0 || target_size == 0) continue;
	if (max_length > 0)
	  if (source_size > max_length || target_size > max_length) continue;
	if (max_fertility > 0)
	  if (double(utils::bithack::max(source_size, target_size)) / double(utils::bithack::min(source_size, target_size)) >= max_fertility) continue;
      
	if (add_bos_eos) {
	  source.insert(source.begin(), bos);
	  source.push_back(eos);
	
	  target.insert(target.begin(), bos);
	  target.push_back(eos);
	
	  alignment_new.clear();
	
	  alignment_new.push_back(std::make_pair(0, 0));
	  alignment_type::const_iterator aiter_end = alignment.end();
	  for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter)
	    alignment_new.push_back(std::make_pair(aiter->source + 1, aiter->target + 1));
	  alignment_new.push_back(std::make_pair(source.size() - 1, target.size() - 1));
	  
	  alignment_new.swap(alignment);
	}

	if (! boost::spirit::karma::generate(oiter_type(os_src), generator, source))
	  throw std::runtime_error("source sentence generation failed at # " + utils::lexical_cast<std::string>(line_no));
	if (! boost::spirit::karma::generate(oiter_type(os_trg), generator, target))
	  throw std::runtime_error("target sentence generation failed at # " + utils::lexical_cast<std::string>(line_no));
	
	if (alignment_mode)
	  *os_align << alignment << '\n';
      }
      
      if (siter != siter_end || titer != titer_end || (alignment_mode && *is_align >> alignment))
	throw std::runtime_error("# of lines do not match: " + source_files[i].string() + " " + target_files[i].string());
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
    ("source",    po::value<path_set_type>(&source_files)->multitoken(),    "source file(s)")
    ("target",    po::value<path_set_type>(&target_files)->multitoken(),    "target file(s)")
    ("alignment", po::value<path_set_type>(&alignment_files)->multitoken(), "alignment file(s)")

    ("list-source",    po::value<path_type>(&list_source_file),    "source list file")
    ("list-target",    po::value<path_type>(&list_target_file),    "target list file")
    ("list-alignment", po::value<path_type>(&list_alignment_file), "alignment list file")

    ("output-source",    po::value<path_type>(&output_source_file)->default_value(output_source_file),       "output source file")
    ("output-target",    po::value<path_type>(&output_target_file)->default_value(output_target_file),       "output target file")
    ("output-alignment", po::value<path_type>(&output_alignment_file)->default_value(output_alignment_file), "output alignment file")
    
    ("max-length",    po::value<int>(&max_length)->default_value(max_length),          "maximum length")
    ("max-fertility", po::value<double>(&max_fertility)->default_value(max_fertility), "maximum fertility")
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

