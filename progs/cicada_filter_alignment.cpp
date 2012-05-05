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
#include <cicada/sentence.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/bithack.hpp"
#include "utils/mathop.hpp"
#include "utils/vector2.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Alignment  alignment_type;
typedef cicada::Dependency permutation_type;
typedef cicada::Sentence   sentence_type;

path_set_type source_files;
path_set_type target_files;
path_set_type alignment_files;
path_set_type permutation_source_files;
path_set_type permutation_target_files;
path_type output_file = "-";

bool inverse_mode = false;
bool visualize_mode = false;

void options(int argc, char** argv);

std::ostream& visualize(std::ostream& os,
			const sentence_type& source,
			const sentence_type& target,
			const alignment_type& alignment);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
        
    if (alignment_files.empty())
      alignment_files.push_back("-");
    
    if (! permutation_source_files.empty())
      if (permutation_source_files.size() != alignment_files.size())
	throw std::runtime_error("# of permutation files does not match");
    
    if (! permutation_target_files.empty())
      if (permutation_target_files.size() != alignment_files.size())
	throw std::runtime_error("# of permutation files does not match");

    if (! source_files.empty())
      if (source_files.size() != alignment_files.size())
	throw std::runtime_error("# of source files does not match");

    if (! target_files.empty())
      if (target_files.size() != alignment_files.size())
	throw std::runtime_error("# of target files does not match");
    
    if (visualize_mode) {
      if (source_files.empty())
	throw std::runtime_error("no source data?");
      if (target_files.empty())
	throw std::runtime_error("no target data?");
    }
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    if (! permutation_source_files.empty() || ! permutation_target_files.empty()) {
      sentence_type    sentence_source;
      sentence_type    sentence_target;
      permutation_type permutation_source;
      permutation_type permutation_target;
      alignment_type   alignment;

      const bool has_source = ! source_files.empty();
      const bool has_target = ! target_files.empty();
      
      const bool has_permutation_source = ! permutation_source_files.empty();
      const bool has_permutation_target = ! permutation_target_files.empty();
      
      for (size_t i = 0; i != alignment_files.size(); ++ i) {
	std::auto_ptr<std::istream> is_source(has_source ? new utils::compress_istream(source_files[i], 1024 * 1024) : 0);
	std::auto_ptr<std::istream> is_target(has_target ? new utils::compress_istream(target_files[i], 1024 * 1024) : 0);
	
	std::auto_ptr<std::istream> ps_source(has_permutation_source ? new utils::compress_istream(permutation_source_files[i], 1024 * 1024) : 0);
	std::auto_ptr<std::istream> ps_target(has_permutation_target ? new utils::compress_istream(permutation_target_files[i], 1024 * 1024) : 0);
	
	utils::compress_istream is(alignment_files[i], 1024 * 1024);
	
	for (;;) {
	  is >> alignment;

	  if (has_source)
	    *is_source >> sentence_source;
	  if (has_target)
	    *is_target >> sentence_target;
	  
	  if (has_permutation_source)
	    *ps_source >> permutation_source;
	  if (has_permutation_target)
	    *ps_target >> permutation_target;
	  
	  if (! is
	      || (ps_source.get() && ! *ps_source)
	      || (ps_target.get() && ! *ps_target)
	      || (is_source.get() && ! *is_source)
	      || (is_target.get() && ! *is_target)) break;
	  
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
	  
	  if (inverse_mode)
	    alignment.inverse();
	  
	  std::sort(alignment.begin(), alignment.end());
	  
	  if (visualize_mode)
	    visualize(os, sentence_source, sentence_target, alignment);
	  else
	    os << alignment << '\n';
	}
	
	if (is
	    || (ps_source.get() && *ps_source)
	    || (ps_target.get() && *ps_target)
	    || (is_source.get() && *is_source)
	    || (is_target.get() && *is_target))
	  throw std::runtime_error("# of samples do not match");
      }
    } else {
      sentence_type  sentence_source;
      sentence_type  sentence_target;
      alignment_type alignment;

      const bool has_source = ! source_files.empty();
      const bool has_target = ! target_files.empty();
      
      for (size_t i = 0; i != alignment_files.size(); ++ i) {
	std::auto_ptr<std::istream> is_source(has_source ? new utils::compress_istream(source_files[i], 1024 * 1024) : 0);
	std::auto_ptr<std::istream> is_target(has_target ? new utils::compress_istream(target_files[i], 1024 * 1024) : 0);
	
	utils::compress_istream is(alignment_files[i], 1024 * 1024);
	
	for (;;) {
	  is >> alignment;
	  
	  if (has_source)
	    *is_source >> sentence_source;
	  if (has_target)
	    *is_target >> sentence_target;
	  
	  if (! is
	      || (is_source.get() && ! *is_source)
	      || (is_target.get() && ! *is_target)) break;
	  
	  if (inverse_mode) {
	    alignment.inverse();
	    std::sort(alignment.begin(), alignment.end());
	  }
	  
	  if (visualize_mode)
	    visualize(os, sentence_source, sentence_target, alignment);
	  else
	    os << alignment << '\n';
	}
	
	if (is
	    || (is_source.get() && *is_source)
	    || (is_target.get() && *is_target))
	  throw std::runtime_error("# of samples do not match");
      }
    }
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}

struct CharacterString
{
  typedef std::vector<UChar32, std::allocator<UChar32> > string_type;
  typedef std::vector<bool, std::allocator<bool> > wide_type;

  CharacterString() : string(), wide() {};
  CharacterString(const std::string& __string)
  {
    icu::UnicodeString ustring = icu::UnicodeString::fromUTF8(__string);
    
    icu::StringCharacterIterator iter(ustring);
    for (iter.setToStart(); iter.hasNext(); /**/) {
      const UChar32 c = iter.next32PostInc();
      const UEastAsianWidth width = static_cast<UEastAsianWidth>(u_getIntPropertyValue(c, UCHAR_EAST_ASIAN_WIDTH));
      
      string.push_back(c);
      wide.push_back(width == U_EA_FULLWIDTH || width == U_EA_WIDE);
    }
  }

  bool empty() const { return string.empty(); }
  size_t size() const { return string.size(); }

  void padding(size_t pad)
  {
    if (pad <= size()) return;

    string_type string_new(pad - size(), ' ');
    wide_type   wide_new(pad - size(), false);

    string_new.insert(string_new.end(), string.begin(), string.end());
    wide_new.insert(wide_new.end(), wide.begin(), wide.end());
    
    string.swap(string_new);
    wide.swap(wide_new);
  }
  
  string_type string;
  wide_type   wide;
};

struct ostream_sink : public ByteSink
{
  
  ostream_sink(std::ostream& _os) : os(_os) {}
  
  virtual void Append(const char* data, int32_t n) 
  {
    os.write((char*) data, n);
  }
  
  void write(char c)
  {
    os.write((char*) & c, sizeof(c));
  }
  
  std::ostream& os;
};

std::ostream& visualize(std::ostream& os,
			const sentence_type& source,
			const sentence_type& target,
			const alignment_type& alignment)
{
  typedef utils::vector2<bool, std::allocator<bool> > matrix_type;
  typedef CharacterString wide_string_type;
  typedef std::vector<wide_string_type, std::allocator<wide_string_type> > wide_string_set_type;
  
  matrix_type matrix(source.size(), target.size(), false);

  alignment_type::const_iterator aiter_end = alignment.end();
  for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
    if (aiter->source >= static_cast<int>(source.size()))
      throw std::runtime_error("invalid alignment");
    if (aiter->target >= static_cast<int>(target.size()))
      throw std::runtime_error("invalid alignment");
    
    matrix(aiter->source, aiter->target) = true;
  }
  
  wide_string_set_type wides(target.size());
  
  size_t size_max = 0;
  for (size_t trg = 0; trg != target.size(); ++ trg) {
    wides[trg] = wide_string_type(target[trg]);
    
    size_max = utils::bithack::max(size_max, wides[trg].size());
  }
  
  for (size_t trg = 0; trg != target.size(); ++ trg)
    wides[trg].padding(size_max);
  
  // start!
  os << '\n';

  // target-side
  ostream_sink sink(os);
  for (size_t i = 0; i != size_max; ++ i) {
    for (size_t trg = 0; trg != target.size(); ++ trg) {
      icu::UnicodeString(wides[trg].string[i]).toUTF8(sink);
      
      if (! wides[trg].wide[i])
	os << ' ';
    }
    os << '\n';
  }
  
  // separator...
  for (size_t trg = 0; trg != target.size(); ++ trg)
    os << '_' << '_';
  os << '\n';
  
  // matrix + source-side
  for (size_t src = 0; src != source.size(); ++ src) {
    const bool bar_horizontal = ! ((src + 1) % 5);
    
    for (size_t trg = 0; trg != target.size(); ++ trg) {
      // blue background: \u001b[44m\u0020\u001b[0m
      if (matrix(src, trg))
	os << char(0x1b) << "[44m" << '*' << char(0x1b) << "[0m";
      else {
	const bool bar_vertical = ! ((trg + 1) % 5);
	
	os << (bar_horizontal && bar_vertical ? '+' : (bar_horizontal ? '-' : (bar_vertical ? ':' : '.')));
      }
      os << ' ';
    }
    
    os << '|' << ' ' << source[src] << '\n';
  }
  
  return os;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("source",             po::value<path_set_type>(&source_files)->multitoken(),             "source file(s)")
    ("target",             po::value<path_set_type>(&target_files)->multitoken(),             "target file(s)")
    ("alignment",          po::value<path_set_type>(&alignment_files)->multitoken(),          "alignment file(s)")
    ("permutation-source", po::value<path_set_type>(&permutation_source_files)->multitoken(), "source side permutation file(s)")
    ("permutation-target", po::value<path_set_type>(&permutation_target_files)->multitoken(), "target side permutation file(s)")
    ("output",             po::value<path_type>(&output_file)->default_value(output_file),    "output file")
    
    ("inverse",   po::bool_switch(&inverse_mode), "inverse alignment")
    ("visualize", po::bool_switch(&visualize_mode), "visualization")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}


