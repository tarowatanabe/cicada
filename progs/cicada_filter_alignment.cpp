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
path_set_type alignment2_files;
path_set_type permutation_source_files;
path_set_type permutation_target_files;
path_type output_file = "-";

bool inverse_mode = false;
bool visualize_mode = false;
bool giza_mode = false;

void options(int argc, char** argv);

std::ostream& visualize(std::ostream& os,
			const sentence_type& source,
			const sentence_type& target,
			const alignment_type& alignment);
std::ostream& visualize(std::ostream& os,
			const sentence_type& source,
			const sentence_type& target,
			const alignment_type& alignment1,
			const alignment_type& alignment2);

std::ostream& giza(std::ostream& os,
		   const sentence_type& source,
		   const sentence_type& target,
		   const alignment_type& alignment);

inline
void permute(alignment_type& alignment, const permutation_type& permutation_source, const permutation_type& permutation_target)
{
  alignment_type::iterator aiter_end = alignment.end();
  for (alignment_type::iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
    if (! permutation_source.empty()) {
      if (aiter->source >= static_cast<int>(permutation_source.size()))
	throw std::runtime_error("invalid source permutation");
      
      aiter->source = permutation_source[aiter->source];
    }
    
    if (! permutation_target.empty()) {
      if (aiter->target >= static_cast<int>(permutation_target.size()))
	throw std::runtime_error("invalid target permutation");
      
      aiter->target = permutation_target[aiter->target];
    }
  }
}

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
    
    if (int(visualize_mode) + giza_mode > 1)
      throw std::runtime_error("either visualize or giza");
    
    if (visualize_mode || giza_mode) {
      if (source_files.empty())
	throw std::runtime_error("no source data?");
      if (target_files.empty())
	throw std::runtime_error("no target data?");
    }
    
    if (! alignment2_files.empty()) {
      if (alignment_files.size() != alignment2_files.size())
	throw std::runtime_error("# of secondary alignment file does not match");
      
      if (! visualize_mode)
	throw std::runtime_error("secondary alignemnt file w/o visualization?");
    }
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
    if (! permutation_source_files.empty() || ! permutation_target_files.empty()) {
      sentence_type    sentence_source;
      sentence_type    sentence_target;
      sentence_type    permuted_source;
      sentence_type    permuted_target;      
      permutation_type permutation_source;
      permutation_type permutation_target;
      alignment_type   alignment;
      alignment_type   alignment2;
      std::vector<bool, std::allocator<bool> > assigned;

      const bool has_source = ! source_files.empty();
      const bool has_target = ! target_files.empty();
      
      const bool has_permutation_source = ! permutation_source_files.empty();
      const bool has_permutation_target = ! permutation_target_files.empty();

      const bool has_alignment2 = ! alignment2_files.empty();
      
      for (size_t i = 0; i != alignment_files.size(); ++ i) {
	std::unique_ptr<std::istream> is_source(has_source ? new utils::compress_istream(source_files[i], 1024 * 1024) : 0);
	std::unique_ptr<std::istream> is_target(has_target ? new utils::compress_istream(target_files[i], 1024 * 1024) : 0);
	
	std::unique_ptr<std::istream> ps_source(has_permutation_source ? new utils::compress_istream(permutation_source_files[i], 1024 * 1024) : 0);
	std::unique_ptr<std::istream> ps_target(has_permutation_target ? new utils::compress_istream(permutation_target_files[i], 1024 * 1024) : 0);
	
	utils::compress_istream is(alignment_files[i], 1024 * 1024);
	std::unique_ptr<std::istream> is2(has_alignment2 ? new utils::compress_istream(alignment2_files[i], 1024 * 1024) : 0);
	
	for (;;) {
	  is >> alignment;
	  if (has_alignment2)
	    *is2 >> alignment2;

	  if (has_source)
	    *is_source >> sentence_source;
	  if (has_target)
	    *is_target >> sentence_target;
	  
	  if (has_permutation_source)
	    *ps_source >> permutation_source;
	  if (has_permutation_target)
	    *ps_target >> permutation_target;
	  
	  if (! is
	      || (is2.get() && ! *is2)
	      || (ps_source.get() && ! *ps_source)
	      || (ps_target.get() && ! *ps_target)
	      || (is_source.get() && ! *is_source)
	      || (is_target.get() && ! *is_target)) break;

	  permute(alignment, permutation_source, permutation_target);
	  permute(alignment2, permutation_source, permutation_target);
	  
	  if (has_permutation_source && has_source) {
	    permuted_source.resize(sentence_source.size());
	    
	    assigned.clear();
	    assigned.resize(sentence_source.size(), false);
	    
	    for (size_t pos = 0; pos != sentence_source.size(); ++ pos) {
	      if (permutation_source[pos] >= static_cast<int>(sentence_source.size()))
		throw std::runtime_error("invalid permutation: out of range");
	      
	      if (assigned[permutation_source[pos]])
		throw std::runtime_error("invalid permutation: duplicates");
	      
	      assigned[permutation_source[pos]] = true;
	      
	      permuted_source[pos] = sentence_source[permutation_source[pos]];
	    }
	    
	    sentence_source.swap(permuted_source);
	  }
	  
	  if (has_permutation_target && has_target) {
	    permuted_target.resize(sentence_target.size());
	    
	    assigned.clear();
	    assigned.resize(sentence_target.size(), false);
	    
	    for (size_t pos = 0; pos != sentence_target.size(); ++ pos) {
	      if (permutation_target[pos] >= static_cast<int>(sentence_target.size()))
		throw std::runtime_error("invalid permutation: out of range");
	      
	      if (assigned[permutation_target[pos]])
		throw std::runtime_error("invalid permutation: duplicates");
	      
	      assigned[permutation_target[pos]] = true;
	      
	      permuted_target[pos] = sentence_target[permutation_target[pos]];
	    }
	    
	    sentence_target.swap(permuted_target);
	  }
	  
	  if (inverse_mode) {
	    alignment.inverse();
	    alignment2.inverse();
	  }
	  
	  std::sort(alignment.begin(), alignment.end());
	  std::sort(alignment2.begin(), alignment2.end());
	  
	  if (visualize_mode) {
	    if (has_alignment2)
	      visualize(os, sentence_source, sentence_target, alignment, alignment2);
	    else
	      visualize(os, sentence_source, sentence_target, alignment);
	  } else if (giza_mode)
	    giza(os, sentence_source, sentence_target, alignment);
	  else
	    os << alignment << '\n';
	}
	
	if (is
	    || (is2.get() && *is2)
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
      alignment_type alignment2;

      const bool has_source = ! source_files.empty();
      const bool has_target = ! target_files.empty();

      const bool has_alignment2 = ! alignment2_files.empty();
      
      for (size_t i = 0; i != alignment_files.size(); ++ i) {
	std::unique_ptr<std::istream> is_source(has_source ? new utils::compress_istream(source_files[i], 1024 * 1024) : 0);
	std::unique_ptr<std::istream> is_target(has_target ? new utils::compress_istream(target_files[i], 1024 * 1024) : 0);
	
	utils::compress_istream is(alignment_files[i], 1024 * 1024);
	std::unique_ptr<std::istream> is2(has_alignment2 ? new utils::compress_istream(alignment2_files[i], 1024 * 1024) : 0);
	
	for (;;) {
	  is >> alignment;
	  if (has_alignment2)
	    *is2 >> alignment2;
	  
	  if (has_source)
	    *is_source >> sentence_source;
	  if (has_target)
	    *is_target >> sentence_target;
	  
	  if (! is
	      || (is2.get() && ! *is2)
	      || (is_source.get() && ! *is_source)
	      || (is_target.get() && ! *is_target)) break;
	  
	  if (inverse_mode) {
	    alignment.inverse();
	    alignment2.inverse();
	  }

	  std::sort(alignment.begin(), alignment.end());
	  std::sort(alignment2.begin(), alignment2.end());
	  
	  if (visualize_mode) {
	    if (has_alignment2)
	      visualize(os, sentence_source, sentence_target, alignment, alignment2);
	    else
	      visualize(os, sentence_source, sentence_target, alignment);  
	  } else if (giza_mode)
	    giza(os, sentence_source, sentence_target, alignment);
	  else
	    os << alignment << '\n';
	}
	
	if (is
	    || (is2.get() && *is2)
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

struct ostream_sink : public icu::ByteSink
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

std::ostream& giza(std::ostream& os,
		   const sentence_type& source,
		   const sentence_type& target,
		   const alignment_type& alignment)
{  
  os << "# Sentence pair (0)"
     << " source length " << source.size() 
     << " target length " << target.size() 
     << " alignment score : " << 0 << '\n';
  os << target << '\n';
  
  if (alignment.empty() || source.empty() || target.empty()) {
    os << "NULL ({ })";
    sentence_type::const_iterator siter_end = source.end();
    for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
      os << ' ' << *siter << " ({ })";
    os << '\n';
  } else {
    typedef std::vector<int, std::allocator<int> > index_set_type;
    typedef std::vector<index_set_type, std::allocator<index_set_type> > align_set_type;
    typedef std::vector<bool, std::allocator<bool> > align_none_type;
    
    align_set_type  aligns(source.size());
    align_none_type aligns_none(target.size(), true);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      aligns[aiter->source].push_back(aiter->target + 1);
      aligns_none[aiter->target] = false;
    }
    
    os << "NULL";
    os << " ({ ";
    for (size_t j = 0; j != target.size(); ++ j)
      if (aligns_none[j])
	os << (j + 1) << ' ';
    os << "})";
    
    for (size_t src = 0; src != source.size(); ++ src) {
      os << ' ' << source[src];
      os << " ({ ";
      std::copy(aligns[src].begin(), aligns[src].end(), std::ostream_iterator<int>(os, " "));
      os << "})";
    }
    os << '\n';
  }
  
  return os;
}


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

std::ostream& visualize(std::ostream& os,
			const sentence_type& source,
			const sentence_type& target,
			const alignment_type& alignment1,
			const alignment_type& alignment2)
{
  typedef utils::vector2<int, std::allocator<int> > matrix_type;
  typedef CharacterString wide_string_type;
  typedef std::vector<wide_string_type, std::allocator<wide_string_type> > wide_string_set_type;
  
  matrix_type matrix(source.size(), target.size(), 0);

  {
    alignment_type::const_iterator aiter_end = alignment1.end();
    for (alignment_type::const_iterator aiter = alignment1.begin(); aiter != aiter_end; ++ aiter) {
      if (aiter->source >= static_cast<int>(source.size()))
	throw std::runtime_error("invalid alignment");
      if (aiter->target >= static_cast<int>(target.size()))
	throw std::runtime_error("invalid alignment");
      
      matrix(aiter->source, aiter->target) |= 1;
    }
  }
  
  {
    alignment_type::const_iterator aiter_end = alignment2.end();
    for (alignment_type::const_iterator aiter = alignment2.begin(); aiter != aiter_end; ++ aiter) {
      if (aiter->source >= static_cast<int>(source.size()))
	throw std::runtime_error("invalid alignment");
      if (aiter->target >= static_cast<int>(target.size()))
	throw std::runtime_error("invalid alignment");
      
      matrix(aiter->source, aiter->target) |= 2;
    }
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
      // intersection: blue background:   \u001b[44m...\u001b[0m
      // align1:       green background:  \u001b[42m...\u001b[0m
      // align2:       yellow background: \u001b[43m...\u001b[0m
      
      if (matrix(src, trg)) {
	if (matrix(src, trg) == 3)
	  os << char(0x1b) << "[44m" << '*' << char(0x1b) << "[0m";
	else if (matrix(src, trg) == 1)
	  os << char(0x1b) << "[42m" << '*' << char(0x1b) << "[0m";
	else if (matrix(src, trg) == 2)
	  os << char(0x1b) << "[43m" << '*' << char(0x1b) << "[0m";
	else
	  throw std::runtime_error("invalid matrix!");
      } else {
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
    ("alignment2",         po::value<path_set_type>(&alignment2_files)->multitoken(),         "secondary alignment file(s)")
    ("permutation-source", po::value<path_set_type>(&permutation_source_files)->multitoken(), "source side permutation file(s)")
    ("permutation-target", po::value<path_set_type>(&permutation_target_files)->multitoken(), "target side permutation file(s)")
    ("output",             po::value<path_type>(&output_file)->default_value(output_file),    "output file")
    
    ("inverse",   po::bool_switch(&inverse_mode),   "inverse alignment")
    ("visualize", po::bool_switch(&visualize_mode), "visualization")
    ("giza",      po::bool_switch(&giza_mode),      "giza-style")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}


