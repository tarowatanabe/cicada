//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// bitext filter

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <boost/iostreams/concepts.hpp>

#include "utils/bithack.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/subprocess.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/utf8.hpp"

#include <unicode/utypes.h>
#include <unicode/unistr.h>
#include <unicode/uchar.h>
#include <unicode/uscript.h>

typedef boost::filesystem::path path_type;

struct color_char_filter : public boost::iostreams::output_filter
{
  color_char_filter() : mode(0) {}
  
  template <typename Sink>
  struct utf8_sink : public icu::ByteSink
  {
    utf8_sink(Sink& __sink) : sink(__sink) {}

    virtual void Append(const char* data, int32_t n) 
    {
      boost::iostreams::write(sink, data, n);
    }
    
    void write(char c)
    {
      boost::iostreams::put(sink, c);
    }

    void put(char c)
    {
      boost::iostreams::put(sink, c);
    }
    
    Sink& sink;
  };

  template<typename Sink>
  bool put(Sink& s, char c)
  { 
    if (mode == 0) {
      mode = utils::utf8_size(c);
      
      if (mode < 1)
	throw std::runtime_error("invalid utf-8 sequence");
      
      switch (mode) {
      case 1: uchar = c; break;
      case 2: uchar = 0x1f & c; break;
      case 3: uchar = 0x0f & c; break;
      case 4: uchar = 0x07 & c; break;
      case 5: uchar = 0x03 & c; break;
      case 6: uchar = 0x01 & c; break;
      default:
	throw std::runtime_error("invalid utf-8 sequence");
      }
    } else
      uchar = (uchar << 6) | (0x3f & c);
    
    -- mode;

    if (! mode) {
      // back into utf-8!
      const UScriptCode script = (UScriptCode) u_getIntPropertyValue(uchar, UCHAR_SCRIPT);
      bool reset = false;
      
      if (script == USCRIPT_KATAKANA) {
	boost::iostreams::put(s, 0x1b);
	boost::iostreams::write(s, "[43m", 4);
	reset = true;
      } else if ((0x3099 <= uchar && uchar <= 0x309c)
		 || uchar == 0x30a0
		 || uchar == 0x30fb
		 || uchar == 0x30fc
		 || uchar == 0xff65
		 || uchar == 0xff70
		 || (0xff9e <= uchar && uchar <= 0xff9f)
		 || uchar == 0x1f201
		 || uchar == 0x1f202
		 || uchar == 0x1f213) {
	// supplemental ranges taken by "grep" of KATAKANA in UnicodeData.txt
	boost::iostreams::put(s, 0x1b);
	boost::iostreams::write(s, "[44m", 4);
	reset = true;
      }
      
      utf8_sink<Sink> sink(s);
      utils::utf8_sink(sink, uchar);
      
      if (reset) {
	boost::iostreams::put(s, 0x1b);
	boost::iostreams::write(s, "[0m", 3);
      }
    }
    
    return true;
  }
  
  UChar32 uchar;
  int     mode;
};

path_type input_file = "-";
path_type output_file = "-";

bool color = false;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  namespace qi = boost::spirit::qi;
  namespace karma = boost::spirit::karma;
  
  try {
    options(argc, argv);

    utils::compress_istream is(input_file, 1024 * 1024);
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    boost::iostreams::filtering_ostream os;
    
    if (color)
      os.push(color_char_filter());
    
    utils::push_compress_ostream(os, output_file, 1024 * 1024 * (! flush_output));
    
    std::string line;
    while (std::getline(is, line))
      os << line << '\n';
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
    ("input",       po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",      po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("color",         po::bool_switch(&color), "colorize output (green for Hiragana, yellow for Katakana, blue for either)")
        
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

