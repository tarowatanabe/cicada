// -*- encoding: utf-8 -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

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

#include "utils/icu_filter.hpp"
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
#include <unicode/translit.h>
#include <unicode/regex.h>
#include <unicode/bytestream.h>
#include <unicode/schriter.h>

typedef boost::filesystem::path path_type;

struct TransLit
{
  TransLit(const std::string& name, 
	   const std::string& pattern) { initialize(name, pattern); }
  
  TransLit(const std::string& name) { initialize(name); }
  
  void operator()(icu::UnicodeString& data) { trans->transliterate(data); }
  
  void initialize(const std::string& name)
  {
    UErrorCode status = U_ZERO_ERROR;
    trans.reset(icu::Transliterator::createInstance(UnicodeString::fromUTF8(name),
						    UTRANS_FORWARD, status));
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("transliterator::create_instance(): ") + u_errorName(status));
  }
  
  void initialize(const std::string& name, const std::string& pattern)
  {
    UErrorCode status = U_ZERO_ERROR;
    UParseError status_parse;
    trans.reset(icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8(name), icu::UnicodeString::fromUTF8(pattern),
						     UTRANS_FORWARD, status_parse, status));
    
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
  }
  
  boost::shared_ptr<Transliterator> trans;
};

class Replace
{
  
public:
  Replace(const char* pattern, const char* subst) : matcher() { initialize(pattern, subst); }
  Replace(const icu::UnicodeString& pattern, const icu::UnicodeString& subst) : matcher() { initialize(pattern, subst); }
  ~Replace() { clear(); }


  const icu::UnicodeString& operator()(icu::UnicodeString& uline)
  {
    UErrorCode status = U_ZERO_ERROR;
    matcher->reset(uline);
    uline = matcher->replaceAll(substitute, status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
    return uline;
  }

  void clear()
  {
    if (matcher) delete matcher;
    matcher = 0;
  }
  
  void initialize(const char* pattern, const char* subst)
  {
    initialize(icu::UnicodeString(pattern, "utf-8"), icu::UnicodeString(subst, "utf-8"));
  }
  
  
  void initialize(const icu::UnicodeString& pattern, const icu::UnicodeString& subst)
  {
    clear();
    
    UErrorCode status = U_ZERO_ERROR;
    matcher = new icu::RegexMatcher(pattern, 0, status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
    substitute = subst;
  }
  

private:
  icu::RegexMatcher* matcher;
  icu::UnicodeString substitute;
};

class ReplaceAll
{
  
public:
  ReplaceAll(const char* pattern, const char* subst) : matcher() { initialize(pattern, subst); }
  ReplaceAll(const icu::UnicodeString& pattern, const icu::UnicodeString& subst) : matcher() { initialize(pattern, subst); }
  ~ReplaceAll() { clear(); }


  const icu::UnicodeString& operator()(icu::UnicodeString& uline)
  {
    
    while (1) {
      matcher->reset(uline);
      if (! matcher->find()) break;
      
      UErrorCode status = U_ZERO_ERROR;
      uline = matcher->replaceAll(substitute, status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
    }
    return uline;
  }

  void clear()
  {
    if (matcher) delete matcher;
    matcher = 0;
  }
  
  void initialize(const char* pattern, const char* subst)
  {
    initialize(icu::UnicodeString(pattern, "utf-8"), icu::UnicodeString(subst, "utf-8"));
  }
  
  
  void initialize(const icu::UnicodeString& pattern, const icu::UnicodeString& subst)
  {
    clear();
    
    UErrorCode status = U_ZERO_ERROR;
    matcher = new icu::RegexMatcher(pattern, 0, status);
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
    substitute = subst;
  }
  

private:
  icu::RegexMatcher* matcher;
  icu::UnicodeString substitute;
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

path_type input_file = "-";
path_type output_file = "-";

bool remove_control = false;

bool sgml_entity = false;

bool simplified = false;
bool traditional = false;

bool fullwidth = false;
bool halfwidth = false;

bool merge_digits = false;
bool split_digits = false;

bool merge_ideographic = false;
bool split_ideographic = false;

bool merge_hiragana = false;
bool split_hiragana = false;

bool merge_katakana = false;
bool split_katakana = false;

bool merge_symbol = false;
bool split_symbol = false;

bool merge_punctuation = false;
bool split_punctuation = false;

bool merge_mark = false;
bool split_mark = false;

bool normalize_nfkc = false;
bool normalize_nfc = false;
bool normalize_space = false;

bool color = false;

std::string codepage_from = "utf-8";
std::string codepage_to   = "utf-8";

Transliterator* initialize();
void options(int argc, char** argv);

int main(int argc, char** argv)
{
  namespace qi = boost::spirit::qi;
  namespace karma = boost::spirit::karma;
  
  try {
    options(argc, argv);

    if (normalize_nfkc && normalize_nfc)
      throw std::runtime_error("You cannot specify both NFKC/NFC normalization");

    if (simplified && traditional)
      throw std::runtime_error("You cannot specify both Simplified/Traditional Hanzi conversion");
    
    if (fullwidth && halfwidth)
      throw std::runtime_error("You cannot specify both Fullwidth and Halfwidth conversion");

    if (split_digits && merge_digits)
      throw std::runtime_error("You cannot split and merge digits");
    
    if (split_ideographic && merge_ideographic)
      throw std::runtime_error("You cannot split and merge ideographic");

    if (split_hiragana && merge_hiragana)
      throw std::runtime_error("You cannot split and merge hiragana");

    if (split_katakana && merge_katakana)
      throw std::runtime_error("You cannot split and merge katakana");

    if (split_symbol && merge_symbol)
      throw std::runtime_error("You cannot split and merge symbol");

    if (split_punctuation && merge_punctuation)
      throw std::runtime_error("You cannot split and merge punctuation");

    if (split_mark && merge_mark)
      throw std::runtime_error("You cannot split and merge mark");

    std::auto_ptr<Transliterator> trans(initialize());

    boost::iostreams::filtering_istream is;
    is.push(utils::icu_filter(codepage_from, "utf-8", utils::icu_filter::stop));
    utils::push_compress_istream(is, input_file, 1024 * 1024);
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    boost::iostreams::filtering_ostream os;
    os.push(utils::icu_filter("utf-8", codepage_to, utils::icu_filter::stop));
    utils::push_compress_ostream(os, output_file, 1024 * 1024 * (! flush_output));
    ostream_sink sink(os);
    
    std::string line;
    while (std::getline(is, line)) {
      icu::UnicodeString uline = UnicodeString::fromUTF8(line);
      
      trans->transliterate(uline);
      
      // digits...
      if (split_digits) {
	static ReplaceAll replace("(?<=[[:White_Space:]])([[:^Numeric_Type=None:]])(?=[[:^Numeric_Type=None:]]+[[:White_Space:]])", "$1 ");
	
	uline = ' ' + uline + ' ';
	replace(uline);
	uline.trim();
      }

      if (merge_digits) {
	static ReplaceAll replace("(?<=[[:White_Space:]])([[:^Numeric_Type=None:]])[[:White_Space:]]+(?=[[:^Numeric_Type=None:]]+[[:White_Space:]])", "$1");
	  
	uline = ' ' + uline + ' ';
	replace(uline);
	uline.trim();
      }
      
      // ideographic...
      if (split_ideographic) {
	static Replace replace("(?<=[[:Ideographic:]])(?=[[:Ideographic:]])", " ");
	replace(uline);
      }
      
      if (merge_ideographic) {
	static Replace replace("(?<=[[:Ideographic:]])[[:White_Space:]]+(?=[[:Ideographic:]])", "");
	replace(uline);
      }

      // hiragana...
      if (split_hiragana) {
	static Replace replace("(?<=[[:Hiragana:][\\u3099-\\u309C][\\u30A0][\\u30FC][\\uFF70]])(?=[[:Hiragana:][\\u3099-\\u309C][\\u30A0][\\u30FC][\\uFF70]])", " ");
	replace(uline);
      }
      
      if (merge_hiragana) {
	static Replace replace("(?<=[[:Hiragana:][\\u3099-\\u309C][\\u30A0][\\u30FC][\\uFF70]])[[:White_Space:]]+(?=[[:Hiragana:][\\u3099-\\u309C][\\u30A0][\\u30FC][\\uFF70]])", "");
	replace(uline);
      }
      
      
      // katakana...
      if (split_katakana) {
	static Replace replace("(?<=[[:Katakana:][\\u3099-\\u309C][\\u30A0][\\u30FB-\\u30FC][\\uFF65][\\uFF70][\\uFF9e-\\uFF9F][\\U0001F201-\\U0001F202][\\U0001F213]])(?=[[:Katakana:][\\u3099-\\u309C][\\u30A0][\\u30FB-\\u30FC][\\uFF65][\\uFF70][\\uFF9e-\\uFF9F][\\U0001F201-\\U0001F202][\\U0001F213]])", " ");
	replace(uline);
      }
      
      if (merge_katakana) {
	static Replace replace("(?<=[[:Katakana:][\\u3099-\\u309C][\\u30A0][\\u30FB-\\u30FC][\\uFF65][\\uFF70][\\uFF9e-\\uFF9F][\\U0001F201-\\U0001F202][\\U0001F213]])[[:White_Space:]]+(?=[[:Katakana:][\\u3099-\\u309C][\\u30A0][\\u30FB-\\u30FC][\\uFF65][\\uFF70][\\uFF9e-\\uFF9F][\\U0001F201-\\U0001F202][\\U0001F213]])", "");
	replace(uline);
      }

      // symbol
      if (split_symbol) {
	static Replace replace("(?<=[[:Symbol:]])(?=[[:Symbol:]])", " ");
	replace(uline);
      }
      
      if (merge_symbol) {
	static Replace replace("(?<=[[:Symbol:]])[[:White_Space:]]+(?=[[:Symbol:]])", "");
	replace(uline);
      }

      // punctuation
      if (split_punctuation) {
	static Replace replace("(?<=[[:Punctuation:]])(?=[[:Punctuation:]])", " ");
	replace(uline);
      }
      
      if (merge_punctuation) {
	static Replace replace("(?<=[[:Punctuation:]])[[:White_Space:]]+(?=[[:Punctuation:]])", "");
	replace(uline);
      }
      
      // mark
      if (split_mark) {
	static Replace replace("(?<=[[:Mark:]])(?=[[:Mark:]])", " ");
	replace(uline);
      }
      
      if (merge_mark) {
	static Replace replace("(?<=[[:Mark:]])[[:White_Space:]]+(?=[[:Mark:]])", "");
	replace(uline);
      }

      
      // space normalization
      if (normalize_space) {
	static TransLit trans("NormalizeSpace", "[:White_Space:]+ > ' ';");
	trans(uline);
	uline.trim();
      }
      
      // color...
      if (color) {
	icu::UnicodeString uline_color;
	  
	icu::StringCharacterIterator iter(uline);
	  
	for (iter.setToStart(); iter.hasNext(); /**/) {
	  const UChar32 uchar = iter.next32PostInc();
	    
	  const UScriptCode   script        = (UScriptCode)   u_getIntPropertyValue(uchar, UCHAR_SCRIPT);
	  const UCharCategory category_mask = (UCharCategory) u_getIntPropertyValue(uchar, UCHAR_GENERAL_CATEGORY_MASK);
	  const UNumericType  numeric       = (UNumericType)  u_getIntPropertyValue(uchar, UCHAR_NUMERIC_TYPE);
	  bool reset = false;
	    
	  if (script == USCRIPT_KATAKANA) {
	    uline_color.append(0x1b);
	    uline_color.append("[43m");
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
	    uline_color.append(0x1b);
	    uline_color.append("[42m");
	    reset = true;
	  } else if (category_mask & (U_GC_P_MASK | U_GC_S_MASK | U_GC_M_MASK)) {
	    uline_color.append(0x1b);
	    uline_color.append("[41m");
	    reset = true;
	  } else if (numeric != U_NT_NONE) {
	    uline_color.append(0x1b);
	    uline_color.append("[45m");
	    reset = true;
	  }
	    
	  uline_color.append(uchar);
	    
	  if (reset) {
	    uline_color.append(0x1b);
	    uline_color.append("[0m");
	  }
	}
	  
	uline = uline_color;
      }
	
      uline.toUTF8(sink);
      sink.write('\n');
    }
  }
  catch(std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return -1;
  }
  return 0;
}

Transliterator* initialize()
{
  icu::UnicodeString rules;
  
  if (remove_control) {
    UErrorCode status = U_ZERO_ERROR;
    UParseError status_parse;
    std::auto_ptr<icu::Transliterator> trans(icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8("RemoveControls"),
										  icu::UnicodeString::fromUTF8("[[:C:]-[:White_Space:]] > ;"),
										  UTRANS_FORWARD, status_parse, status));
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
    
    icu::Transliterator::registerInstance(trans.release());
    
    rules += icu::UnicodeString::fromUTF8(":: RemoveControls ;\n");
  }
  
    
  if (sgml_entity) {
    static const char* table_sgml2entity[] = {
#include "../utils/sgml_entity_table.hpp"
    };
    
    const int sgml_table_size = sizeof(table_sgml2entity) / sizeof(char*);

    icu::UnicodeString rules_entity;
    for (int i = 0; i < sgml_table_size; ++ i) {
      rules_entity += icu::UnicodeString::fromUTF8(table_sgml2entity[i]);
      rules_entity += '\n';
    }
    
    UErrorCode status = U_ZERO_ERROR;
    UParseError status_parse;
    std::auto_ptr<icu::Transliterator> trans(icu::Transliterator::createFromRules(icu::UnicodeString::fromUTF8("SGMLEntities"),
										  rules_entity,
										  UTRANS_FORWARD, status_parse, status));
    if (U_FAILURE(status))
      throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
    
    icu::Transliterator::registerInstance(trans.release());
    
    rules += icu::UnicodeString::fromUTF8(":: SGMLEntities ;\n");
    rules += icu::UnicodeString::fromUTF8(":: Hex-Any;\n");
  }
  
  if (fullwidth)
    rules += icu::UnicodeString::fromUTF8(":: Halfwidth-Fullwidth ;\n");
  if (halfwidth)
    rules += icu::UnicodeString::fromUTF8(":: Fullwidth-Halfwidth ;\n");

  if (simplified)
    rules += icu::UnicodeString::fromUTF8(":: Traditional-Simplified; \n");
  if (traditional)
    rules += icu::UnicodeString::fromUTF8(":: Simplified-Traditional; \n");
  
  if (normalize_nfc)
    rules += icu::UnicodeString::fromUTF8(":: NFC; \n");
  if (normalize_nfkc)
    rules += icu::UnicodeString::fromUTF8(":: NFKC; \n");
  
  UErrorCode status = U_ZERO_ERROR;
  UParseError status_parse;
  std::auto_ptr<icu::Transliterator> trans(icu::Transliterator::createFromRules("FullNormalizer", rules,
										UTRANS_FORWARD, status_parse, status));
  
  if (U_FAILURE(status))
    throw std::runtime_error(std::string("transliterator::create_from_rules(): ") + u_errorName(status));
  
  return trans.release();
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",       po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",      po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("remove-control", po::bool_switch(&remove_control), "remove non white-space controls")
    
    ("entity", po::bool_switch(&sgml_entity), "convert SGML entities")
    
    ("fullwidth", po::bool_switch(&fullwidth),  "halfwidth to fullwidth conversion")
    ("halfwidth", po::bool_switch(&halfwidth),  "fullwidth to halfwidth conversion")
    
    ("simplified",  po::bool_switch(&simplified),  "traditional to simplified Hanzi conversion")
    ("traditional", po::bool_switch(&traditional), "simplied to traditional Hanzi conversion")

    ("merge-digits", po::bool_switch(&merge_digits), "merge digits")
    ("split-digits", po::bool_switch(&split_digits), "split digits")
    
    ("merge-ideographic", po::bool_switch(&merge_ideographic), "merge ideographic")
    ("split-ideographic", po::bool_switch(&split_ideographic), "split ideographic")

    ("merge-hiragana", po::bool_switch(&merge_hiragana), "merge hiragana")
    ("split-hiragana", po::bool_switch(&split_hiragana), "split hiragana")

    ("merge-katakana", po::bool_switch(&merge_katakana), "merge katakana")
    ("split-katakana", po::bool_switch(&split_katakana), "split katakana")

    ("merge-symbol", po::bool_switch(&merge_symbol), "merge symbol")
    ("split-symbol", po::bool_switch(&split_symbol), "split symbol")

    ("merge-punctuation", po::bool_switch(&merge_punctuation), "merge punctuation")
    ("split-punctuation", po::bool_switch(&split_punctuation), "split punctuation")

    ("merge-mark", po::bool_switch(&merge_mark), "merge mark")
    ("split-mark", po::bool_switch(&split_mark), "split mark")
    
    ("normalize-nfkc",  po::bool_switch(&normalize_nfkc),  "normalize NFKC")
    ("normalize-nfc",   po::bool_switch(&normalize_nfc),   "normalize NFC")
    ("normalize-space", po::bool_switch(&normalize_space), "normalize space")
    
    ("color",         po::bool_switch(&color), "colorize output (yellow for Katakana, green for maybe-symbol Katakana, red for punctuation/symbol/mark)")
    
    ("codepage-from", po::value<std::string>(&codepage_from)->default_value(codepage_from), "input codepage")
    ("codepage-to",   po::value<std::string>(&codepage_to)->default_value(codepage_to),     "output codepage")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

