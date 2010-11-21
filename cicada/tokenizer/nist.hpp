// -*- mode: c++ -*-

#ifndef __CICADA__TOKENIZER__NIST__HPP__
#define __CICADA__TOKENIZER__NIST__HPP__ 1

#include <vector>
#include <string>

#include <cicada/tokenizer.hpp>

#include <boost/regex.hpp>

// NIST mteval script style tokenizer...
/*
  $norm_text =~ s/([\{-\~\[-\` -\&\(-\+\:-\@\/])/ $1 /g;   # tokenize punctuation
  $norm_text =~ s/([^0-9])([\.,])/$1 $2 /g; # tokenize period and comma unless preceded by a digit
  $norm_text =~ s/([\.,])([^0-9])/ $1 $2/g; # tokenize period and comma unless followed by a digit
  $norm_text =~ s/([0-9])(-)/$1 $2 /g; # tokenize dash when preceded by a digit
  $norm_text =~ s/\s+/ /g; # one space only between words
  $norm_text =~ s/^\s+//;  # no leading space
  $norm_text =~ s/\s+$//;  # no trailing space
 */

namespace cicada
{
  namespace tokenizer
  {
    class Nist : public cicada::Tokenizer
    {
    private:
      typedef boost::regex regex_type;
      typedef std::string  replace_type;
      typedef std::pair<regex_type, replace_type> pattern_type;

    public:
      Nist()
	: pattern1(regex_type("([\\x7B-\\x7E\\x5B-\\x60\\x21-\\x26\\x28-\\x2B\\x3A-\\x40\\x2F])"), " $1 "),
	  pattern2(regex_type("([^0-9])([\\x2E\\x2C])"), "$1 $2 "),
	  pattern3(regex_type("([\\x2E\\x2C])([^0-9])"), " $1 $2"),
	  pattern4(regex_type("([0-9])(-)"), "$1 $2 ")
      {}
    protected:
      virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	tokenized.clear();
	
	if (source.empty()) return;
	
	std::string tokenized_string;
	sentence_type::const_iterator siter_end = source.end() - 1;
	for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter) {
	  tokenized_string += static_cast<const std::string&>(*siter);
	  tokenized_string += ' ';
	}
	tokenized_string += static_cast<const std::string&>(source.back());

	std::string result;
	boost::regex_replace(std::back_inserter(result), tokenized_string.begin(), tokenized_string.end(), pattern1.first, pattern1.second);
	
	tokenized_string.swap(result);
	result.clear();
	boost::regex_replace(std::back_inserter(result), tokenized_string.begin(), tokenized_string.end(), pattern2.first, pattern2.second);
	
	tokenized_string.swap(result);
	result.clear();
	boost::regex_replace(std::back_inserter(result), tokenized_string.begin(), tokenized_string.end(), pattern3.first, pattern3.second);
	
	tokenized_string.swap(result);
	result.clear();
	boost::regex_replace(std::back_inserter(result), tokenized_string.begin(), tokenized_string.end(), pattern4.first, pattern4.second);
	
	tokenized.assign(result);
      }
      
    private:
      pattern_type pattern1;
      pattern_type pattern2;
      pattern_type pattern3;
      pattern_type pattern4;
    };
  };
};

#endif
