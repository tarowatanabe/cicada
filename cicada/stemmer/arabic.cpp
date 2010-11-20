// -*- encoding: utf-8 -*-

#include "stemmer/arabic.hpp"

#include <boost/thread.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/translit.h>
#include <unicode/bytestream.h>
#include <unicode/regex.h>

#include "utils/config.hpp"

/*
#!/usr/bin/perl -w

use strict;

# This is Kareem Darwish's stem_cp1256.pl modified by
# Leah Larkey, Alexander Fraser and Xiaoyi Ma

my $atb_stems = "$ENV{CTK}/lib/atb_stem_tbl.txt";
my %stem;

open S, "<$atb_stems" || die "$0: Cannot open $atb_stems\n";
while (<S>) {
    chomp;
    if (/^(.+)\s+(.+)$/) {
	$stem{$1} = $2;
    }
}
close S;

while (<>) {
    chomp;
    # split on spaces and punctuation
    s/،/\,/g;
    s/؟/\?/g;

    # split on spaces since tokenization was done by atoken.pl
    my @tokens = split ' ', $_;
    for my $token (@tokens) {
      # remove all non-letters (diacritics, punctuation)
      my $newtoken = "";
      while ($token =~ /\G.*?((ء|آ|أ|ؤ|إ|ئ|ا|ب|ة|ت|ث|ج|ح|خ|د|ذ|ر|ز|س|ش|ص|ض|ط|ظ|ع|غ|ف|ق|ك|ل|م|ن|ه|و|ي|ى|[\x21-\x7E])+)/g) {
	$newtoken .= $1;
      }
      $token = $newtoken;

      # normalize ya and Alef Maqsoura
      $token =~ s/ى/ي/g;

      # normalizing different alef-maad, alef-hamza-top,
      # alef-hamza-bottom, bare-alef you can choose between light and
      # aggressive normalization.  The default is aggressive.

      # light normalization
      # $token =~ s/(آ|أ|إ)/ا/g;
      # aggressive normalization
      $token =~ s/(ء|آ|أ|ؤ|إ|ئ)/ا/g;

      if (exists $stem{$token}) {
	  print "$stem{$token} ";
	  next;
      }# else {
#	  print STDERR "$token\n";
#      }


      # this regexp will match every string. It tries to take the longest
      # possible prefix and suffix. $2 will always be defined but can be empty.
      if ($token =~ /^(وال|فال|بال|بت|يت|لت|مت|تت|وت|ست|نت|بم|لم|وم|كم|فم|ال|لل|وي|لي|سي|في|وا|فا|لا|با)(.+)$/) {
	  $token = $2;
      }
      while ($token =~ /^(.+)(ات|وا|تا|ون|وه|ان|تي|ته|تم|كم|هن|هم|ها|ية|تك|نا|ين|يه|ة|ه|ي|ا)$/) {
	  $token = $1;
      }
#    if ($token =~ /^(وال|فال|بال|بت|يت|لت|مت|تت|وت|ست|نت|بم|لم|وم|كم|فم|ال|لل|وي|لي|سي|في|وا|فا|لا|با|)(.+?)(ات|وا|تا|ون|وه|ان|تي|ته|تم|كم|هن|هم|ها|ية|تك|نا|ين|يه|ة|ه|ي|ا)$/) {
#	  print "$2 ";
#      } else {
#	  print "$token ";
#      }
      print "$token ";
  }
    print "\n";
}

## Saved for possible future use
## remove diacritics and kashida
#s/(ً|ٌ|ٍ|َ|ُ|ِ|ّ|ْ|ـ)//g;
*/

namespace cicada
{
  namespace stemmer
  {
    struct ArabicImpl
    {
      struct Replace
      {
      public:
	Replace(const char* pattern, const char* subst) : matcher(0) { initialize(pattern, subst); }
	~Replace() { std::auto_ptr<RegexMatcher> tmp(matcher); }
      
	const UnicodeString& operator()(UnicodeString& uline)
	{
	  if (! matcher)
	    throw std::runtime_error("no matcher...");
	  
	  UErrorCode status = U_ZERO_ERROR;
	  matcher->reset(uline);
	  uline = matcher->replaceAll(substitute, status);
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
	  return uline;
	}

      private:
	void initialize(const char* pattern, const char* subst)
	{
	  UErrorCode status = U_ZERO_ERROR;
	  matcher = new RegexMatcher(UnicodeString::fromUTF8(pattern), 0, status);
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
	  
	  substitute = UnicodeString::fromUTF8(subst);
	}
	
      private:
	RegexMatcher* matcher;
	UnicodeString substitute;
      };
      
      struct ReplaceAll
      {
  
      public:
	ReplaceAll(const char* pattern, const char* subst) : matcher(0) { initialize(pattern, subst); }
	~ReplaceAll() { std::auto_ptr<RegexMatcher> tmp(matcher); }
	
	const UnicodeString& operator()(UnicodeString& uline)
	{
	  if (! matcher)
	    throw std::runtime_error("no matcher...");
	
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
      
      private:
	void initialize(const char* pattern, const char* subst)
	{
	  UErrorCode status = U_ZERO_ERROR;
	  matcher = new RegexMatcher(UnicodeString::fromUTF8(pattern), 0, status);
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RegexMatcher: ") + u_errorName(status));
	  
	  substitute = UnicodeString::fromUTF8(subst);
	}

      private:
	RegexMatcher* matcher;
	UnicodeString substitute;
      };
      
      ArabicImpl()
	: normalize_ya_alef_maqsoura("ى", "ي"),
	  normalize_light("(آ|أ|إ)", "ا"),
	  normalize_aggressive("(ء|آ|أ|ؤ|إ|ئ)", "ا"),
	  remove_prefix("^(وال|فال|بال|بت|يت|لت|مت|تت|وت|ست|نت|بم|لم|وم|كم|فم|ال|لل|وي|لي|سي|في|وا|فا|لا|با)(.+)$", "$2"),
	  remove_suffix("^(.+)(ات|وا|تا|ون|وه|ان|تي|ته|تم|كم|هن|هم|ها|ية|تك|نا|ين|يه|ة|ه|ي|ا)$", "$1"),
	  remove_diacritics_kashida("(ً|ٌ|ٍ|َ|ُ|ِ|ّ|ْ|ـ)", "") {}
      
      
      std::string operator()(const std::string& word)
      {
	UnicodeString uword = UnicodeString::fromUTF8(word);
	
	normalize_ya_alef_maqsoura(uword);
	
	normalize_aggressive(uword);
	
	remove_prefix(uword);
	
	remove_suffix(uword);
	
	std::string word_stemmed;
	StringByteSink<std::string> __sink(&word_stemmed);
	uword.toUTF8(__sink);
	
	return word_stemmed;
      }
      
      
    private:
      Replace    normalize_ya_alef_maqsoura;
      Replace    normalize_light;
      Replace    normalize_aggressive;
      Replace    remove_prefix;
      ReplaceAll remove_suffix;
      Replace    remove_diacritics_kashida;
    };

    Arabic::Arabic() : pimpl(new impl_type()) {}
    Arabic::~Arabic() { std::auto_ptr<impl_type>(pimpl); }
    
    Stemmer::symbol_type Arabic::operator[](const symbol_type& word) const
    {
      if (word == vocab_type::EMPTY || word.is_non_terminal()) return word;
    
      const size_type word_size = word.size();
      
      // SGML-like symbols are not arabiced...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      symbol_set_type& __cache = const_cast<symbol_set_type&>(cache);
      
      if (word.id() >= __cache.size())
	__cache.resize(word.id() + 1, vocab_type::EMPTY);
    
      if (__cache[word.id()] == vocab_type::EMPTY)
	__cache[word.id()] = pimpl->operator()(word);
    
      return __cache[word.id()];
    }

  };
};
