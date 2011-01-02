// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TOKENIZER__PENNTREEBANK__HPP__
#define __CICADA__TOKENIZER__PENNTREEBANK__HPP__ 1

#include <vector>
#include <string>

#include <cicada/tokenizer.hpp>

#include <boost/regex.hpp>

// Penntreebank style tokenizer
/*
#!/bin/sed -f

# Sed script to produce Penn Treebank tokenization on arbitrary raw text.
# Yeah, sure.

# expected input: raw text with ONE SENTENCE TOKEN PER LINE

# by Robert MacIntyre, University of Pennsylvania, late 1995.

# If this wasn't such a trivial program, I'd include all that stuff about
# no warrantee, free use, etc. from the GNU General Public License.  If you
# want to be picky, assume that all of its terms apply.  Okay?

# attempt to get correct directional quotes
s=^"=`` =g
s=\([ ([{<]\)"=\1 `` =g
# close quotes handled at end

s=\.\.\.= ... =g
s=[,;:@#$%&]= & =g

# Assume sentence tokenization has been done first, so split FINAL periods
# only. 
s=\([^.]\)\([.]\)\([])}>"']*\)[ 	]*$=\1 \2\3 =g
# however, we may as well split ALL question marks and exclamation points,
# since they shouldn't have the abbrev.-marker ambiguity problem
s=[?!]= & =g

# parentheses, brackets, etc.
s=[][(){}<>]= & =g
# Some taggers, such as Adwait Ratnaparkhi's MXPOST, use the parsed-file
# version of these symbols.
# UNCOMMENT THE FOLLOWING 6 LINES if you're using MXPOST.
# s/(/-LRB-/g
# s/)/-RRB-/g
# s/\[/-LSB-/g
# s/\]/-RSB-/g
# s/{/-LCB-/g
# s/}/-RCB-/g

s=--= -- =g

# NOTE THAT SPLIT WORDS ARE NOT MARKED.  Obviously this isn't great, since
# you might someday want to know how the words originally fit together --
# but it's too late to make a better system now, given the millions of
# words we've already done "wrong".

# First off, add a space to the beginning and end of each line, to reduce
# necessary number of regexps.
s=$= =
s=^= =

s="= '' =g
# possessive or close-single-quote
s=\([^']\)' =\1 ' =g
# as in it's, I'm, we'd
s='\([sSmMdD]\) = '\1 =g
s='ll = 'll =g
s='re = 're =g
s='ve = 've =g
s=n't = n't =g
s='LL = 'LL =g
s='RE = 'RE =g
s='VE = 'VE =g
s=N'T = N'T =g

s= \([Cc]\)annot = \1an not =g
s= \([Dd]\)'ye = \1' ye =g
s= \([Gg]\)imme = \1im me =g
s= \([Gg]\)onna = \1on na =g
s= \([Gg]\)otta = \1ot ta =g
s= \([Ll]\)emme = \1em me =g
s= \([Mm]\)ore'n = \1ore 'n =g
s= '\([Tt]\)is = '\1 is =g
s= '\([Tt]\)was = '\1 was =g
s= \([Ww]\)anna = \1an na =g
# s= \([Ww]\)haddya = \1ha dd ya =g
# s= \([Ww]\)hatcha = \1ha t cha =g

# clean out extra spaces
s=  *= =g
s=^ *==g
 */

namespace cicada
{
  namespace tokenizer
  {
    class Penntreebank : public cicada::Tokenizer
    {
    private:
      typedef boost::regex regex_type;
      typedef const char*  replace_type;
      typedef std::pair<regex_type, replace_type> pattern_type;

      typedef std::vector<pattern_type, std::allocator<pattern_type> > pattern_set_type;

      
    public:
      Penntreebank()
	: patterns() 
      {
	patterns.push_back(pattern_type(regex_type("^[[:space:]]*\""), " `` "));
	
	patterns.push_back(pattern_type(regex_type("([\\x20\\x28\\x5B\\x7B\\x3C])\""), "$1 `` "));
	
	patterns.push_back(pattern_type(regex_type("([,;:@#$%&]|\\.\\.\\.\\.)"), " $1 "));
	
	patterns.push_back(pattern_type(regex_type("([^.])([.])([\\x5B\\x5D\\x29\\x7D\\x3E\\x22\\x27]*)(?=[[:space:]]*$)"), "$1 $2$3"));
	patterns.push_back(pattern_type(regex_type("([?!\\x5B\\x5D\\x28\\x29\\x7B\\x7D\\x3C\\x3E]|--)"), " $1 "));
	
	patterns.push_back(pattern_type(regex_type("\""), " \'\' "));
	
	patterns.push_back(pattern_type(regex_type("([^\'])\' "), "$1 \' ")); // possessive or close-single-quote
	patterns.push_back(pattern_type(regex_type("\'([sSmMdD]|ll|LL|re|RE|ve|VE) "), " \'$1 ")); // as in it's I'm we'd, 'll 're
	patterns.push_back(pattern_type(regex_type("(n\'t|N\'T) "), " $1 "));
	
	patterns.push_back(pattern_type(regex_type(" ([Cc])annot "), " $1an not"));
	patterns.push_back(pattern_type(regex_type(" ([Dd])\'ye "), " $1\' ye "));
	patterns.push_back(pattern_type(regex_type(" ([Gg])imme "), " $1im me "));
	patterns.push_back(pattern_type(regex_type(" ([Gg])onna "), " $1on na "));
	patterns.push_back(pattern_type(regex_type(" ([Gg])otta "), " $1ot ta "));
	patterns.push_back(pattern_type(regex_type(" ([Ll])emme "), " $1em me "));
	patterns.push_back(pattern_type(regex_type(" ([Mm])ore\'n "), " $1ore \'n "));
	patterns.push_back(pattern_type(regex_type(" \'([Tt])is "), " \'$1 is "));
	patterns.push_back(pattern_type(regex_type(" \'([Tt])was "), " \'$1 was "));
	patterns.push_back(pattern_type(regex_type(" ([Ww])anna "), " $1an na "));
      }
      
    protected:
      virtual void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	tokenized.clear();
	
	if (source.empty()) return;
	
	std::string result = " ";
	sentence_type::const_iterator siter_end = source.end();
	for (sentence_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter) {
	  result += static_cast<const std::string&>(*siter);
	  result += ' ';
	}
	
	std::string source_string;
	for (size_t i = 0; i != patterns.size(); ++ i) {
	  source_string.swap(result);
	  result.clear();
	  
	  boost::regex_replace(std::back_inserter(result),
			       source_string.begin(), source_string.end(),
			       patterns[i].first, patterns[i].second);
	}
	
	tokenized.assign(result);
      }
      
    private:
      pattern_set_type patterns;
    };
  };
};

#endif
