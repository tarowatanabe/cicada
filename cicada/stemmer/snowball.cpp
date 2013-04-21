//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <vector>
#include <stdexcept>

#include "stemmer/snowball.hpp"

#include "libstemmer_c/include/libstemmer.h"

#include "utils/piece.hpp"

namespace cicada
{
  namespace stemmer
  {
    struct SnowballImpl
    {
      typedef sb_stemmer stemmer_type;
      typedef std::vector<sb_symbol, std::allocator<sb_symbol> > buffer_type;

      SnowballImpl(const std::string& lang) : pimpl(sb_stemmer_new(lang.c_str(), 0)) {}
      ~SnowballImpl() 
      {
	if (pimpl)
	  sb_stemmer_delete(pimpl);
      }
      
      std::string operator()(const utils::piece& word) const
      {
	if (! pimpl) return word;
	
	buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	
	buffer.clear();
	buffer.insert(buffer.end(), word.begin(), word.end());
	
	const sb_symbol* stemmed = sb_stemmer_stem(pimpl, static_cast<const sb_symbol*>(&(*buffer.begin())), buffer.size());
	const size_t length = sb_stemmer_length(pimpl);
	
	return std::string(stemmed, stemmed + length);
      }
      
      stemmer_type* pimpl;
      buffer_type buffer_impl;
    };

    Snowball::Snowball(const std::string& language)
      : pimpl(new impl_type(language))
    {
      if (! pimpl->pimpl) {
	std::string message = "We do not support the stemming algorithm: " + language + ". Supported languages are:";
	
	const char **names = sb_stemmer_list();
	for (const char** first = names; *first; ++ first)
	  message += std::string(first != names ? ", " : " ") + *first;
	
	throw std::runtime_error(message);
      }
    }
    
    Snowball::~Snowball() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    std::string Snowball::operator()(const utils::piece& word) const
    {
      if (word.empty()) return word;
      
      const size_type word_size = word.size();
    
      // SGML-like symbols are not stemmed...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;
      
      return pimpl->operator()(word);
    }

  };
};
