// -*- encoding: utf-8 -*-

#include <sstream>
#include <pthread.h>

#include "symbol.hpp"
#include "vocab.hpp"

#include "utils/filesystem.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/atomicop.hpp"
#include "utils/spinlock.hpp"
#include "utils/config.hpp"

#include <boost/thread.hpp>

namespace cicada
{
  // constants...
  const Vocab::symbol_type Vocab::EMPTY   = Vocab::symbol_type("");
  const Vocab::symbol_type Vocab::NONE    = Vocab::symbol_type("<none>");
  const Vocab::symbol_type Vocab::EPSILON = Vocab::symbol_type("<epsilon>");
  
  const Vocab::symbol_type Vocab::UNK     = Vocab::symbol_type("<unk>");
  const Vocab::symbol_type Vocab::STAR    = Vocab::symbol_type("<star>");
  const Vocab::symbol_type Vocab::BOS     = Vocab::symbol_type("<s>");
  const Vocab::symbol_type Vocab::EOS     = Vocab::symbol_type("</s>");

  // X_{n} indicated subscript...
  const Vocab::symbol_type Vocab::GOAL  = Vocab::symbol_type("[goal]");

  const Vocab::symbol_type Vocab::S  = Vocab::symbol_type("[s]");
  const Vocab::symbol_type Vocab::S1 = Vocab::symbol_type("[s,1]");
  const Vocab::symbol_type Vocab::S2 = Vocab::symbol_type("[s,2]");

  const Vocab::symbol_type Vocab::X  = Vocab::symbol_type("[x]");
  const Vocab::symbol_type Vocab::X1 = Vocab::symbol_type("[x,1]");
  const Vocab::symbol_type Vocab::X2 = Vocab::symbol_type("[x,2]");
  

  // normal vocab services...
  
  Symbol::id_type Vocab::insert(const std::string& word)
  {
    if (__succinct_hash_stream)
      return __succinct_hash_stream->insert(word.c_str(), word.size(), __hasher(word.begin(), word.end(), 0));

    if (! __succinct_hash)
      __succinct_hash.reset(new succinct_hash_type(1024 * 1024));

    if (__succinct_hash_mapped) {
      const hash_value_type hash_value = __hasher(word.begin(), word.end(), 0);
      
      symbol_type::id_type word_id_mapped = __succinct_hash_mapped->find(word.c_str(), word.size(), hash_value);
      if (word_id_mapped != succinct_hash_mapped_type::npos())
	return word_id_mapped;
      else
	return __succinct_hash->insert(word.c_str(), word.size(), hash_value) + __succinct_hash_mapped->size();
      
    } else
      return __succinct_hash->insert(word.c_str(), word.size(), __hasher(word.begin(), word.end(), 0));
  }
  
  void Vocab::write(const path_type& path) const
  {
    //if both of dynamic/static hash are open,

    if (__succinct_hash_mapped) {
      
      if (__succinct_hash && ! __succinct_hash->empty()) {
	succinct_hash_stream_type succinct_hash(path, (__succinct_hash_mapped->size() + __succinct_hash->size()) >> 1);
	
	{
	  // insert data from mapped file
	  succinct_hash_mapped_type::const_iterator iter_end = __succinct_hash_mapped->end();
	  for (succinct_hash_mapped_type::const_iterator iter = __succinct_hash_mapped->begin(); iter != iter_end; ++ iter) {
	    const std::string word(iter.begin(), iter.end());
	    succinct_hash.insert(word.c_str(), word.size(), __hasher(word.begin(), word.end(), 0));
	  }
	}
	
	{
	  // insert data from raw storage...
	  succinct_hash_type::const_iterator iter_end = __succinct_hash->end();
	  for (succinct_hash_type::const_iterator iter = __succinct_hash->begin(); iter != iter_end; ++ iter)
	    succinct_hash.insert(&(*iter.begin()), iter.size(), __hasher(iter.begin(), iter.end(), 0));
	}
	
	// finally, dump!
	succinct_hash.close();
	
      } else
	__succinct_hash_mapped->write(path);
      
    } else if (__succinct_hash) {
      // we have only dynamic db... dump!
      __succinct_hash->write(path);
    } 
  }
  
};
