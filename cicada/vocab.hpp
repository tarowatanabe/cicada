// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__VOCAB__HPP__
#define __CICADA__VOCAB__HPP__ 1

#include <string>
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/shared_ptr.hpp>

#include <cicada/symbol.hpp>

#include <succinct_db/succinct_hash.hpp>

#include <utils/atomicop.hpp>
#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/array_power2.hpp>
#include <utils/piece.hpp>

namespace cicada
{  
  // We will share this vocabulary structure with symbol
  
  
  class Vocab
  {
  public:
    typedef Symbol    symbol_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint64_t  hash_value_type;
    
    typedef boost::filesystem::path path_type;

    typedef uint32_t id_type;
    
  public:
    // constants...
    static const symbol_type EMPTY;
    static const symbol_type NONE;
    static const symbol_type EPSILON;
        
    // for language models
    static const symbol_type UNK;
    static const symbol_type STAR;
    static const symbol_type BOS;
    static const symbol_type EOS;

    // constants for non-terminals..
    static const symbol_type GOAL;

    static const symbol_type S;
    static const symbol_type S1;
    static const symbol_type S2;
    
    static const symbol_type X;
    static const symbol_type X1;
    static const symbol_type X2;
    
  private:
    typedef succinctdb::succinct_hash<char, std::allocator<char> >        succinct_hash_type;
    typedef succinctdb::succinct_hash_mapped<char, std::allocator<char> > succinct_hash_mapped_type;
    typedef succinctdb::succinct_hash_stream<char, std::allocator<char> > succinct_hash_stream_type;
    typedef utils::hashmurmur<hash_value_type>                            hasher_type;
    
    struct Cache
    {
      typedef int64_t value_type;
      
      Cache() : value(value_type(-1)) {}
      
      volatile value_type value;
    };
    typedef Cache cache_type;
    typedef std::vector<cache_type, std::allocator<cache_type> > cache_set_type;
    
  public:
    Vocab(size_type size_hint=1024 * 1024 * 4) : __succinct_hash(new succinct_hash_type(size_hint)) {}
    Vocab(const path_type& path, size_type bin_size=0) : __succinct_hash_mapped() { open(path, bin_size); }

    uint64_t size_bytes() const 
    {
      return (__succinct_hash_mapped ? __succinct_hash_mapped->size_bytes() : uint64_t(0));
    }
    uint64_t size_compressed() const
    {
      return (__succinct_hash_mapped ? __succinct_hash_mapped->size_compressed() : uint64_t(0));
    }
    uint64_t size_cache() const
    {
      return ((__succinct_hash_mapped ? __succinct_hash_mapped->size_cache() : uint64_t(0))
	      + __cache_word.size() * sizeof(cache_type)
	      + __cache_id.size() * sizeof(cache_type));
    }
    
    void clear()
    {
      __succinct_hash.reset();
      __succinct_hash_mapped.reset();
      __succinct_hash_stream.reset();
      
      __cache_word.clear();
      __cache_id.clear();
    }
    void close() { clear(); }
    
    // insert a word as a vocabulary
    // we will support dynamic updating even if we write into a file!
    id_type insert(const utils::piece& word);
    
    // dump the vocabularies content, inserted so-forth...
    // we will perform vocabulary merging...
    
    static bool exists(const path_type& path)
    {
      return succinct_hash_mapped_type::exists(path);
    }

    void open(const path_type& path, size_type bin_size = 0)
    {
      clear();

      if (bin_size > 0)
	__succinct_hash_stream.reset(new succinct_hash_stream_type(path, bin_size));
      else {
	__succinct_hash_mapped.reset(new succinct_hash_mapped_type(path));
	
	const size_type cache_size = std::max(size_type(utils::bithack::next_largest_power2(__succinct_hash_mapped->size() >> 5)),
					      size_type(1024 * 16));
	
	__cache_word.reserve(cache_size);
	__cache_word.resize(cache_size, cache_type());
	
	__cache_id.reserve(cache_size);
	__cache_id.resize(cache_size, cache_type());
      }
    }
    void write(const path_type& path) const;
    
    id_type operator[](const symbol_type& word) const
    {
      // we will try twice to extract id for UNK
      const id_type id = find(word);
      return (id == id_type(-1) ? find(UNK) : id);
    }
    
    symbol_type operator[](const id_type& id) const
    {
      return find(id);
    }
    
    bool exists(const utils::piece& word) const
    {
      const hash_value_type hash_value = __hasher(word.begin(), word.end(), 0);
      
      if (__succinct_hash_mapped)
	if (__succinct_hash_mapped->find(word.c_str(), word.size(), hash_value) != succinct_hash_mapped_type::npos())
	  return true;
      
      if (__succinct_hash)
	if (__succinct_hash->find(word.c_str(), word.size(), hash_value) != succinct_hash_type::npos())
	  return true;
      
      return false;
    }
    
    
    id_type find(const symbol_type& word) const
    {
      if (__succinct_hash_mapped) {
	const id_type id = __find(word);
	if (id != succinct_hash_mapped_type::npos())
	  return id_type(id);
      }
      
      if (__succinct_hash) {
	const std::string& word_str = static_cast<const std::string&>(word);
	const hash_value_type hash_value = __hasher(word_str.begin(), word_str.end(), 0);
	
	const id_type id = __succinct_hash->find(word_str.c_str(), word_str.size(), hash_value);
	const size_type offset = (__succinct_hash_mapped ? __succinct_hash_mapped->size() : size_type(0));
	
	return (id != succinct_hash_type::npos() ? id + offset : id_type(-1));
      }
      
      return id_type(-1);
    }
    
    symbol_type find(const id_type& id) const
    {
      if (id == id_type(-1))
	return UNK;
      
      if (__succinct_hash_mapped) {
	if (id < __succinct_hash_mapped->size()) {
	  return __find(id);
	} else if (__succinct_hash && id - __succinct_hash_mapped->size() < __succinct_hash->size()) {
	  succinct_hash_type::const_iterator iter = (__succinct_hash->begin() + id - __succinct_hash_mapped->size());
	  return symbol_type(iter.begin(), iter.end());
	} else
	  return UNK;
      } else if (__succinct_hash && id < __succinct_hash->size()) {
	succinct_hash_type::const_iterator iter = (__succinct_hash->begin() + id);
	return symbol_type(iter.begin(), iter.end());
      } else
	return UNK;
    }

  private:
    id_type __find(const symbol_type& word) const
    {
      cache_set_type& caches = const_cast<cache_set_type&>(__cache_word);
	
      cache_type cache;
      cache_type cache_new;
      
      cache.value = utils::atomicop::fetch_and_add(caches[word.id() & (caches.size() - 1)].value, int64_t(0));
      
      uint32_t __word = (cache.value >> 32) & 0xffffffff;
      uint32_t __id   = (cache.value) & 0xffffffff;
      
      if (__word == word.id())
	return id_type(__id);
      
      const std::string& word_str = static_cast<const std::string&>(word);
      const hash_value_type hash_value = __hasher(word_str.begin(), word_str.end(), 0);
      __id = __succinct_hash_mapped->find(word_str.c_str(), word_str.size(), hash_value);
      
      cache_new.value = (uint64_t(word.id()) << 32) | (uint64_t(__id) & 0xffffffff);
      
      utils::atomicop::compare_and_swap(caches[word.id() & (caches.size() - 1)].value, cache.value, cache_new.value);
      
      return id_type(__id);
    }

    symbol_type __find(const id_type& id) const
    {
      cache_set_type& caches = const_cast<cache_set_type&>(__cache_id);
      
      cache_type cache;
      cache_type cache_new;
      
      cache.value = utils::atomicop::fetch_and_add(caches[id & (caches.size() - 1)].value, int64_t(0));
      
      uint32_t __word = (cache.value >> 32) & 0xffffffff;
      uint32_t __id   = (cache.value) & 0xffffffff;
      
      if (__id == id)
	return symbol_type(__word);
      
      succinct_hash_mapped_type::const_iterator iter = (__succinct_hash_mapped->begin() + id);
      __word = symbol_type(iter.begin(), iter.end()).id();
      
      cache_new.value = (uint64_t(__word) << 32) | (uint64_t(id) & 0xffffffff);
      
      utils::atomicop::compare_and_swap(caches[id & (caches.size() - 1)].value, cache.value, cache_new.value);
      
      return symbol_type(__word);
    }
    
  private:
    boost::shared_ptr<succinct_hash_type>        __succinct_hash;
    boost::shared_ptr<succinct_hash_mapped_type> __succinct_hash_mapped;
    boost::shared_ptr<succinct_hash_stream_type> __succinct_hash_stream;
    hasher_type                                  __hasher;
    
    cache_set_type __cache_word;
    cache_set_type __cache_id;
  };
};


#endif

