//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include "utils/spinlock.hpp"

#include "wordnet.hpp"
#include "wn.h"

namespace wn
{
  // this is the global guard!
  typedef utils::spinlock              mutex_type;
  typedef utils::spinlock::scoped_lock lock_type;

  static mutex_type __wordnet_mutex;

  struct EscapeIterator : public std::string::const_iterator
  {
    typedef std::string::const_iterator iterator;
    
    EscapeIterator(iterator x) : iterator(x) {}
    
    std::string::value_type operator*() const { return (*static_cast<iterator>(*this) == ' ' ? '_' : *static_cast<iterator>(*this)); }
  };
  
  // todo: implement convertor with '_' into ' '
  
  struct UnescapeIterator
  {
    typedef const char* iter_type;
    typedef std::iterator_traits<iter_type>::iterator_category iterator_category;
    typedef std::iterator_traits<iter_type>::value_type  value_type;
    typedef std::iterator_traits<iter_type>::difference_type difference_type;
    typedef std::iterator_traits<iter_type>::reference reference;
    typedef std::iterator_traits<iter_type>::pointer   pointer;

    UnescapeIterator(iter_type __iter) : iter(__iter) {}
    
    char operator*() const{ return *iter == '_' ? ' ' : *iter; }
    
    UnescapeIterator&
    operator++()
    {
      ++ iter;
      return *this;
    }
    
    UnescapeIterator
    operator++(int)
    { return UnescapeIterator(iter++); }

    UnescapeIterator
    operator+(const difference_type& __n) const
    { return UnescapeIterator(iter + __n); }

    iter_type iter;
  };

  inline
  bool operator!=(const UnescapeIterator& x, const UnescapeIterator& y)
  {
    return x.iter != y.iter;
  }

  inline
  bool operator==(const UnescapeIterator& x, const UnescapeIterator& y)
  {
    return x.iter == y.iter;
  }
  
  inline
  ptrdiff_t operator-(const UnescapeIterator& x, const UnescapeIterator& y)
  {
    return x.iter - y.iter;
  }
  
  template <typename _Ptr>
  struct __wordnet_autoptr
  {
    __wordnet_autoptr(_Ptr __ptr) : ptr(__ptr) {}
    ~__wordnet_autoptr() { free_syns(ptr); }
    
    _Ptr operator->() { return ptr; }
    _Ptr get() { return ptr; }
    
    _Ptr ptr;
  };


  template <typename Ptr, typename Synsets>
  inline
  void __wordnet_synset(Ptr ptr, Synsets& synsets)
  {
    for (Ptr current = ptr; current; current = current->nextss) {
      const std::string pos(UnescapeIterator(current->pos), UnescapeIterator(current->pos + std::strlen(current->pos)));
      
      for (int i = 0; i != current->wcount; ++ i) {
	synsets.resize(synsets.size() + 1);
	
	synsets.back().pos  = pos;
	synsets.back().word = std::string(UnescapeIterator(current->words[i]),
					  UnescapeIterator(current->words[i] + std::strlen(current->words[i])));
	synsets.back().sense = current->wnsns[i];
      }
      
      // we do not deep copy...
      for (int i = 0; i != current->ptrcount; ++ i) 
	if (current->ptrtyp[i] == HYPERPTR) {
	  __wordnet_autoptr<Ptr> curr(read_synset(current->ppos[i], current->ptroff[i], ""));
	  
	  const std::string pos(UnescapeIterator(curr->pos), UnescapeIterator(curr->pos + std::strlen(curr->pos)));
	  
	  for (int j = 0; j < curr->wcount; ++ j) {
	    synsets.resize(synsets.size() + 1);
	    
	    synsets.back().pos  = pos;
	    synsets.back().word = std::string(UnescapeIterator(curr->words[j]),
					      UnescapeIterator(curr->words[j] + std::strlen(curr->words[j])));
	    synsets.back().sense = curr->wnsns[j];
	  }
	}
    }
  }

  void WordNet::operator()(const std::string& word, morph_set_type& morphs) const
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    morphs.clear();
    buffer_type buffer(word.size() + 1, 0);
    std::copy(EscapeIterator(word.begin()), EscapeIterator(word.end()), buffer.begin());
    
    for (int pos = 1; pos <= NUMPARTS; ++ pos) {
      lock_type lock(__wordnet_mutex);
      
      char* morphword = 0;
      if (morphword = morphstr(&(*buffer.begin()), pos))
	do {
	  morphs.push_back(std::string(UnescapeIterator(morphword), UnescapeIterator(morphword + std::strlen(morphword))));
	} while (morphword = morphstr(0, pos));
    }
  }
  
  void WordNet::operator()(const std::string& word, synset_set_type& synsets) const
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    synsets.clear();
    buffer_type buffer(word.size() + 1, 0);
    std::copy(EscapeIterator(word.begin()), EscapeIterator(word.end()), buffer.begin());
    
    for (int pos = 1; pos <= NUMPARTS; ++ pos) {
      lock_type lock(__wordnet_mutex);
      
      __wordnet_autoptr<SynsetPtr> synset_ptr(findtheinfo_ds(&(*buffer.begin()), pos, 0, ALLSENSES));
      __wordnet_synset(synset_ptr.get(), synsets);
      
      char* morphword = 0;
      if (morphword = morphstr(&(*buffer.begin()), pos))
	do {
	  __wordnet_autoptr<SynsetPtr> synset_ptr(findtheinfo_ds(morphword, pos, 0, ALLSENSES));
	  __wordnet_synset(synset_ptr.get(), synsets);
	  
	} while (morphword = morphstr(0, pos));
    }
  }
  
  static boost::once_flag __wordnet_init_once = BOOST_ONCE_INIT;
  static std::string      __wordnet_path;
  
  static void __wordnet_init()
  {
    if (! __wordnet_path.empty()) {
      if (! boost::filesystem::exists(__wordnet_path))
	throw std::runtime_error("no path? " + __wordnet_path);
      
      setenv("WNHOME", __wordnet_path.c_str(), 1);
      
      if (wninit()) {
	setenv("WNSEARCHDIR", __wordnet_path.c_str(), 1);
	
	if (wninit())
	  throw std::runtime_error("no db at: " + __wordnet_path);
      }
    } if (wninit())
	throw std::runtime_error("no wordnet database? check WNHOME, WNSEARCHDIR or supply path with \"path to db\"");
  }

  
  void WordNet::initialize(const std::string& path)
  {
    // we will assure locking + once to make sure this will be called only once...
    lock_type lock(__wordnet_mutex);
    
    __wordnet_path = path;
    
    boost::call_once(__wordnet_init_once, __wordnet_init);
  }
};
