#include <cstdlib>
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
  
  ptrdiff_t operator-(const UnescapeIterator& x, const UnescapeIterator& y)
  {
    return x.iter - y.iter;
  }
  
  void WordNet::operator()(const std::string& word, synset_set_type& synsets)
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    synsets.clear();
    buffer_type buffer(word.size() + 1, 0);
    std::copy(EscapeIterator(word.begin()), EscapeIterator(word.end()), buffer.begin());
    
    for (int pos = 1; pos <= NUMPARTS; ++ pos) {
      SynsetPtr synset_ptr;
      {
	lock_type lock(__wordnet_mutex);
	
	// retrieven all senses, but do not traverse edges...
	synset_ptr = findtheinfo_ds(&(*buffer.begin()), pos, 0, ALLSENSES);
      }
      
      for (SynsetPtr current = synset_ptr; current; current = current->nextss) {
	const std::string pos = current->pos;

	for (int i = 0; i != current->wcount; ++ i) {
	  synsets.resize(synsets.size() + 1);
	  
	  synsets.back().pos = std::string(UnescapeIterator(pos.c_str()), UnescapeIterator(pos.c_str() + pos.size()));
	  synsets.back().word = std::string(UnescapeIterator(current->words[i]), UnescapeIterator(current->words[i] + std::strlen(current->words[i])));
	  synsets.back().sense = current->wnsns[i];
	}
      }
      
      free_syns(synset_ptr);
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
