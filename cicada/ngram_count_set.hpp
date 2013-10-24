// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_COUNT_SET__HPP__
#define __CICADA__NGRAM_COUNT_SET__HPP__ 1

#include <iostream>
#include <string>

#include <cicada/symbol.hpp>
#include <cicada/symbol_vector.hpp>

#include <utils/unordered_map.hpp>
#include <utils/piece.hpp>

namespace cicada
{
  
  class NGramCountSet
  {
  public:
    typedef Symbol       word_type;
    typedef SymbolVector ngram_type;
    typedef double       count_type;

  private:
    typedef utils::unordered_map<ngram_type, count_type, boost::hash<ngram_type>, std::equal_to<ngram_type>,
				 std::allocator<std::pair<const ngram_type, count_type> > >::type ngram_set_type;
    
  public:
    typedef ngram_set_type::size_type        size_type;
    typedef ngram_set_type::difference_type  difference_type;
    
    typedef ngram_set_type::value_type  value_type;
    typedef ngram_set_type::key_type    key_type;
    typedef ngram_set_type::mapped_type mapped_type;

    typedef ngram_set_type::const_iterator  const_iterator;
    typedef ngram_set_type::iterator        iterator;

  public:
    NGramCountSet() : ngrams() {}
    NGramCountSet(size_type hint) : ngrams(hint) {}
    NGramCountSet(const utils::piece& x) : ngrams() { assign(x); }

  public:
    size_type size() const { return ngrams.size(); }
    bool empty() const { return ngrams.empty(); }

    void swap(NGramCountSet& x) { ngrams.swap(x.ngrams); }
    
    count_type sum_unigram() const
    {
      count_type sum = 0.0;
      const_iterator iter_end = ngrams.end();
      for (const_iterator iter = ngrams.begin(); iter != iter_end; ++ iter)
	if (iter->first.size() == 1)
	  sum += iter->second;
      return sum;
    }
    
  public:
    count_type& operator[](const ngram_type& ngram) { return ngrams[ngram]; }
    
    void clear() { ngrams.clear(); }

    template <typename Iterator>
    Iterator insert(const value_type& x) { return ngrams.insert(x); }
    template <typename Iterator>
    Iterator insert(Iterator iter, const value_type& x) { return ngrams.insert(iter, x); }
    template <typename Iterator>
    void insert(Iterator first, Iterator last) { ngrams.insert(first, last); }
    
    size_type erase(const key_type& key) { return ngrams.erase(key); }
    template <typename Iterator>
    Iterator erase(Iterator iter) { return ngrams.erase(iter); }
    template <typename Iterator>
    Iterator erase(Iterator first, Iterator last) { return ngrams.erase(first, last); }
    
  public:
    const_iterator begin() const { return ngrams.begin(); }
    iterator begin() { return ngrams.begin(); }

    const_iterator end() const { return ngrams.end(); }
    iterator end() { return ngrams.end(); }

  public:
    void assign(const utils::piece& x);
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    bool assign(utils::piece::const_iterator& iter, utils::piece::const_iterator end);
    
    friend
    std::ostream& operator<<(std::ostream& os, const NGramCountSet& x);
    friend
    std::istream& operator>>(std::istream& is, NGramCountSet& x);

  public:    
    friend
    bool operator==(const NGramCountSet& x, const NGramCountSet& y)
    {
      if (x.size() != y.size()) return false;
      
      NGramCountSet::const_iterator iter_end = y.end();
      for (NGramCountSet::const_iterator iter = y.begin(); iter != iter_end; ++ iter) {
	NGramCountSet::const_iterator fiter = x.ngrams.find(iter->first);
	
	if (fiter == x.end() || *fiter != *iter)
	  return false;
      }
      return true;
    }
    
    friend
    bool operator!=(const NGramCountSet& x, const NGramCountSet& y)
    {
      return !(x == y);
    }
    
private:
    ngram_set_type ngrams;
  };
  
};

namespace std
{
  inline
  void swap(cicada::NGramCountSet& x, cicada::NGramCountSet& y)
  {
    x.swap(y);
  }
};

#endif
