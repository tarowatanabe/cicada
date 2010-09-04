// -*- mode: c++ -*-

#ifndef __CICADA__SENTENCE_VECTOR__HPP__
#define __CICADA__SENTENCE_VECTOR__HPP__ 1

#include <vector>
#include <iostream>

#include <cicada/sentence.hpp>

namespace cicada
{
  class SentenceVector
  {
  public:
    typedef Sentence sentence_type;
    typedef sentence_type::symbol_type symbol_type;
    typedef sentence_type::word_type   word_type;

  private:
    typedef std::vector<sentence_type, std::allocator<sentence_type> > sent_set_type;
    
  public:
    typedef sent_set_type::size_type              size_type;
    typedef sent_set_type::difference_type        difference_type;
    typedef sent_set_type::value_type             value_type;
    
    typedef sent_set_type::iterator               iterator;
    typedef sent_set_type::const_iterator         const_iterator;
    typedef sent_set_type::reverse_iterator       reverse_iterator;
    typedef sent_set_type::const_reverse_iterator const_reverse_iterator;
    typedef sent_set_type::reference              reference;
    typedef sent_set_type::const_reference        const_reference;
    
  public:
    // constructor etc...
    SentenceVector() {}
    SentenceVector(const std::string& x) { assign(x); }
    SentenceVector(const SentenceVector& x) : __sents(x.__sents) {}

  public:
    void assign(const std::string& x);
    bool assign(std::string::const_iterator& iter, std::string::const_iterator end);
    
  public:
    void push_back(const sentence_type& sentence) { __sents.push_back(sentence); }
    void pop_back() { __sents.pop_back(); }
    
    void swap(SentenceVector& __x) { __sents.swap(__x.__sents); }
    
    void clear() { __sents.clear(); }
    void reserve(size_type __n) { __sents.reserve(__n); }
    void resize(size_type __n) { __sents.resize(__n); }
    
    size_type size() const { return __sents.size(); }
    bool empty() const { return __sents.empty(); }

  public:
    inline const_iterator begin() const { return __sents.begin(); }
    inline       iterator begin()       { return __sents.begin(); }
    
    inline const_iterator end() const { return __sents.end(); }
    inline       iterator end()       { return __sents.end(); }
    
    inline const_reverse_iterator rbegin() const { return __sents.rbegin(); }
    inline       reverse_iterator rbegin()       { return __sents.rbegin(); }
    
    inline const_reverse_iterator rend() const { return __sents.rend(); }
    inline       reverse_iterator rend()       { return __sents.rend(); }
    
    inline const_reference operator[](size_type pos) const { return __sents[pos]; }
    inline       reference operator[](size_type pos)       { return __sents[pos]; }

    inline const_reference front() const { return __sents.front(); }
    inline       reference front()       { return __sents.front(); }
    
    inline const_reference back() const { return __sents.back(); }
    inline       reference back()       { return __sents.back(); }

  public:
    friend
    std::ostream& operator<<(std::ostream& os, const SentenceVector& x);
    friend
    std::istream& operator>>(std::istream& is, SentenceVector& x);
    
  private:
    sent_set_type __sents;
  };
};

namespace std
{
  inline
  void swap(cicada::SentenceVector& x, cicada::SentenceVector& y)
  {
    x.swap(y);
  }
};

#endif
