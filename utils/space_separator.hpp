// -*- mode: c++ -*-

#ifndef __UTILS__SPACE_SEPARATOR__HPP__
#define __UTILS__SPACE_SEPARATOR__HPP__ 1

#include <cctype>
#include <iterator>

namespace utils
{

  template <typename IteratorTag>
  struct __space_separator_assigner
  {
    template<typename Iterator, typename Token>
    static void assign(Iterator first, Iterator last, Token& tok) { tok.assign(first, last); }
    
    template <typename Token>
    static void clear(Token& tok) { }
    
    template <typename Token, typename Value>
    static void plus_equal(Token& tok, const Value& value) { }
  };
  
  template <>
  struct __space_separator_assigner<std::input_iterator_tag >
  {
    template<typename Iterator, typename Token>
    static void assign(Iterator first, Iterator last, Token& tok) { }
    
    template <typename Token>
    static void clear(Token& tok) { tok = Token(); }
    
    template <typename Token, typename Value>
    static void plus_equal(Token& tok, const Value& value) { tok += value; }
  };
  
  struct space_separator
  {
    void reset() {}
    
    template <typename Iterator, typename Token>
    bool operator()(Iterator& next, Iterator end, Token& tok)
    {
      typedef typename std::iterator_traits<Iterator>::iterator_category __category;
      typedef __space_separator_assigner<__category> token_assigner;
      
      token_assigner::clear(tok);
      
      // skip past all dropped_delims
      for (/**/; next != end && __is_dropped(*next); ++ next) {}
      
      if (next == end)
	return false;
      
      Iterator start(next);
      for (/**/; next != end && ! __is_dropped(*next); ++ next)
	token_assigner::plus_equal(tok, *next);
      
      token_assigner::assign(start, next, tok);
      
      return true; 
    }
    
    template <typename Char>
    bool __is_dropped(Char E) const
    {
      return std::isspace(E) != 0;
    }
  };
  
};

#endif
