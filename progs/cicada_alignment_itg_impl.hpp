// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_ALIGNMENT_ITG_IMPL__HPP__
#define __CICADA_ALIGNMENT_ITG_IMPL__HPP__ 1

#include <numeric>
#include <limits>

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <cicada/sentence.hpp>
#include <cicada/alignment.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>

struct ITG
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef uint32_t  id_type;

  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

  struct word_pair_type
  {
    word_type source;
    word_type target;
    
    word_pair_type() : source(), target() {}
    word_pair_type(const word_type& __source, const word_type& __target)
      : source(__source), target(__target) {}
    
    friend
    bool operator==(const word_pair_type& x, const word_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    size_t hash_value(word_pair_type const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.source.id(), x.target.id());
    }
  };

  struct word_pair_unassigned : public utils::unassigned<word_type>
  {
    typedef utils::unassigned<word_type> unassigned_type;

    word_pair_type operator()() const
    {
      return word_pair_type(unassigned_type::operator()(),
			    unassigned_type::operator()());
    }
  };

  struct span_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    size_type first;
    size_type last;
    
    span_type() : first(0), last(0) {}
    span_type(const size_type& __first, const size_type& __last)
      : first(__first), last(__last) { }

    bool empty() const { return first == last; }
    size_type size() const { return last - first; }
    
    friend
    bool operator==(const span_type& x, const span_type& y)
    {
      return x.first == y.first && x.last == y.last;
    }
    
    friend
    size_t hash_value(span_type const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.first, x.last);
    }

    friend
    bool disjoint(const span_type& x, const span_type& y)
    {
      return x.last <= y.first || y.last <= x.first;
    }

    friend
    bool adjacent(const span_type& x, const span_type& y)
    {
      return x.last == y.first || y.last == x.first;
    }
    
  };

  struct span_pair_type
  {
    typedef span_type::size_type      size_type;
    typedef span_type::difference_type difference_type;

    span_type source;
    span_type target;
    
    span_pair_type() : source(), target() {}
    span_pair_type(const span_type& __source, const span_type& __target)
      : source(__source), target(__target) {}
    span_pair_type(size_type s1, size_type s2, size_type t1, size_type t2)
      : source(s1, s2), target(t1, t2) {}
    
    bool empty() const { return source.empty() && target.empty(); }
    size_type size() const { return source.size() + target.size(); }
    
    friend
    bool operator==(const span_pair_type& x, const span_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    size_t hash_value(span_pair_type const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x);
    }
    
    friend
    bool disjoint(const span_pair_type& x, const span_pair_type& y)
    {
      return disjoint(x.source, y.source) && disjoint(x.target, y.target);
    }

    friend
    bool adjacent(const span_pair_type& x, const span_pair_type& y)
    {
      return adjacent(x.source, y.source) && adjacent(x.target, y.target);
    }
    
  };

  // How to differentiate base or generative?
  struct rule_type
  {
    span_pair_type span;
    span_pair_type left;
    span_pair_type right;
    
    rule_type()
      : span(), left(), right() {}
    rule_type(const span_pair_type& __span)
      : span(__span), left(), right() {}
    rule_type(const span_pair_type& __span, const span_pair_type& __left, const span_pair_type& __right)
      : span(__span), left(__left), right(__right) {}
    
    bool is_terminal() const { return left.empty() && right.empty(); }
    bool is_straight() const { return ! is_terminal() && left.target.last  == right.target.first; }
    bool is_inverted() const { return ! is_terminal() && left.target.first == right.target.last; }
  };
  
  
  
  logprob_type prob_straight;
  logprob_type prob_inverted;
  logprob_type prob_terminal;
  
  // fallback probabilities...
  logprob_type p0_source;
  logprob_type p0_target;
  logprob_type p0_epsilon;
  logprob_type p0_word;
};

#endif
