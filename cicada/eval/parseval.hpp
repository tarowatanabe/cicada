// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__PARSEVAL__HPP__
#define __CICADA__EVAL__PARSEVAL__HPP__ 1

#include <set>

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class ParsevalScorer;
    
    class Parseval : public Score
    {
    private:
      friend class ParsevalScorer;
      
    public:
      typedef double count_type;
      
    public:
      Parseval() : matched(0), test(0), reference(0) {}
      
      double score() const
      {
	return 2.0 * matched / (test + reference);
      }

      bool equal(const score_type& score) const
      {
	const Parseval* rhs = dynamic_cast<const Parseval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Parseval score");

	return (matched == rhs->matched
		&& test == rhs->test
		&& reference == rhs->reference);
      }

      void assign(const score_type& score)
      {
	const Parseval* rhs = dynamic_cast<const Parseval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Parseval score");

	matched   = rhs->matched;
	test      = rhs->test;
	reference = rhs->reference;
      }

      void plus_equal(const score_type& score)
      {
	const Parseval* rhs = dynamic_cast<const Parseval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Parseval score");

	matched   += rhs->matched;
	test      += rhs->test;
	reference += rhs->reference;
      }
      
      void minus_equal(const score_type& score)
      {
	const Parseval* rhs = dynamic_cast<const Parseval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Parseval score");

	matched   -= rhs->matched;
	test      -= rhs->test;
	reference -= rhs->reference;
      }

      void multiplies_equal(const double& scale)
      {
	matched   *= scale;
	test      *= scale;
	reference *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	matched   /= scale;
	test      /= scale;
	reference /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new Parseval());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new Parseval(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);

    private:
      count_type matched;
      count_type test;
      count_type reference;
    };

    class ParsevalScorerImpl;
    
    class ParsevalScorer : public Scorer
    {
    private:
      friend class ParsevalScorerImpl;
      
    public:
      typedef double count_type;
      
    private:
      typedef ParsevalScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
      typedef std::set<word_type, std::less<word_type>, std::allocator<word_type> > ignored_type;
      
    public:
      ParsevalScorer() : impl(), ignored() { }
      template <typename Iterator>
      ParsevalScorer(Iterator first, Iterator last) : impl(), ignored(first, last) { }
      
      ParsevalScorer(const ParsevalScorer& x);
      ~ParsevalScorer();
      ParsevalScorer& operator=(const ParsevalScorer& x);
      
      bool error_metric() const { return false; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new ParsevalScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
      ignored_type  ignored;
    };
    
  };
};


#endif
