// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__DEPEVAL__HPP__
#define __CICADA__EVAL__DEPEVAL__HPP__ 1

#include <set>

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class DepevalScorer;
    
    class Depeval : public Score
    {
    private:
      friend class DepevalScorer;
      
    public:
      typedef double count_type;
      
    public:
      Depeval() : matched(0), total(0) {}
      
      double score() const
      {
	return matched / total;
      }

      bool equal(const score_type& score) const
      {
	const Depeval* rhs = dynamic_cast<const Depeval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Depeval score");

	return (matched == rhs->matched && total == rhs->total);
      }

      void assign(const score_type& score)
      {
	const Depeval* rhs = dynamic_cast<const Depeval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Depeval score");

	matched   = rhs->matched;
	total     = rhs->total;
      }

      void plus_equal(const score_type& score)
      {
	const Depeval* rhs = dynamic_cast<const Depeval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Depeval score");

	matched   += rhs->matched;
	total     += rhs->total;
      }
      
      void minus_equal(const score_type& score)
      {
	const Depeval* rhs = dynamic_cast<const Depeval*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Depeval score");

	matched   -= rhs->matched;
	total     -= rhs->total;
      }

      void multiplies_equal(const double& scale)
      {
	matched   *= scale;
	total     *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	matched   /= scale;
	total     /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new Depeval());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new Depeval(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);

    private:
      count_type matched;
      count_type total;
    };

    class DepevalScorerImpl;
    
    class DepevalScorer : public Scorer
    {
    private:
      friend class DepevalScorerImpl;
      
    public:
      typedef double count_type;
      
    private:
      typedef DepevalScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      DepevalScorer() : impl() { }
      
      DepevalScorer(const DepevalScorer& x);
      ~DepevalScorer();
      DepevalScorer& operator=(const DepevalScorer& x);
      
      bool error_metric() const { return false; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new DepevalScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
    };
    
  };
};


#endif
