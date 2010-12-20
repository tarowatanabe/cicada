// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__PER__HPP__
#define __CICADA__EVAL__PER__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class PERScorer;

    class PER : public Score
    {
    private:
      friend class PERScorer;
      
    public:
      typedef double count_type;
      
      
    public:
      PER() : insertion(0), deletion(0), substitution(0), references(0) {}
      
      double score() const
      {
	const count_type edits = insertion + deletion + substitution;
	return edits / references;
      }

      void assign(const score_type& score)
      {
	const PER* rhs = dynamic_cast<const PER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid PER score");

	insertion    = rhs->insertion;
	deletion     = rhs->deletion;
	substitution = rhs->substitution;
	references   = rhs->references;
      }

      void plus_equal(const score_type& score)
      {
	const PER* rhs = dynamic_cast<const PER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid PER score");

	insertion    += rhs->insertion;
	deletion     += rhs->deletion;
	substitution += rhs->substitution;
	references   += rhs->references;
      }
      
      void minus_equal(const score_type& score)
      {
	const PER* rhs = dynamic_cast<const PER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid PER score");

	insertion    -= rhs->insertion;
	deletion     -= rhs->deletion;
	substitution -= rhs->substitution;
	references   -= rhs->references;
      }

      void multiplies_equal(const double& scale)
      {
	insertion    *= scale;
	deletion     *= scale;
	substitution *= scale;
	references   *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	insertion    /= scale;
	deletion     /= scale;
	substitution /= scale;
	references   /= scale;
      }

      score_ptr_type zero() const
      {
	return score_ptr_type(new PER());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new PER(*this));
      }
      
      std::string description() const;

    private:
      count_type insertion;
      count_type deletion;
      count_type substitution;
      count_type references;
    };

    class PERScorerImpl;
    
    class PERScorer : public Scorer
    {
    public:
      typedef double count_type;
      
    private:
      typedef PERScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      PERScorer() : impl() { }
      PERScorer(const PERScorer& x);
      ~PERScorer();
      PERScorer& operator=(const PERScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new PERScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
    };
  };
};

#endif
