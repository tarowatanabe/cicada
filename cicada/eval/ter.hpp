// -*- mode: c++ -*-

#ifndef __CICADA__EVAL__TER__HPP__
#define __CICADA__EVAL__TER__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class TERScorer;

    class TER : public Score
    {
    private:
      friend class TERScorer;
      
    public:
      typedef double count_type;
      
      
    public:
      TER() : insertion(0), deletion(0), substitution(0), shift(0), references(0) {}
      
      std::pair<double, double> score() const
      {
	const count_type edits = insertion + deletion + substitution + shift;
	return std::make_pair(edits / references, 0.0);
      }

      void assign(const score_type& score)
      {
	const TER* rhs = dynamic_cast<const TER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid TER score");

	insertion    = rhs->insertion;
	deletion     = rhs->deletion;
	substitution = rhs->substitution;
	shift        = rhs->shift;
	references   = rhs->references;
      }

      void plus_equal(const score_type& score)
      {
	const TER* rhs = dynamic_cast<const TER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid TER score");

	insertion    += rhs->insertion;
	deletion     += rhs->deletion;
	substitution += rhs->substitution;
	shift        += rhs->shift;
	references   += rhs->references;
      }
      
      void minus_equal(const score_type& score)
      {
	const TER* rhs = dynamic_cast<const TER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid TER score");

	insertion    -= rhs->insertion;
	deletion     -= rhs->deletion;
	substitution -= rhs->substitution;
	shift        -= rhs->shift;
	references   -= rhs->references;
      }

      void multiplies_equal(const double& scale)
      {
	insertion    *= scale;
	deletion     *= scale;
	substitution *= scale;
	shift        *= scale;
	references   *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	insertion    /= scale;
	deletion     /= scale;
	substitution /= scale;
	shift        /= scale;
	references   /= scale;
      }

      score_ptr_type zero() const
      {
	return score_ptr_type(new TER());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new TER(*this));
      }
      
    private:
      count_type insertion;
      count_type deletion;
      count_type substitution;
      count_type shift;
      count_type references;
    };

    class TERScorerImpl;
    
    class TERScorer : public Scorer
    {
    public:
      typedef double count_type;
      
    private:
      typedef TERScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      TERScorer() : impl() { }
      TERScorer(const TERScorer& x);
      ~TERScorer();
      TERScorer& operator=(const TERScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new TERScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
    };
  };
};

#endif
