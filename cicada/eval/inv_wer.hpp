// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__INV_WER__HPP__
#define __CICADA__EVAL__INV_WER__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>
#include <cicada/matcher.hpp>

namespace cicada
{
  namespace eval
  {
    class InvWERScorer;

    class InvWER : public Score
    {
    private:
      friend class InvWERScorer;
      
    public:
      typedef double count_type;
      
      
    public:
      InvWER() : insertion(0), deletion(0), substitution(0), inversion(0), references(0) {}
      
      double score() const
      {
	const count_type edits = insertion + deletion + substitution + inversion;
	return edits / references;
      }

      double loss() const { return score(); }

      double reward() const { return 1.0 - score(); }
      
      bool error_metric() const { return true; }

      bool equal(const score_type& score) const
      {
	const InvWER* rhs = dynamic_cast<const InvWER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid InvWER score");

	return (insertion == rhs->insertion
		&& deletion == rhs->deletion
		&& substitution == rhs->substitution
		&& inversion == rhs->inversion
		&& references == rhs->references);
      }

      void assign(const score_type& score)
      {
	const InvWER* rhs = dynamic_cast<const InvWER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid InvWER score");

	insertion    = rhs->insertion;
	deletion     = rhs->deletion;
	substitution = rhs->substitution;
	inversion    = rhs->inversion;
	references   = rhs->references;
      }

      void plus_equal(const score_type& score)
      {
	const InvWER* rhs = dynamic_cast<const InvWER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid InvWER score");

	insertion    += rhs->insertion;
	deletion     += rhs->deletion;
	substitution += rhs->substitution;
	inversion    += rhs->inversion;
	references   += rhs->references;
      }
      
      void minus_equal(const score_type& score)
      {
	const InvWER* rhs = dynamic_cast<const InvWER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid InvWER score");

	insertion    -= rhs->insertion;
	deletion     -= rhs->deletion;
	substitution -= rhs->substitution;
	inversion    -= rhs->inversion;
	references   -= rhs->references;
      }

      void multiplies_equal(const double& scale)
      {
	insertion    *= scale;
	deletion     *= scale;
	substitution *= scale;
	inversion    *= scale;
	references   *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	insertion    /= scale;
	deletion     /= scale;
	substitution /= scale;
	inversion    /= scale;
	references   /= scale;
      }

      score_ptr_type zero() const
      {
	return score_ptr_type(new InvWER());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new InvWER(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(utils::piece::const_iterator& iter, utils::piece::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);
      
    private:
      count_type insertion;
      count_type deletion;
      count_type substitution;
      count_type inversion;
      count_type references;
    };

    class InvWERScorerImpl;
    
    class InvWERScorer : public Scorer
    {
    public:
      typedef double count_type;
      typedef double weight_type;
      typedef cicada::Matcher matcher_type;

      struct weights_type
      {
	weight_type match;
	weight_type substitution;
	weight_type insertion;
	weight_type deletion;
	weight_type inversion;
	
	weights_type() : match(0.2), substitution(1.0), insertion(1.0), deletion(1.0), inversion(1.0) {}
      };

    private:
      typedef InvWERScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      InvWERScorer() : impl(), weights(), matcher(0) { }
      InvWERScorer(const weights_type& __weights) : impl(), weights(__weights), matcher(0) {}
      InvWERScorer(const weights_type& __weights, const matcher_type* __matcher) : impl(), weights(__weights), matcher(__matcher) {}
      
      InvWERScorer(const InvWERScorer& x);
      ~InvWERScorer();
      InvWERScorer& operator=(const InvWERScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new InvWERScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
      weights_type  weights;
      
      const matcher_type* matcher;
    };
  };
};

#endif
