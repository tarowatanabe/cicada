// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__WER__HPP__
#define __CICADA__EVAL__WER__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>
#include <cicada/matcher.hpp>

namespace cicada
{
  namespace eval
  {
    class WERScorer;

    class WER : public Score
    {
    private:
      friend class WERScorer;
      
    public:
      typedef double count_type;
      
      
    public:
      WER() : insertion(0), deletion(0), substitution(0), references(0) {}
      
      double score() const
      {
	const count_type edits = insertion + deletion + substitution;
	return edits / references;
      }

      bool equal(const score_type& score) const
      {
	const WER* rhs = dynamic_cast<const WER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid WER score");

	return (insertion == rhs->insertion
		&& deletion == rhs->deletion
		&& substitution == rhs->substitution
		&& references == rhs->references);
      }

      void assign(const score_type& score)
      {
	const WER* rhs = dynamic_cast<const WER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid WER score");

	insertion    = rhs->insertion;
	deletion     = rhs->deletion;
	substitution = rhs->substitution;
	references   = rhs->references;
      }

      void plus_equal(const score_type& score)
      {
	const WER* rhs = dynamic_cast<const WER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid WER score");

	insertion    += rhs->insertion;
	deletion     += rhs->deletion;
	substitution += rhs->substitution;
	references   += rhs->references;
      }
      
      void minus_equal(const score_type& score)
      {
	const WER* rhs = dynamic_cast<const WER*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid WER score");

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
	return score_ptr_type(new WER());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new WER(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const std::string& encoded);
      
    private:
      count_type insertion;
      count_type deletion;
      count_type substitution;
      count_type references;
    };

    class WERScorerImpl;
    
    class WERScorer : public Scorer
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
	
	weights_type() : match(0.2), substitution(1.0), insertion(1.0), deletion(1.0) {}
      };

    private:
      typedef WERScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      WERScorer() : impl(), weights(), matcher(0) { }
      WERScorer(const weights_type& __weights) : impl(), weights(__weights), matcher(0) {}
      WERScorer(const weights_type& __weights, const matcher_type* __matcher) : impl(), weights(__weights), matcher(__matcher) {}
      
      WERScorer(const WERScorer& x);
      ~WERScorer();
      WERScorer& operator=(const WERScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new WERScorer(*this)); }
      
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
