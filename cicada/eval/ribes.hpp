// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__RIBES__HPP__
#define __CICADA__EVAL__RIBES__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>
#include <cicada/matcher.hpp>

namespace cicada
{
  namespace eval
  {
    class RIBESScorer;

    class RIBES : public Score
    {
    private:
      friend class RIBESScorer;
      
    public:
      typedef double count_type;
      
      
    public:
      RIBES() : score(0), norm(0) {}
      
      double score() const
      {
	return score / norm;
      }

      bool equal(const score_type& score) const
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");

	return (score == rhs->score && norm == rhs->norm);
	
      }
      
      void assign(const score_type& score)
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");

	score = rhs->score;
	norm  = rhs->norm;
      }

      void plus_equal(const score_type& score)
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");

	score += rhs->score;
	norm  += rhs->norm;
      }
      
      void minus_equal(const score_type& score)
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");

	score -= rhs->score;
	norm  -= rhs->norm;
      }

      void multiplies_equal(const double& scale)
      {
	score *= scale;
	norm  *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	score /= scale;
	norm  /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new RIBES());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new RIBES(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);
      
    private:
      count_type score;
      count_type norm;
    };

    class RIBESScorerImpl;
    
    class RIBESScorer : public Scorer
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
      typedef RIBESScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      RIBESScorer() : impl(), weights(), matcher(0) { }
      RIBESScorer(const weights_type& __weights) : impl(), weights(__weights), matcher(0) {}
      RIBESScorer(const weights_type& __weights, const matcher_type* __matcher) : impl(), weights(__weights), matcher(__matcher) {}
      
      RIBESScorer(const RIBESScorer& x);
      ~RIBESScorer();
      RIBESScorer& operator=(const RIBESScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new RIBESScorer(*this)); }
      
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
