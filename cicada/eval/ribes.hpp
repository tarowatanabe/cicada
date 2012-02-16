// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__RIBES__HPP__
#define __CICADA__EVAL__RIBES__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class RibesScorer;

    class Ribes : public Score
    {
    private:
      friend class RibesScorer;
      
    public:
      typedef double count_type;
      
    public:
      Ribes() : distance(0), penalty(0) {}
      
      double score() const
      {
	return (penalty != 0.0 ? distance / penalty : 0.0);
      }
      
      bool equal(const score_type& score) const
      {
	const Ribes* rhs = dynamic_cast<const Ribes*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Ribes score");
	
	return (distance == rhs->distance && penalty == rhs->penalty);
      }
      
      void assign(const score_type& score)
      {
	const Ribes* rhs = dynamic_cast<const Ribes*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Ribes score");
	
	distance  = rhs->distance;
	penalty   = rhs->penalty;
      }

      void plus_equal(const score_type& score)
      {
	const Ribes* rhs = dynamic_cast<const Ribes*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Ribes score");

	distance  += rhs->distance;
	penalty   += rhs->penalty;
      }
      
      void minus_equal(const score_type& score)
      {
	const Ribes* rhs = dynamic_cast<const Ribes*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid Ribes score");

	distance  -= rhs->distance;
	penalty   -= rhs->penalty;
      }

      void multiplies_equal(const double& scale)
      {
	distance  *= scale;
	penalty   *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	distance  /= scale;
	penalty   /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new Ribes());
      }
      
      score_ptr_type clone() const
      {
	return score_ptr_type(new Ribes(*this));
      }
      
      std::string description() const;
      std::string encode() const;
      
      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);
      
    public:
      count_type distance;
      count_type penalty;
    };

    class RibesScorerImpl;
    
    class RibesScorer : public Scorer
    {
    public:
      typedef double count_type;
      typedef double weight_type;
      
    private:
      typedef RibesScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      RibesScorer() : impl(), alpha(), beta(), spearman(false), kendall(true) { }
      RibesScorer(const weight_type& __alpha, const weight_type& __beta) : impl(), alpha(__alpha), beta(__beta), spearman(false), kendall(true) {}
      RibesScorer(const weight_type& __alpha, const weight_type& __beta, const bool __spearman)
	: impl(), alpha(__alpha), beta(__beta), spearman(__spearman), kendall(! __spearman) {}
      
      RibesScorer(const RibesScorer& x);
      ~RibesScorer();
      RibesScorer& operator=(const RibesScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new RibesScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
      weight_type  alpha;
      weight_type  beta;
      bool spearman;
      bool kendall;
    };
  };
};

#endif
