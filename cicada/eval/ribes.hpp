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
    class RIBESScorer;

    class RIBES : public Score
    {
    private:
      friend class RIBESScorer;
      
    public:
      typedef double count_type;
      
    public:
      RIBES() : distance(0), penalty(0) {}
      
      double score() const
      {
	return (penalty != 0.0 ? distance / penalty : 0.0);
      }
      
      bool equal(const score_type& score) const
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");
	
	return (distance == rhs->distance
		&& penalty == rhs->penalty);
	
      }
      
      void assign(const score_type& score)
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");
	
	distance  = rhs->distance;
	penalty   = rhs->penalty;
      }

      void plus_equal(const score_type& score)
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");

	distance  += rhs->distance;
	penalty   += rhs->penalty;
      }
      
      void minus_equal(const score_type& score)
      {
	const RIBES* rhs = dynamic_cast<const RIBES*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid RIBES score");

	distance -= rhs->distance;
	penalty  -= rhs->penalty;
      }

      void multiplies_equal(const double& scale)
      {
	distance *= scale;
	penalty  *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	distance /= scale;
	penalty  /= scale;
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
      count_type distance;
      count_type penalty;
    };

    class RIBESScorerImpl;
    
    class RIBESScorer : public Scorer
    {
    public:
      typedef double count_type;
      typedef double weight_type;
      
    private:
      typedef RIBESScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      RIBESScorer() : impl(), weight() { }
      RIBESScorer(const weight_type& __weight) : impl(), weight(__weight) {}
      
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
      weight_type  weight;
    };
  };
};

#endif
