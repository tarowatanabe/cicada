// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__COMBINED__HPP__
#define __CICADA__EVAL__COMBINED__HPP__ 1

#include <vector>

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class CombinedScorer;

    class Combined : public Score
    {
    private:
      friend class CombinedScorer;
      
    private:
      typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;
      typedef std::vector<double, std::allocator<double> > weight_set_type;
      
    public:
      Combined() : scores(), weights() {}
      
      std::pair<double, double> score() const
      {
	std::pair<double, double> result(0.0, 0.0);
	
	weight_set_type::const_iterator witer = weights.begin();
	score_ptr_set_type::const_iterator siter_end = scores.end();
	for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter, ++ witer) {
	  const std::pair<double, double> tmp = (*siter)->score();
	  
	  result.first  += (*witer) * tmp.first;
	  result.second += (*witer) * tmp.second;
	}
	
	return result;
      }

      void assign(const score_type& score)
      {
	const Combined* rhs = dynamic_cast<const Combined*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid combined score");
	
	scores.resize(rhs->scores.size());
	
	score_ptr_set_type::iterator iter = scores.begin();
	score_ptr_set_type::const_iterator siter_end = rhs->scores.end();
	for (score_ptr_set_type::const_iterator siter = rhs->scores.begin(); siter != siter_end; ++ siter, ++ iter)
	  (*iter)->assign(*(*ster));
	
	weights = rhs->weights;
      }

      void plus_equal(const score_type& score)
      {
	const Combined* rhs = dynamic_cast<const Combined*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid combined score");
	
	scores.resize(utils::bithack::max(scores.size(), rhs->scores.size()));
	
	score_ptr_set_type::iterator iter = scores.begin();
	score_ptr_set_type::const_iterator siter_end = rhs->scores.end();
	for (score_ptr_set_type::const_iterator siter = rhs->scores.begin(); siter != siter_end; ++ siter, ++ iter)
	  (*iter)->plus_equal(*(*siter));
	
	if (! rhs->weights.empty())
	  weights = rhs->weights;
      }
      
      void minus_equal(const score_type& score)
      {
	const Combined* rhs = dynamic_cast<const Combined*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid combined score");	

	scores.resize(utils::bithack::max(scores.size(), rhs->scores.size()));
	
	score_ptr_set_type::iterator iter = scores.begin();
	score_ptr_set_type::const_iterator siter_end = rhs->scores.end();
	for (score_ptr_set_type::const_iterator siter = rhs->scores.begin(); siter != siter_end; ++ siter, ++ iter)
	  (*iter)->minus_equal(*(*siter));
	
	if (! rhs->weights.empty())
	  weights = rhs->weights;
      }
      
      void multiplies_equal(const double& scale)
      {
	score_ptr_set_type::iterator siter_end = scores.end();
	for (score_ptr_set_type::iterator siter = scores.begin(); siter != siter_end; ++ siter, ++ iter)
	  (*siter)->multiplies_equal(scale);
      }
      
      void divides_equal(const double& scale)
      {
	score_ptr_set_type::iterator siter_end = scores.end();
	for (score_ptr_set_type::iterator siter = scores.begin(); siter != siter_end; ++ siter, ++ iter)
	  (*siter)->divides_equal(scale);
      }

      score_ptr_type zero() const
      {
	return score_ptr_type(new Combined());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new Combined(*this));
      }
      
    private:
      score_set_type  scores;
      weight_set_type weights;
    };
    
    class CombinedScorerImpl;
    
    class CombinedScorer : public Scorer
    {
    public:
      typedef double count_type;
      
    private:
      typedef CombinedScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
      
    public:
      CombinedScorer() : impl() { }
      CombinedScorer(const CombinedScorer& x);
      ~CombinedScorer();
      CombinedScorer& operator=(const CombinedScorer& x);
      
      bool error_metric() const { return true; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new CombinedScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
    };
  };
};

#endif
