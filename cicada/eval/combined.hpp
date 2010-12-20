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
      
      double score() const
      {
	double result = 0.0;
	weight_set_type::const_iterator witer = weights.begin();
	score_ptr_set_type::const_iterator siter_end = scores.end();
	for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter, ++ witer)
	  result += (*witer) * (*siter)->score();
	
	return result;
      }

      void assign(const score_type& score)
      {
	const Combined* rhs = dynamic_cast<const Combined*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid combined score");
	
	scores.reserve(utils::bithack::max(scores.size(), rhs->scores.size()));
	scores.resize(rhs->scores.size());
	
	score_ptr_set_type::iterator iter = scores.begin();
	score_ptr_set_type::const_iterator siter_end = rhs->scores.end();
	for (score_ptr_set_type::const_iterator siter = rhs->scores.begin(); siter != siter_end; ++ siter, ++ iter)
	  (*iter)->assign(*(*siter));
	
	weights = rhs->weights;
      }

      void plus_equal(const score_type& score)
      {
	const Combined* rhs = dynamic_cast<const Combined*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid combined score");
	
	scores.reserve(utils::bithack::max(scores.size(), rhs->scores.size()));
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

	scores.reserve(utils::bithack::max(scores.size(), rhs->scores.size()));
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
	for (score_ptr_set_type::iterator siter = scores.begin(); siter != siter_end; ++ siter)
	  (*siter)->multiplies_equal(scale);
      }
      
      void divides_equal(const double& scale)
      {
	score_ptr_set_type::iterator siter_end = scores.end();
	for (score_ptr_set_type::iterator siter = scores.begin(); siter != siter_end; ++ siter)
	  (*siter)->divides_equal(scale);
      }

      score_ptr_type zero() const
      {
	return score_ptr_type(new Combined());
      }

      score_ptr_type clone() const
      {
	std::auto_ptr<Combined> combined(new Combined());
	
	combined->scores.reserve(scores.size());
	score_ptr_set_type::const_iterator siter_end = scores.end();
	for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter)
	  combined->scores.push_back((*siter)->clone());
	
	combined->weights = weights;
	
	return score_ptr_type(combined.release());
      }

      std::string description() const;
      
    private:
      score_ptr_set_type scores;
      weight_set_type    weights;
    };
    
    class CombinedScorer : public Scorer
    {
    public:
      friend class Scorer;
      
    private:
      typedef std::vector<scorer_ptr_type, std::allocator<scorer_ptr_type> > scorer_ptr_set_type;
      typedef std::vector<double, std::allocator<double> > weight_set_type;
      
    public:
      CombinedScorer() : scorers() { }
      CombinedScorer(const CombinedScorer& x) : scorers(x.scorers), weights(x.weights) {}
      CombinedScorer& operator=(const CombinedScorer& x)
      {
	scorers = x.scorers;
	weights = x.weights;
	return *this;
      }
      
      bool error_metric() const { return error; }
      
      scorer_ptr_type clone() const
      {
	std::auto_ptr<CombinedScorer> scorer(new CombinedScorer());
	
	scorer->scorers.reserve(scorers.size());
	scorer_ptr_set_type::const_iterator siter_end = scorers.end();
	for (scorer_ptr_set_type::const_iterator siter = scorers.begin(); siter != siter_end; ++ siter)
	  scorer->scorers.push_back((*siter)->clone());

	scorer->weights = weights;

	return scorer_ptr_type(scorer.release());
      }
      
      void clear()
      {
	scorers.clear();
	weights.clear();
      }
      
      void insert(const sentence_type& sentence)
      {
	scorer_ptr_set_type::iterator siter_end = scorers.end();
	for (scorer_ptr_set_type::iterator siter = scorers.begin(); siter != siter_end; ++ siter)
	  (*siter)->insert(sentence);
      }
      
      score_ptr_type score(const sentence_type& __sentence) const
      {
	std::auto_ptr<Combined> combined(new Combined());
	
	combined->scores.reserve(scorers.size());
	scorer_ptr_set_type::const_iterator siter_end = scorers.end();
	for (scorer_ptr_set_type::const_iterator siter = scorers.begin(); siter != siter_end; ++ siter)
	  combined->scores.push_back((*siter)->score(__sentence));
	
	combined->weights = weights;
	
	return score_ptr_type(combined.release());
      }
      
    private:
      scorer_ptr_set_type scorers;
      weight_set_type weights;
      bool error;
    };
  };
};

#endif
