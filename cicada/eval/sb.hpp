// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__SB__HPP__
#define __CICADA__EVAL__SB__HPP__ 1

// skip-bigram

#include <cicada/eval/f.hpp>

namespace cicada
{
  namespace eval
  {
    class SBScorer;

    class SB : public F
    {
      friend class SBScorer;
    public:
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new SB());
      }
      
      score_ptr_type clone() const
      {
	return score_ptr_type(new SB(*this));
      }

      
    protected:
      const char* __description() const;
    };

    class SBScorerImpl;
    
    class SBScorer : public Scorer
    {
    public:
      typedef double count_type;
      
    private:
      typedef SBScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
  
    public:
      SBScorer(const int& __window) : impl(), window(__window) { }
      SBScorer(const SBScorer& x);
      ~SBScorer();
      SBScorer& operator=(const SBScorer& x);

    public:      
      bool error_metric() const { return false; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new SBScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
      int window;
    };
    
  };
};


#endif
