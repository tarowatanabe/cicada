// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__SK__HPP__
#define __CICADA__EVAL__SK__HPP__ 1

// string-kernel

#include <cicada/eval/f.hpp>

namespace cicada
{
  namespace eval
  {
    class SKScorer;

    class SK : public F
    {
      friend class SKScorer;

    public:      
      score_ptr_type zero() const
      {
	return score_ptr_type(new SK());
      }
      
      score_ptr_type clone() const
      {
	return score_ptr_type(new SK(*this));
      }

    protected:
      const char* __description() const;
    };

    class SKScorerImpl;
    
    class SKScorer : public Scorer
    {
    public:
      typedef double count_type;
      
    private:
      typedef SKScorerImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> >  impl_set_type;
  
    public:
      SKScorer(const int __p, const double& __decay) : impl(), p(__p), decay(__decay) { }
      SKScorer(const SKScorer& x);
      ~SKScorer();
      SKScorer& operator=(const SKScorer& x);

    public:      
      bool error_metric() const { return false; }
      
      scorer_ptr_type clone() const { return scorer_ptr_type(new SKScorer(*this)); }
      
      void clear();
      
      void insert(const sentence_type& sentence);
      score_ptr_type score(const sentence_type& __sentence) const;
      
    private:
      impl_set_type impl;
      int p;
      double decay;
    };
    
  };
};


#endif
