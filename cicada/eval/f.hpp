// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__F__HPP__
#define __CICADA__EVAL__F__HPP__ 1

#include <cicada/eval/score.hpp>

namespace cicada
{
  namespace eval
  {
    class F : public Score
    {
    public:
      typedef double count_type;
      
    public:
      F() : match_ref(0), match_hyp(0), norm_ref(0), norm_hyp(0) {}
      
    public:
      double score() const
      {
	const double precision = (norm_hyp != 0.0 ? match_hyp / norm_hyp : 0.0);
	const double recall    = (norm_ref != 0.0 ? match_ref / norm_ref : 0.0);
	
	return precision * recall / (0.5 * recall + 0.5 * precision);
      }

      void assign(const score_type& score)
      {
	const F* rhs = dynamic_cast<const F*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid F measure");

	match_ref = rhs->match_ref;
	match_hyp = rhs->match_hyp;
	norm_ref  = rhs->norm_ref;
	norm_hyp  = rhs->norm_hyp;
      }

      void plus_equal(const score_type& score)
      {
	const F* rhs = dynamic_cast<const F*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid F measure");

	match_ref += rhs->match_ref;
	match_hyp += rhs->match_hyp;
	norm_ref  += rhs->norm_ref;
	norm_hyp  += rhs->norm_hyp;
      }
      
      void minus_equal(const score_type& score)
      {
	const F* rhs = dynamic_cast<const F*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid F measure");

	match_ref -= rhs->match_ref;
	match_hyp -= rhs->match_hyp;
	norm_ref  -= rhs->norm_ref;
	norm_hyp  -= rhs->norm_hyp;
      }

      void multiplies_equal(const double& scale)
      {
	match_ref *= scale;
	match_hyp *= scale;
	norm_ref  *= scale;
	norm_hyp  *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	match_ref /= scale;
	match_hyp /= scale;
	norm_ref  /= scale;
	norm_hyp  /= scale;
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const std::string& encoded);
      
    protected:
      virtual const char* __description() const = 0;
      
    protected:
      count_type match_ref;
      count_type match_hyp;
      count_type norm_ref;
      count_type norm_hyp;
    };
    
  };
};

#endif
