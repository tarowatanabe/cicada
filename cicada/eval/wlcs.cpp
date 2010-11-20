
#include <numeric>
#include <memory>

#include "wlcs.hpp"

#include <utils/vector2.hpp>

// @inproceedings{lin-och:2004:ACL,
//  author    = {Lin, Chin-Yew  and  Och, Franz Josef},
//  title     = {Automatic Evaluation of Machine Translation Quality Using Longest Common Subsequence and Skip-Bigram Statistics},
//  booktitle = {Proceedings of the 42nd Meeting of the Association for Computational Linguistics (ACL'04), Main Volume},
//  year      = 2004,
//  month     = {July},
//  address   = {Barcelona, Spain},
//  pages     = {605--612}
//}


namespace cicada
{
  namespace eval
  {
    class WLCSScorerImpl
    {
    public:
      typedef Sentence sentence_type;
      typedef double count_type;
      
      typedef utils::vector2<count_type, std::allocator<count_type> > c_matrix_type;
      typedef utils::vector2<int, std::allocator<int> > w_matrix_type;

      struct Score
      {
	Score() : match(0), norm_ref(0), norm_hyp(0), precision(0), recall(0) {}

	count_type match;
	count_type norm_ref;
	count_type norm_hyp;
	
	count_type precision;
	count_type recall;
      };
      typedef Score value_type;
      
    public:
      WLCSScorerImpl(const sentence_type& __ref)
	: ref(__ref) {}
      
      value_type operator()(const sentence_type& hyp, const double& alpha)
      {
	value_type value;
	
	value.match = length_weight(distance(hyp, ref, alpha), true);
	value.norm_hyp = hyp.size();
	value.norm_ref = ref.size();
	
	value.precision = (value.norm_hyp > 0.0 ? value.match / value.norm_hyp : 0.0);
	value.recall    = (value.norm_ref > 0.0 ? value.match / value.norm_ref : 0.0);

	return value;
      }


      double length_weight(const double length, const double& alpha, bool reverse=false)
      {
	return (alpha == 1.0 ? length : pow(length, reverse ? 1.0 / alpha : alpha));
      }

      
      double distance(const sentence_type& x,
		      const sentence_type& y,
		      const double& alpha)
      {
	const int m = x.size();
	const int n = y.size();

	c.clear();
	w.clear();
	
	c.resize(m + 1, n + 1, 0.0);
	w.resize(m + 1, n + 1, 0);
	
	for (int i = 1; i <= m; ++ i)
	  for (int j = 1; j <= n; ++ j) {
	    
	    if (x[i - 1] == y[j - 1]) {
	      const int k = w(i - 1, j - 1);
	      
	      c(i, j) = c(i - 1, j - 1) + length_weight(k + 1, alpha) - length_weight(k, alpha);
	      w(i, j) = k + 1;
	    } else if (c(i - 1, j) > c(i, j - 1)) {
	      c(i, j) = c(i - 1, j);
	      w(i, j) = 0;
	    } else {
	      c(i, j) = c(i, j - 1);
	      w(i, j) = 0;
	    }
	  }
	
	return c(m, n);
      }
      
      sentence_type ref;
      
      c_matrix_type c;
      w_matrix_type w;
    };
    
    WLCSScorer::WLCSScorer(const WLCSScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	alpha(x.alpha)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    WLCSScorer::~WLCSScorer()
    {
      clear();
    }

    WLCSScorer& WLCSScorer::operator=(const WLCSScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));

      alpha = x.alpha;
      
      return *this;
    }

    void WLCSScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void WLCSScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    

    WLCSScorer::score_ptr_type WLCSScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double precision_max = 0.0;
      double recall_max = 0.0;
      
      std::auto_ptr<WLCS> wlcs(new WLCS());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	const impl_type::value_type value = evaluator(sentence, alpha);
	
	if (value.precision > precision_max) {
	  precision_max = value.precision;
	  wlcs->match_hyp = value.match;
	  wlcs->norm_hyp  = value.norm_hyp;
	}
	if (value.recall > recall_max) {
	  recall_max = value.recall;
	  wlcs->match_ref = value.match;
	  wlcs->norm_ref  = value.norm_ref;
	}
      }
      
      return score_ptr_type(wlcs.release());
    }
    
  };
};

