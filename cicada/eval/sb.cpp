
#include <numeric>
#include <memory>

#include "sb.hpp"

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
    class SBScorerImpl
    {
    public:
      typedef Sentence sentence_type;
      typedef double count_type;
      
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
      SBScorerImpl(const sentence_type& __ref)
	: ref(__ref) {}
      
      value_type operator()(const sentence_type& hyp, const int& window)
      {
	value_type value;
	
	value.match    = distance(hyp, ref, window);
	value.norm_hyp = distance(hyp, hyp, window);
	value.norm_ref = distance(ref, ref, window);
	
	value.precision = (value.norm_hyp > 0.0 ? value.match / value.norm_hyp : 0.0);
	value.recall    = (value.norm_ref > 0.0 ? value.match / value.norm_ref : 0.0);

	return value;
      }
      
      
      double distance(const sentence_type& sent1,
		      const sentence_type& sent2,
		      const int& window)
      {
	double count = 0;
	
	const int size1 = sent1.size();
	const int size2 = sent2.size();
	
	for (int i = 0; i < size1; ++ i)
	  for (int j = 0; j < size2; ++ j)
	    if (sent1[i] == sent2[j])
	      for (int k = i + 1; k < size1 && (window < 0 || k - i - 1 <= window); ++ k)
		for (int l = j + 1; l < size2 && (window < 0 || l - j - 1 <= window); ++ l)
		  if (sent1[k] == sent2[l])
		    count ++;
	return count;
      }
      
      sentence_type ref;
    };
    
    SBScorer::SBScorer(const SBScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	window(x.window)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    SBScorer::~SBScorer()
    {
      clear();
    }

    SBScorer& SBScorer::operator=(const SBScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));

      window = x.window;
      
      return *this;
    }

    void SBScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void SBScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence_split;
      sentence_type sentence_lower;
      
      if (split)
	split_non_ascii_characters(__sentence, sentence_split);
      const sentence_type& __sentence_split = (split ? sentence_split : __sentence);
      
      if (lower)
	lower_case(__sentence_split, sentence_lower);
      const sentence_type& sentence = (lower ? sentence_lower : __sentence_split);

      impl.push_back(new impl_type(sentence));
    }
    

    SBScorer::score_ptr_type SBScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence_split;
      sentence_type sentence_lower;
      
      if (split)
	split_non_ascii_characters(__sentence, sentence_split);
      const sentence_type& __sentence_split = (split ? sentence_split : __sentence);
      
      if (lower)
	lower_case(__sentence_split, sentence_lower);
      const sentence_type& sentence = (lower ? sentence_lower : __sentence_split);
      
      
      double precision_max = 0.0;
      double recall_max = 0.0;
      
      std::auto_ptr<SB> sb(new SB());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	const impl_type::value_type value = evaluator(sentence, window);
	
	if (value.precision > precision_max) {
	  precision_max = value.precision;
	  sb->match_hyp = value.match;
	  sb->norm_hyp  = value.norm_hyp;
	}
	if (value.recall > recall_max) {
	  recall_max = value.recall;
	  sb->match_ref = value.match;
	  sb->norm_ref  = value.norm_ref;
	}
      }
      
      return score_ptr_type(sb.release());
    }
    
  };
};

