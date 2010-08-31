#include "ter.hpp"

#include <google/dense_hash_set>

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace eval
  {

    struct TERScorerConstant
    {
      static const double cost_ins;
      static const double cost_del;
      static const double cost_sub;
      static const double cost_shi;
      
      static const int max_shift_size = 10;
      static const int max_shift_dist = 50;
    };
    
    const double TERScorerConstant::cost_ins = 1.0;
    const double TERScorerConstant::cost_del = 1.0;
    const double TERScorerConstant::cost_sub = 1.0;
    const double TERScorerConstant::cost_shi = 1.0;

    
    
    // Do we really implement this...???
    class TERScorerImpl : public TERScorerConstant
    {
    private:
      friend class TERScorer;
      
    public:
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      
      typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;

      struct Score
      {
	double score;
	int insertion;
	int deletion;
	int substitution;
	int shift;
	
	Score() : score(), insertion(0), deletion(0), substitution(0), shift(0) {}
      };

      typedef Score value_type;
      
      TERScorerImpl() { words_unique.set_empty_key(word_type()); }
      TERScorerImpl(const sentence_type& __ref)
	: ref(__ref) { words_unique.set_empty_key(word_type()); words_unique.insert(__ref.begin(), __ref.end()); }

      
      value_type operator()(const sentence_type& sentence)
      {
	
	
	
	
      }
      
    private:
      sentence_type ref;
      word_set_type words_unique;
    };
   
    TERScorer::TERScorer(const TERScorer& x)
      : Scorer(static_cast<const Scorer&>(*this))
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    TERScorer::~TERScorer()
    {
      clear();
    }
    
    TERScorer& TERScorer::operator=(const TERScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);

      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      return *this;
    }
    
    void TERScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void TERScorer::insert(const sentence_type& sentence)
    {
      if (split) {
	sentence_type sentence_split;
	split_non_ascii_characters(sentence, sentence_split);
	impl.push_back(new impl_type(sentence_split));
      } else
	impl.push_back(new impl_type(sentence));
    }
    
    TERScorer::score_ptr_type TERScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence_split;
      if (split)
	split_non_ascii_characters(__sentence, sentence_split);
      
      const sentence_type& sentence = (split ? sentence_split : __sentence);

      double score_best = std::numeric_limits<double>::infinity();

      std::auto_ptr<TER> ter(new TER());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence);
	
	if (value.score < score_best) {
	  score_best = value.score;
	  
	  ter->insertion    = value.insertion;
	  ter->deletion     = value.deletion;
	  ter->substitution = value.substitution;
	  ter->shift        = value.shift;
	}

	ter->references += evaluator.ref.size();
      }
      
      if (! impl.empty())
	ter->references /= impl.size();

      return score_ptr_type(ter.release());
    }
  };
};

