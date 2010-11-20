#include <set>
#include <algorithm>

#include "per.hpp"

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace eval
  {

    struct PERScorerConstant
    {
      struct COSTS
      {
	static const double insertion;
	static const double deletion;
	static const double substitution;
      };
    };
    
    const double PERScorerConstant::COSTS::insertion    = 1.0;
    const double PERScorerConstant::COSTS::deletion     = 1.0;
    const double PERScorerConstant::COSTS::substitution = 1.0;

    class PERScorerImpl : public PERScorerConstant
    {
    private:
      friend class PERScorer;
      
    public:
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      
      struct Score
      {
	double score;
	double insertion;
	double deletion;
	double substitution;
	
	Score() : score(), insertion(0), deletion(0), substitution(0) {}

	Score(const double& __score,
	      const double& __insertion,
	      const double& __deletion,
	      const double& __substitution)
	  : score(__score), insertion(__insertion), deletion(__deletion), substitution(__substitution) {}
      };

      typedef Score value_type;
      
      PERScorerImpl() { }
      PERScorerImpl(const sentence_type& __ref) : ref(__ref) {}

      
      value_type operator()(const sentence_type& sentence)
      {
	typedef std::multiset<word_type, std::less<word_type>, std::allocator<word_type> > word_set_type;

	word_set_type seq1(sentence.begin(), sentence.end());
	word_set_type seq2(ref.begin(), ref.end());

	sentence_type ins_words;
	sentence_type del_words;
	
	std::set_difference(seq1.begin(), seq1.end(), seq2.begin(), seq2.end(), std::back_inserter(ins_words));
	std::set_difference(seq2.begin(), seq2.end(), seq1.begin(), seq1.end(), std::back_inserter(del_words));
	
	const size_t i_size = ins_words.size();
	const size_t d_size = del_words.size();

	value_type value;
	
	if (i_size > d_size) {
	  value.insertion    = 0;
	  value.deletion     = COSTS::deletion * d_size;
	  value.substitution = COSTS::substitution * (i_size - d_size);
	} else {
	  value.insertion    = COSTS::insertion * i_size;
	  value.deletion     = 0;
	  value.substitution = COSTS::substitution * (d_size - i_size);
	}
	value.score = value.insertion + value.deletion + value.substitution;
	
	return value;
      }
      
    private:
      sentence_type ref;
    };
   
    PERScorer::PERScorer(const PERScorer& x)
      : Scorer(static_cast<const Scorer&>(*this))
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    PERScorer::~PERScorer()
    {
      clear();
    }
    
    PERScorer& PERScorer::operator=(const PERScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);

      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      return *this;
    }
    
    void PERScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void PERScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    
    PERScorer::score_ptr_type PERScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double score_best = std::numeric_limits<double>::infinity();
      std::auto_ptr<PER> per(new PER());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence);
	
	if (value.score < score_best) {
	  score_best = value.score;
	  
	  per->insertion    = value.insertion;
	  per->deletion     = value.deletion;
	  per->substitution = value.substitution;
	}
	
	per->references += evaluator.ref.size();
      }
      
      if (! impl.empty())
	per->references /= impl.size();

      return score_ptr_type(per.release());
    }
  };
};

