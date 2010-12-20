//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <sstream>
#include <iterator>

#include "wer.hpp"

#include "utils/base64.hpp"

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace eval
  {

    std::string WER::description() const
    {
      std::ostringstream stream;
      stream << "wer: " << score()
	     << " " << insertion << '|' << deletion << '|' << substitution
	     << " length: " << references;
      
      return stream.str();
    }

    inline
    std::string escape_base64(const std::string& x)
    {
      std::string result;
      
      std::string::const_iterator iter_end = x.end();
      for (std::string::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	if (*iter == '/')
	  result += "\\/";
	else
	  result += *iter;
      
      return result;
    }

    std::string WER::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"wer\",";
      stream << "\"insertion\":\"" <<  escape_base64(utils::encode_base64(insertion)) << "\",";
      stream << "\"deletion\":\"" <<  escape_base64(utils::encode_base64(deletion)) << "\",";
      stream << "\"substitution\":\"" <<  escape_base64(utils::encode_base64(substitution)) << "\",";
      stream << "\"reference\":\"" <<  escape_base64(utils::encode_base64(references)) << "\"";
      stream << '}';
      return stream.str();
    }



    class WERScorerImpl
    {
    private:
      friend class WERScorer;
      
    public:
      typedef WERScorer wer_scorer_type;
      typedef wer_scorer_type::weights_type weights_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      typedef cicada::Matcher  matcher_type;
      
      struct Score
      {
	double score;
	
	double insertion;
	double deletion;
	double substitution;
	
	Score() : score(0), insertion(0), deletion(0), substitution(0) {}
	
	Score(const double& __score,
	      const double& __insertion,
	      const double& __deletion,
	      const double& __substitution)
	  : score(__score), insertion(__insertion), deletion(__deletion), substitution(__substitution) {}
      };

      typedef Score value_type;
      
      WERScorerImpl() { }
      WERScorerImpl(const sentence_type& __ref) : ref(__ref) {}
      
      value_type operator()(const sentence_type& sentence, const weights_type& weights, const matcher_type* matcher)
      {
	const sentence_type& a = ref;
	const sentence_type& b = sentence;

	const size_t a_size = a.size();
	const size_t b_size = b.size();
	
	if (a_size == 0) return value_type(b_size * weights.insertion, b_size * weights.insertion, 0, 0);
	if (b_size == 0) return value_type(a_size * weights.deletion, 0, a_size * weights.deletion, 0);
	
	std::vector<value_type, std::allocator<value_type> > curr(b_size + 1, value_type());
	std::vector<value_type, std::allocator<value_type> > prev(b_size + 1, value_type());
	
	for (size_t j = 1; j <= b_size; ++ j) {
	  prev[j].score     = prev[j - 1].score + weights.insertion;
	  prev[j].insertion = prev[j - 1].insertion + weights.insertion;
	}
	
	for (size_t i = 1; i <= a_size; ++ i) {
	  curr[0].score    = prev[0].score + weights.deletion;
	  curr[0].deletion = prev[0].deletion + weights.deletion;
	  
	  for (size_t j = 1; j <= b_size; ++ j) {
	    const double subst = (a[i - 1] == b[j - 1] ? 0.0 : (matcher && matcher->operator()(a[i - 1], b[j - 1])
								? weights.match
								: weights.substitution));
	    
	    const double score_del = prev[j].score     + weights.deletion;
	    const double score_ins = curr[j - 1].score + weights.insertion;
	    //const double score_sub = prev[j - 1].score + subst;
	    
	    curr[j] = prev[j - 1];
	    curr[j].score        += subst;
	    curr[j].substitution += subst;
	    
	    if (score_ins < curr[j].score) {
	      curr[j] = curr[j - 1];
	      curr[j].score     += weights.insertion;
	      curr[j].insertion += weights.insertion;
	    }
	    
	    if (score_del < curr[j].score) {
	      curr[j] = prev[j];
	      curr[j].score    += weights.deletion;
	      curr[j].deletion += weights.deletion;
	    }
	  }
	  std::swap(prev, curr);
	}
	return prev[b_size];
      }
      
    private:
      sentence_type ref;
    };
   
    WERScorer::WERScorer(const WERScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	weights(x.weights),
	matcher(0)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      if (x.matcher)
	matcher = &matcher_type::create(x.matcher->algorithm());
    }
    
    WERScorer::~WERScorer()
    {
      clear();
    }
    
    WERScorer& WERScorer::operator=(const WERScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      weights = x.weights;
      matcher = 0;
      if (x.matcher)
	matcher = &matcher_type::create(x.matcher->algorithm());
      
      return *this;
    }
    
    void WERScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void WERScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    
    WERScorer::score_ptr_type WERScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double score_best = std::numeric_limits<double>::infinity();

      std::auto_ptr<WER> wer(new WER());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence, weights, matcher);
	
	if (value.score < score_best) {
	  score_best = value.score;
	  
	  wer->insertion    = value.insertion;
	  wer->deletion     = value.deletion;
	  wer->substitution = value.substitution;
	}
	
	wer->references += evaluator.ref.size();
      }
      
      if (! impl.empty())
	wer->references /= impl.size();

      return score_ptr_type(wer.release());
    }
  };
};

