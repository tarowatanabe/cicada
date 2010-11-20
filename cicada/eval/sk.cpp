
#include <numeric>
#include <memory>

#include "sk.hpp"

#include <utils/vector2.hpp>

namespace cicada
{
  namespace eval
  {
    class SKScorerImpl
    {
    public:
      typedef Sentence sentence_type;
      typedef double count_type;
      typedef std::vector<count_type, std::allocator<count_type> > spectrum_type;
      typedef utils::vector2<count_type, std::allocator<count_type> > table_type;

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
      SKScorerImpl(const sentence_type& __ref)
	: ref(__ref) {}
      
      value_type operator()(const sentence_type& hyp, const int p, const double& decay)
      {
	value_type value;

	distance(hyp, ref, p, decay, spectrum);
	value.match = std::accumulate(spectrum.begin(), spectrum.end(), 0.0);
	
	distance(ref, ref, p, decay, spectrum_ref);
	value.norm_ref = std::accumulate(spectrum_ref.begin(), spectrum_ref.end(), 0.0);
	
	distance(hyp, hyp, p, decay, spectrum_hyp);
	value.norm_hyp = std::accumulate(spectrum_hyp.begin(), spectrum_hyp.end(), 0.0);
	
	value.precision = (value.norm_hyp > 0.0 ? value.match / value.norm_hyp : 0.0);
	value.recall    = (value.norm_ref > 0.0 ? value.match / value.norm_ref : 0.0);

	return value;
      }
      
      void distance(const sentence_type& sent1,
		    const sentence_type& sent2,
		    const int p,
		    const double& decay,
		    spectrum_type& spectrum)
      {
	const int n = sent1.size();
	const int m = sent2.size();
	
	spectrum.clear();
	spectrum.resize(p + 1, 0.0);
	
	dps.clear();
	dps.resize(n + 1, m + 1, 0.0);
	
	dp.clear();
	dp.resize(n + 1, m + 1, 0.0);

	// setup dps for p == 1
	for (int i = 1; i <= n; ++ i)
	  for (int j = 1; j <= m; ++ j)
	    if (sent1[i - 1] == sent2[j - 1]) {
	      dps(i, j) = decay * decay;
	      spectrum[1] += decay * decay;
	    }
	
	for (int l = 2; l <= p; ++ l)
	  for (int i = 1; i <= n; ++ i)
	    for (int j = 1; j <= m; ++ j) {
	      dp(i, j) = dps(i, j) + decay * dp(i - 1, j) + decay * dp(i, j - 1) - decay * decay * dp(i - 1, j - 1);
	      
	      if (sent1[i - 1] == sent2[j - 1]) {
		dps(i, j) = decay * decay * dp(i - 1, j - 1);
		spectrum[l] += dps(i, j);
	      }
	    }
      }

      
      sentence_type ref;
      spectrum_type spectrum;
      spectrum_type spectrum_ref;
      spectrum_type spectrum_hyp;
      table_type dps;
      table_type dp;
    };
    
    SKScorer::SKScorer(const SKScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	p(x.p),
	decay(x.decay)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    SKScorer::~SKScorer()
    {
      clear();
    }

    SKScorer& SKScorer::operator=(const SKScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));

      p     = x.p;
      decay = x.decay;
      
      return *this;
    }

    void SKScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void SKScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    

    SKScorer::score_ptr_type SKScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double precision_max = 0.0;
      double recall_max = 0.0;
      
      std::auto_ptr<SK> sk(new SK());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	const impl_type::value_type value = evaluator(sentence, p, decay);
	
	if (value.precision > precision_max) {
	  precision_max = value.precision;
	  sk->match_hyp = value.match;
	  sk->norm_hyp  = value.norm_hyp;
	}
	if (value.recall > recall_max) {
	  recall_max = value.recall;
	  sk->match_ref = value.match;
	  sk->norm_ref  = value.norm_ref;
	}
      }
      
      return score_ptr_type(sk.release());
    }
    
  };
};
