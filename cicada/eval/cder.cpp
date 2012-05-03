//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>

#include "cder.hpp"

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace eval
  {

    std::string CDER::description() const
    {
      std::ostringstream stream;
      stream << "cder: " << score()
	     << ' ' << insertion << '|' << deletion << '|' << substitution << '|' << jump
	     << " length: " << references;
      
      return stream.str();
    }
    
    std::string CDER::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"cder\",";
      stream << "\"edits\":[";
      stream << escaper(insertion)
	     << ',' << escaper(deletion)
	     << ',' << escaper(substitution)
	     << ',' << escaper(jump)
	     << ',' << escaper(references);
      stream << "]}";
      return stream.str();
    }

    typedef boost::fusion::tuple<double, double, double, double, double> cder_parsed_type;
    
    template <typename Iterator>
    struct cder_parser : boost::spirit::qi::grammar<Iterator, cder_parsed_type(), boost::spirit::standard::space_type>
    {
      cder_parser() : cder_parser::base_type(cder_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	cder_parsed %= (qi::lit('{')
			>> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"cder\"") >> qi::lit(',')
			>> qi::lit("\"edits\"") >> qi::lit(':')
			>> qi::lit('[')
			>> double_value >> qi::lit(',')
			>> double_value >> qi::lit(',')
			>> double_value >> qi::lit(',')
			>> double_value >> qi::lit(',')
			>> double_value
			>> qi::lit(']')
			>> qi::lit('}'));
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, cder_parsed_type(), space_type> cder_parsed;
    };

    Score::score_ptr_type CDER::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      cder_parser<iterator_type> parser;
      cder_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<CDER> cder(new CDER());
      cder->insertion    = boost::fusion::get<0>(parsed);
      cder->deletion     = boost::fusion::get<1>(parsed);
      cder->substitution = boost::fusion::get<2>(parsed);
      cder->jump         = boost::fusion::get<3>(parsed);
      cder->references   = boost::fusion::get<4>(parsed);
      
      return score_ptr_type(cder.release());
    }
    
    Score::score_ptr_type CDER::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

    class CDERScorerImpl
    {
    private:
      friend class CDERScorer;
      
    public:
      typedef CDERScorer cder_scorer_type;
      typedef cder_scorer_type::weights_type weights_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      typedef cicada::Matcher  matcher_type;
      
      struct Score
      {
	double score;
	
	double insertion;
	double deletion;
	double substitution;
	double jump;
	
	Score() : score(0), insertion(0), deletion(0), substitution(0), jump(0) {}
	
	Score(const double& __score,
	      const double& __insertion,
	      const double& __deletion,
	      const double& __substitution,
	      const double& __jump)
	  : score(__score), insertion(__insertion), deletion(__deletion), substitution(__substitution), jump(__jump) {}
      };

      typedef Score value_type;
      
      CDERScorerImpl() { }
      CDERScorerImpl(const sentence_type& __ref) : ref(__ref) {}
      
      value_type operator()(const sentence_type& sentence, const weights_type& weights, const matcher_type* matcher)
      {
	const sentence_type& a = ref;
	const sentence_type& b = sentence;

	const size_t a_size = a.size();
	const size_t b_size = b.size();
	
	if (a_size == 0) return value_type(b_size * weights.insertion, b_size * weights.insertion, 0, 0, 0);
	if (b_size == 0) return value_type(a_size * weights.deletion, 0, a_size * weights.deletion, 0, 0);
	
	std::vector<value_type, std::allocator<value_type> > curr(b_size + 1, value_type());
	std::vector<value_type, std::allocator<value_type> > prev(b_size + 1, value_type());
	
	for (size_t j = 1; j <= b_size; ++ j) {
	  prev[j].score     = prev[j - 1].score + weights.insertion;
	  prev[j].insertion = prev[j - 1].insertion + weights.insertion;
	}
	
	for (size_t i = 1; i <= a_size; ++ i) {
	  curr[0].score    = prev[0].score + weights.deletion;
	  curr[0].deletion = prev[0].deletion + weights.deletion;
	  
	  size_t pos_min = 0;
	  double score_min = curr[0].score;
	  
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
	    
	    if (curr[j].score < score_min) {
	      pos_min = j;
	      score_min = curr[j].score;
	    }
	  }
	  
	  const double score_jump = score_min + weights.jump;
	  
	  for (size_t j = 0; j <= b_size; ++ j)
	    if (score_jump < curr[j].score) {
	      curr[j] = curr[pos_min];
	      curr[j].score += weights.jump;
	      curr[j].jump  += weights.jump;
	    }
	  
	  std::swap(prev, curr);
	}
	
	return prev[b_size];
      }
      
    private:
      sentence_type ref;
    };
   
    CDERScorer::CDERScorer(const CDERScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	weights(x.weights),
	matcher(0)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      if (x.matcher)
	matcher = &matcher_type::create(x.matcher->algorithm());
    }
    
    CDERScorer::~CDERScorer()
    {
      clear();
    }
    
    CDERScorer& CDERScorer::operator=(const CDERScorer& x)
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
    
    void CDERScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void CDERScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    
    CDERScorer::score_ptr_type CDERScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double score_best = std::numeric_limits<double>::infinity();

      std::auto_ptr<CDER> cder(new CDER());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence, weights, matcher);
	
	if (value.score < score_best) {
	  score_best = value.score;
	  
	  cder->insertion    = value.insertion;
	  cder->deletion     = value.deletion;
	  cder->substitution = value.substitution;
	  cder->jump         = value.jump;
	}
	
	cder->references += evaluator.ref.size();
      }
      
      if (! impl.empty())
	cder->references /= impl.size();

      return score_ptr_type(cder.release());
    }
  };
};
