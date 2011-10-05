//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>
#include <set>
#include <algorithm>

#include <cicada/span_vector.hpp>

#include "depeval.hpp"

#include <boost/functional/hash.hpp>

#include "utils/bithack.hpp"

namespace cicada
{
  namespace eval
  {
    std::string Depeval::description() const
    {
      std::ostringstream stream;
      stream << "depeval: " << score()
	     << ' ' << matched << '|' << total;
      return stream.str();
    }
    
    std::string Depeval::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"depeval\",";
      stream << "\"dependencies\":[";
      stream << escaper(matched) << ',' << escaper(total);
      stream << "]}";
      return stream.str();
    }
    
    typedef boost::fusion::tuple<double, double> depeval_parsed_type;

    template <typename Iterator>
    struct depeval_parser : boost::spirit::qi::grammar<Iterator, depeval_parsed_type(), boost::spirit::standard::space_type>
    {
      depeval_parser() : depeval_parser::base_type(depeval_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	depeval_parsed %= (qi::lit('{')
			   >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"depeval\"") >> qi::lit(',')
			   >> qi::lit("\"dependencies\"") >> qi::lit(':')
			   >> qi::lit('[')
			   >> double_value >> qi::lit(',')
			   >> double_value
			   >> qi::lit(']')
			   >> qi::lit('}'));
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, depeval_parsed_type(), space_type> depeval_parsed;
    };
    
    Score::score_ptr_type Depeval::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      depeval_parser<iterator_type> parser;
      depeval_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<Depeval> depeval(new Depeval());
      depeval->matched   = boost::fusion::get<0>(parsed);
      depeval->total     = boost::fusion::get<1>(parsed);
      
      return score_ptr_type(depeval.release());
    }

    Score::score_ptr_type Depeval::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

    class DepevalScorerImpl
    {
    private:
      friend class DepevalScorer;
    public:
      typedef DepevalScorer depeval_scorer_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;

      struct Score
      {
	double score;
	double matched;
	double total;
	
	Score() : score(0), matched(0), total(0) {}
	
	Score(const double& __score,
	      const double& __matched,
	      const double& __total)
	  : score(__score), matched(__matched), total(__total) {}
      };

      typedef Score value_type;
      
      DepevalScorerImpl() {}
      DepevalScorerImpl(const sentence_type& __ref) : ref(__ref) { }
      
      value_type operator()(const sentence_type& sentence)
      {
	const size_t size = utils::bithack::min(sentence.size(), ref.size());
	size_t matched = 0;
	for (size_t i = 0; i != size; ++ i)
	  matched += sentence[i] == ref[i];
	
	return value_type(double(matched) / double(sentence.size()), matched, sentence.size());
      }
      
    private:
      sentence_type ref;
    };

    
    DepevalScorer::DepevalScorer(const DepevalScorer& x)
      : Scorer(static_cast<const Scorer&>(*this))
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    DepevalScorer::~DepevalScorer()
    {
      clear();
    }
    
    DepevalScorer& DepevalScorer::operator=(const DepevalScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      return *this;
    }
    
    void DepevalScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void DepevalScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      impl.push_back(new impl_type(sentence));
    }
    
    DepevalScorer::score_ptr_type DepevalScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double score_best = - std::numeric_limits<double>::infinity();
      
      std::auto_ptr<Depeval> depeval(new Depeval());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence);
	
	if (value.score > score_best) {
	  score_best = value.score;
	  
	  depeval->matched = value.matched;
	  depeval->total   = value.total;
	}
      }
      
      return score_ptr_type(depeval.release());
    }
  };
};
