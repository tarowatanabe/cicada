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

#include "parseval.hpp"

#include <boost/functional/hash.hpp>

namespace cicada
{
  namespace eval
  {
    std::string Parseval::description() const
    {
      std::ostringstream stream;
      stream << "parseval: " << score()
	     << ' ' << matched << '|' << test << '|' << reference;
      return stream.str();
    }
    
    std::string Parseval::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"parseval\",";
      stream << "\"brackets\":[";
      stream << '\"' << escaper(matched) << "\",";
      stream << '\"' << escaper(test) << "\",";
      stream << '\"' << escaper(reference) << '\"';
      stream << "]}";
      return stream.str();
    }
    
    typedef boost::fusion::tuple<double, double, double> parseval_parsed_type;

    template <typename Iterator>
    struct parseval_parser : boost::spirit::qi::grammar<Iterator, parseval_parsed_type(), boost::spirit::standard::space_type>
    {
      parseval_parser() : parseval_parser::base_type(parseval_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	parseval_parsed %= (qi::lit('{')
			    >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"parseval\"") >> qi::lit(',')
			    >> qi::lit("\"brackets\"") >> qi::lit(':')
			    >> qi::lit('[')
			    >> double_value >> qi::lit(',')
			    >> double_value >> qi::lit(',')
			    >> double_value
			    >> qi::lit(']')
			    >> qi::lit('}'));
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, parseval_parsed_type(), space_type> parseval_parsed;
    };
    
    Score::score_ptr_type Parseval::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      parseval_parser<iterator_type> parser;
      parseval_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<Parseval> parseval(new Parseval());
      parseval->matched   = boost::fusion::get<0>(parsed);
      parseval->test      = boost::fusion::get<1>(parsed);
      parseval->reference = boost::fusion::get<2>(parsed);
      
      return score_ptr_type(parseval.release());
    }

    Score::score_ptr_type Parseval::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

    class ParsevalScorerImpl
    {
    private:
      friend class ParsevalScorer;
    public:
      typedef ParsevalScorer parseval_scorer_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;

      typedef parseval_scorer_type::ignored_type ignored_type;

      typedef std::multiset<word_type, std::less<word_type>, std::allocator<word_type> > word_set_type;
      
      struct Score
      {
	double score;
	
	double matched;
	double test;
	double reference;
	
	Score() : score(0), matched(0), test(0), reference(0) {}
	
	Score(const double& __score,
	      const double& __matched,
	      const double& __test,
	      const double& __reference)
	  : score(__score), matched(__matched), test(__test), reference(__reference) {}
      };

      typedef Score value_type;
      
      ParsevalScorerImpl() {}
      ParsevalScorerImpl(const sentence_type& __ref) : ref(__ref) { }
      
      value_type operator()(const sentence_type& sentence, const ignored_type& ignored)
      {
	value_type value;
	

	if (ignored.empty()) {
	  word_set_type test(sentence.begin(), sentence.end());
	  word_set_type oracle(ref.begin(), ref.end());
	  
	  evaluate(test, oracle, value);
	} else {
	  // it is specially handled for parseval!
	  word_set_type test;
	  word_set_type oracle;

	  sentence_type::const_iterator siter_end = sentence.end();
	  for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	    cicada::SpanVector::span_type span(*siter);
	    if (ignored.find(span.label) == ignored.end())
	      test.insert(*siter);
	  }
	  
	  sentence_type::const_iterator riter_end = ref.end();
	  for (sentence_type::const_iterator riter = ref.begin(); riter != riter_end; ++ riter) {
	    cicada::SpanVector::span_type span(*riter);
	    if (ignored.find(span.label) == ignored.end())
	      oracle.insert(*riter);
	  }

	  evaluate(test, oracle, value);
	}

	return value;
      }

      
      struct counter_iterator
      {
	counter_iterator(int& __counter) : counter(__counter) {}

	counter_iterator& operator*() { return *this; }

	template <typename Tp>
	counter_iterator& operator=(const Tp& x)
	{
	  ++ counter;
	  return *this;
	}
	
	counter_iterator& operator++() { return *this; }
	counter_iterator operator++(int) { return *this; }

	int& counter;
      };
      
      void evaluate(const word_set_type& test, const word_set_type& oracle, value_type& value)
      {
	int counter = 0;
	
	std::set_intersection(test.begin(), test.end(), oracle.begin(), oracle.end(), counter_iterator(counter));
	
	value.matched   = counter;
	value.test      = test.size();
	value.reference = oracle.size();
	
	// compute f-measure!
	value.score = 2.0 * value.matched / (value.test + value.reference);
      }
      
    private:
      sentence_type ref;
    };

    
    ParsevalScorer::ParsevalScorer(const ParsevalScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	ignored(x.ignored)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    ParsevalScorer::~ParsevalScorer()
    {
      clear();
    }
    
    ParsevalScorer& ParsevalScorer::operator=(const ParsevalScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      ignored = x.ignored;
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      return *this;
    }
    
    void ParsevalScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
      ignored.clear();
    }
    
    void ParsevalScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      impl.push_back(new impl_type(sentence));
    }
    
    ParsevalScorer::score_ptr_type ParsevalScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      double score_best = - std::numeric_limits<double>::infinity();

      std::auto_ptr<Parseval> parseval(new Parseval());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	impl_type::value_type value = evaluator(sentence, ignored);
	
	if (value.score > score_best) {
	  score_best = value.score;
	  
	  parseval->matched   = value.matched;
	  parseval->test      = value.test;
	  parseval->reference = value.reference;
	}
      }
      
      return score_ptr_type(parseval.release());
    }
  };
};
