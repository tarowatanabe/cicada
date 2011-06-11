//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>

#include "ribes.hpp"

#include "utils/mathop.hpp"
#include "utils/bit_vector.hpp"

#include <boost/functional/hash.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace cicada
{
  namespace eval
  {

    std::string RIBES::description() const
    {
      std::ostringstream stream;
      stream << "ribes: " << score() << ' ' << distance << " penalty: " << penalty;
      return stream.str();
    }
    
    std::string RIBES::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"ribes\",";
      stream << "\"distance\":";
      stream << escaper(distance);
      stream << ',';
      stream << "\"penalty\":";
      stream << escaper(penalty);
      stream << '}';
      return stream.str();
    }

    typedef boost::fusion::tuple<double, double> ribes_parsed_type;
    
    template <typename Iterator>
    struct ribes_parser : boost::spirit::qi::grammar<Iterator, ribes_parsed_type(), boost::spirit::standard::space_type>
    {
      ribes_parser() : ribes_parser::base_type(ribes_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	ribes_parsed %= (qi::lit('{')
			 >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"ribes\"") >> qi::lit(',')
			 >> qi::lit("\"distance\"") >> qi::lit(':') >> double_value >> ','
			 >> qi::lit("\"penalty\"") >> qi::lit(':') >> double_value
			 >> qi::lit('}'));
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, ribes_parsed_type(), space_type> ribes_parsed;
    };

    Score::score_ptr_type RIBES::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      ribes_parser<iterator_type> parser;
      ribes_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<RIBES> ribes(new RIBES());
      ribes->distance = boost::fusion::get<0>(parsed);
      ribes->penalty  = boost::fusion::get<1>(parsed);
      
      return score_ptr_type(ribes.release());
    }
    
    Score::score_ptr_type RIBES::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

    class RIBESScorerImpl
    {
    private:
      friend class RIBESScorer;
      
    public:
      typedef RIBESScorer ribes_scorer_type;
      
      typedef ribes_scorer_type::count_type count_type;
      typedef ribes_scorer_type::weight_type weight_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      
      struct Score
      {
	double distance;
	double penalty;
	
	Score() : distance(0.0), penalty(0.0) {}
	
	Score(const double& __distance,
	      const double& __penalty)
	  : distance(__distance), penalty(__penalty) {}
      };
      
      typedef Score value_type;
      
      RIBESScorerImpl() { }
      RIBESScorerImpl(const sentence_type& __ref) : ref(__ref) {}
      
      typedef std::vector<int, std::allocator<int> > alignment_type;
      typedef utils::bit_vector<4096> aligned_type;
      
      value_type operator()(const sentence_type& hyp, const weight_type& weight, const bool spearman) const
      {
	alignment_type& align = const_cast<alignment_type&>(align_impl);
	aligned_type&   aligned = const_cast<aligned_type&>(aligned_impl);
	
	align.clear();
	aligned.clear();
	
	for (size_t i = 0; i != hyp.size(); ++ i)
	  for (size_t j = 0; j != ref.size(); ++ j)
	    if (! aligned[j] && ref[j] == hyp[i])
	      if ((j + 1 < ref.size() && i + 1 < hyp.size() && ref[j + 1] == hyp[i + 1])
		  || (j != 0 && i != 0 && ref[j - 1] == hyp[i - 1])) {
		// how to handle boundary condition...?
		// see ribes code...
		aligned.set(j, true);
		align.push_back(j);
	      }

	if (align.size() <= 1)
	  return value_type(0.0, utils::mathop::pow(static_cast<double>(align.size()) / hyp.size(), weight));
	
	//
	// after filling aligned, then, we need to re-number indicated by the "aligned" vector...
	//
	// we simply check "rank" of "alinged" as our new index...
	//
	alignment_type::iterator aiter_end = align.end();
	for (alignment_type::iterator aiter = align.begin(); aiter != aiter_end; ++ aiter)
	  *aiter = aligned.rank(*aiter, true) - 1;
	
	value_type value;
	value.penalty = utils::mathop::pow(static_cast<double>(align.size()) / hyp.size(), weight);
	
	// we use binomial coefficient found in boost.math
	if (spearman) {
	  // Spearman
	  double distance = 0.0;
	  for (size_t i = 0; i != align.size(); ++ i)
	    distance += (align[i] - i) * (align[i] - i);
	  
	  const double rho = 1.0 - distance / boost::math::binomial_coefficient<double>(align.size() + 1, 3);
	  
	  value.distance = (rho + 1.0) * 0.5;
	} else {
	  // Kendall
	  size_t num_increasing = 0;
	  for (size_t i = 0; i != align.size() - 1; ++ i)
	    for (size_t j = i + 1; j != align.size(); ++ j)
	      num_increasing += (align[j] > align[i]);
	  
	  const double tau = 2.0 * num_increasing / boost::math::binomial_coefficient<double>(align.size(), 2) - 1.0;
	  
	  value.distance = (tau + 1.0) * 0.5;
	}
	
	return value;
      }
      
    private:
      sentence_type ref;
      
      alignment_type align_impl;
      aligned_type   aligned_impl;
    };
   
    RIBESScorer::RIBESScorer(const RIBESScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	weight(x.weight)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    RIBESScorer::~RIBESScorer()
    {
      clear();
    }
    
    RIBESScorer& RIBESScorer::operator=(const RIBESScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      weight = x.weight;
      
      return *this;
    }
    
    void RIBESScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void RIBESScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);

      impl.push_back(new impl_type(sentence));
    }
    
    RIBESScorer::score_ptr_type RIBESScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      std::auto_ptr<RIBES> ribes(new RIBES());
      
      if (! impl.empty()) {
	ribes->distance = - std::numeric_limits<double>::infinity();
	
	for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	  const impl_type& evaluator = *(*iter);
	  
	  const impl_type::value_type value = evaluator(sentence, weight, spearman);
	  
	  ribes->distance = std::max(ribes->distance, value.distance * value.penalty);
	}
	
	ribes->penalty = 1;
      }
      
      return score_ptr_type(ribes.release());
    }
  };
};

