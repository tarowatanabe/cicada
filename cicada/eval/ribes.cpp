//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>

#include "ribes.hpp"

#include "utils/mathop.hpp"
#include "utils/bit_vector.hpp"
#include "utils/unordered_map.hpp"
#include "utils/hashmurmur3.hpp"

#include <boost/functional/hash.hpp>
#include <boost/math/special_functions/binomial.hpp>

namespace cicada
{
  namespace eval
  {

    std::string Ribes::description() const
    {
      std::ostringstream stream;
      stream << "ribes: " << score() << ' ' << distance << " norm: " << norm;
      return stream.str();
    }
    
    std::string Ribes::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"ribes\",";
      stream << "\"distance\":";
      stream << escaper(distance);
      stream << ',';
      stream << "\"norm\":";
      stream << escaper(norm);
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
			 >> qi::lit("\"norm\"") >> qi::lit(':') >> double_value
			 >> qi::lit('}'));
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, ribes_parsed_type(), space_type> ribes_parsed;
    };

    Score::score_ptr_type Ribes::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      ribes_parser<iterator_type> parser;
      ribes_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<Ribes> ribes(new Ribes());
      ribes->distance = boost::fusion::get<0>(parsed);
      ribes->norm     = boost::fusion::get<1>(parsed);
      
      return score_ptr_type(ribes.release());
    }
    
    Score::score_ptr_type Ribes::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

    class RibesScorerImpl
    {
    private:
      friend class RibesScorer;
      
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef RibesScorer ribes_scorer_type;
      
      typedef ribes_scorer_type::count_type count_type;
      typedef ribes_scorer_type::weight_type weight_type;
      
      typedef cicada::Sentence                sentence_type;
      typedef cicada::Symbol                  word_type;
      typedef word_type                       unigram_type;
      typedef std::pair<word_type, word_type> bigram_type;
      
      struct Score
      {
	double distance;
	double precision;
	double penalty;
	
	Score() : distance(0.0), precision(0.0), penalty(0.0) {}
	Score(const double& __distance,
	      const double& __precision,
	      const double& __penalty)
	  : distance(__distance), precision(__precision), penalty(__penalty) {}
      };
      
      typedef Score value_type;
      
      RibesScorerImpl(const sentence_type& __ref) : ref(__ref) { collect_stats(ref, ref_unigrams); }
      
      typedef std::vector<int, std::allocator<int> > alignment_type;
      typedef utils::bit_vector<4096> aligned_type;

      typedef std::vector<size_type, std::allocator<size_type> > matched_set_type;
      
      typedef utils::unordered_map<unigram_type, matched_set_type, boost::hash<unigram_type>, std::equal_to<unigram_type>,
				   std::allocator<std::pair<const unigram_type, matched_set_type> > >::type unigram_count_type;
      typedef utils::unordered_map<bigram_type, matched_set_type, utils::hashmurmur3<size_t>, std::equal_to<bigram_type>,
				   std::allocator<std::pair<const bigram_type, matched_set_type> > >::type bigram_count_type;
      
      void collect_stats(const sentence_type& sentence, unigram_count_type& unigrams, bigram_count_type& bigrams) const
      {
	unigrams.clear();
	bigrams.clear();
	
	sentence_type::const_iterator siter_begin = sentence.begin();
	sentence_type::const_iterator siter_end  = sentence.end();
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  unigrams[*siter].push_back(siter - siter_begin);
	
	for (sentence_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter)
	  bigrams[bigram_type(*(siter - 1), *siter)].push_back(siter - siter_begin);
      }
 
      void collect_stats(const sentence_type& sentence, unigram_count_type& unigrams) const
      {
	unigrams.clear();
	
	sentence_type::const_iterator siter_begin = sentence.begin();
	sentence_type::const_iterator siter_end  = sentence.end();
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  unigrams[*siter].push_back(siter - siter_begin);
      }
      
      value_type operator()(const sentence_type& hyp, const bool spearman) const
      {
	// no score
	if (hyp.empty())
	  return value_type();
	
	alignment_type& align = const_cast<alignment_type&>(align_impl);
	aligned_type&   aligned = const_cast<aligned_type&>(aligned_impl);
	
	unigram_count_type& hyp_unigrams = const_cast<unigram_count_type&>(hyp_unigrams_impl);
	//bigram_count_type&  hyp_bigrams  = const_cast<bigram_count_type&>(hyp_bigrams_impl);

	collect_stats(hyp, hyp_unigrams);
	
	align.clear();
	aligned.clear();
	
	const size_type hyp_size = hyp.size();
	const size_type ref_size = ref.size();
	
	const double bp = std::min(1.0, std::exp(1.0 - double(ref_size) / hyp_size));

	matched_set_type ref_matched_left;
	matched_set_type hyp_matched_left;
	matched_set_type ref_matched_right;
	matched_set_type hyp_matched_right;
	matched_set_type ref_matched_next;
	matched_set_type hyp_matched_next;
	
	for (size_type i = 0; i != hyp_size; ++ i) {
	  unigram_count_type::const_iterator riter = ref_unigrams.find(hyp[i]);

	  if (riter == ref_unigrams.end()) continue;
	  
	  const matched_set_type& ref_matched = riter->second;
	  const matched_set_type& hyp_matched = hyp_unigrams.find(hyp[i])->second;
	  
	  if (ref_matched.size() == 1 && hyp_matched.size() == 1) {
	    aligned.set(ref_matched.front(), true);
	    align.push_back(ref_matched.front());

	    //std::cerr << "matched: " << hyp[i] << " i = " << i << " j = " << align.back() << std::endl;
	  } else {
	    // we will try matching ngrams from lower order
	    
	    ref_matched_left  = ref_matched;
	    ref_matched_right = ref_matched;
	    hyp_matched_left  = hyp_matched;
	    hyp_matched_right = hyp_matched;
	    
	    for (size_type order = 1; order <= utils::bithack::max(i, hyp_size - i); ++ order) {
	      // try matching with the ngram to the left
	      if (i >= order && ! ref_matched_left.empty() && ! hyp_matched_left.empty()) {
		ref_matched_next.clear();
		hyp_matched_next.clear();
		
		matched_set_type::const_iterator riter_end = ref_matched_left.end();
		for (matched_set_type::const_iterator riter = ref_matched_left.begin(); riter != riter_end; ++ riter)
		  if (*riter >= order && ref[*riter - order] == hyp[i - order])
		    ref_matched_next.push_back(*riter);
		
		matched_set_type::const_iterator hiter_end = hyp_matched_left.end();
		for (matched_set_type::const_iterator hiter = hyp_matched_left.begin(); hiter != hiter_end; ++ hiter)
		  if (*hiter >= order && hyp[*hiter - order] == hyp[i - order])
		    hyp_matched_next.push_back(*hiter);

		if (ref_matched_next.size() == 1 && hyp_matched_next.size() == 1) {
		  aligned.set(ref_matched_next.front(), true);
		  align.push_back(ref_matched_next.front());
		  
		  //std::cerr << "matched: " << hyp[i] << " i = " << i << " j = " << align.back() << " left-order = " << order << std::endl;
		  break;
		}
		
		ref_matched_left.swap(ref_matched_next);
		hyp_matched_left.swap(hyp_matched_next);
	      }
	      
	      // try matching with the ngram to the right
	      if (i + order < hyp.size() && ! ref_matched_right.empty() && ! hyp_matched_right.empty()) {
		ref_matched_next.clear();
		hyp_matched_next.clear();
		
		matched_set_type::const_iterator riter_end = ref_matched_right.end();
		for (matched_set_type::const_iterator riter = ref_matched_right.begin(); riter != riter_end; ++ riter)
		  if (*riter + order < ref.size() && ref[*riter + order] == hyp[i + order])
		    ref_matched_next.push_back(*riter);
		
		matched_set_type::const_iterator hiter_end = hyp_matched_right.end();
		for (matched_set_type::const_iterator hiter = hyp_matched_right.begin(); hiter != hiter_end; ++ hiter)
		  if (*hiter + order < hyp.size() && hyp[*hiter + order] == hyp[i + order])
		    hyp_matched_next.push_back(*hiter);
		
		if (ref_matched_next.size() == 1 && hyp_matched_next.size() == 1) {
		  aligned.set(ref_matched_next.front(), true);
		  align.push_back(ref_matched_next.front());
		  
		  //std::cerr << "matched: " << hyp[i] << " i = " << i << " j = " << align.back() << " right-order = " << order << std::endl;
		  break;
		}
		
		ref_matched_right.swap(ref_matched_next);
		hyp_matched_right.swap(hyp_matched_next);
	      }
	      
	      // early termination...
	      if ((ref_matched_left.empty() || hyp_matched_left.empty())
		  && (ref_matched_right.empty() || hyp_matched_right.empty())) break;
	    }
	  }
	}
	
#if 0
	// this is an old RIBES which can match upto bigrams
	for (size_type i = 0; i != hyp_size; ++ i)
	  for (size_type j = 0; j != ref_size; ++ j)
	    if (hyp[i] == ref[j]) {
	      
	      // check unique unigram
	      if (ref_unigrams.find(ref[j])->second.size() == 1
		  && hyp_unigrams.find(hyp[i])->second.size() == 1) {
		aligned.set(j, true);
		align.push_back(j);
		continue;
	      }
	      
	      // bigram with the next word
	      if (i + 1 != hyp_size && j + 1 != ref_size && hyp[i + 1] == ref[j + 1]
		  && ref_bigrams.find(bigram_type(ref[j], ref[j + 1]))->second.size() == 1
		  && hyp_bigrams.find(bigram_type(hyp[i], hyp[i + 1]))->second.size() == 1) {
		aligned.set(j, true);
		align.push_back(j);
		continue;
	      }
	      
	      // bigram with the previos word
	      if (i != 0 && j != 0 && hyp[i - 1] == ref[j - 1]
		  && ref_bigrams.find(bigram_type(ref[j - 1], ref[j]))->second.size() == 1
		  && hyp_bigrams.find(bigram_type(hyp[i - 1], hyp[i]))->second.size() == 1) {
		aligned.set(j, true);
		align.push_back(j);
		continue;
	      }
	    }
#endif
	
	if (align.size() == 1 && ref_size == 1)
	  return value_type(1.0, 1.0 / hyp.size(), bp);
	else if (align.size() <= 1)
	  return value_type(0.0, 0.0, bp);
	
	value_type value(0.0, static_cast<double>(align.size()) / hyp.size(), bp);
	
	// we use binomial coefficient found in boost.math
	if (spearman) {
	  //
	  // after filling aligned, then, we need to re-number indicated by the "aligned" vector...
	  //
	  // we simply check "rank" of "alinged" as our new index...
	  //
	  
	  alignment_type::iterator aiter_end = align.end();
	  for (alignment_type::iterator aiter = align.begin(); aiter != aiter_end; ++ aiter)
	    *aiter = aligned.rank(*aiter, true) - 1;

	  // Spearman
	  double distance = 0.0;
	  for (size_type i = 0; i != align.size(); ++ i)
	    distance += (align[i] - i) * (align[i] - i);
	  
	  const double rho = 1.0 - distance / boost::math::binomial_coefficient<double>(align.size() + 1, 3);
	  
	  value.distance = (rho + 1.0) * 0.5;
	} else {
	  // Kendall
	  size_type num_increasing = 0;
	  for (size_type i = 0; i != align.size() - 1; ++ i)
	    for (size_type j = i + 1; j != align.size(); ++ j)
	      num_increasing += (align[j] > align[i]);
	  
	  const double tau = 2.0 * num_increasing / boost::math::binomial_coefficient<double>(align.size(), 2) - 1.0;
	  
	  value.distance = (tau + 1.0) * 0.5;
	}
	
	return value;
      }
      
    private:
      sentence_type ref;
      unigram_count_type ref_unigrams;
      //bigram_count_type  ref_bigrams;

      unigram_count_type hyp_unigrams_impl;
      //bigram_count_type  hyp_bigrams_impl;
      
      alignment_type align_impl;
      aligned_type   aligned_impl;
    };
   
    RibesScorer::RibesScorer(const RibesScorer& x)
      : Scorer(static_cast<const Scorer&>(*this)),
	alpha(x.alpha),
	beta(x.beta),
	spearman(x.spearman)
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    RibesScorer::~RibesScorer()
    {
      clear();
    }
    
    RibesScorer& RibesScorer::operator=(const RibesScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);
      
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      alpha = x.alpha;
      beta = x.beta;
      spearman = x.spearman;
      
      return *this;
    }
    
    void RibesScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void RibesScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      impl.push_back(new impl_type(sentence));
    }
    
    RibesScorer::score_ptr_type RibesScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      std::auto_ptr<Ribes> ribes(new Ribes());
      
      if (! impl.empty()) {
	ribes->distance = - std::numeric_limits<double>::infinity();
	
	for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	  const impl_type& evaluator = *(*iter);
	  
	  const impl_type::value_type value = evaluator(sentence, spearman);
	  
	  ribes->distance = std::max(ribes->distance, (value.distance
						       * utils::mathop::pow(value.precision, alpha)
						       * utils::mathop::pow(value.penalty, beta)));
	}
	
	ribes->norm = 1;
      }
      
      return score_ptr_type(ribes.release());
    }
  };
};

