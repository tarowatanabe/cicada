// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__BLEU__HPP__
#define __CICADA__EVAL__BLEU__HPP__ 1

#include <cmath>
#include <cfloat>

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <utility>

#include <cicada/eval/score.hpp>
#include <cicada/symbol_vector.hpp>

#include <utils/trie_compact.hpp>
#include <utils/simple_vector.hpp>
#include <utils/bithack.hpp>

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace eval
  {
    class Bleu : public Score
    {
    public:
      typedef double count_type;
      
      typedef utils::simple_vector<count_type, std::allocator<count_type> > ngram_counts_type;

    public:      
      Bleu(const int order)
	: ngrams_hypothesis(order, 0),  ngrams_matched(order, 0),
	  length_reference(0), length_hypothesis(0) {}
      
      double score() const
      {
	if (ngrams_matched.empty() || ngrams_matched[0] == 0.0) return 0.0;
	
	const double penalty = std::min(1.0 - length_reference / length_hypothesis, 0.0);
	
	double smooth = 0.5;
	double score = 0.0;
	int norm = 0;
	
	for (size_t n = 0; n != ngrams_hypothesis.size(); ++ n) {
	  if (ngrams_hypothesis[n] > 0) {
	    score += std::log((n < ngrams_matched.size() && ngrams_matched[n] > 0 ? ngrams_matched[n] : smooth) / ngrams_hypothesis[n]);
	    ++ norm;
	  }
	  
	  smooth *= 0.5;
	}
	
	return std::exp(score / norm + penalty);
      }
      
      double loss() const { return 1.0 - score(); }
      
      double reward() const { return score(); }
      
      bool error_metric() const { return false; }

      bool equal(const score_type& score) const
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	return (ngrams_hypothesis == rhs->ngrams_hypothesis
		&& ngrams_matched == rhs->ngrams_matched
		&& length_reference == rhs->length_reference
		&& length_hypothesis == rhs->length_hypothesis);
      }

      void assign(const score_type& score)
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	ngrams_hypothesis  = rhs->ngrams_hypothesis;
	ngrams_matched = rhs->ngrams_matched;
	
	length_reference  = rhs->length_reference;
	length_hypothesis = rhs->length_hypothesis;
      }
      
      void plus_equal(const score_type& score)
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	std::transform(ngrams_hypothesis.begin(),  ngrams_hypothesis.end(),  rhs->ngrams_hypothesis.begin(),  ngrams_hypothesis.begin(),  std::plus<count_type>());
	std::transform(ngrams_matched.begin(), ngrams_matched.end(), rhs->ngrams_matched.begin(), ngrams_matched.begin(), std::plus<count_type>());
	
	length_reference  += rhs->length_reference;
	length_hypothesis += rhs->length_hypothesis;
      }
      
      void minus_equal(const score_type& score)
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	std::transform(ngrams_hypothesis.begin(),  ngrams_hypothesis.end(),  rhs->ngrams_hypothesis.begin(),  ngrams_hypothesis.begin(),  std::minus<count_type>());
	std::transform(ngrams_matched.begin(), ngrams_matched.end(), rhs->ngrams_matched.begin(), ngrams_matched.begin(), std::minus<count_type>());
	
	length_reference  -= rhs->length_reference;
	length_hypothesis -= rhs->length_hypothesis;
      }
      
      void multiplies_equal(const double& scale)
      {
	std::transform(ngrams_hypothesis.begin(),  ngrams_hypothesis.end(),  ngrams_hypothesis.begin(),  std::bind2nd(std::multiplies<count_type>(), scale));
	std::transform(ngrams_matched.begin(), ngrams_matched.end(), ngrams_matched.begin(), std::bind2nd(std::multiplies<count_type>(), scale));
	
	length_reference  *= scale;
	length_hypothesis *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	std::transform(ngrams_hypothesis.begin(),  ngrams_hypothesis.end(),  ngrams_hypothesis.begin(),  std::bind2nd(std::divides<count_type>(), scale));
	std::transform(ngrams_matched.begin(), ngrams_matched.end(), ngrams_matched.begin(), std::bind2nd(std::divides<count_type>(), scale));
	
	length_reference  /= scale;
	length_hypothesis /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new Bleu(ngrams_matched.size()));
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new Bleu(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);

      ngram_counts_type ngrams_hypothesis;
      ngram_counts_type ngrams_matched;
      count_type        length_reference;
      count_type        length_hypothesis;
    };
    

    class BleuScorer : public Scorer
    {
    private:
      typedef double count_type;
      
      typedef std::allocator<std::pair<const word_type, count_type> > ngram_allocator_type;
      typedef utils::trie_compact<word_type, count_type,
				  utils::unassigned<word_type>, 
				  boost::hash<word_type>, std::equal_to<word_type>,
				  ngram_allocator_type> ngram_set_type;
      
      typedef std::vector<int, std::allocator<int> > size_set_type;

    public:
      BleuScorer(int __order = 4) : ngrams(), sizes(), order(__order) { }

      bool error_metric() const { return false; }

      scorer_ptr_type clone() const { return scorer_ptr_type(new BleuScorer(*this)); }

      void clear()
      {
	ngrams.clear();
	sizes.clear();
      }
      
      void insert(const sentence_type& __sentence);
      
      score_ptr_type score(const sentence_type& __sentence) const;
      
      // bleu specific interface...
      count_type reference_length(const double& length) const;
      count_type find(const SymbolVector& ngram) const;

    private:
      ngram_set_type ngrams;
      size_set_type  sizes;
      
      int            order;
    };
  };
};

#endif
